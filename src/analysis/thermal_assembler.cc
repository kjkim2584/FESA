#include "analysis/thermal_assembler.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <new>

#include "adapters/openmp/thread_runtime.h"

namespace fesa {
namespace {
StatusCode Error(const char* text, Diagnostic* d,
                 StatusCode code = StatusCode::kInvalidArgument) {
  return ReportFailure(code,
                       code == StatusCode::kAllocationFailure
                           ? DiagnosticStage::kAllocation
                           : DiagnosticStage::kComputation,
                       text, d);
}
StatusCode OptionalVector(int count, Vector** out, Diagnostic* d) {
  *out = nullptr;
  return count > 0 ? Vector::Create(count, out, d) : StatusCode::kOk;
}
StatusCode MakeVector(int count, std::unique_ptr<Vector>* out, Diagnostic* d) {
  Vector* raw = nullptr;
  const auto status = Vector::Create(count, &raw, d);
  out->reset(raw);
  return status;
}
StatusCode MakeMatrix(int count, std::unique_ptr<Matrix>* out, Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(count, count, &raw, d);
  out->reset(raw);
  return status;
}
}  // namespace
ThermalAssembly::~ThermalAssembly() {
  delete rows_;
  delete columns_;
  delete matrix_indices_;
  delete globals_;
  delete element_offsets_;
  delete coo_offsets_;
  delete load_flux_offsets_;
  delete load_starts_;
  delete load_indices_;
  delete external_globals_;
  delete free_globals_;
  delete constrained_globals_;
  delete participating_globals_;
  delete values_;
  delete residual_;
  delete reaction_;
  delete element_flux_;
  delete external_flux_;
  delete free_residual_;
  delete increment_;
  delete global_increment_;
  delete constrained_scratch_;
  delete zero_;
  delete step_start_;
}
StatusCode ThermalAssembly::Create(const ThermalModel* model,
                                   const ThermalDofMap* map,
                                   const ThermalConditions* conditions,
                                   ThermalAssembly** out, Diagnostic* d) {
  return CreateForStep(model, map, conditions, nullptr, out, d);
}
StatusCode ThermalAssembly::CreateForStep(const ThermalModel* model,
                                          const ThermalDofMap* map,
                                          const ThermalConditions* conditions,
                                          const Vector* step_start,
                                          ThermalAssembly** out,
                                          Diagnostic* d) {
  if (!out) return Error("Null thermal assembly output", d);
  *out = nullptr;
  if (!model || !map || !conditions || model->DofMap() != map ||
      model->Conditions() != conditions ||
      model->SourceDomain() != conditions->SourceDomain())
    return Error("Thermal assembly binding mismatch", d);
  try {
    std::unique_ptr<ThermalAssembly> a(new ThermalAssembly);
    a->model_ = model;
    a->map_ = map;
    a->conditions_ = conditions;
    bool film = false;
    for (const auto* row : *conditions->ElementLoads())
      film = film || (row->origin->specification()->kind ==
                          ThermalConditionKind::kFilm ||
                      row->origin->specification()->kind ==
                          ThermalConditionKind::kRadiation);
    if (film) {
      if (!step_start || step_start->Size() != map->GlobalCount())
        return Error("Film assembly requires frozen step-start temperature", d);
      auto status = Vector::Create(map->GlobalCount(), &a->step_start_, d);
      if (status == StatusCode::kOk)
        status = a->step_start_->CopyFrom(step_start, d);
      if (status != StatusCode::kOk) return status;
    }
    std::vector<int>** lists[] = {&a->rows_,
                                  &a->columns_,
                                  &a->matrix_indices_,
                                  &a->globals_,
                                  &a->element_offsets_,
                                  &a->coo_offsets_,
                                  &a->load_flux_offsets_,
                                  &a->load_starts_,
                                  &a->load_indices_,
                                  &a->external_globals_,
                                  &a->free_globals_,
                                  &a->constrained_globals_,
                                  &a->participating_globals_};
    for (auto* list : lists) *list = new std::vector<int>(0);
    a->element_offsets_->push_back(0);
    a->coo_offsets_->push_back(0);
    const int limit = std::numeric_limits<int>::max();
    for (const auto* element : *model->Elements()) {
      const int nodes = element->node_count(), ns = element->dofs_per_node();
      if (nodes <= 0 || ns <= 0 || nodes > limit / ns)
        return Error("Thermal element layout overflow", d,
                     StatusCode::kInputError);
      const int n = nodes * ns;
      if (n > limit / n ||
          a->globals_->size() > static_cast<std::size_t>(limit - n))
        return Error("Thermal element index range", d, StatusCode::kInputError);
      std::vector<int> equations(0);
      for (int local = 0; local < n; ++local) {
        const int g =
            map->GlobalDof(element->GetNode(local / ns)->Id(), 11 + local % ns);
        if (g < 0)
          return Error("Thermal assembly element mapping", d,
                       StatusCode::kInputError);
        a->globals_->push_back(g);
        equations.push_back(map->Equation(g));
      }
      for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j) {
          if (equations[i] < 0 || equations[j] < 0) continue;
          if (a->rows_->size() == static_cast<std::size_t>(limit))
            return Error("Thermal COO index range", d, StatusCode::kInputError);
          a->rows_->push_back(std::min(equations[i], equations[j]));
          a->columns_->push_back(std::max(equations[i], equations[j]));
          a->matrix_indices_->push_back(i * n + j);
        }
      a->element_offsets_->push_back(static_cast<int>(a->globals_->size()));
      a->coo_offsets_->push_back(static_cast<int>(a->rows_->size()));
    }
    // Cache immutable row/element routing once; workers then visit only
    // their own rows, writing disjoint packed origin regions.
    std::map<int, int> owners;
    for (std::size_t i = 0; i < model->Elements()->size(); ++i)
      owners.emplace((*model->Elements())[i]->Id(), static_cast<int>(i));
    a->load_starts_->resize(model->Elements()->size() + 1, 0);
    std::vector<int> load_owners(0);
    a->load_flux_offsets_->push_back(static_cast<int>(a->globals_->size()));
    if (conditions->ElementLoads()->size() > static_cast<std::size_t>(limit))
      return Error("Thermal distributed row count range", d,
                   StatusCode::kInputError);
    for (const auto* row : *conditions->ElementLoads()) {
      const auto found = owners.find(row->element_id);
      if (found == owners.end())
        return Error("Distributed thermal element not active", d);
      const auto* element = (*model->Elements())[found->second];
      if (a->globals_->size() >
          static_cast<std::size_t>(limit - element->dof_count()))
        return Error("Thermal distributed origin index range", d,
                     StatusCode::kInputError);
      for (int local = 0; local < element->node_count(); ++local)
        for (int label = 11; label < 11 + element->dofs_per_node(); ++label) {
          const int global =
              map->GlobalDof(element->GetNode(local)->Id(), label);
          if (global < 0) return Error("Distributed thermal DOF not active", d);
          a->globals_->push_back(global);
        }
      a->load_flux_offsets_->push_back(static_cast<int>(a->globals_->size()));
      ++(*a->load_starts_)[found->second + 1];
      load_owners.push_back(found->second);
    }
    for (std::size_t i = 1; i < a->load_starts_->size(); ++i)
      (*a->load_starts_)[i] += (*a->load_starts_)[i - 1];
    a->load_indices_->resize(load_owners.size());
    std::vector<int> cursors(a->load_starts_->begin(), a->load_starts_->end());
    for (std::size_t i = 0; i < load_owners.size(); ++i)
      (*a->load_indices_)[cursors[load_owners[i]]++] = static_cast<int>(i);
    if (conditions->Fluxes()->size() > static_cast<std::size_t>(limit))
      return Error("Thermal external origin index range", d,
                   StatusCode::kInputError);
    for (const auto* row : *conditions->Fluxes())
      a->external_globals_->push_back(map->GlobalDof(row->node_id, row->label));
    for (int eq = 0; eq < map->FreeCount(); ++eq)
      a->free_globals_->push_back(map->GlobalDofForEquation(eq));
    for (int row = 0; row < map->NodeCount(); ++row) {
      const int id = map->NodeId(row);
      if (!map->Participates(id)) continue;
      for (int label = 11; label < 11 + map->NodeDofCount(id); ++label) {
        const int g = map->GlobalDof(id, label);
        a->participating_globals_->push_back(g);
        if (map->Equation(g) == -1) a->constrained_globals_->push_back(g);
      }
    }
    Vector** globals[] = {&a->residual_, &a->reaction_, &a->global_increment_,
                          &a->zero_};
    for (auto* target : globals) {
      const auto status = Vector::Create(map->GlobalCount(), target, d);
      if (status != StatusCode::kOk) return status;
    }
    auto status =
        OptionalVector(static_cast<int>(a->rows_->size()), &a->values_, d);
    if (status == StatusCode::kOk)
      status = OptionalVector(static_cast<int>(a->globals_->size()),
                              &a->element_flux_, d);
    if (status == StatusCode::kOk)
      status = OptionalVector(static_cast<int>(a->external_globals_->size()),
                              &a->external_flux_, d);
    if (status == StatusCode::kOk)
      status = OptionalVector(map->FreeCount(), &a->free_residual_, d);
    if (status == StatusCode::kOk)
      status =
          OptionalVector(static_cast<int>(a->participating_globals_->size()),
                         &a->increment_, d);
    if (status == StatusCode::kOk)
      status = OptionalVector(static_cast<int>(a->constrained_globals_->size()),
                              &a->constrained_scratch_, d);
    if (status != StatusCode::kOk) return status;
    *out = a.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return Error("Thermal assembly allocation", d,
                 StatusCode::kAllocationFailure);
  }
}
bool ThermalAssembly::IsValid() const { return valid_; }
const ThermalDofMap* ThermalAssembly::DofMap() const { return map_; }
const std::vector<int>* ThermalAssembly::Rows() const { return rows_; }
const std::vector<int>* ThermalAssembly::Columns() const { return columns_; }
const Vector* ThermalAssembly::Values() const { return values_; }
const Vector* ThermalAssembly::Residual() const { return residual_; }
const Vector* ThermalAssembly::Reaction() const { return reaction_; }
const Vector* ThermalAssembly::ElementFluxes() const { return element_flux_; }
const Vector* ThermalAssembly::ExternalFluxes() const { return external_flux_; }
const Vector* ThermalAssembly::FreeResidual() const { return free_residual_; }
const Vector* ThermalAssembly::ParticipatingIncrement() const {
  return increment_;
}
struct ThermalAssembler::Context {
  const ThermalState* state = nullptr;  ///< Borrowed until join.
  ThermalAssembly* out = nullptr;       ///< Exclusive disjoint element regions.
  double coefficient = 0;
  double start = 0, end = 0;
  bool extrapolated_start = false;  ///< D13 evaluation, never committed state.
};
StatusCode ThermalAssembler::ElementTask(int item, int, void* data,
                                         Diagnostic* d) {
  auto* context = static_cast<Context*>(data);
  auto* a = context->out;
  const auto* element = (*a->model_->Elements())[item];
  if (d) d->element_id = element->Id();
  try {
    const int offset = (*a->element_offsets_)[item];
    const int n = (*a->element_offsets_)[item + 1] - offset;
    const int* globals = a->globals_->data() + offset;
    std::vector<int> predicted_free(0);
    if (context->extrapolated_start) {
      predicted_free.resize(n);
      for (int i = 0; i < n; ++i)
        predicted_free[i] = a->map_->Equation(globals[i]) >= 0 ? 1 : 0;
    }
    std::vector<int> identity(n);
    for (int i = 0; i < n; ++i) identity[i] = i;
    std::unique_ptr<Matrix> k;
    std::unique_ptr<Vector> theta, internal;
    auto status = MakeMatrix(n, &k, d);
    if (status == StatusCode::kOk) status = MakeVector(n, &theta, d);
    if (status == StatusCode::kOk) status = MakeVector(n, &internal, d);
    if (status != StatusCode::kOk) return status;
    status = theta->GatherFrom(context->state->CandidateTemperature(), globals,
                               n, 0, d);
    if (status == StatusCode::kOk)
      status = element->ComputeConduction(theta.get(), k.get(), d);
    if (status == StatusCode::kOk)
      status = Matrix::MultiplyVector(k.get(), false, theta.get(),
                                      internal.get(), d);
    if (status != StatusCode::kOk) return status;
    if (context->coefficient > 0) {
      std::unique_ptr<Matrix> capacity;
      std::unique_ptr<Vector> old, rate, weights, capacity_flux;
      status = MakeMatrix(n, &capacity, d);
      if (status == StatusCode::kOk) status = MakeVector(n, &old, d);
      if (status == StatusCode::kOk) status = MakeVector(n, &rate, d);
      if (status == StatusCode::kOk) status = MakeVector(n, &weights, d);
      if (status == StatusCode::kOk) status = MakeVector(n, &capacity_flux, d);
      if (status != StatusCode::kOk) return status;
      status = old->GatherFrom(context->state->CommittedTemperature(), globals,
                               n, 0, d);
      if (status == StatusCode::kOk) status = rate->CopyFrom(theta.get(), d);
      if (status == StatusCode::kOk) status = rate->AddScaled(-1, old.get(), d);
      if (status == StatusCode::kOk)
        status = rate->Scale(context->coefficient, d);
      if (status == StatusCode::kOk)
        status = element->ComputeCapacityWeights(weights.get(), d);
      // Task16 G-M2-C: define each diagonal constitutive coefficient at
      // the current nodal/section temperature. Deliberately omit c' * dT.
      for (int i = 0; i < n && status == StatusCode::kOk; ++i) {
        double rho_c = 0;
        status = element->material()->Capacity(theta->Data()[i], &rho_c, d);
        if (status == StatusCode::kOk)
          status = weights->Set(i, weights->Data()[i] * rho_c, d);
      }
      if (status == StatusCode::kOk)
        status = Matrix::Diagonal(weights.get(), capacity.get(), d);
      if (status == StatusCode::kOk)
        status = Matrix::MultiplyVector(capacity.get(), false, rate.get(),
                                        capacity_flux.get(), d);
      if (status == StatusCode::kOk)
        status = internal->AddScaled(1, capacity_flux.get(), d);
      if (status == StatusCode::kOk)
        status = k->AddScaled(context->coefficient, capacity.get(), d);
      if (status != StatusCode::kOk) return status;
    }
    if ((*a->load_starts_)[item] != (*a->load_starts_)[item + 1]) {
      std::unique_ptr<Vector> contribution, face_flux, step_start, prior_flux,
          first_accepted;
      std::unique_ptr<Matrix> film_jacobian;
      if (a->step_start_) {
        status = MakeVector(n, &step_start, d);
        if (status == StatusCode::kOk) status = MakeVector(n, &prior_flux, d);
        if (status == StatusCode::kOk)
          status = MakeMatrix(n, &film_jacobian, d);
        if (status == StatusCode::kOk)
          status = step_start->GatherFrom(a->step_start_, globals, n, 0, d);
        if (status == StatusCode::kOk &&
            context->state->FirstAcceptedTemperature()) {
          status = MakeVector(n, &first_accepted, d);
          if (status == StatusCode::kOk)
            status = first_accepted->GatherFrom(
                context->state->FirstAcceptedTemperature(), globals, n, 0, d);
        }
        if (status != StatusCode::kOk) return status;
      }
      status = MakeVector(n, &contribution, d);
      if (status == StatusCode::kOk) status = MakeVector(n, &face_flux, d);
      if (status != StatusCode::kOk) return status;
      for (int slot = (*a->load_starts_)[item];
           slot < (*a->load_starts_)[item + 1]; ++slot) {
        const int index = (*a->load_indices_)[slot];
        const auto* row = (*a->conditions_->ElementLoads())[index];
        const auto* spec = row->origin->specification();
        if (d) {
          d->source_file = spec->source_file;
          d->source_line = spec->source_line;
          d->step = a->conditions_->StepNumber();
        }
        double magnitude = 0, endpoint = 0;
        const bool exchange = spec->kind == ThermalConditionKind::kFilm ||
                              spec->kind == ThermalConditionKind::kRadiation;
        if (!exchange)
          status = a->conditions_->EvaluateElementFlux(
              row, context->start, context->end, &magnitude, &endpoint, d);
        if (status != StatusCode::kOk) return status;
        bool first = true;
        for (int face = 0; face < spec->expanded_count; ++face) {
          if (spec->expanded_element_ids[face] != element->Id()) continue;
          ThermalLoadContext load;
          load.kind = spec->kind;
          load.extrapolated_free_dofs =
              predicted_free.empty() ? nullptr : predicted_free.data();
          load.face_label =
              NormalizedThermalFaceLabel(spec->expanded_face_codes[face]);
          load.magnitude = magnitude;
          if (exchange) {
            status = prior_flux->Scale(0, d);
            for (int old = 0; old < row->previous_exchange_count &&
                              status == StatusCode::kOk;
                 ++old) {
              const auto* record = a->conditions_->PreviousExchangeAt(
                  row->previous_exchange_begin + old);
              if (!record) return Error("Film prior record range", d);
              ThermalLoadContext previous_load;
              previous_load.kind = record->origin->specification()->kind;
              previous_load.face_label = load.face_label;
              previous_load.coefficient = record->coefficient;
              previous_load.film_property =
                  record->origin->specification()->film_property;
              previous_load.sink_temperature = record->reference_temperature;
              if (previous_load.kind == ThermalConditionKind::kRadiation) {
                const auto* constants =
                    a->conditions_->SourceDomain()->ThermalConstants();
                if (!constants || !constants->has_absolute_zero ||
                    !constants->has_stefan_boltzmann)
                  return Error("Radiation prior physical constants missing", d);
                previous_load.emissivity = record->coefficient;
                previous_load.ambient_temperature =
                    record->reference_temperature;
                previous_load.ambient_baseline = record->reference_temperature;
                previous_load.sigma = constants->stefan_boltzmann;
                previous_load.absolute_zero = constants->absolute_zero;
              }
              status = element->EvaluateLoadContribution(
                  &previous_load, step_start.get(), face_flux.get(), nullptr,
                  d);
              if (status == StatusCode::kOk)
                status = prior_flux->AddScaled(-1, face_flux.get(), d);
            }
            if (status == StatusCode::kOk)
              status = a->conditions_->PopulateExchangeContext(
                  row, context->end, &load, d);
            if (status != StatusCode::kOk) return status;
            load.step_start_temperature = step_start.get();
            load.previous_outward_flux = prior_flux.get();
            load.first_accepted_temperature = first_accepted.get();
          }
          status = element->EvaluateLoadContribution(
              &load, theta.get(), face_flux.get(),
              exchange ? film_jacobian.get() : nullptr, d);
          if (status == StatusCode::kOk && exchange)
            status = k->AddScaled(1, film_jacobian.get(), d);
          if (status == StatusCode::kOk)
            status = first ? contribution->CopyFrom(face_flux.get(), d)
                           : contribution->AddScaled(1, face_flux.get(), d);
          if (status != StatusCode::kOk) return status;
          first = false;
        }
        if (first) return Error("Empty thermal row/element expansion", d);
        status = contribution->Scale(-1, d);
        if (status == StatusCode::kOk)
          status = a->element_flux_->GatherFrom(
              contribution.get(), identity.data(), n,
              (*a->load_flux_offsets_)[index], d);
        if (status != StatusCode::kOk) return status;
      }
    }
    status = a->element_flux_->GatherFrom(internal.get(), identity.data(), n,
                                          offset, d);
    const int coo = (*a->coo_offsets_)[item],
              count = (*a->coo_offsets_)[item + 1] - coo;
    if (status == StatusCode::kOk && count > 0)
      status = k->GatherTo(a->matrix_indices_->data() + coo, count, a->values_,
                           coo, d);
    return status;
  } catch (const std::bad_alloc&) {
    return Error("Thermal element workspace allocation", d,
                 StatusCode::kAllocationFailure);
  }
}
StatusCode ThermalAssembler::Assemble(
    const ThermalModel* model, const ThermalDofMap* map,
    const ThermalConditions* conditions, const ThermalState* state, double time,
    double period, const FirstOrderTimeIntegrator* integrator, double dt,
    ThermalAssembly* out, Diagnostic* d, bool extrapolated_start) {
  if (!model || !map || !conditions || !state || !out || out->model_ != model ||
      out->map_ != map || out->conditions_ != conditions ||
      model->DofMap() != map || model->Conditions() != conditions ||
      !state->HasValidCandidate())
    return Error("Thermal assembly argument/binding mismatch", d);
  const StepDefinition* step = nullptr;
  if (conditions->SourceDomain()->FindStepDefinition(
          conditions->StepNumber(), &step, nullptr) != StatusCode::kOk ||
      !step->IsHeatTransfer())
    return Error("Thermal assembly step missing", d);
  const bool steady = step->HeatTransfer()->steady;
  if (!std::isfinite(time) || time <= 0 || !std::isfinite(period) ||
      period != step->Period() || time > period ||
      (steady && (integrator || dt != 0)) ||
      (!steady && (!integrator || !std::isfinite(dt) || dt <= 0 || dt > time)))
    return Error("Thermal assembly time/integrator mode mismatch", d);
  auto status = state->ValidateMapping(map, d);
  if (status != StatusCode::kOk) return status;
  double coefficient = 0;
  if (!steady) {
    status = integrator->CapacityCoefficient(dt, &coefficient, d);
    if (status != StatusCode::kOk) return status;
    if (!std::isfinite(coefficient) || coefficient <= 0)
      return Error("Invalid thermal capacity coefficient", d,
                   StatusCode::kLibraryFailure);
  }
  out->valid_ = false;
  try {
    const int workers = ConfiguredThreadCount();
    std::unique_ptr<Diagnostic[]> diagnostics(new Diagnostic[workers]);
    Context context;
    context.state = state;
    context.out = out;
    context.coefficient = coefficient;
    context.extrapolated_start = extrapolated_start;
    context.start = steady ? 0 : time - dt;
    context.end = time;
    int observed = 0;
    status = ExecuteElementTasks(static_cast<int>(model->Elements()->size()),
                                 workers, ElementTask, &context,
                                 diagnostics.get(), &observed, d);
    if (status != StatusCode::kOk) return status;
    status = out->residual_->CopyFrom(out->zero_, d);
    if (status == StatusCode::kOk)
      status = out->reaction_->CopyFrom(out->zero_, d);
    if (status == StatusCode::kOk)
      status =
          out->residual_->ScatterAdd(out->element_flux_, out->globals_->data(),
                                     out->element_flux_->Size(), d);
    if (status == StatusCode::kOk) status = out->residual_->Scale(-1, d);
    if (status != StatusCode::kOk) return status;
    const int nx = static_cast<int>(conditions->Fluxes()->size());
    for (int i = 0; i < nx; ++i) {
      double effective = 0, endpoint = 0;
      status = conditions->EvaluateFlux((*conditions->Fluxes())[i],
                                        steady ? 0 : time - dt, time,
                                        &effective, &endpoint, d);
      if (status == StatusCode::kOk)
        status = out->external_flux_->Set(i, effective, d);
      if (status != StatusCode::kOk) return status;
    }
    if (nx > 0)
      status = out->residual_->ScatterAdd(
          out->external_flux_, out->external_globals_->data(), nx, d);
    if (status == StatusCode::kOk && out->free_residual_)
      status = out->free_residual_->GatherFrom(
          out->residual_, out->free_globals_->data(), map->FreeCount(), 0, d);
    const int constrained = static_cast<int>(out->constrained_globals_->size());
    if (status == StatusCode::kOk && constrained > 0) {
      status = out->constrained_scratch_->GatherFrom(
          out->residual_, out->constrained_globals_->data(), constrained, 0, d);
      if (status == StatusCode::kOk)
        status = out->constrained_scratch_->Scale(-1, d);
      if (status == StatusCode::kOk)
        status = out->reaction_->ScatterAdd(out->constrained_scratch_,
                                            out->constrained_globals_->data(),
                                            constrained, d);
    }
    if (status == StatusCode::kOk)
      status =
          out->global_increment_->CopyFrom(state->CandidateTemperature(), d);
    if (status == StatusCode::kOk)
      status = out->global_increment_->AddScaled(
          -1, state->CommittedTemperature(), d);
    if (status == StatusCode::kOk)
      status = out->increment_->GatherFrom(out->global_increment_,
                                           out->participating_globals_->data(),
                                           out->increment_->Size(), 0, d);
    if (status != StatusCode::kOk) return status;
    out->valid_ = true;
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return Error("Thermal assembly dispatch allocation", d,
                 StatusCode::kAllocationFailure);
  }
}
}  // namespace fesa
