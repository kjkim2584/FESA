#include "fesa/analysis/effective_conditions.h"

#include <map>
#include <memory>
#include <new>
#include <utility>

#include "model/model_storage.h"

namespace fesa {
namespace {

using BoundaryMap =
    std::map<std::pair<int, int>, const BoundaryDefinition*>;

bool SameName(const char* a, const char* b) {
  for (;;) {
    while (*a == ' ' || *a == '\t') ++a;
    while (*b == ' ' || *b == '\t') ++b;
    const char ac = *a >= 'A' && *a <= 'Z' ? *a - 'A' + 'a' : *a;
    const char bc = *b >= 'A' && *b <= 'Z' ? *b - 'A' + 'a' : *b;
    if (ac != bc) return false;
    if (ac == '\0') return true;
    ++a;
    ++b;
  }
}

template <typename T>
bool SameTarget(const T* a, const T* b) {
  return a->Target() == b->Target() && a->Dof() == b->Dof() &&
         SameName(a->Scope(), b->Scope()) &&
         (a->Target() == ConditionTarget::kNode
              ? a->NodeId() == b->NodeId()
              : SameName(a->Name(), b->Name()));
}

template <typename T>
bool HasTarget(const std::vector<const T*>* definitions, const T* target) {
  for (const auto* definition : *definitions)
    if (SameTarget(definition, target)) return true;
  return false;
}

template <typename T>
void SourceContext(const T* definition, int node_id, int step, Diagnostic* d) {
  if (d == nullptr) return;
  d->source_file = definition->SourceFile();
  d->duplicate_boundary = false;
  d->path_kind = DiagnosticPathKind::kLiteralFile;
  d->source_line = definition->SourceLine();
  d->first_source_line = 0;
  d->step = step;
  d->node_id = node_id;
  d->element_id = 0;
  d->condition_dof = definition->Dof();
}

template <typename T>
StatusCode Targets(const Domain* domain, const T* definition, int step,
                   std::vector<int>* nodes, Diagnostic* d) {
  nodes->clear();
  if (definition->Target() == ConditionTarget::kNode) {
    nodes->push_back(definition->NodeId());
  } else {
    const IndexSet* set = nullptr;
    StatusCode status = domain->FindIndexSetByName(
        definition->Scope(), definition->Name(), false, &set, d);
    if (status != StatusCode::kOk) {
      SourceContext(definition, 0, step, d);
      if (status == StatusCode::kAllocationFailure) return status;
      return ReportFailure(StatusCode::kInputError,
                           DiagnosticStage::kValidation,
                           "BuildEffectiveConditions missing set", d);
    }
    for (int i = 0; i < set->Size(); ++i) nodes->push_back(set->Ids()[i]);
  }
  for (int id : *nodes) {
    const Node* node = nullptr;
    if (domain->FindNode(id, &node, nullptr) != StatusCode::kOk) {
      SourceContext(definition, id, step, d);
      return ReportFailure(StatusCode::kInputError,
                           DiagnosticStage::kValidation,
                           "BuildEffectiveConditions missing node", d);
    }
  }
  return StatusCode::kOk;
}

// Reuse the same two common vectors for differences and assembled totals.
StatusCode Add(double a, double scale, double b, Vector* addend, Vector* sum,
               double* out, Diagnostic* d) {
  StatusCode status = addend->Set(0, b, d);
  if (status == StatusCode::kOk) status = sum->Set(0, a, d);
  if (status == StatusCode::kOk) status = sum->AddScaled(scale, addend, d);
  if (status == StatusCode::kOk) *out = sum->Data()[0];
  return status;
}

void AppendBoundary(const BoundaryDefinition* definition,
                    const BoundaryDefinition* origin,
                    const BoundaryDefinition* previous, int node,
                    std::vector<const EffectiveBoundary*>* rows,
                    const BoundaryDefinition* physical_previous = nullptr) {
  std::unique_ptr<EffectiveBoundary> entry(new EffectiveBoundary);
  entry->definition = definition;
  entry->origin_definition = origin;
  entry->previous_definition = previous;
  entry->physical_previous_definition = physical_previous;
  entry->node_id = node;
  rows->push_back(entry.get());
  entry.release();
}

StatusCode DuplicateBoundary(const BoundaryDefinition* definition,
                             const BoundaryDefinition* first, int node,
                             int step, Diagnostic* d) {
  SourceContext(definition, node, step, d);
  if (d) {
    d->first_source_line = first->SourceLine();
    d->duplicate_boundary = true;
    d->boundary_dof = definition->Dof();
    d->first_boundary_value = first->Value();
    d->boundary_value = definition->Value();
  }
  return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                       "BuildEffectiveConditions duplicate boundary", d);
}

// Physical history borrows Domain definitions, never previous snapshot rows.
StatusCode BuildPhysicalBoundaries(
    const Domain* domain, int step,
    const std::vector<const EffectiveBoundary*>* previous,
    std::vector<int>* targets, BoundaryMap* physical, Diagnostic* d) {
  if (previous) {
    for (const auto* boundary : *previous)
      physical->emplace(
          std::make_pair(boundary->node_id, boundary->definition->Dof()),
          boundary->definition);
    return StatusCode::kOk;
  }
  // Initial constraints are prior physical conditions even for OP=NEW.
  for (const auto& item : *domain->Boundaries()) {
    const auto* definition = item.second;
    if (definition->Step() != 0) continue;
    StatusCode status = Targets(domain, definition, step, targets, d);
    if (status != StatusCode::kOk) return status;
    for (int id : *targets) {
      const auto inserted =
          physical->emplace(std::make_pair(id, definition->Dof()), definition);
      if (!inserted.second)
        return DuplicateBoundary(definition, inserted.first->second, id, step,
                                 d);
    }
  }
  return StatusCode::kOk;
}

void ResolveBoundaries(const Domain* domain, const StepDefinition* step,
                       const std::vector<const EffectiveBoundary*>* previous,
                       std::vector<const EffectiveBoundary*>* rows) {
  if (step->BoundaryOperation() != ConditionOperation::kNew) {
    if (previous) {
      for (const auto* row : *previous)
        if (!HasTarget(step->Boundaries(), row->definition))
          AppendBoundary(row->definition, row->origin_definition,
                         row->definition, 0, rows);
    } else {
      for (const auto& item : *domain->Boundaries()) {
        const auto* initial = item.second;
        if (initial->Step() == 0 && !HasTarget(step->Boundaries(), initial))
          AppendBoundary(initial, initial, nullptr, 0, rows);
      }
    }
  }
  for (const auto* definition : *step->Boundaries()) {
    const BoundaryDefinition* origin = definition;
    const BoundaryDefinition* prior = nullptr;
    if (step->BoundaryOperation() != ConditionOperation::kNew) {
      if (previous) {
        for (const auto* row : *previous)
          if (SameTarget(row->definition, definition)) {
            origin = row->origin_definition;
            prior = row->definition;
            break;
          }
      } else {
        for (const auto& item : *domain->Boundaries())
          if (item.second->Step() == 0 && SameTarget(item.second, definition)) {
            origin = prior = item.second;
            break;
          }
      }
    }
    AppendBoundary(definition, origin, prior, 0, rows);
  }
}

StatusCode ExpandBoundaries(const Domain* domain, int step,
                            const std::vector<const EffectiveBoundary*>* rows,
                            BoundaryMap* previous_physical,
                            std::vector<const EffectiveBoundary*>* boundaries,
                            std::vector<const EffectiveBoundary*>* released,
                            std::vector<int>* targets, Diagnostic* d) {
  BoundaryMap seen;
  for (const auto* row : *rows) {
    const auto* definition = row->definition;
    StatusCode status = Targets(domain, definition, step, targets, d);
    if (status != StatusCode::kOk) return status;
    for (int id : *targets) {
      const auto key = std::make_pair(id, definition->Dof());
      const auto found = seen.find(key);
      if (found != seen.end())
        return DuplicateBoundary(definition, found->second, id, step, d);
      const BoundaryDefinition* physical_previous = nullptr;
      if (previous_physical) {
        const auto prior = previous_physical->find(key);
        if (prior != previous_physical->end()) {
          physical_previous = prior->second;
          previous_physical->erase(prior);
        }
      }
      AppendBoundary(definition, row->origin_definition,
                     row->previous_definition, id, boundaries,
                     physical_previous);
      seen.emplace(key, definition);
    }
  }
  if (previous_physical)
    for (const auto& item : *previous_physical)
      AppendBoundary(item.second, item.second, nullptr, item.first.first,
                     released);
  return StatusCode::kOk;
}

StatusCode AppendLoad(const LoadDefinition* definition,
                      const LoadDefinition* origin,
                      const LoadDefinition* previous, bool active, double value,
                      double prior_value, int step, Vector* addend, Vector* sum,
                      std::vector<const EffectiveLoad*>* rows, Diagnostic* d) {
  std::unique_ptr<EffectiveLoad> row(new EffectiveLoad);
  row->definition = definition;
  row->origin_definition = origin;
  row->previous_definition = previous;
  row->active = active;
  row->value = value;
  row->contributes_flux = value != 0;
  StatusCode status = Add(value, -1, prior_value, addend, sum, &row->change, d);
  if (status != StatusCode::kOk) {
    SourceContext(definition, definition->NodeId(), step, d);
    return status;
  }
  rows->push_back(row.get());
  row.release();
  return StatusCode::kOk;
}

StatusCode ResolveLoads(const StepDefinition* step,
                        const std::vector<const EffectiveLoad*>* previous,
                        Vector* addend, Vector* sum,
                        std::vector<const EffectiveLoad*>* rows,
                        Diagnostic* d) {
  const bool reset = step->LoadOperation() == ConditionOperation::kNew;
  if (previous) {
    for (const auto* row : *previous) {
      if (!row->active) continue;
      if (!reset && HasTarget(step->Loads(), row->definition)) continue;
      const StatusCode status = AppendLoad(
          row->definition, row->origin_definition, row->definition, !reset,
          reset ? 0 : row->value, row->value, step->Id(), addend, sum, rows, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  for (auto current = step->Loads()->begin(); current != step->Loads()->end();
       ++current) {
    const auto* definition = *current;
    const EffectiveLoad* prior = nullptr;
    if (!reset && previous) {
      for (const auto* row : *previous) {
        if (!row->active || !SameTarget(row->definition, definition)) continue;
        if (prior) {
          SourceContext(definition, definition->NodeId(), step->Id(), d);
          if (d) d->first_source_line = prior->definition->SourceLine();
          return ReportFailure(
              StatusCode::kInputError, DiagnosticStage::kValidation,
              "Ambiguous CLOAD original target/DOF; use NEW", d);
        }
        prior = row;
      }
    }
    // Only the first current row replaces the prior singleton. Other rows
    // create separate origins, so their deltas do not subtract it again.
    for (auto earlier = step->Loads()->begin(); earlier != current; ++earlier)
      if (SameTarget(*earlier, definition)) {
        prior = nullptr;
        break;
      }
    const StatusCode status =
        AppendLoad(definition, prior ? prior->origin_definition : definition,
                   prior ? prior->definition : nullptr,
                   prior || definition->Value() != 0, definition->Value(),
                   prior ? prior->value : 0, step->Id(), addend, sum, rows, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

StatusCode ExpandLoads(const Domain* domain, int step,
                       const std::vector<const EffectiveLoad*>* rows,
                       Vector* addend, Vector* sum,
                       std::vector<const EffectiveLoad*>* loads,
                       std::vector<const EffectiveLoadTotal*>* output_totals,
                       std::vector<int>* targets, Diagnostic* d) {
  std::map<std::pair<int, int>, EffectiveLoadTotal*> totals;
  for (const auto* row : *rows) {
    const auto* definition = row->definition;
    StatusCode status = Targets(domain, definition, step, targets, d);
    if (status != StatusCode::kOk) return status;
    for (int id : *targets) {
      std::unique_ptr<EffectiveLoad> entry(new EffectiveLoad);
      entry->node_id = id;
      entry->definition = definition;
      entry->origin_definition = row->origin_definition;
      entry->previous_definition = row->previous_definition;
      entry->active = row->active;
      entry->contributes_flux = row->contributes_flux;
      entry->value = row->value;
      entry->change = row->change;
      loads->push_back(entry.get());
      entry.release();
      const auto key = std::make_pair(id, definition->Dof());
      auto found = totals.find(key);
      if (found == totals.end()) {
        std::unique_ptr<EffectiveLoadTotal> total(new EffectiveLoadTotal);
        total->node_id = id;
        total->dof = definition->Dof();
        output_totals->push_back(total.get());
        auto* stored = total.release();
        found = totals.emplace(key, stored).first;
      }
      status = Add(found->second->value, 1, row->value, addend, sum,
                   &found->second->value, d);
      if (status != StatusCode::kOk) {
        SourceContext(definition, id, step, d);
        return status;
      }
    }
  }
  return StatusCode::kOk;
}

}  // namespace

EffectiveConditions::~EffectiveConditions() {
  delete line_load_state_;
  if (loads_)
    for (const auto* entry : *loads_) delete entry;
  if (boundaries_)
    for (const auto* entry : *boundaries_) delete entry;
  if (released_boundaries_)
    for (const auto* entry : *released_boundaries_) delete entry;
  if (totals_)
    for (const auto* entry : *totals_) delete entry;
  if (load_rows_)
    for (const auto* entry : *load_rows_) delete entry;
  if (boundary_rows_)
    for (const auto* entry : *boundary_rows_) delete entry;
  delete loads_;
  delete boundaries_;
  delete released_boundaries_;
  delete totals_;
  delete load_rows_;
  delete boundary_rows_;
}

StatusCode BuildEffectiveConditions(const Domain* domain, int step_number,
                                    const EffectiveConditions* previous,
                                    EffectiveConditions** out, Diagnostic* d) {
  if (out == nullptr)
    return model_storage::Invalid("BuildEffectiveConditions", d);
  *out = nullptr;
  if (domain == nullptr || domain->Nodes()->empty() || step_number <= 0 ||
      (step_number == 1 && previous != nullptr) ||
      (step_number > 1 &&
       (previous == nullptr || previous->SourceDomain() != domain ||
        previous->StepNumber() != step_number - 1)))
    return model_storage::Invalid("BuildEffectiveConditions history", d);
  const StepDefinition* step = nullptr;
  StatusCode status = domain->FindStepDefinition(step_number, &step, d);
  if (status != StatusCode::kOk) return status;
  if (step->IsHeatTransfer()) {
    if (d) d->step = step_number;
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Heat transfer step requires thermal conditions", d);
  }
  try {
    std::unique_ptr<EffectiveConditions> result(new EffectiveConditions);
    result->domain_ = domain;
    result->step_number_ = step_number;
    result->loads_ = new std::vector<const EffectiveLoad*>(0);
    result->boundaries_ = new std::vector<const EffectiveBoundary*>(0);
    result->totals_ = new std::vector<const EffectiveLoadTotal*>(0);
    result->load_rows_ = new std::vector<const EffectiveLoad*>(0);
    result->boundary_rows_ = new std::vector<const EffectiveBoundary*>(0);
    std::vector<int> targets(0);
    std::unique_ptr<BoundaryMap> previous_physical;
    if (step->HasFiniteGeometry()) {
      result->released_boundaries_ =
          new std::vector<const EffectiveBoundary*>(0);
      previous_physical.reset(new BoundaryMap);
      status = BuildPhysicalBoundaries(
          domain, step_number, previous ? previous->boundaries_ : nullptr,
          &targets, previous_physical.get(), d);
      if (status != StatusCode::kOk) return status;
    }
    ResolveBoundaries(domain, step,
                      previous ? previous->boundary_rows_ : nullptr,
                      result->boundary_rows_);
    status = ExpandBoundaries(domain, step_number, result->boundary_rows_,
                              previous_physical.get(), result->boundaries_,
                              result->released_boundaries_, &targets, d);
    if (status != StatusCode::kOk) return status;
    std::unique_ptr<Vector> addend, sum;
    const auto* prior_rows = previous ? previous->load_rows_ : nullptr;
    if (!step->Loads()->empty() || (prior_rows && !prior_rows->empty())) {
      Vector* raw = nullptr;
      status = Vector::Create(1, &raw, d);
      if (status != StatusCode::kOk) return status;
      addend.reset(raw);
      status = Vector::Create(1, &raw, d);
      if (status != StatusCode::kOk) return status;
      sum.reset(raw);
    }
    status = ResolveLoads(step, prior_rows, addend.get(), sum.get(),
                          result->load_rows_, d);
    if (status == StatusCode::kOk)
      status =
          ExpandLoads(domain, step_number, result->load_rows_, addend.get(),
                      sum.get(), result->loads_, result->totals_, &targets, d);
    if (status != StatusCode::kOk) return status;
    if (domain->LineLoadProgram()) {
      status = domain->LineLoadProgram()->ResolveStep(
          step_number, previous ? previous->line_load_state_ : nullptr,
          &result->line_load_state_, d);
      if (status != StatusCode::kOk) return status;
    }
    *out = result.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("BuildEffectiveConditions", d);
  }
}

const std::vector<const EffectiveLoad*>* EffectiveConditions::Loads() const {
  return loads_;
}
const BeamLineLoadState* EffectiveConditions::LineLoadState() const {
  return line_load_state_;
}
const std::vector<const EffectiveBoundary*>* EffectiveConditions::Boundaries()
    const {
  return boundaries_;
}
const std::vector<const EffectiveBoundary*>*
EffectiveConditions::ReleasedBoundaries() const {
  return released_boundaries_;
}
const std::vector<const EffectiveLoadTotal*>* EffectiveConditions::LoadTotals()
    const {
  return totals_;
}
const Domain* EffectiveConditions::SourceDomain() const { return domain_; }
int EffectiveConditions::StepNumber() const { return step_number_; }

}  // namespace fesa
