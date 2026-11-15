#include "fesa/analysis/thermal_conditions.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <new>
#include <set>
#include <utility>

namespace fesa {
namespace {
using Definition = ThermalConditionDefinition;
using Flux = EffectiveConcentratedFlux;
using Temperature = EffectiveTemperature;
using ElementLoad = EffectiveElementLoad;
StatusCode Failure(const ThermalConditionSpec* s, int step, const char* message,
                   Diagnostic* d, StatusCode code = StatusCode::kInputError) {
  if (d && s) {
    d->source_file = s->source_file;
    d->source_line = s->source_line;
    d->step = step;
    d->node_id = s->target_id;
    d->condition_dof = s->first;
    d->path_kind = DiagnosticPathKind::kLiteralFile;
  }
  return ReportFailure(code, DiagnosticStage::kValidation, message, d);
}
int CompareName(const char* a, const char* b) {
  for (;;) {
    const char ac = *a >= 'a' && *a <= 'z' ? *a - 'a' + 'A' : *a;
    const char bc = *b >= 'a' && *b <= 'z' ? *b - 'a' + 'A' : *b;
    if (ac != bc) return ac < bc ? -1 : 1;
    if (!ac) return 0;
    ++a;
    ++b;
  }
}
// Pointer keys borrow immutable Domain rows. Comparing their original target
// avoids allocating/copying a string key and retains equal-key input order.
struct KeyLess {
  bool operator()(const Definition* a, const Definition* b) const {
    const auto* x = a->specification();
    const auto* y = b->specification();
    if (x->kind != y->kind) return x->kind < y->kind;
    if (x->family != y->family) return x->family < y->family;
    const int face = CompareName(x->face_label, y->face_label);
    if (face) return face < 0;
    if (x->target != y->target) return x->target < y->target;
    const int scope = CompareName(x->scope, y->scope);
    if (scope) return scope < 0;
    if (x->target_id != y->target_id) return x->target_id < y->target_id;
    const int name = CompareName(x->name, y->name);
    if (name) return name < 0;
    return x->first < y->first;
  }
};
void CopyFlux(const Flux* from, Flux* to) {
  to->node_id = from->node_id;
  to->label = from->label;
  to->end_value = from->end_value;
  to->previous_end_value = from->previous_end_value;
  to->amplitude = from->amplitude;
  to->ramp = from->ramp;
  to->definition_step = from->definition_step;
  to->definition_period = from->definition_period;
  to->evaluation = from->evaluation;
  to->history = from->history;
  to->origin = from->origin;
}
void CopyElementLoad(const ElementLoad* from, ElementLoad* to) {
  to->element_id = from->element_id;
  to->previous_exchange_begin = from->previous_exchange_begin;
  to->previous_exchange_count = from->previous_exchange_count;
  to->end_value = from->end_value;
  to->previous_end_value = from->previous_end_value;
  to->amplitude = from->amplitude;
  to->ramp = from->ramp;
  to->definition_step = from->definition_step;
  to->definition_period = from->definition_period;
  to->history = from->history;
  to->origin = from->origin;
}
void CopyTemperature(const Temperature* from, Temperature* to) {
  to->node_id = from->node_id;
  to->label = from->label;
  to->end_value = from->end_value;
  to->amplitude = from->amplitude;
  to->ramp = from->ramp;
  to->start_from_committed_temperature = from->start_from_committed_temperature;
  to->definition_step = from->definition_step;
  to->definition_period = from->definition_period;
  to->origin = from->origin;
}
StatusCode Targets(const Domain* domain, const Definition* definition,
                   const std::vector<int>* participating, int step,
                   std::vector<int>* ids, Diagnostic* d) {
  ids->clear();
  const auto* s = definition->specification();
  if (s->target == ThermalConditionTarget::kNode) {
    if (!std::binary_search(participating->begin(), participating->end(),
                            s->target_id))
      return Failure(s, step, "Thermal target node is not participating", d);
    ids->push_back(s->target_id);
  } else if (s->target == ThermalConditionTarget::kNodeSet) {
    const IndexSet* set = nullptr;
    const auto lookup =
        domain->FindIndexSetByName(s->scope, s->name, false, &set, nullptr);
    if (lookup == StatusCode::kAllocationFailure)
      return Failure(s, step, "Thermal target set lookup allocation", d,
                     lookup);
    if (lookup != StatusCode::kOk)
      return Failure(s, step, "Undefined thermal target set", d);
    for (int i = 0; i < set->Size(); ++i)
      if (std::binary_search(participating->begin(), participating->end(),
                             set->Ids()[i]))
        ids->push_back(set->Ids()[i]);
    if (ids->empty())
      return Failure(s, step, "Thermal target set has no participating nodes",
                     d);
  } else {
    return Failure(s, step, "Thermal target kind unsupported", d,
                   StatusCode::kUnsupported);
  }
  for (int id : *ids) {
    const Node* node = nullptr;
    if (domain->FindNode(id, &node, nullptr) != StatusCode::kOk ||
        s->first < 11 || s->last >= 11 + node->temperature_dof_count())
      return Failure(s, step, "Thermal target label outside node components",
                     d);
  }
  return StatusCode::kOk;
}
StatusCode Coefficient(const Amplitude* amplitude, int definition_step,
                       double definition_period, int current_step,
                       double total_start, double time, double* out,
                       Diagnostic* d) {
  AmplitudeEvaluationContext context;
  context.definition_step = definition_step;
  context.definition_period = definition_period;
  context.current_step = current_step;
  context.step_time = time;
  context.total_time = total_start + time;
  return amplitude->EvaluateInStep(&context, out, d);
}
template <typename T>
StatusCode RowValue(const T* row, int step, double period, double total_start,
                    double time, double* out, Diagnostic* d) {
  double value = 0;
  if (row->amplitude) {
    double coefficient = 0;
    const auto status = Coefficient(row->amplitude, row->definition_step,
                                    row->definition_period, step, total_start,
                                    time, &coefficient, d);
    if (status != StatusCode::kOk) return status;
    value = coefficient * row->origin->specification()->magnitude;
  } else if (row->history == ThermalFluxHistory::kCarried) {
    value = row->previous_end_value;
  } else {
    const double fraction = row->ramp ? time / period : 1;
    const double end = row->history == ThermalFluxHistory::kRemoved
                           ? 0
                           : row->origin->specification()->magnitude;
    value = (1 - fraction) * row->previous_end_value + fraction * end;
  }
  if (!std::isfinite(value))
    return Failure(row->origin->specification(), step,
                   "Nonfinite thermal power", d, StatusCode::kLibraryFailure);
  *out = value;
  return StatusCode::kOk;
}
template <typename T>
void DeleteEntries(std::vector<const T*>* entries) {
  if (!entries) return;
  for (const auto* entry : *entries) delete entry;
  delete entries;
}
}  // namespace
ThermalConditions::~ThermalConditions() {
  DeleteEntries(temperatures_);
  DeleteEntries(fluxes_);
  DeleteEntries(rows_);
  DeleteEntries(element_rows_);
  DeleteEntries(element_loads_);
  DeleteEntries(previous_exchanges_);
}
const std::vector<const EffectiveTemperature*>*
ThermalConditions::Temperatures() const {
  return temperatures_;
}
const std::vector<const EffectiveConcentratedFlux*>* ThermalConditions::Fluxes()
    const {
  return fluxes_;
}
const std::vector<const EffectiveElementLoad*>*
ThermalConditions::ElementLoads() const {
  return element_loads_;
}
const Domain* ThermalConditions::SourceDomain() const { return domain_; }
int ThermalConditions::StepNumber() const { return step_number_; }
StatusCode ThermalConditions::EvaluateTemperature(const Temperature* entry,
                                                  double time,
                                                  double committed_start,
                                                  double* out,
                                                  Diagnostic* d) const {
  if (!entry || !out || !std::isfinite(time) || time < 0 || time > period_ ||
      !std::isfinite(committed_start))
    return Failure(nullptr, step_number_,
                   "Invalid thermal temperature evaluation", d,
                   StatusCode::kInvalidArgument);
  double value = entry->end_value;
  if (entry->amplitude) {
    double coefficient = 0;
    const auto status = Coefficient(entry->amplitude, entry->definition_step,
                                    entry->definition_period, step_number_,
                                    total_start_, time, &coefficient, d);
    if (status != StatusCode::kOk) return status;
    value = entry->origin->specification()->magnitude * coefficient;
  } else if (entry->start_from_committed_temperature) {
    const double fraction = time / period_;
    value = (1 - fraction) * committed_start + fraction * entry->end_value;
  }
  if (!std::isfinite(value))
    return Failure(entry->origin->specification(), step_number_,
                   "Nonfinite prescribed temperature", d,
                   StatusCode::kLibraryFailure);
  *out = value;
  return StatusCode::kOk;
}
StatusCode ThermalConditions::EvaluateFlux(const Flux* entry, double start,
                                           double end, double* effective,
                                           double* endpoint,
                                           Diagnostic* d) const {
  if (!entry || !effective || !endpoint || effective == endpoint ||
      !std::isfinite(start) || !std::isfinite(end) || start < 0 ||
      end <= start || end > period_)
    return Failure(nullptr, step_number_, "Invalid thermal power evaluation", d,
                   StatusCode::kInvalidArgument);
  double last = 0;
  auto status =
      RowValue(entry, step_number_, period_, total_start_, end, &last, d);
  if (status != StatusCode::kOk) return status;
  double value = last;
  if (!steady_ &&
      entry->evaluation == ThermalFluxEvaluation::kEndpointAverage) {
    double first = 0;
    status =
        RowValue(entry, step_number_, period_, total_start_, start, &first, d);
    if (status != StatusCode::kOk) return status;
    value = .5 * first + .5 * last;
  }
  *effective = value;
  *endpoint = last;
  return StatusCode::kOk;
}
StatusCode ThermalConditions::EvaluateElementFlux(const ElementLoad* entry,
                                                  double start, double end,
                                                  double* effective,
                                                  double* endpoint,
                                                  Diagnostic* d) const {
  if (!entry || !effective || !endpoint || effective == endpoint ||
      !std::isfinite(start) || !std::isfinite(end) || start < 0 ||
      end <= start || end > period_)
    return Failure(nullptr, step_number_, "Invalid thermal power evaluation", d,
                   StatusCode::kInvalidArgument);
  double last = 0;
  auto status =
      RowValue(entry, step_number_, period_, total_start_, end, &last, d);
  if (status != StatusCode::kOk) return status;
  double value = last;
  if (!steady_) {
    double first = 0;
    status =
        RowValue(entry, step_number_, period_, total_start_, start, &first, d);
    if (status != StatusCode::kOk) return status;
    value = .5 * first + .5 * last;
  }
  *effective = value;
  *endpoint = last;
  return StatusCode::kOk;
}
const PreviousThermalExchange* ThermalConditions::PreviousExchangeAt(
    int index) const {
  return previous_exchanges_ && index >= 0 &&
                 static_cast<std::size_t>(index) < previous_exchanges_->size()
             ? (*previous_exchanges_)[index]
             : nullptr;
}
StatusCode ThermalConditions::PopulateExchangeContext(const ElementLoad* entry,
                                                      double time,
                                                      ThermalLoadContext* out,
                                                      Diagnostic* d) const {
  if (entry && entry->origin && out && std::isfinite(time) && time >= 0 &&
      time <= period_ &&
      entry->origin->specification()->kind ==
          ThermalConditionKind::kRadiation) {
    const auto* spec = entry->origin->specification();
    const auto* constants = domain_->ThermalConstants();
    if (!constants || !constants->has_absolute_zero ||
        !constants->has_stefan_boltzmann)
      return Failure(spec, step_number_,
                     "Radiation requires physical constants", d);
    if (spec->amplitude) {
      double start = total_start_, end = total_start_ + period_;
      if (spec->amplitude->TimeBasis() == AmplitudeTimeBasis::kStepTime) {
        start = entry->definition_step < step_number_ ? entry->definition_period
                                                      : 0;
        end = entry->definition_step < step_number_ ? entry->definition_period
                                                    : period_;
      }
      double minimum = 0, maximum = 0;
      const auto range_status =
          spec->amplitude->EvaluateRange(start, end, &minimum, &maximum, d);
      if (range_status != StatusCode::kOk) return range_status;
      const double low = minimum * spec->ambient_temperature;
      const double high = maximum * spec->ambient_temperature;
      if (!std::isfinite(low) || !std::isfinite(high) ||
          std::min(low, high) < constants->absolute_zero)
        return Failure(spec, step_number_,
                       "Radiation ambient interval below absolute zero", d);
    }
    double ambient = 1, baseline = 1;
    auto status = StatusCode::kOk;
    if (spec->amplitude)
      status = Coefficient(spec->amplitude, entry->definition_step,
                           entry->definition_period, step_number_, total_start_,
                           time, &ambient, d);
    if (status == StatusCode::kOk && spec->amplitude)
      status = Coefficient(spec->amplitude, entry->definition_step,
                           entry->definition_period, step_number_, total_start_,
                           0, &baseline, d);
    if (status != StatusCode::kOk) return status;
    ambient *= spec->ambient_temperature;
    baseline *= spec->ambient_temperature;
    if (!std::isfinite(ambient) || !std::isfinite(baseline) ||
        ambient < constants->absolute_zero ||
        baseline < constants->absolute_zero)
      return Failure(spec, step_number_,
                     "Radiation ambient below absolute zero", d);
    out->kind = ThermalConditionKind::kRadiation;
    out->ambient_temperature = ambient;
    out->ambient_baseline = baseline;
    out->emissivity = spec->emissivity;
    out->sigma = constants->stefan_boltzmann;
    out->absolute_zero = constants->absolute_zero;
    out->step_scale = entry->history == ThermalFluxHistory::kCarried
                          ? 1
                          : (entry->ramp ? time / period_ : 1);
    return StatusCode::kOk;
  }
  if (!entry || !entry->origin || !out || !std::isfinite(time) || time < 0 ||
      time > period_ ||
      entry->origin->specification()->kind != ThermalConditionKind::kFilm)
    return Failure(nullptr, step_number_, "Invalid film time evaluation", d,
                   StatusCode::kInvalidArgument);
  const auto* spec = entry->origin->specification();
  double multiplier = 1, sink = 1, baseline = 1;
  auto status = StatusCode::kOk;
  if (spec->film_amplitude && !spec->film_property)
    status = Coefficient(spec->film_amplitude, entry->definition_step,
                         entry->definition_period, step_number_, total_start_,
                         time, &multiplier, d);
  if (status == StatusCode::kOk && spec->amplitude)
    status = Coefficient(spec->amplitude, entry->definition_step,
                         entry->definition_period, step_number_, total_start_,
                         time, &sink, d);
  if (status == StatusCode::kOk && spec->amplitude)
    status = Coefficient(spec->amplitude, entry->definition_step,
                         entry->definition_period, step_number_, total_start_,
                         0, &baseline, d);
  if (status != StatusCode::kOk) return status;
  const double h = spec->coefficient * multiplier;
  sink *= spec->sink_temperature;
  baseline *= spec->sink_temperature;
  if (!std::isfinite(h) || h < 0 || !std::isfinite(sink) ||
      !std::isfinite(baseline))
    return Failure(spec, step_number_, "Film time coefficient range", d,
                   StatusCode::kLibraryFailure);
  const bool explicit_h = spec->film_amplitude && !spec->film_property;
  out->kind = ThermalConditionKind::kFilm;
  out->coefficient = h;
  out->film_property = spec->film_property;
  out->sink_temperature = sink;
  out->sink_baseline = baseline;
  out->step_scale = entry->history == ThermalFluxHistory::kCarried || explicit_h
                        ? 1
                        : (entry->ramp ? time / period_ : 1);
  return StatusCode::kOk;
}
StatusCode BuildThermalConditions(const Domain* domain, int step,
                                  const ThermalConditions* previous,
                                  ThermalConditions** out, Diagnostic* d) {
  if (!out)
    return Failure(nullptr, step, "Null thermal snapshot output", d,
                   StatusCode::kInvalidArgument);
  *out = nullptr;
  if (!domain || step < 1 || (step == 1 && previous) ||
      (step > 1 && (!previous || previous->domain_ != domain ||
                    previous->step_number_ != step - 1)))
    return Failure(nullptr, step, "Invalid thermal snapshot history", d,
                   StatusCode::kInvalidArgument);
  const StepDefinition* definition = nullptr;
  if (domain->FindStepDefinition(step, &definition, nullptr) != StatusCode::kOk)
    return Failure(nullptr, step, "Missing thermal step", d);
  if (!definition->IsHeatTransfer())
    return Failure(nullptr, step, "Mechanical step in thermal resolver", d,
                   StatusCode::kUnsupported);
  try {
    std::unique_ptr<ThermalConditions> result(new ThermalConditions);
    result->domain_ = domain;
    result->step_number_ = step;
    result->period_ = definition->Period();
    result->steady_ = definition->HeatTransfer()->steady;
    const bool ramp = definition->HeatTransfer()->DefaultAmplitude() ==
                      ThermalStepAmplitude::kRamp;
    for (int i = 1; i < step; ++i) {
      const StepDefinition* old = nullptr;
      if (domain->FindStepDefinition(i, &old, nullptr) != StatusCode::kOk ||
          !old->IsHeatTransfer())
        return Failure(nullptr, step, "Missing preceding thermal step", d);
      result->total_start_ += old->Period();
    }
    if (!std::isfinite(result->total_start_ + result->period_))
      return Failure(nullptr, step, "Nonfinite thermal total time", d,
                     StatusCode::kLibraryFailure);
    result->temperatures_ = new std::vector<const Temperature*>(0);
    result->fluxes_ = new std::vector<const Flux*>(0);
    result->rows_ = new std::vector<const Flux*>(0);
    result->element_rows_ = new std::vector<const ElementLoad*>(0);
    result->element_loads_ = new std::vector<const ElementLoad*>(0);
    std::vector<int> participating(0), ids(0);
    auto status = CollectThermalNodeIds(domain, &participating, d);
    if (status != StatusCode::kOk) return status;
    using PhysicalKey = std::pair<int, int>;
    std::map<PhysicalKey, std::unique_ptr<Temperature>> temperatures;
    const bool new_temperature =
        definition->HasBoundaryBlock() &&
        definition->BoundaryOperation() == ConditionOperation::kNew;
    if (previous && !new_temperature) {
      for (const auto* old : *previous->temperatures_) {
        std::unique_ptr<Temperature> copy(new Temperature);
        CopyTemperature(old, copy.get());
        copy->ramp = false;
        copy->start_from_committed_temperature = false;
        if (copy->amplitude) {
          double coefficient = 0;
          status = Coefficient(
              copy->amplitude, copy->definition_step, copy->definition_period,
              step, result->total_start_, result->period_, &coefficient, d);
          if (status != StatusCode::kOk) return status;
          copy->end_value =
              coefficient * copy->origin->specification()->magnitude;
          if (!std::isfinite(copy->end_value))
            return Failure(copy->origin->specification(), step,
                           "Nonfinite carried temperature", d,
                           StatusCode::kLibraryFailure);
        }
        temperatures.emplace(PhysicalKey(copy->node_id, copy->label),
                             std::move(copy));
      }
    }
    std::set<PhysicalKey> current_temperatures;
    std::multimap<const Definition*, const Definition*, KeyLess> current_flux;
    std::multimap<const Definition*, const Definition*, KeyLess>
        current_element;
    // Model rows seed step1. StepDefinition retains actual input order; IDs
    // identify definitions but do not decide which current row inherits.
    std::vector<const Definition*> inputs(0);
    if (step == 1)
      for (const auto entry : *domain->ThermalConditions())
        if (entry.second->specification()->step == 0)
          inputs.push_back(entry.second);
    for (const auto* original : *definition->ThermalConditions())
      inputs.push_back(original);
    for (const auto* original : inputs) {
      const auto* s = original->specification();
      const bool model_temperature = step == 1 && s->step == 0;
      if (s->kind == ThermalConditionKind::kBodyFlux ||
          s->kind == ThermalConditionKind::kSurfaceFlux ||
          s->kind == ThermalConditionKind::kFilm ||
          s->kind == ThermalConditionKind::kRadiation) {
        if (model_temperature || s->expanded_count <= 0)
          return Failure(s, step, "Invalid distributed thermal row", d);
        current_element.emplace(original, original);
        continue;
      }
      if (s->kind != ThermalConditionKind::kTemperature &&
          s->kind != ThermalConditionKind::kConcentratedFlux)
        return Failure(s, step, "M2 thermal condition gate pending", d,
                       StatusCode::kUnsupported);
      status = Targets(domain, original, &participating, step, &ids, d);
      if (status != StatusCode::kOk) return status;
      if (s->kind == ThermalConditionKind::kConcentratedFlux) {
        if (s->first != s->last || model_temperature)
          return Failure(s, step, "Invalid concentrated thermal component", d);
        current_flux.emplace(original, original);
        continue;
      }
      if (model_temperature && s->magnitude != 0)
        return Failure(s, step, "Nonzero model temperature boundary", d);
      if (model_temperature && new_temperature) continue;
      for (int id : ids) {
        for (int label = s->first; label <= s->last; ++label) {
          const PhysicalKey key(id, label);
          if (!model_temperature && !current_temperatures.insert(key).second)
            return Failure(s, step, "Duplicate thermal temperature component",
                           d);
          std::unique_ptr<Temperature> row(new Temperature);
          row->node_id = id;
          row->label = label;
          row->origin = original;
          row->amplitude = s->amplitude;
          row->definition_step = step;
          row->definition_period = result->period_;
          row->ramp = ramp;
          row->start_from_committed_temperature = ramp && !s->amplitude;
          row->end_value = s->magnitude;
          if (s->amplitude) {
            double coefficient = 0;
            status = Coefficient(s->amplitude, step, result->period_, step,
                                 result->total_start_, result->period_,
                                 &coefficient, d);
            if (status != StatusCode::kOk) return status;
            row->end_value *= coefficient;
            if (!std::isfinite(row->end_value))
              return Failure(s, step, "Nonfinite prescribed temperature", d,
                             StatusCode::kLibraryFailure);
          }
          temperatures[key] = std::move(row);
        }
      }
    }
    for (auto& entry : temperatures) {
      result->temperatures_->push_back(entry.second.get());
      entry.second.release();
    }
    std::multimap<const Definition*, const Flux*, KeyLess> prior;
    if (previous)
      for (const auto* row : *previous->rows_)
        if (row->history != ThermalFluxHistory::kRemoved)
          prior.emplace(row->origin, row);
    const bool new_flux =
        definition->HasLoadBlock() &&
        definition->LoadOperation() == ConditionOperation::kNew;
    for (auto group = current_flux.begin(); group != current_flux.end();) {
      const auto range = current_flux.equal_range(group->first);
      const auto old = prior.equal_range(group->first);
      const bool redefined = old.first != old.second;
      if (redefined && !new_flux && std::distance(old.first, old.second) > 1)
        return Failure(group->first->specification(), step,
                       "Cannot MOD multiple thermal power rows; use NEW", d);
      double previous_end = 0;
      for (auto i = old.first; i != old.second; ++i)
        previous_end += i->second->end_value;
      if (!std::isfinite(previous_end))
        return Failure(group->first->specification(), step,
                       "Nonfinite prior thermal power", d,
                       StatusCode::kLibraryFailure);
      for (auto i = range.first; i != range.second; ++i) {
        const auto* s = i->second->specification();
        std::unique_ptr<Flux> row(new Flux);
        row->label = s->first;
        row->origin = i->second;
        row->amplitude = s->amplitude;
        row->ramp = ramp;
        row->definition_step = step;
        row->definition_period = result->period_;
        if (redefined && i == std::prev(range.second)) {
          row->history = ThermalFluxHistory::kRedefined;
          row->previous_end_value = previous_end;
          if (!s->amplitude) row->evaluation = ThermalFluxEvaluation::kEndpoint;
        }
        status =
            RowValue(row.get(), step, result->period_, result->total_start_,
                     result->period_, &row->end_value, d);
        if (status != StatusCode::kOk) return status;
        result->rows_->push_back(row.get());
        row.release();
      }
      group = range.second;
    }
    for (const auto entry : prior) {
      if (current_flux.find(entry.first) != current_flux.end()) continue;
      std::unique_ptr<Flux> row(new Flux);
      CopyFlux(entry.second, row.get());
      row->previous_end_value = entry.second->end_value;
      row->history = new_flux ? ThermalFluxHistory::kRemoved
                              : ThermalFluxHistory::kCarried;
      if (new_flux) {
        row->amplitude = nullptr;
        row->ramp = ramp;
        row->evaluation = ThermalFluxEvaluation::kEndpoint;
      }
      status = RowValue(row.get(), step, result->period_, result->total_start_,
                        result->period_, &row->end_value, d);
      if (status != StatusCode::kOk) return status;
      result->rows_->push_back(row.get());
      row.release();
    }
    for (const auto* row : *result->rows_) {
      status = Targets(domain, row->origin, &participating, step, &ids, d);
      if (status != StatusCode::kOk) return status;
      for (int id : ids) {
        std::unique_ptr<Flux> copy(new Flux);
        CopyFlux(row, copy.get());
        copy->node_id = id;
        result->fluxes_->push_back(copy.get());
        copy.release();
      }
    }
    // Distributed histories preserve original targets/labels. Expanded faces
    // never become temporal keys, and their multiple faces do not multiply
    // the number of original-row/element norm contributions.
    std::multimap<const Definition*, const ElementLoad*, KeyLess> prior_element;
    if (previous)
      for (const auto* row : *previous->element_rows_)
        if (row->history != ThermalFluxHistory::kRemoved)
          prior_element.emplace(row->origin, row);
    for (auto group = current_element.begin();
         group != current_element.end();) {
      const auto range = current_element.equal_range(group->first);
      const auto old = prior_element.equal_range(group->first);
      const bool redefined = old.first != old.second;
      const auto* spec = group->first->specification();
      const bool new_family = definition->ThermalLoadOperation(spec->family) ==
                              ConditionOperation::kNew;
      if (redefined && !new_family && std::distance(old.first, old.second) > 1)
        return Failure(spec, step,
                       "Cannot MOD multiple distributed rows; use NEW", d);
      double previous_end = 0;
      for (auto i = old.first; i != old.second; ++i)
        previous_end += i->second->end_value;
      if (!std::isfinite(previous_end))
        return Failure(spec, step, "Nonfinite prior distributed magnitude", d,
                       StatusCode::kLibraryFailure);
      int prior_begin = 0, prior_count = 0;
      const bool exchange = spec->kind == ThermalConditionKind::kFilm ||
                            spec->kind == ThermalConditionKind::kRadiation;
      if (exchange && redefined) {
        if (!result->previous_exchanges_)
          result->previous_exchanges_ =
              new std::vector<const PreviousThermalExchange*>(0);
        if (result->previous_exchanges_->size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max()))
          return Failure(spec, step, "Film prior index range", d);
        prior_begin = static_cast<int>(result->previous_exchanges_->size());
        for (auto previous_row = old.first; previous_row != old.second;
             ++previous_row) {
          if (result->previous_exchanges_->size() ==
              static_cast<std::size_t>(std::numeric_limits<int>::max()))
            return Failure(spec, step, "Film prior count range", d);
          ThermalLoadContext values;
          status = previous->PopulateExchangeContext(
              previous_row->second, previous->period_, &values, d);
          if (status != StatusCode::kOk) return status;
          std::unique_ptr<PreviousThermalExchange> record(
              new PreviousThermalExchange);
          record->origin = previous_row->second->origin;
          record->coefficient = values.kind == ThermalConditionKind::kFilm
                                    ? values.coefficient
                                    : values.emissivity;
          record->reference_temperature =
              values.kind == ThermalConditionKind::kFilm
                  ? values.sink_temperature
                  : values.ambient_temperature;
          result->previous_exchanges_->push_back(record.get());
          record.release();
          ++prior_count;
        }
      }
      for (auto i = range.first; i != range.second; ++i) {
        std::unique_ptr<ElementLoad> row(new ElementLoad);
        row->origin = i->second;
        row->amplitude = i->second->specification()->amplitude;
        row->ramp = ramp;
        row->definition_step = step;
        row->definition_period = result->period_;
        // Flux FIRST; exchange all_p uses the full prior key list on every row.
        if (redefined && (exchange || i == range.first)) {
          row->history = ThermalFluxHistory::kRedefined;
          row->previous_end_value = previous_end;
          row->previous_exchange_begin = prior_begin;
          row->previous_exchange_count = prior_count;
        }
        status =
            RowValue(row.get(), step, result->period_, result->total_start_,
                     result->period_, &row->end_value, d);
        if (status != StatusCode::kOk) return status;
        result->element_rows_->push_back(row.get());
        row.release();
      }
      group = range.second;
    }
    for (const auto& entry : prior_element) {
      if (current_element.find(entry.first) != current_element.end()) continue;
      std::unique_ptr<ElementLoad> row(new ElementLoad);
      CopyElementLoad(entry.second, row.get());
      row->previous_end_value = entry.second->end_value;
      const bool removed =
          definition->ThermalLoadOperation(
              entry.first->specification()->family) == ConditionOperation::kNew;
      if ((entry.first->specification()->kind == ThermalConditionKind::kFilm ||
           entry.first->specification()->kind ==
               ThermalConditionKind::kRadiation) &&
          removed)
        continue;  // Exchange removal is immediate, including RAMP steps.
      row->previous_exchange_begin = 0;
      row->previous_exchange_count =
          0;  // Carried rows have A=1; no old snapshot.
      row->history =
          removed ? ThermalFluxHistory::kRemoved : ThermalFluxHistory::kCarried;
      if (removed) {
        row->amplitude = nullptr;
        row->ramp = ramp;
      }
      status = RowValue(row.get(), step, result->period_, result->total_start_,
                        result->period_, &row->end_value, d);
      if (status != StatusCode::kOk) return status;
      result->element_rows_->push_back(row.get());
      row.release();
    }
    for (const auto* row : *result->element_rows_) {
      if (row->history == ThermalFluxHistory::kRemoved && !row->ramp) continue;
      const auto* spec = row->origin->specification();
      std::set<int> element_ids;
      for (int i = 0; i < spec->expanded_count; ++i)
        element_ids.insert(spec->expanded_element_ids[i]);
      for (int id : element_ids) {
        std::unique_ptr<ElementLoad> copy(new ElementLoad);
        CopyElementLoad(row, copy.get());
        copy->element_id = id;
        result->element_loads_->push_back(copy.get());
        copy.release();
      }
    }
    *out = result.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return Failure(nullptr, step, "Thermal snapshot allocation", d,
                   StatusCode::kAllocationFailure);
  }
}
}  // namespace fesa
