#include "fesa/analysis/thermal_state.h"

#include <algorithm>
#include <cassert>
#include <limits>
#include <memory>
#include <new>
#include <utility>
#include <vector>

#include "fesa/analysis/thermal_dof_map.h"
#include "fesa/model/domain.h"

namespace fesa {
namespace {
StatusCode Invalid(const char* message, Diagnostic* d,
                   StatusCode code = StatusCode::kInvalidArgument) {
  return ReportFailure(code, DiagnosticStage::kValidation, message, d);
}
}  // namespace
ThermalState::~ThermalState() {
  delete[] nodes_;
  delete[] offsets_;
  delete committed_;
  delete candidate_;
  delete committed_reaction_;
  delete candidate_reaction_;
  delete zero_;
  delete first_accepted_;
  delete flux_;
  delete rates_;
}
StatusCode ThermalState::CreateForDomain(const Domain* domain,
                                         ThermalState** out, Diagnostic* d) {
  if (!out) return Invalid("Null thermal state output", d);
  *out = nullptr;
  if (!domain || domain->Nodes()->empty())
    return Invalid("Thermal state requires nodes", d);
  const auto maximum = std::numeric_limits<int>::max();
  if (domain->Nodes()->size() >= static_cast<std::size_t>(maximum))
    return Invalid("Thermal state node count exceeds index range", d,
                   StatusCode::kInputError);
  try {
    std::unique_ptr<ThermalState> state(new ThermalState);
    state->node_count_ = static_cast<int>(domain->Nodes()->size());
    state->nodes_ = new int[state->node_count_];
    state->offsets_ = new int[state->node_count_ + 1];
    int row = 0;
    for (const auto entry : *domain->Nodes()) {
      const int count = entry.second->temperature_dof_count();
      if (count < 1 || count > 20 || count > maximum - state->global_count_)
        return Invalid("Invalid thermal state node layout", d,
                       StatusCode::kInputError);
      state->nodes_[row] = entry.first;
      state->offsets_[row++] = state->global_count_;
      state->global_count_ += count;
    }
    state->offsets_[row] = state->global_count_;
    Vector** buffers[] = {&state->committed_, &state->candidate_,
                          &state->committed_reaction_,
                          &state->candidate_reaction_, &state->zero_};
    for (auto* buffer : buffers) {
      const auto status = Vector::Create(state->global_count_, buffer, d);
      if (status != StatusCode::kOk) return status;
    }
    for (const auto& condition : *domain->ThermalConditions()) {
      const auto* spec = condition.second->specification();
      if (spec->kind == ThermalConditionKind::kFilm && spec->film_property) {
        const auto status =
            Vector::Create(state->global_count_, &state->first_accepted_, d);
        if (status != StatusCode::kOk) return status;
        break;
      }
    }
    std::vector<int> active(0);
    auto status = CollectThermalNodeIds(domain, &active, d);
    if (status != StatusCode::kOk) return status;
    int indices[20];
    for (int i = 0; i < 20; ++i) indices[i] = i;
    for (const auto entry : *domain->ThermalInitialTemperatures()) {
      const auto* initial = entry.second;
      if (!std::binary_search(active.begin(), active.end(), entry.first)) {
        if (d) {
          d->source_file = initial->source_file();
          d->source_line = initial->source_line();
          d->node_id = entry.first;
          d->path_kind = DiagnosticPathKind::kLiteralFile;
        }
        return Invalid("Thermal initial condition node is not participating", d,
                       StatusCode::kInputError);
      }
      const auto* found = std::lower_bound(
          state->nodes_, state->nodes_ + state->node_count_, entry.first);
      const int n = static_cast<int>(found - state->nodes_);
      if (found == state->nodes_ + state->node_count_ ||
          *found != entry.first ||
          initial->values()->Size() !=
              state->offsets_[n + 1] - state->offsets_[n])
        return Invalid("Thermal initial condition layout mismatch", d,
                       StatusCode::kInputError);
      status = state->committed_->GatherFrom(initial->values(), indices,
                                             initial->values()->Size(),
                                             state->offsets_[n], d);
      if (status != StatusCode::kOk) return status;
    }
    status = HeatFluxHistory::Create(&state->flux_, d);
    if (status == StatusCode::kOk)
      status = RateHistory::Create(&state->rates_, d);
    if (status != StatusCode::kOk) return status;
    *out = state.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return Invalid("Thermal state allocation", d,
                   StatusCode::kAllocationFailure);
  }
}
StatusCode ThermalState::ValidateMapping(const ThermalDofMap* map,
                                         Diagnostic* d) const {
  if (!map || map->NodeCount() != node_count_ ||
      map->GlobalCount() != global_count_)
    return Invalid("Thermal state mapping size mismatch", d);
  for (int row = 0; row < node_count_; ++row)
    if (map->NodeId(row) != nodes_[row] ||
        map->NodeDofCount(nodes_[row]) != offsets_[row + 1] - offsets_[row])
      return Invalid("Thermal state mapping layout mismatch", d);
  return StatusCode::kOk;
}
StatusCode ThermalState::ValidateDomainLayout(const Domain* domain,
                                              Diagnostic* d) const {
  if (!domain ||
      domain->Nodes()->size() != static_cast<std::size_t>(node_count_))
    return Invalid("Thermal state domain node count mismatch", d);
  int row = 0;
  for (const auto entry : *domain->Nodes()) {
    if (entry.first != nodes_[row] || entry.second->temperature_dof_count() !=
                                          offsets_[row + 1] - offsets_[row])
      return Invalid("Thermal state domain layout mismatch", d);
    ++row;
  }
  return StatusCode::kOk;
}
int ThermalState::GlobalCount() const { return global_count_; }
StatusCode ThermalState::PrepareCandidate(Diagnostic* d) {
  DiscardCandidate();
  auto status = candidate_->CopyFrom(committed_, d);
  if (status == StatusCode::kOk)
    status = candidate_reaction_->CopyFrom(zero_, d);
  if (status == StatusCode::kOk) candidate_ready_ = true;
  return status;
}
bool ThermalState::HasValidCandidate() const { return candidate_ready_; }
Vector* ThermalState::CandidateTemperature() {
  return candidate_ready_ ? candidate_ : nullptr;
}
const Vector* ThermalState::CandidateTemperature() const {
  return candidate_ready_ ? candidate_ : nullptr;
}
Vector* ThermalState::CandidateReactionFlux() {
  return candidate_ready_ ? candidate_reaction_ : nullptr;
}
const Vector* ThermalState::CandidateReactionFlux() const {
  return candidate_ready_ ? candidate_reaction_ : nullptr;
}
const Vector* ThermalState::CommittedTemperature() const { return committed_; }
const Vector* ThermalState::CommittedReactionFlux() const {
  return committed_reaction_;
}
StatusCode ThermalState::StageFirstAcceptedTemperature(Diagnostic* d) {
  if (!candidate_ready_)
    return Invalid("Film freeze needs a prepared candidate", d);
  if (!first_accepted_ || first_accepted_set_) return StatusCode::kOk;
  first_accepted_staged_ = false;
  const auto status = first_accepted_->CopyFrom(candidate_, d);
  if (status == StatusCode::kOk) first_accepted_staged_ = true;
  return status;
}
const Vector* ThermalState::FirstAcceptedTemperature() const {
  return first_accepted_set_ ? first_accepted_ : nullptr;
}
HeatFluxHistory* ThermalState::FluxHistory() { return flux_; }
const HeatFluxHistory* ThermalState::FluxHistory() const { return flux_; }
RateHistory* ThermalState::IncrementRates() { return rates_; }
const RateHistory* ThermalState::IncrementRates() const { return rates_; }
StatusCode ThermalState::BeginStep(int step, const FieldControls* c,
                                   Diagnostic* d) {
  if (candidate_ready_)
    return Invalid("Thermal step has a prepared attempt", d);
  const auto status = flux_->BeginStep(step, c, d);
  if (status != StatusCode::kOk) return status;
  rates_->BeginStep();
  step_ = step;
  accepted_increment_ = 0;
  first_accepted_set_ = false;
  first_accepted_staged_ = false;
  return StatusCode::kOk;
}
void ThermalState::DiscardCandidate() noexcept {
  candidate_ready_ = false;
  first_accepted_staged_ = false;
  flux_->DiscardCandidate();
  rates_->DiscardCandidate();
}
void ThermalState::MarkResultsUnconfirmed(int step, int increment) {
  assert(candidate_ready_ && step == step_ && increment > 0 &&
         increment == accepted_increment_ + 1 &&
         (!first_accepted_ || first_accepted_set_ || first_accepted_staged_));
  result_state_ = ResultState::kUnconfirmed;
  result_step_ = step;
  result_increment_ = increment;
}
void ThermalState::Commit(int step, int increment) noexcept {
  assert(candidate_ready_ && flux_->CanCommit() && step == step_ &&
         increment > 0 && increment == accepted_increment_ + 1 &&
         result_state_ == ResultState::kUnconfirmed && result_step_ == step &&
         result_increment_ == increment);
  assert(!first_accepted_ || first_accepted_set_ || first_accepted_staged_);
  if (first_accepted_ && !first_accepted_set_) first_accepted_set_ = true;
  first_accepted_staged_ = false;
  std::swap(candidate_, committed_);
  std::swap(candidate_reaction_, committed_reaction_);
  flux_->Commit();
  if (rates_->CanCommit()) rates_->Commit();
  accepted_increment_ = increment;
  candidate_ready_ = false;
  result_state_ = ResultState::kConfirmed;
}
ResultState ThermalState::Results() const { return result_state_; }
int ThermalState::ResultStepNumber() const { return result_step_; }
int ThermalState::ResultIncrementNumber() const { return result_increment_; }
}  // namespace fesa
