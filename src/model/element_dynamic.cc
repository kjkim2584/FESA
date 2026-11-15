#include "fesa/model/element_dynamic.h"

#include <cmath>
#include <cstddef>
#include <limits>

namespace fesa {
namespace {
StatusCode Invalid(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation, operation, d);
}
bool InputVector(const Vector* value, int size) {
  if (!value || value->Size() != size) return false;
  for (int i = 0; i < size; ++i)
    if (!std::isfinite(value->Data()[i])) return false;
  return true;
}
bool InputMatrix(const Matrix* value, int rows, int columns) {
  if (!value || value->Rows() != rows || value->Columns() != columns)
    return false;
  for (int row = 0; row < rows; ++row)
    for (int column = 0; column < columns; ++column)
      if (!std::isfinite(
          value->Data()[static_cast<std::size_t>(row) * columns + column]))
        return false;
  return true;
}
StatusCode ValidateDynamicEvaluation(int node_count, int dofs_per_node,
                                     int internal_count,
                                     const ElementMotionInput* input,
                                     DynamicEvaluationRequest request,
                                     const ElementDynamicResult* out,
                                     Diagnostic* d) {
  if (!input || !out || node_count <= 0 ||
      (dofs_per_node != 3 && dofs_per_node != 6) ||
      node_count >
          (std::numeric_limits<int>::max() - internal_count) / dofs_per_node ||
      (request != DynamicEvaluationRequest::kForces &&
       request != DynamicEvaluationRequest::kForcesAndLinearization) ||
      (input->kinematics != DynamicKinematics::kLinear &&
       input->kinematics != DynamicKinematics::kFiniteRotation))
    return Invalid("Dynamic element metadata", d);
  const int physical_size = node_count * dofs_per_node;
  const int size = physical_size + internal_count;
  const auto* trial = input->finite_state;
  if (internal_count == 0 &&
      (input->internal_velocity || input->internal_acceleration ||
       (trial && trial->beam_internal_position)))
    return Invalid("Dynamic unexpected Beam input", d);
  if (input->kinematics == DynamicKinematics::kLinear) {
    if (internal_count != 0 || trial ||
        !InputVector(input->linear_displacement, physical_size))
      return Invalid("Dynamic linear pose", d);
  } else {
    if (dofs_per_node != 6)
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Finite solid dynamics unsupported", d);
    if (input->linear_displacement || !trial ||
        !InputVector(trial->translations, 3 * node_count) ||
        !InputMatrix(trial->orientations, node_count, 4) ||
        (trial->committed_orientations &&
         !InputMatrix(trial->committed_orientations, node_count, 4)) ||
        (trial->body_spin && !InputVector(trial->body_spin, 3 * node_count)) ||
        (internal_count != 0 &&
         (trial->incremental_solid ||
          !InputVector(trial->beam_internal_position, internal_count))))
      return Invalid("Dynamic finite pose", d);
  }
  if (!InputVector(input->velocity, physical_size) ||
      !InputVector(input->acceleration, physical_size) ||
      (internal_count != 0 &&
       (!InputVector(input->internal_velocity, internal_count) ||
        !InputVector(input->internal_acceleration, internal_count))))
    return Invalid("Dynamic velocity/acceleration", d);
  const Vector* sources[] = {input->linear_displacement,
                             input->velocity,
                             input->acceleration,
                             trial ? trial->translations : nullptr,
                             trial ? trial->body_spin : nullptr,
                             trial ? trial->beam_internal_position : nullptr,
                             input->internal_velocity,
                             input->internal_acceleration};
  const Vector* forces[] = {out->inertia_force, out->damping_force};
  for (const auto* force : forces) {
    if (!force || force->Size() != size)
      return Invalid("Dynamic force shape", d);
    for (const auto* source : sources)
      if (force == source) return Invalid("Dynamic force/input alias", d);
  }
  if (forces[0] == forces[1]) return Invalid("Dynamic force/output alias", d);
  const Matrix* derivatives[] = {out->inertia_displacement,
      out->inertia_velocity, out->inertia_acceleration,
      out->damping_displacement, out->damping_velocity};
  for (int i = 0; i < 5; ++i) {
    const auto* derivative = derivatives[i];
    if (request == DynamicEvaluationRequest::kForces) {
      if (derivative) return Invalid("Unrequested dynamic derivative", d);
      continue;
    }
    if (!derivative || derivative->Rows() != size ||
        derivative->Columns() != size)
      return Invalid("Dynamic derivative shape", d);
    for (int j = 0; j < i; ++j)
      if (derivative == derivatives[j])
        return Invalid("Dynamic derivative/output alias", d);
    if (trial && (derivative == trial->orientations ||
                  derivative == trial->committed_orientations))
      return Invalid("Dynamic derivative/input alias", d);
  }
  return StatusCode::kOk;
}
}  // namespace
StatusCode ValidateElementDynamicEvaluation(int node_count, int dofs_per_node,
                                            const ElementMotionInput* input,
                                            DynamicEvaluationRequest request,
                                            const ElementDynamicResult* out,
                                            Diagnostic* d) {
  return ValidateDynamicEvaluation(node_count, dofs_per_node, 0, input, request,
                                   out, d);
}
StatusCode ValidateBeamElementDynamicEvaluation(
    const ElementMotionInput* input, DynamicEvaluationRequest request,
    const ElementDynamicResult* out, Diagnostic* d) {
  return ValidateDynamicEvaluation(2, 6, 2, input, request, out, d);
}
}  // namespace fesa
