#include "fesa/analysis/rotation_prescription.h"

#include <cmath>
#include <limits>
#include <memory>

#include "fesa/math/rotation.h"
#include "model/model_storage.h"

namespace fesa {

StatusCode EvaluateRotationPrescription(const Vector* q_step_start,
                                        const Vector* command_start,
                                        const Vector* command_end,
                                        double lambda, double committed_lambda,
                                        Vector* out, bool* branch_exceeded,
                                        Diagnostic* d) {
  if (branch_exceeded) *branch_exceeded = false;
  if (!q_step_start || !command_start || !command_end || !out ||
      !branch_exceeded || q_step_start->Size() != 4 ||
      command_start->Size() != 3 || command_end->Size() != 3 ||
      out->Size() != 4 || out == q_step_start || out == command_start ||
      out == command_end || !std::isfinite(lambda) ||
      !std::isfinite(committed_lambda) || committed_lambda < 0.0 ||
      lambda < committed_lambda || lambda > 1.0)
    return model_storage::Invalid("Rotation prescription inputs/range/alias",
                                  d);
  const Vector* inputs[3] = {q_step_start, command_start, command_end};
  for (const Vector* input : inputs)
    for (int component = 0; component < input->Size(); ++component)
      if (!std::isfinite(input->Data()[component]))
        return model_storage::Invalid("Rotation prescription nonfinite input",
                                      d);
  std::unique_ptr<Vector> owners[4];
  const int sizes[4] = {4, 3, 3, 4};
  for (int i = 0; i < 4; ++i) {
    Vector* raw = nullptr;
    const auto status = Vector::Create(sizes[i], &raw, d);
    owners[i].reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  auto status = NormalizeQuaternion(q_step_start, owners[0].get(), d);
  if (status != StatusCode::kOk) return status;
  Vector* delta = owners[1].get();
  Vector* path = owners[2].get();
  status = delta->CopyFrom(command_end, d);
  if (status != StatusCode::kOk) return status;
  status = delta->AddScaled(-1.0, command_start, d);
  if (status != StatusCode::kOk) return status;
  status = path->CopyFrom(delta, d);
  if (status != StatusCode::kOk) return status;
  status = path->Scale(lambda - committed_lambda, d);
  if (status != StatusCode::kOk) return status;
  double angle = 0.0;
  status = path->Norm2(&angle, d);
  if (status != StatusCode::kOk) return status;
  constexpr double kPi = 3.141592653589793238462643383279502884;
  if (angle >= kPi ||
      std::cos(angle / 2.0) <= 256 * std::numeric_limits<double>::epsilon()) {
    *branch_exceeded = true;
    return ReportFailure(StatusCode::kEquilibriumFailure,
                         DiagnosticStage::kComputation,
                         "Rotation prescription path crosses pi", d);
  }
  // B uses the fixed step's total command difference, not a new increment
  // composed onto the previous accepted orientation. Cutback changes the path
  // guard above but cannot change the target at the same lambda.
  status = path->CopyFrom(delta, d);
  if (status != StatusCode::kOk) return status;
  status = path->Scale(lambda, d);
  if (status != StatusCode::kOk) return status;
  status = RotationExp(path, owners[3].get(), d);
  if (status != StatusCode::kOk) return status;
  return ComposeRotations(owners[3].get(), owners[0].get(), out, d);
}

}  // namespace fesa
