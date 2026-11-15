#include "fesa/math/vector.h"

#include <cmath>
#include <new>

#include "adapters/mkl/dense_ops.h"

namespace fesa {
namespace {

StatusCode Invalid(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation, operation, d);
}

StatusCode BackendFailure(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kLibraryFailure,
                       DiagnosticStage::kComputation, operation, d);
}

}  // namespace

Vector::Vector(int size, double* data) : size_(size), data_(data) {}
Vector::~Vector() { dense_ops::Release(data_); }
bool Vector::IsFinite() const {
  return finite_verified_ || dense_ops::AllFinite(size_, data_);
}
bool Vector::VerifyFinite() {
  finite_verified_ = dense_ops::AllFinite(size_, data_);
  return finite_verified_;
}

StatusCode Vector::Create(int size, Vector** out, Diagnostic* d) {
  if (out == nullptr) return Invalid("Vector::Create", d);
  *out = nullptr;
  if (!dense_ops::ValidDimensions(size, 1)) return Invalid("Vector::Create", d);
  double* data = dense_ops::Allocate(static_cast<std::size_t>(size));
  if (data != nullptr) {
    if (dense_ops::AllowObjectAllocation()) {
      *out = new (std::nothrow) Vector(size, data);
    }
    if (*out != nullptr) return StatusCode::kOk;
    dense_ops::Release(data);
  }
  return ReportFailure(StatusCode::kAllocationFailure,
                       DiagnosticStage::kAllocation, "Vector::Create", d);
}

StatusCode Vector::Get(int index, double* out, Diagnostic* d) const {
  if (index < 0 || index >= size_ || out == nullptr) {
    return Invalid("Vector::Get", d);
  }
  *out = data_[index];
  return StatusCode::kOk;
}

StatusCode Vector::Set(int index, double value, Diagnostic* d) {
  if (index < 0 || index >= size_ || !std::isfinite(value)) {
    return Invalid("Vector::Set", d);
  }
  data_[index] = value;
  return StatusCode::kOk;
}

StatusCode Vector::CopyFrom(const Vector* input, Diagnostic* d) {
  if (input == nullptr || input->size_ != size_) {
    return Invalid("Vector::CopyFrom", d);
  }
  if (input == this) return StatusCode::kOk;
  if (!dense_ops::BeginOperation()) {
    return BackendFailure("Vector::CopyFrom", d);
  }
  finite_verified_ = false;
  dense_ops::Copy(size_, input->data_, data_);
  finite_verified_ = input->finite_verified_;
  return StatusCode::kOk;
}

StatusCode Vector::AddScaled(double alpha, const Vector* x, Diagnostic* d) {
  if (x == nullptr || x == this || x->size_ != size_ || !std::isfinite(alpha)) {
    return Invalid("Vector::AddScaled", d);
  }
  if (!dense_ops::BeginOperation()) {
    return BackendFailure("Vector::AddScaled", d);
  }
  finite_verified_ = false;
  dense_ops::AddScaled(size_, alpha, x->data_, data_);
  return VerifyFinite() ? StatusCode::kOk
                        : BackendFailure("Vector::AddScaled", d);
}

StatusCode Vector::Scale(double alpha, Diagnostic* d) {
  if (!std::isfinite(alpha)) return Invalid("Vector::Scale", d);
  if (!dense_ops::BeginOperation()) return BackendFailure("Vector::Scale", d);
  finite_verified_ = false;
  dense_ops::Scale(size_, alpha, data_);
  return VerifyFinite() ? StatusCode::kOk : BackendFailure("Vector::Scale", d);
}

StatusCode Vector::Dot(const Vector* x, double* out, Diagnostic* d) const {
  if (x == nullptr || x->size_ != size_ || out == nullptr) {
    return Invalid("Vector::Dot", d);
  }
  if (!dense_ops::BeginOperation()) return BackendFailure("Vector::Dot", d);
  const double value = dense_ops::Dot(size_, data_, x->data_);
  if (!std::isfinite(value)) return BackendFailure("Vector::Dot", d);
  *out = value;
  return StatusCode::kOk;
}

StatusCode Vector::Norm2(double* out, Diagnostic* d) const {
  if (out == nullptr) return Invalid("Vector::Norm2", d);
  if (!dense_ops::BeginOperation()) return BackendFailure("Vector::Norm2", d);
  const double value = dense_ops::Norm2(size_, data_);
  if (!std::isfinite(value)) return BackendFailure("Vector::Norm2", d);
  *out = value;
  return StatusCode::kOk;
}

const double* Vector::Data() const { return data_; }
int Vector::Size() const { return size_; }

StatusCode Vector::GatherFrom(const Vector* source, const int* indices,
                              int count, int offset, Diagnostic* d) {
  if (!source || source == this || count < 0 || offset < 0 || offset > size_ ||
      count > size_ - offset || (count > 0 && !indices))
    return Invalid("Vector::GatherFrom", d);
  for (int i = 0; i < count; ++i)
    if (indices[i] < 0 || indices[i] >= source->size_)
      return Invalid("Vector::GatherFrom", d);
  if (count == 0) return StatusCode::kOk;
  if (!dense_ops::BeginOperation())
    return BackendFailure("Vector::GatherFrom", d);
  const bool verified = finite_verified_ && source->finite_verified_;
  finite_verified_ = false;
  dense_ops::Gather(count, source->data_, indices, data_ + offset);
  finite_verified_ = verified;
  return StatusCode::kOk;
}
StatusCode Vector::ScatterAdd(const Vector* source, const int* indices,
                              int count, Diagnostic* d) {
  if (!source || source == this || count < 0 || count > source->size_ ||
      (count > 0 && !indices))
    return Invalid("Vector::ScatterAdd", d);
  for (int i = 0; i < count; ++i)
    if (indices[i] < 0 || indices[i] >= size_)
      return Invalid("Vector::ScatterAdd", d);
  if (count == 0) return StatusCode::kOk;
  if (!dense_ops::BeginOperation())
    return BackendFailure("Vector::ScatterAdd", d);
  finite_verified_ = false;
  const int status =
      dense_ops::ScatterAdd(size_, count, indices, source->data_, data_);
  if (status == -1 || status == 2)
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation, "Vector::ScatterAdd", d);
  if (status != 0 || !VerifyFinite())
    return BackendFailure("Vector::ScatterAdd", d);
  return StatusCode::kOk;
}
StatusCode Vector::AbsoluteSum(double* out, Diagnostic* d) const {
  if (!out) return Invalid("Vector::AbsoluteSum", d);
  if (!dense_ops::BeginOperation())
    return BackendFailure("Vector::AbsoluteSum", d);
  const double value = dense_ops::AbsoluteSum(size_, data_);
  if (!std::isfinite(value)) return BackendFailure("Vector::AbsoluteSum", d);
  *out = value;
  return StatusCode::kOk;
}
StatusCode Vector::MaximumAbs(double* out, int* index, Diagnostic* d) const {
  if (!out || !index) return Invalid("Vector::MaximumAbs", d);
  if (!dense_ops::BeginOperation())
    return BackendFailure("Vector::MaximumAbs", d);
  *index = dense_ops::MaximumAbsIndex(size_, data_);
  *out = std::abs(data_[*index]);
  return StatusCode::kOk;
}

}  // namespace fesa
