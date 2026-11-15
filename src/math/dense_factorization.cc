#include "fesa/math/dense_factorization.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <new>

#include "adapters/mkl/dense_ops.h"

namespace fesa {
namespace {
StatusCode Invalid(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation, operation, d);
}
StatusCode Allocation(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kAllocationFailure,
                       DiagnosticStage::kAllocation, operation, d);
}
StatusCode Backend(const char* operation, int native, Diagnostic* d,
                   bool singular = false) {
  const bool allocation = dense_ops::LapackAllocationFailure(native);
  const auto code = allocation               ? StatusCode::kAllocationFailure
                    : singular && native > 0 ? StatusCode::kSingular
                                             : StatusCode::kLibraryFailure;
  if (d)
    d->Record(code,
              allocation ? DiagnosticStage::kAllocation
                         : DiagnosticStage::kComputation,
              operation, native);
  return code;
}
StatusCode RankFailure(Diagnostic* d) {
  return ReportFailure(StatusCode::kSingular, DiagnosticStage::kComputation,
                       "DenseQrFactorization numerical column rank", d);
}
bool Finite(const Matrix* input) {
  return input &&
         dense_ops::AllFinite(input->Rows() * input->Columns(), input->Data());
}
}  // namespace

DenseLuFactorization::DenseLuFactorization(Matrix* factors, int* pivots)
    : factors_(factors), pivots_(pivots) {}
DenseLuFactorization::~DenseLuFactorization() {
  delete factors_;
  delete[] pivots_;
}
StatusCode DenseLuFactorization::Create(const Matrix* input,
                                        DenseLuFactorization** out,
                                        Diagnostic* d) {
  if (!out) return Invalid("DenseLuFactorization output", d);
  *out = nullptr;
  if (!Finite(input) || input->Rows() != input->Columns())
    return Invalid("DenseLuFactorization finite square input", d);
  Matrix* raw = nullptr;
  auto status = Matrix::Create(input->Rows(), input->Columns(), &raw, d);
  std::unique_ptr<Matrix> factors(raw);
  if (status != StatusCode::kOk) return status;
  status = factors->CopyFrom(input, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<int[]> pivots(new (std::nothrow) int[input->Rows()]);
  if (!pivots) return Allocation("DenseLuFactorization pivots", d);
  if (!dense_ops::BeginOperation()) return Backend("Dense LU preflight", 0, d);
  factors->finite_verified_ = false;
  const int native =
      dense_ops::LuFactor(input->Rows(), factors->data_, pivots.get());
  if (native != 0) return Backend("Dense LU factor", native, d, true);
  if (!Finite(factors.get()))
    return Backend("Dense LU nonfinite factors", 0, d);
  *out = new (std::nothrow) DenseLuFactorization(factors.get(), pivots.get());
  if (!*out) return Allocation("DenseLuFactorization object", d);
  factors.release();
  pivots.release();
  return StatusCode::kOk;
}
StatusCode DenseLuFactorization::Solve(const Matrix* rhs, bool transpose,
                                       Matrix* out, Diagnostic* d) const {
  if (!Finite(rhs) || !out || out == rhs || rhs->Rows() != factors_->Rows() ||
      out->Rows() != rhs->Rows() || out->Columns() != rhs->Columns())
    return Invalid("Dense LU solve shape/alias/finiteness", d);
  if (!dense_ops::BeginOperation())
    return Backend("Dense LU solve preflight", 0, d);
  auto status = out->CopyFrom(rhs, d);
  if (status != StatusCode::kOk) return status;
  out->finite_verified_ = false;
  const int native =
      dense_ops::LuSolve(factors_->Rows(), rhs->Columns(), transpose,
                         factors_->Data(), pivots_, out->data_);
  if (native != 0) return Backend("Dense LU solve", native, d);
  return Finite(out) ? StatusCode::kOk
                     : Backend("Dense LU nonfinite solution", 0, d);
}
StatusCode DenseLuFactorization::Determinant(double* out, Diagnostic* d) const {
  if (!out) return Invalid("Dense LU determinant output", d);
  const int size = factors_->Rows();
  Vector *diagonal_raw = nullptr, *logs_raw = nullptr, *ones_raw = nullptr;
  auto status = Vector::Create(size, &diagonal_raw, d);
  std::unique_ptr<Vector> diagonal(diagonal_raw);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(size, &logs_raw, d);
  std::unique_ptr<Vector> logs(logs_raw);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(size, &ones_raw, d);
  std::unique_ptr<Vector> ones(ones_raw);
  if (status != StatusCode::kOk) return status;
  int parity = 0;
  for (int i = 0; i < size; ++i) {
    // Integer pivot parity and scalar sign classification, not a diagonal
    // product.
    if (pivots_[i] != i + 1) parity ^= 1;
    if (factors_->Data()[i * size + i] < 0) parity ^= 1;
    status = ones->Set(i, 1, d);  // Constant reduction coefficients.
    if (status != StatusCode::kOk) return status;
  }
  if (!dense_ops::BeginOperation())
    return Backend("Dense determinant preflight", 0, d);
  diagonal->finite_verified_ = false;
  logs->finite_verified_ = false;
  dense_ops::CopyStrided(size, factors_->Data(), size + 1, diagonal->data_, 1);
  dense_ops::LogAbs(size, diagonal->Data(), logs->data_);
  if (!dense_ops::AllFinite(size, logs->Data()))
    return Backend("Dense determinant logarithm", 0, d);
  double exponent = 0;
  status = logs->Dot(ones.get(), &exponent, d);
  if (status != StatusCode::kOk) return status;
  if (exponent < std::log(std::numeric_limits<double>::denorm_min()) ||
      exponent > std::log(std::numeric_limits<double>::max()))
    return Backend("Dense determinant representable range", 0, d);
  const double magnitude = std::exp(exponent);
  if (!std::isfinite(magnitude) || magnitude == 0)
    return Backend("Dense determinant overflow/underflow", 0, d);
  *out = parity ? -magnitude : magnitude;
  return StatusCode::kOk;
}

DenseQrFactorization::DenseQrFactorization(Matrix* factors, Vector* tau)
    : factors_(factors), tau_(tau) {}
DenseQrFactorization::~DenseQrFactorization() {
  delete factors_;
  delete tau_;
}
StatusCode DenseQrFactorization::Create(const Matrix* input,
                                        DenseQrFactorization** out,
                                        Diagnostic* d) {
  if (!out) return Invalid("DenseQrFactorization output", d);
  *out = nullptr;
  if (!Finite(input) || input->Rows() < input->Columns())
    return Invalid("Dense QR finite tall input", d);
  Matrix* raw = nullptr;
  auto status = Matrix::Create(input->Rows(), input->Columns(), &raw, d);
  std::unique_ptr<Matrix> factors(raw);
  if (status != StatusCode::kOk) return status;
  status = factors->CopyFrom(input, d);
  if (status != StatusCode::kOk) return status;
  Vector* tau_raw = nullptr;
  status = Vector::Create(input->Columns(), &tau_raw, d);
  std::unique_ptr<Vector> tau(tau_raw);
  if (status != StatusCode::kOk) return status;
  if (!dense_ops::BeginOperation()) return Backend("Dense QR preflight", 0, d);
  const double norm =
      dense_ops::Norm2(input->Rows() * input->Columns(), input->Data());
  if (!std::isfinite(norm)) return Backend("Dense QR input norm", 0, d);
  if (norm == 0) return RankFailure(d);
  const double threshold = std::numeric_limits<double>::epsilon() *
                           std::max(input->Rows(), input->Columns()) * norm;
  if (!std::isfinite(threshold) || threshold == 0)
    return Backend("Dense QR rank threshold range", 0, d);
  factors->finite_verified_ = false;
  tau->finite_verified_ = false;
  const int native = dense_ops::QrFactor(input->Rows(), input->Columns(),
                                         factors->data_, tau->data_);
  if (native != 0) return Backend("Dense QR factor", native, d);
  if (!Finite(factors.get()) || !dense_ops::AllFinite(tau->Size(), tau->Data()))
    return Backend("Dense QR nonfinite factors", 0, d);
  for (int i = 0; i < input->Columns(); ++i)
    if (std::abs(factors->Data()[i * input->Columns() + i]) <= threshold)
      return RankFailure(d);
  *out = new (std::nothrow) DenseQrFactorization(factors.get(), tau.get());
  if (!*out) return Allocation("DenseQrFactorization object", d);
  factors.release();
  tau.release();
  return StatusCode::kOk;
}
StatusCode DenseQrFactorization::SolveLeastSquares(const Matrix* rhs,
                                                   Matrix* out,
                                                   Diagnostic* d) const {
  if (!Finite(rhs) || !out || rhs == out || rhs->Rows() != factors_->Rows() ||
      out->Rows() != factors_->Columns() || out->Columns() != rhs->Columns())
    return Invalid("Dense QR RHS/output shape/alias/finiteness", d);
  Matrix* raw = nullptr;
  auto status = Matrix::Create(rhs->Rows(), rhs->Columns(), &raw, d);
  std::unique_ptr<Matrix> work(raw);
  if (status != StatusCode::kOk) return status;
  if (!dense_ops::BeginOperation())
    return Backend("Dense QR solve preflight", 0, d);
  status = work->CopyFrom(rhs, d);
  if (status != StatusCode::kOk) return status;
  work->finite_verified_ = false;
  int native =
      dense_ops::ApplyQt(factors_->Rows(), factors_->Columns(), rhs->Columns(),
                         factors_->Data(), tau_->Data(), work->data_);
  if (native != 0) return Backend("Dense QR apply Q transpose", native, d);
  if (!Finite(work.get()))
    return Backend("Dense QR nonfinite transformed RHS", 0, d);
  status = out->CopyBlockFrom(work.get(), 0, 0, out->Rows(), out->Columns(), 0,
                              0, d);
  if (status != StatusCode::kOk) return status;
  out->finite_verified_ = false;
  native = dense_ops::SolveUpper(out->Rows(), out->Columns(), factors_->Data(),
                                 factors_->Columns(), out->data_);
  if (native != 0) return Backend("Dense QR triangular solve", native, d, true);
  return Finite(out) ? StatusCode::kOk
                     : Backend("Dense QR nonfinite solution", 0, d);
}

}  // namespace fesa
