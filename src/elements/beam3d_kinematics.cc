#include "elements/beam3d_kinematics.h"

#include <cassert>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <memory>
#include <new>

#include "elements/beam3d_kinematic_graph.h"

namespace fesa {
namespace {
// Four exact blocks, not a coordinate Hessian or a retained AD expression tape.
enum Block { kValue, kVirtual, kDirection, kMixed, kBlockCount };
constexpr int kVariables = 14;
constexpr int kFieldCount = 6;
constexpr int kSlots = 64;
using beam3d_graph::Geometry;
using beam3d_graph::Invalid;
using beam3d_graph::Valid;
StatusCode Allocation(Diagnostic* d) {
  return ReportFailure(StatusCode::kAllocationFailure,
                       DiagnosticStage::kAllocation,
                       "BeamKinematics workspace allocation", d);
}
// Owns only numerical matrices; every copy is explicit and library-backed.
struct Differential {
  Matrix* block[kBlockCount] = {};
  Differential() = default;
  Differential(const Differential&) = delete;
  Differential& operator=(const Differential&) = delete;
  Differential(Differential&&) = delete;
  Differential& operator=(Differential&&) = delete;
  ~Differential() {
    for (Matrix* matrix : block) delete matrix;
  }
  StatusCode Initialize(int rows, Diagnostic* d) {
    for (int i = 0; i < kBlockCount; ++i) {
      const StatusCode status = Matrix::Create(
          rows, i == kVirtual || i == kMixed ? kVariables : 1, &block[i], d);
      if (status != StatusCode::kOk) return status;
    }
    return StatusCode::kOk;
  }
};

// Fixed scratch capacity for this single expression, reused at every point.
// Slots hold numerical values only; no graph or model state is retained.
struct Arena {
  Differential* scalar[kSlots] = {};
  Differential* vector[kSlots] = {};
  Differential* zero_scalar = nullptr;
  Differential* zero_vector = nullptr;
  Differential* strain = nullptr;
  Matrix* skew[4] = {};
  Matrix* selector[4] = {};  // xa, wa, xb, wb selectors, all3x14.
  Matrix* node_selector[2] = {};
  Matrix* direction_row = nullptr;
  int scalar_count = 0;
  int vector_count = 0;
  Arena() = default;
  Arena(const Arena&) = delete;
  Arena& operator=(const Arena&) = delete;
  Arena(Arena&&) = delete;
  Arena& operator=(Arena&&) = delete;
  ~Arena() {
    for (int i = 0; i < kSlots; ++i) {
      delete scalar[i];
      delete vector[i];
    }
    delete zero_scalar;
    delete zero_vector;
    delete strain;
    for (Matrix* matrix : skew) delete matrix;
    for (Matrix* matrix : selector) delete matrix;
    for (Matrix* matrix : node_selector) delete matrix;
    delete direction_row;
  }
  StatusCode CreateDifferential(int rows, Differential** out, Diagnostic* d) {
    *out = new (std::nothrow) Differential;
    if (!*out) return Allocation(d);
    return (*out)->Initialize(rows, d);
  }
  StatusCode Initialize(Diagnostic* d) {
    StatusCode status = StatusCode::kOk;
    for (int i = 0; i < kSlots && status == StatusCode::kOk; ++i) {
      status = CreateDifferential(1, &scalar[i], d);
      if (status == StatusCode::kOk)
        status = CreateDifferential(3, &vector[i], d);
    }
    if (status == StatusCode::kOk)
      status = CreateDifferential(1, &zero_scalar, d);
    if (status == StatusCode::kOk)
      status = CreateDifferential(3, &zero_vector, d);
    if (status == StatusCode::kOk) status = CreateDifferential(4, &strain, d);
    for (int i = 0; i < 4 && status == StatusCode::kOk; ++i) {
      status = Matrix::Create(3, 3, &skew[i], d);
      if (status == StatusCode::kOk)
        status = Matrix::Create(3, kVariables, &selector[i], d);
      for (int j = 0; j < 3 && status == StatusCode::kOk; ++j)
        status = selector[i]->Set(j, 3 * i + j, 1, d);
    }
    for (int i = 0; i < 2 && status == StatusCode::kOk; ++i) {
      status = Matrix::Create(2, 1, &node_selector[i], d);
      if (status == StatusCode::kOk) status = node_selector[i]->Set(i, 0, 1, d);
    }
    if (status == StatusCode::kOk)
      status = Matrix::Create(1, kVariables, &direction_row, d);
    return status;
  }
};

// A call-local evaluator with sticky error propagation. The low-level wrappers
// skip all later operations after failure, preserving the original diagnostic.
class DifferentialOps final {
 public:
  DifferentialOps(Arena* arena, int blocks, bool reuse, Diagnostic* d)
      : arena_(arena), d_(d), blocks_(blocks), first_(reuse ? 2 : 0) {
    arena_->scalar_count = 0;
    arena_->vector_count = 0;
  }
  DifferentialOps(const DifferentialOps&) = delete;
  DifferentialOps& operator=(const DifferentialOps&) = delete;
  DifferentialOps(DifferentialOps&&) = delete;
  DifferentialOps& operator=(DifferentialOps&&) = delete;
  StatusCode status() const { return status_; }
  bool Writes(int block) const { return block >= first_ && block < blocks_; }
  double Value(const Differential* a) const {
    return a->block[kValue]->Data()[0];
  }
  void Copy(Matrix* out, const Matrix* a) {
    if (status_ == StatusCode::kOk) status_ = out->CopyFrom(a, d_);
  }
  void ScaleMatrix(Matrix* out, double coefficient) {
    if (status_ == StatusCode::kOk) status_ = out->Scale(coefficient, d_);
  }
  void AddMatrix(Matrix* out, double coefficient, const Matrix* a) {
    if (status_ == StatusCode::kOk)
      status_ = out->AddScaled(coefficient, a, d_);
  }
  void Set(Matrix* out, int row, int column, double value) {
    if (status_ == StatusCode::kOk) status_ = out->Set(row, column, value, d_);
  }
  void Multiply(Matrix* out, const Matrix* a, bool ta, const Matrix* b, bool tb,
                double alpha = 1, double beta = 0) {
    if (status_ == StatusCode::kOk)
      status_ = Matrix::MultiplyAdd(a, ta, b, tb, alpha, beta, out, d_);
  }
  void Direction(const Vector* direction) {
    if (status_ == StatusCode::kOk)
      status_ = arena_->direction_row->SetRow(0, direction, d_);
  }
  Differential* Next(int rows) {
    assert(rows == 1 || rows == 3);
    if (rows == 1) {
      assert(arena_->scalar_count < kSlots);
      return arena_->scalar[arena_->scalar_count++];
    }
    assert(arena_->vector_count < kSlots);
    return arena_->vector[arena_->vector_count++];
  }
  void Clear(Differential* out) {
    const Differential* zero =
        out->block[0]->Rows() == 1 ? arena_->zero_scalar : arena_->zero_vector;
    for (int i = first_; i < blocks_; ++i) Copy(out->block[i], zero->block[i]);
  }
  void Skew(const Matrix* a, Matrix* out) {
    if (status_ != StatusCode::kOk) return;
    // Coefficient assembly of [a]x; its constant zero diagonal is never
    // changed.
    const double* v = a->Data();
    Set(out, 0, 1, -v[2]);
    Set(out, 0, 2, v[1]);
    Set(out, 1, 0, v[2]);
    Set(out, 1, 2, -v[0]);
    Set(out, 2, 0, -v[1]);
    Set(out, 2, 1, v[0]);
  }
  Differential* Seed(const Matrix* rows, int node, bool angular) {
    Differential* out = Next(3);
    Clear(out);
    if (Writes(kValue))
      Multiply(out->block[kValue], rows, true, arena_->node_selector[node], false);
    const Matrix* selector = arena_->selector[2 * node + (angular ? 1 : 0)];
    if (Writes(kVirtual) && angular) {
      Skew(out->block[kValue], arena_->skew[0]);
      Multiply(out->block[kVirtual], arena_->skew[0], false, selector, false,
               -1);
    } else if (Writes(kVirtual)) {
      Copy(out->block[kVirtual], selector);
    }
    if (Writes(kDirection))
      Multiply(out->block[kDirection], out->block[kVirtual], false,
               arena_->direction_row, true);
    if (Writes(kMixed) && angular) {
      // D(delta a) = w_virtual x (w_direction x a), in this order.
      Skew(out->block[kDirection], arena_->skew[0]);
      Multiply(out->block[kMixed], arena_->skew[0], false, selector, false, -1);
    }
    return out;
  }
  Differential* LengthSeed(double value, int node) {
    Differential* out = Next(1);
    Clear(out);
    if (Writes(kValue)) Set(out->block[kValue], 0, 0, value);
    if (Writes(kVirtual)) Set(out->block[kVirtual], 0, 12 + node, 1);
    if (Writes(kDirection))
      Set(out->block[kDirection], 0, 0, arena_->direction_row->Data()[12 + node]);
    return out;
  }
  Differential* Scale(const Differential* a, double coefficient) {
    Differential* out = Next(a->block[0]->Rows());
    for (int i = first_; i < blocks_; ++i) {
      Copy(out->block[i], a->block[i]);
      ScaleMatrix(out->block[i], coefficient);
    }
    return out;
  }
  Differential* Add(const Differential* a, const Differential* b,
                    double coefficient = 1) {
    Differential* out = Next(a->block[0]->Rows());
    for (int i = first_; i < blocks_; ++i) {
      Copy(out->block[i], a->block[i]);
      AddMatrix(out->block[i], coefficient, b->block[i]);
    }
    return out;
  }
  Differential* Offset(const Differential* a, double constant) {
    Differential* out = Scale(a, 1);
    if (Writes(kValue) && status_ == StatusCode::kOk)
      Set(out->block[kValue], 0, 0, Value(a) + constant);
    return out;
  }
  Differential* Combination(const Differential* a, const Differential* b,
                            const Differential* c, const Differential* e,
                            double ca, double cb, double cc, double ce) {
    Differential* out = Scale(a, ca);
    for (int i = first_; i < blocks_; ++i) {
      AddMatrix(out->block[i], cb, b->block[i]);
      AddMatrix(out->block[i], cc, c->block[i]);
      AddMatrix(out->block[i], ce, e->block[i]);
    }
    return out;
  }
  Differential* Product(const Differential* a, const Differential* scalar) {
    Differential* out = Next(a->block[0]->Rows());
    if (status_ != StatusCode::kOk) return out;
    const double value = Value(scalar);
    const double direction =
        Writes(kDirection) ? scalar->block[kDirection]->Data()[0] : 0;
    for (int i = first_; i < blocks_; ++i) {
      Copy(out->block[i], a->block[i]);
      ScaleMatrix(out->block[i], value);
    }
    if (Writes(kVirtual))
      Multiply(out->block[kVirtual], a->block[kValue], false,
               scalar->block[kVirtual], false, 1, 1);
    if (!Writes(kDirection)) return out;
    AddMatrix(out->block[kDirection], direction, a->block[kValue]);
    AddMatrix(out->block[kMixed], direction, a->block[kVirtual]);
    Multiply(out->block[kMixed], a->block[kDirection], false,
             scalar->block[kVirtual], false, 1, 1);
    Multiply(out->block[kMixed], a->block[kValue], false, scalar->block[kMixed],
             false, 1, 1);
    return out;
  }
  Differential* Dot(const Differential* a, const Differential* b) {
    Differential* out = Next(1);
    if (Writes(kValue))
      Multiply(out->block[kValue], a->block[kValue], true, b->block[kValue],
               false);
    if (Writes(kVirtual)) {
      Multiply(out->block[kVirtual], a->block[kValue], true, b->block[kVirtual],
               false);
      Multiply(out->block[kVirtual], b->block[kValue], true, a->block[kVirtual],
               false, 1, 1);
    }
    if (!Writes(kDirection)) return out;
    Multiply(out->block[kDirection], a->block[kValue], true,
             b->block[kDirection], false);
    Multiply(out->block[kDirection], b->block[kValue], true,
             a->block[kDirection], false, 1, 1);
    Multiply(out->block[kMixed], a->block[kDirection], true, b->block[kVirtual],
             false);
    Multiply(out->block[kMixed], a->block[kValue], true, b->block[kMixed],
             false, 1, 1);
    Multiply(out->block[kMixed], b->block[kDirection], true, a->block[kVirtual],
             false, 1, 1);
    Multiply(out->block[kMixed], b->block[kValue], true, a->block[kMixed],
             false, 1, 1);
    return out;
  }
  Differential* Cross(const Differential* a, const Differential* b) {
    Differential* out = Next(3);
    Skew(a->block[kValue], arena_->skew[0]);
    Skew(b->block[kValue], arena_->skew[1]);
    if (Writes(kMixed)) {
      Skew(a->block[kDirection], arena_->skew[2]);
      Skew(b->block[kDirection], arena_->skew[3]);
    }
    if (Writes(kValue))
      Multiply(out->block[kValue], arena_->skew[0], false, b->block[kValue],
               false);
    for (int i : {kVirtual, kDirection, kMixed}) {
      if (!Writes(i)) continue;
      Multiply(out->block[i], arena_->skew[0], false, b->block[i], false);
      Multiply(out->block[i], arena_->skew[1], false, a->block[i], false, -1,
               1);
    }
    if (!Writes(kMixed)) return out;
    Multiply(out->block[kMixed], arena_->skew[2], false, b->block[kVirtual],
             false, 1, 1);
    Multiply(out->block[kMixed], arena_->skew[3], false, a->block[kVirtual],
             false, -1, 1);
    return out;
  }
  Differential* Chain(const Differential* a, double value, double first,
                      double second) {
    Differential* out = Next(1);
    if (status_ != StatusCode::kOk) return out;
    const double direction =
        Writes(kDirection) ? a->block[kDirection]->Data()[0] : 0;
    if (Writes(kValue)) Set(out->block[kValue], 0, 0, value);
    if (Writes(kVirtual)) {
      Copy(out->block[kVirtual], a->block[kVirtual]);
      ScaleMatrix(out->block[kVirtual], first);
    }
    if (!Writes(kDirection)) return out;
    Set(out->block[kDirection], 0, 0, first * direction);
    Copy(out->block[kMixed], a->block[kMixed]);
    ScaleMatrix(out->block[kMixed], first);
    AddMatrix(out->block[kMixed], second * direction, a->block[kVirtual]);
    return out;
  }
  Differential* Reciprocal(const Differential* a) {
    const double inverse = status_ == StatusCode::kOk ? 1 / Value(a) : 0;
    return Chain(a, inverse, -inverse * inverse,
                 2 * inverse * inverse * inverse);
  }
  Differential* SquareRoot(const Differential* a) {
    const double root = status_ == StatusCode::kOk ? std::sqrt(Value(a)) : 1;
    return Chain(a, root, .5 / root, -.25 / (root * root * root));
  }
  Differential* Sin(const Differential* a) {
    const double value = status_ == StatusCode::kOk ? Value(a) : 0;
    return Chain(a, std::sin(value), std::cos(value), -std::sin(value));
  }
  Differential* Cos(const Differential* a) {
    const double value = status_ == StatusCode::kOk ? Value(a) : 0;
    return Chain(a, std::cos(value), -std::sin(value), -std::cos(value));
  }
  Differential* Angle(const Differential* y, const Differential* x) {
    Differential* out = Next(1);
    if (status_ != StatusCode::kOk) return out;
    const double xv = Value(x), yv = Value(y);
    const double xd = Writes(kDirection) ? x->block[kDirection]->Data()[0] : 0;
    const double yd = Writes(kDirection) ? y->block[kDirection]->Data()[0] : 0;
    const double r2 = xv * xv + yv * yv;
    const double fx = -yv / r2, fy = xv / r2;
    const double fxx = 2 * xv * yv / (r2 * r2);
    const double fxy = (yv * yv - xv * xv) / (r2 * r2);
    if (Writes(kValue)) Set(out->block[kValue], 0, 0, std::atan2(yv, xv));
    if (Writes(kDirection)) Set(out->block[kDirection], 0, 0, fx * xd + fy * yd);
    for (int i : {kVirtual, kMixed}) {
      if (!Writes(i)) continue;
      Copy(out->block[i], x->block[i]);
      ScaleMatrix(out->block[i], fx);
      AddMatrix(out->block[i], fy, y->block[i]);
    }
    if (!Writes(kMixed)) return out;
    AddMatrix(out->block[kMixed], fxx * xd + fxy * yd, x->block[kVirtual]);
    AddMatrix(out->block[kMixed], fxy * xd - fxx * yd, y->block[kVirtual]);
    return out;
  }
  Differential* Strain() const { return arena_->strain; }
  void Pack(const Differential* a, int row) {
    for (int i = first_; i < blocks_ && status_ == StatusCode::kOk; ++i)
      status_ = arena_->strain->block[i]->CopyBlockFrom(
          a->block[i], 0, 0, 1, a->block[i]->Columns(), row, 0, d_);
  }

 private:
  Arena* arena_;  // Borrowed for this call; never stored by any model object.
  Diagnostic* d_;
  int blocks_;  // Requested contiguous prefix: value1, first2, full4.
  int first_;   // Reused primal/virtual blocks are read-only this call.
  StatusCode status_ = StatusCode::kOk;
};

}  // namespace

struct BeamKinematicWorkspace::Storage {
  Arena* arena = nullptr;
  Differential* fields[kFieldCount] = {};
  bool valid = false;
  int blocks = 0;
  bool base_valid = false;
  bool has_relative = false;
  int base_blocks = 0;
  double length = 0, xi = 0;
  Matrix* key[4] = {};  // Owned positions, tangents, normals, relative positions.
  Vector* eta = nullptr;  // Owned numeric snapshot; no borrowed input retained.
  Storage() = default;
  Storage(const Storage&) = delete;
  Storage& operator=(const Storage&) = delete;
  Storage(Storage&&) = delete;
  Storage& operator=(Storage&&) = delete;
  ~Storage() {
    delete arena;
    for (auto* matrix : key) delete matrix;
    delete eta;
  }
  bool Matches(const BeamKinematicInput* input, double point) const {
    if (!base_valid || has_relative != (input->relative_positions != nullptr) ||
        std::memcmp(&length, &input->length, sizeof(double)) ||
        std::memcmp(&xi, &point, sizeof(double))) return false;
    const Matrix* values[] = {input->positions, input->tangents, input->normals,
                              input->relative_positions};
    for (int i = 0; i < 4; ++i)
      if (values[i] && std::memcmp(key[i]->Data(), values[i]->Data(),
                                    6 * sizeof(double))) return false;
    return std::memcmp(eta->Data(), input->internal_position->Data(),
                        2 * sizeof(double)) == 0;
  }
  StatusCode Save(const BeamKinematicInput* input, double point, Diagnostic* d) {
    const Matrix* values[] = {input->positions, input->tangents, input->normals,
                              input->relative_positions};
    for (int i = 0; i < 4; ++i) {
      if (!values[i]) continue;
      const auto status = key[i]->CopyFrom(values[i], d);
      if (status != StatusCode::kOk) return status;
    }
    const auto status = eta->CopyFrom(input->internal_position, d);
    if (status != StatusCode::kOk) return status;
    length = input->length;
    xi = point;
    has_relative = input->relative_positions != nullptr;
    return StatusCode::kOk;
  }
};

BeamKinematicWorkspace::BeamKinematicWorkspace() = default;
BeamKinematicWorkspace::~BeamKinematicWorkspace() { delete storage_; }

StatusCode BeamKinematicWorkspace::Create(BeamKinematicWorkspace** out,
                                          Diagnostic* d) {
  if (!out) return Invalid(d);
  *out = nullptr;
  std::unique_ptr<BeamKinematicWorkspace> work(new (std::nothrow)
                                                   BeamKinematicWorkspace);
  if (!work) return Allocation(d);
  work->storage_ = new (std::nothrow) Storage;
  if (!work->storage_) return Allocation(d);
  work->storage_->arena = new (std::nothrow) Arena;
  if (!work->storage_->arena) return Allocation(d);
  StatusCode status = work->storage_->arena->Initialize(d);
  if (status != StatusCode::kOk) return status;
  for (auto*& matrix : work->storage_->key) {
    status = Matrix::Create(2, 3, &matrix, d);
    if (status != StatusCode::kOk) return status;
  }
  status = Vector::Create(2, &work->storage_->eta, d);
  if (status != StatusCode::kOk) return status;
  *out = work.release();
  return StatusCode::kOk;
}

StatusCode BeamKinematicWorkspace::Evaluate(const BeamKinematicInput* input,
                                            double xi,
                                            TrialEvaluationFailure* failure,
                                            Diagnostic* d,
                                            BeamKinematicDerivative derivative) {
  storage_->valid = false;
  const bool previous_base = storage_->base_valid;
  storage_->base_valid = false;
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (!failure || !Valid(input, xi) ||
      (derivative != BeamKinematicDerivative::kValue &&
       derivative != BeamKinematicDerivative::kFirst &&
       derivative != BeamKinematicDerivative::kFull)) return Invalid(d);
  const int blocks = derivative == BeamKinematicDerivative::kValue ? 1
                   : derivative == BeamKinematicDerivative::kFirst ? 2 : 4;
  const double length = input->length;
  const double la = length + input->internal_position->Data()[0];
  const double lb = length + input->internal_position->Data()[1];
  if (!(la > 0) || !(lb > 0) || !std::isfinite(la) || !std::isfinite(lb))
    return Geometry(TrialEvaluationFailure::kGeometry,
                    "BeamKinematics nonpositive endpoint stretch", failure, d);

  // Match only after the same full input/geometry preflight as a fresh call.
  storage_->base_valid = previous_base;
  const bool reuse = storage_->Matches(input, xi) &&
                     storage_->base_blocks >= (blocks == 1 ? 1 : 2);
  storage_->base_valid = false;
  DifferentialOps ops(storage_->arena, blocks, reuse, d);
  ops.Direction(input->direction);
  const Differential* xa = ops.Seed(input->positions, 0, false);
  const Differential* xb = ops.Seed(input->positions, 1, false);
  const Differential* relative_xa = input->relative_positions
      ? ops.Seed(input->relative_positions, 0, false) : xa;
  const Differential* relative_xb = input->relative_positions
      ? ops.Seed(input->relative_positions, 1, false) : xb;
  const Differential* ta = ops.Seed(input->tangents, 0, true);
  const Differential* tb = ops.Seed(input->tangents, 1, true);
  const Differential* na = ops.Seed(input->normals, 0, true);
  const Differential* nb = ops.Seed(input->normals, 1, true);
  const Differential* length_a = ops.LengthSeed(la, 0);
  const Differential* length_b = ops.LengthSeed(lb, 1);
  const StatusCode status =
      beam3d_graph::Evaluate(&ops, xa, xb, relative_xa, relative_xb,
                             ta, tb, na, nb, length_a, length_b,
                             length, xi, storage_->fields, failure, d);
  if (status != StatusCode::kOk) return status;
  if (!reuse) {
    const auto saved = storage_->Save(input, xi, d);
    if (saved != StatusCode::kOk) return saved;
    storage_->base_blocks = blocks == 1 ? 1 : 2;
  }
  storage_->base_valid = true;
  storage_->blocks = blocks;
  storage_->valid = true;
  return StatusCode::kOk;
}

void BeamKinematicWorkspace::Invalidate() {
  storage_->valid = false;
  storage_->base_valid = false;
}

const Matrix* BeamKinematicWorkspace::Value(BeamKinematicField field) const {
  const int index = static_cast<int>(field);
  return storage_->valid && index >= 0 && index < kFieldCount
             ? storage_->fields[index]->block[kValue]
             : nullptr;
}
const Matrix* BeamKinematicWorkspace::Variation(
    BeamKinematicField field) const {
  const int index = static_cast<int>(field);
  return storage_->valid && storage_->blocks >= 2 && index >= 0 &&
                 index < kFieldCount
             ? storage_->fields[index]->block[kVirtual]
             : nullptr;
}
const Matrix* BeamKinematicWorkspace::Directional(
    BeamKinematicField field) const {
  const int index = static_cast<int>(field);
  return storage_->valid && storage_->blocks == 4 && index >= 0 &&
                 index < kFieldCount
             ? storage_->fields[index]->block[kDirection]
             : nullptr;
}
const Matrix* BeamKinematicWorkspace::Mixed(BeamKinematicField field) const {
  const int index = static_cast<int>(field);
  return storage_->valid && storage_->blocks == 4 && index >= 0 &&
                 index < kFieldCount
             ? storage_->fields[index]->block[kMixed]
             : nullptr;
}
}  // namespace fesa
