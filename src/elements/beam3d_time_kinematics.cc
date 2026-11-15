#include "elements/beam3d_time_kinematics.h"

#include <cassert>
#include <cmath>
#include <memory>
#include <new>

#include "elements/beam3d_kinematic_graph.h"

namespace fesa {
namespace {
constexpr int kOrders = 3;
constexpr int kVariables = 14;
constexpr int kFields = 6;
constexpr int kSlots = 64;

StatusCode Allocation(Diagnostic* d) {
  return ReportFailure(StatusCode::kAllocationFailure,
                       DiagnosticStage::kAllocation,
                       "BeamTimeKinematics workspace allocation", d);
}

// Numerical coefficient/D pairs, not an expression tape. Coefficient2=.5ddot.
struct TimePair {
  Matrix* value[kOrders] = {};
  Matrix* derivative[kOrders] = {};
  TimePair() = default;
  TimePair(const TimePair&) = delete;
  TimePair& operator=(const TimePair&) = delete;
  TimePair(TimePair&&) = delete;
  TimePair& operator=(TimePair&&) = delete;
  ~TimePair() {
    for (int i = 0; i < kOrders; ++i) {
      delete value[i];
      delete derivative[i];
    }
  }
  StatusCode Initialize(int rows, Diagnostic* d) {
    for (int i = 0; i < kOrders; ++i) {
      StatusCode status = Matrix::Create(rows, 1, &value[i], d);
      if (status != StatusCode::kOk) return status;
      status = Matrix::Create(rows, kVariables, &derivative[i], d);
      if (status != StatusCode::kOk) return status;
    }
    return StatusCode::kOk;
  }
};

// Owns only fixed-size math buffers. The shared graph fits64 scalar/vector
// slots.
struct TimeArena {
  TimePair* scalar[kSlots] = {};
  TimePair* vector[kSlots] = {};
  TimePair* zero_scalar = nullptr;
  TimePair* zero_vector = nullptr;
  TimePair* strain = nullptr;
  Matrix* selector[4] = {};
  Matrix* node_selector[2] = {};
  Matrix* motion_row[2] = {};  // V and A, 1x14.
  Matrix* angular[2] = {};     // omega and alpha, 3x1.
  Matrix* skew[4] = {};
  Matrix* angle_derivative[4] = {};  // DR, DN, DP, DQ, each1x14.
  int scalar_count = 0;
  int vector_count = 0;
  TimeArena() = default;
  TimeArena(const TimeArena&) = delete;
  TimeArena& operator=(const TimeArena&) = delete;
  TimeArena(TimeArena&&) = delete;
  TimeArena& operator=(TimeArena&&) = delete;
  ~TimeArena() {
    for (int i = 0; i < kSlots; ++i) {
      delete scalar[i];
      delete vector[i];
    }
    delete zero_scalar;
    delete zero_vector;
    delete strain;
    for (int i = 0; i < 4; ++i) {
      delete selector[i];
      delete skew[i];
      delete angle_derivative[i];
    }
    for (int i = 0; i < 2; ++i) {
      delete node_selector[i];
      delete motion_row[i];
      delete angular[i];
    }
  }
  StatusCode CreatePair(int rows, TimePair** out, Diagnostic* d) {
    *out = new (std::nothrow) TimePair;
    if (!*out) return Allocation(d);
    return (*out)->Initialize(rows, d);
  }
  StatusCode Initialize(Diagnostic* d) {
    StatusCode status = StatusCode::kOk;
    for (int i = 0; i < kSlots && status == StatusCode::kOk; ++i) {
      status = CreatePair(1, &scalar[i], d);
      if (status == StatusCode::kOk) status = CreatePair(3, &vector[i], d);
    }
    if (status == StatusCode::kOk) status = CreatePair(1, &zero_scalar, d);
    if (status == StatusCode::kOk) status = CreatePair(3, &zero_vector, d);
    if (status == StatusCode::kOk) status = CreatePair(4, &strain, d);
    for (int i = 0; i < 4 && status == StatusCode::kOk; ++i) {
      status = Matrix::Create(3, kVariables, &selector[i], d);
      if (status == StatusCode::kOk) status = Matrix::Create(3, 3, &skew[i], d);
      if (status == StatusCode::kOk)
        status = Matrix::Create(1, kVariables, &angle_derivative[i], d);
      for (int j = 0; j < 3 && status == StatusCode::kOk; ++j)
        status = selector[i]->Set(j, 3 * i + j, 1, d);
    }
    for (int i = 0; i < 2 && status == StatusCode::kOk; ++i) {
      status = Matrix::Create(2, 1, &node_selector[i], d);
      if (status == StatusCode::kOk) status = node_selector[i]->Set(i, 0, 1, d);
      if (status == StatusCode::kOk)
        status = Matrix::Create(1, kVariables, &motion_row[i], d);
      if (status == StatusCode::kOk)
        status = Matrix::Create(3, 1, &angular[i], d);
    }
    return status;
  }
};

// Call-local exact algebra. Sticky status preserves the first backend failure.
class TimePairOps final {
 public:
  TimePairOps(TimeArena* arena, BeamTimeDerivative mode, Diagnostic* d)
      : arena_(arena), mode_(mode), d_(d) {
    arena_->scalar_count = 0;
    arena_->vector_count = 0;
  }
  TimePairOps(const TimePairOps&) = delete;
  TimePairOps& operator=(const TimePairOps&) = delete;
  TimePairOps(TimePairOps&&) = delete;
  TimePairOps& operator=(TimePairOps&&) = delete;
  StatusCode status() const { return status_; }
  double Value(const TimePair* a) const { return a->value[0]->Data()[0]; }
  TimePair* Strain() const { return arena_->strain; }
  void Copy(Matrix* out, const Matrix* a) {
    if (status_ == StatusCode::kOk) status_ = out->CopyFrom(a, d_);
  }
  void ScaleMatrix(Matrix* out, double factor) {
    if (status_ == StatusCode::kOk) status_ = out->Scale(factor, d_);
  }
  void AddMatrix(Matrix* out, double factor, const Matrix* a) {
    if (status_ == StatusCode::kOk) status_ = out->AddScaled(factor, a, d_);
  }
  void Set(Matrix* out, int row, int column, double value) {
    if (status_ == StatusCode::kOk) status_ = out->Set(row, column, value, d_);
  }
  void Multiply(Matrix* out, const Matrix* a, bool ta, const Matrix* b, bool tb,
                double alpha = 1, double beta = 0) {
    if (status_ == StatusCode::kOk)
      status_ = Matrix::MultiplyAdd(a, ta, b, tb, alpha, beta, out, d_);
  }
  void Motion(const Vector* velocity, const Vector* acceleration) {
    if (status_ == StatusCode::kOk)
      status_ = arena_->motion_row[0]->SetRow(0, velocity, d_);
    if (status_ == StatusCode::kOk)
      status_ = arena_->motion_row[1]->SetRow(0, acceleration, d_);
  }
  TimePair* Next(int rows) {
    assert(rows == 1 || rows == 3);
    if (rows == 1) {
      assert(arena_->scalar_count < kSlots);
      return arena_->scalar[arena_->scalar_count++];
    }
    assert(arena_->vector_count < kSlots);
    return arena_->vector[arena_->vector_count++];
  }
  void Clear(TimePair* out) {
    const TimePair* zero =
        out->value[0]->Rows() == 1 ? arena_->zero_scalar : arena_->zero_vector;
    for (int order = 0; order < kOrders; ++order) {
      Copy(out->value[order], zero->value[order]);
      if (Differentiate())
        Copy(out->derivative[order], zero->derivative[order]);
    }
  }
  void Skew(const Matrix* a, Matrix* out) {
    if (status_ != StatusCode::kOk) return;
    const double* v = a->Data();
    // Fixed scalar cross-product coefficient assembly; diagonal stays zero.
    Set(out, 0, 1, -v[2]);
    Set(out, 0, 2, v[1]);
    Set(out, 1, 0, v[2]);
    Set(out, 1, 2, -v[0]);
    Set(out, 2, 0, -v[1]);
    Set(out, 2, 1, v[0]);
  }
  TimePair* Seed(const Matrix* rows, int node, bool angular) {
    TimePair* out = Next(3);
    Clear(out);
    Multiply(out->value[0], rows, true, arena_->node_selector[node], false);
    const Matrix* selector = arena_->selector[2 * node + (angular ? 1 : 0)];
    if (!angular) {
      Multiply(out->value[1], selector, false, arena_->motion_row[0], true);
      Multiply(out->value[2], selector, false, arena_->motion_row[1], true, .5);
      if (mode_ == BeamTimeDerivative::kDisplacement)
        Copy(out->derivative[0], selector);
      if (mode_ == BeamTimeDerivative::kVelocity)
        Copy(out->derivative[1], selector);
      if (mode_ == BeamTimeDerivative::kAcceleration) {
        Copy(out->derivative[2], selector);
        ScaleMatrix(out->derivative[2], .5);
      }
      return out;
    }
    // Inertia contract endpoint seed: a1=omega x a0,
    // a2=(alpha x a0+omega x a1)/2. q/V/A remain independent.
    for (int i = 0; i < 2; ++i)
      Multiply(arena_->angular[i], selector, false, arena_->motion_row[i],
               true);
    Skew(arena_->angular[0], arena_->skew[0]);
    Skew(arena_->angular[1], arena_->skew[1]);
    Multiply(out->value[1], arena_->skew[0], false, out->value[0], false);
    Multiply(out->value[2], arena_->skew[1], false, out->value[0], false, .5);
    Multiply(out->value[2], arena_->skew[0], false, out->value[1], false, .5,
             1);
    if (!Differentiate()) return out;
    Skew(out->value[0], arena_->skew[2]);
    Skew(out->value[1], arena_->skew[3]);
    if (mode_ == BeamTimeDerivative::kDisplacement)
      Multiply(out->derivative[0], arena_->skew[2], false, selector, false, -1);
    // Da1=-[a0]x Domega+[omega]x Da0.
    Multiply(out->derivative[1], arena_->skew[0], false, out->derivative[0],
             false);
    if (mode_ == BeamTimeDerivative::kVelocity)
      Multiply(out->derivative[1], arena_->skew[2], false, selector, false, -1,
               1);
    // Da2=.5(-[a0]x Dalpha+[alpha]x Da0-[a1]x Domega+[omega]x Da1).
    Multiply(out->derivative[2], arena_->skew[1], false, out->derivative[0],
             false, .5);
    Multiply(out->derivative[2], arena_->skew[0], false, out->derivative[1],
             false, .5, 1);
    if (mode_ == BeamTimeDerivative::kVelocity)
      Multiply(out->derivative[2], arena_->skew[3], false, selector, false, -.5,
               1);
    if (mode_ == BeamTimeDerivative::kAcceleration)
      Multiply(out->derivative[2], arena_->skew[2], false, selector, false, -.5,
               1);
    return out;
  }
  TimePair* LengthSeed(double value, int node) {
    TimePair* out = Next(1);
    Clear(out);
    Set(out->value[0], 0, 0, value);
    Set(out->value[1], 0, 0, arena_->motion_row[0]->Data()[12 + node]);
    Set(out->value[2], 0, 0, .5 * arena_->motion_row[1]->Data()[12 + node]);
    const int order = mode_ == BeamTimeDerivative::kDisplacement ? 0
                      : mode_ == BeamTimeDerivative::kVelocity   ? 1
                                                                 : 2;
    if (Differentiate())
      Set(out->derivative[order], 0, 12 + node, order == 2 ? .5 : 1);
    return out;
  }
  TimePair* Scale(const TimePair* a, double factor) {
    TimePair* out = Next(a->value[0]->Rows());
    for (int i = 0; i < kOrders; ++i) {
      Copy(out->value[i], a->value[i]);
      ScaleMatrix(out->value[i], factor);
      if (Differentiate()) {
        Copy(out->derivative[i], a->derivative[i]);
        ScaleMatrix(out->derivative[i], factor);
      }
    }
    return out;
  }
  TimePair* Add(const TimePair* a, const TimePair* b, double factor = 1) {
    TimePair* out = Next(a->value[0]->Rows());
    for (int i = 0; i < kOrders; ++i) {
      Copy(out->value[i], a->value[i]);
      AddMatrix(out->value[i], factor, b->value[i]);
      if (Differentiate()) {
        Copy(out->derivative[i], a->derivative[i]);
        AddMatrix(out->derivative[i], factor, b->derivative[i]);
      }
    }
    return out;
  }
  TimePair* Offset(const TimePair* a, double constant) {
    TimePair* out = Scale(a, 1);
    if (status_ == StatusCode::kOk)
      Set(out->value[0], 0, 0, Value(a) + constant);
    return out;
  }
  TimePair* Combination(const TimePair* a, const TimePair* b, const TimePair* c,
                        const TimePair* e, double ca, double cb, double cc,
                        double ce) {
    TimePair* out = Scale(a, ca);
    for (int i = 0; i < kOrders; ++i) {
      AddMatrix(out->value[i], cb, b->value[i]);
      AddMatrix(out->value[i], cc, c->value[i]);
      AddMatrix(out->value[i], ce, e->value[i]);
      if (Differentiate()) {
        AddMatrix(out->derivative[i], cb, b->derivative[i]);
        AddMatrix(out->derivative[i], cc, c->derivative[i]);
        AddMatrix(out->derivative[i], ce, e->derivative[i]);
      }
    }
    return out;
  }
  TimePair* Product(const TimePair* a, const TimePair* scalar) {
    TimePair* out = Next(a->value[0]->Rows());
    if (status_ != StatusCode::kOk) return out;
    // Time convolution, each term differentiated by the Pair product rule.
    for (int order = 0; order < kOrders; ++order) {
      Copy(out->value[order], a->value[0]);
      ScaleMatrix(out->value[order], scalar->value[order]->Data()[0]);
      if (Differentiate()) {
        Copy(out->derivative[order], a->derivative[0]);
        ScaleMatrix(out->derivative[order], scalar->value[order]->Data()[0]);
        Multiply(out->derivative[order], a->value[0], false,
                 scalar->derivative[order], false, 1, 1);
      }
      for (int j = 1; j <= order; ++j) {
        const double factor = scalar->value[order - j]->Data()[0];
        AddMatrix(out->value[order], factor, a->value[j]);
        if (Differentiate()) {
          AddMatrix(out->derivative[order], factor, a->derivative[j]);
          Multiply(out->derivative[order], a->value[j], false,
                   scalar->derivative[order - j], false, 1, 1);
        }
      }
    }
    return out;
  }
  TimePair* Dot(const TimePair* a, const TimePair* b) {
    TimePair* out = Next(1);
    for (int order = 0; order < kOrders; ++order)
      for (int j = 0; j <= order; ++j) {
        Multiply(out->value[order], a->value[j], true, b->value[order - j],
                 false, 1, j == 0 ? 0 : 1);
        if (Differentiate()) {
          Multiply(out->derivative[order], a->value[j], true,
                   b->derivative[order - j], false, 1, j == 0 ? 0 : 1);
          Multiply(out->derivative[order], b->value[order - j], true,
                   a->derivative[j], false, 1, 1);
        }
      }
    return out;
  }
  TimePair* Cross(const TimePair* a, const TimePair* b) {
    TimePair* out = Next(3);
    for (int order = 0; order < kOrders; ++order)
      for (int j = 0; j <= order; ++j) {
        Skew(a->value[j], arena_->skew[0]);
        Multiply(out->value[order], arena_->skew[0], false, b->value[order - j],
                 false, 1, j == 0 ? 0 : 1);
        if (Differentiate()) {
          Skew(b->value[order - j], arena_->skew[1]);
          Multiply(out->derivative[order], arena_->skew[0], false,
                   b->derivative[order - j], false, 1, j == 0 ? 0 : 1);
          Multiply(out->derivative[order], arena_->skew[1], false,
                   a->derivative[j], false, -1, 1);
        }
      }
    return out;
  }
  TimePair* Chain(const TimePair* a, double value, double first, double second,
                  double third) {
    TimePair* out = Next(1);
    if (status_ != StatusCode::kOk) return out;
    const double a1 = a->value[1]->Data()[0], a2 = a->value[2]->Data()[0];
    // f2=f' a2+.5 f'' a1². Its D requires f''' as well as f''.
    Set(out->value[0], 0, 0, value);
    Set(out->value[1], 0, 0, first * a1);
    Set(out->value[2], 0, 0, first * a2 + .5 * second * a1 * a1);
    if (Differentiate()) {
      Copy(out->derivative[0], a->derivative[0]);
      ScaleMatrix(out->derivative[0], first);
      Copy(out->derivative[1], a->derivative[1]);
      ScaleMatrix(out->derivative[1], first);
      AddMatrix(out->derivative[1], second * a1, a->derivative[0]);
      Copy(out->derivative[2], a->derivative[2]);
      ScaleMatrix(out->derivative[2], first);
      AddMatrix(out->derivative[2], second * a1, a->derivative[1]);
      AddMatrix(out->derivative[2], second * a2 + .5 * third * a1 * a1,
                a->derivative[0]);
    }
    return out;
  }
  TimePair* Reciprocal(const TimePair* a) {
    const double v = status_ == StatusCode::kOk ? 1 / Value(a) : 0;
    return Chain(a, v, -v * v, 2 * v * v * v, -6 * v * v * v * v);
  }
  TimePair* SquareRoot(const TimePair* a) {
    const double r = status_ == StatusCode::kOk ? std::sqrt(Value(a)) : 1;
    return Chain(a, r, .5 / r, -.25 / (r * r * r), .375 / (r * r * r * r * r));
  }
  TimePair* Sin(const TimePair* a) {
    const double v = status_ == StatusCode::kOk ? Value(a) : 0;
    return Chain(a, std::sin(v), std::cos(v), -std::sin(v), -std::cos(v));
  }
  TimePair* Cos(const TimePair* a) {
    const double v = status_ == StatusCode::kOk ? Value(a) : 0;
    return Chain(a, std::cos(v), -std::sin(v), -std::cos(v), std::sin(v));
  }
  TimePair* Angle(const TimePair* y, const TimePair* x) {
    TimePair* out = Next(1);
    if (status_ != StatusCode::kOk) return out;
    const double x0 = Value(x), y0 = Value(y);
    const double x1 = x->value[1]->Data()[0], y1 = y->value[1]->Data()[0];
    const double x2 = x->value[2]->Data()[0], y2 = y->value[2]->Data()[0];
    const double r = x0 * x0 + y0 * y0;
    const double n = x0 * y1 - y0 * x1;
    const double p = x0 * y2 - y0 * x2;
    const double q = x0 * x1 + y0 * y1;
    // Inertia contract exact atan2 coefficients; D below is its Pair algebra.
    Set(out->value[0], 0, 0, std::atan2(y0, x0));
    Set(out->value[1], 0, 0, n / r);
    Set(out->value[2], 0, 0, p / r - n * q / (r * r));
    if (!Differentiate()) return out;
    Matrix* dr = arena_->angle_derivative[0];
    Matrix* dn = arena_->angle_derivative[1];
    Matrix* dp = arena_->angle_derivative[2];
    Matrix* dq = arena_->angle_derivative[3];
    Copy(dr, x->derivative[0]);
    ScaleMatrix(dr, 2 * x0);
    AddMatrix(dr, 2 * y0, y->derivative[0]);
    Copy(dn, x->derivative[0]);
    ScaleMatrix(dn, y1);
    AddMatrix(dn, x0, y->derivative[1]);
    AddMatrix(dn, -x1, y->derivative[0]);
    AddMatrix(dn, -y0, x->derivative[1]);
    Copy(dp, x->derivative[0]);
    ScaleMatrix(dp, y2);
    AddMatrix(dp, x0, y->derivative[2]);
    AddMatrix(dp, -x2, y->derivative[0]);
    AddMatrix(dp, -y0, x->derivative[2]);
    Copy(dq, x->derivative[0]);
    ScaleMatrix(dq, x1);
    AddMatrix(dq, x0, x->derivative[1]);
    AddMatrix(dq, y1, y->derivative[0]);
    AddMatrix(dq, y0, y->derivative[1]);
    Copy(out->derivative[0], x->derivative[0]);
    ScaleMatrix(out->derivative[0], -y0 / r);
    AddMatrix(out->derivative[0], x0 / r, y->derivative[0]);
    Copy(out->derivative[1], dn);
    ScaleMatrix(out->derivative[1], 1 / r);
    AddMatrix(out->derivative[1], -n / (r * r), dr);
    Copy(out->derivative[2], dp);
    ScaleMatrix(out->derivative[2], 1 / r);
    AddMatrix(out->derivative[2], -p / (r * r), dr);
    AddMatrix(out->derivative[2], -q / (r * r), dn);
    AddMatrix(out->derivative[2], -n / (r * r), dq);
    AddMatrix(out->derivative[2], 2 * n * q / (r * r * r), dr);
    return out;
  }
  void Pack(const TimePair* a, int row) {
    for (int order = 0; order < kOrders; ++order) {
      if (status_ == StatusCode::kOk)
        status_ = arena_->strain->value[order]->CopyBlockFrom(
            a->value[order], 0, 0, 1, 1, row, 0, d_);
      if (Differentiate() && status_ == StatusCode::kOk)
        status_ = arena_->strain->derivative[order]->CopyBlockFrom(
            a->derivative[order], 0, 0, 1, kVariables, row, 0, d_);
    }
  }

 private:
  bool Differentiate() const { return mode_ != BeamTimeDerivative::kNone; }
  TimeArena* arena_;  // Borrowed only during Evaluate.
  BeamTimeDerivative mode_;
  Diagnostic* d_;
  StatusCode status_ = StatusCode::kOk;
};

bool MotionValid(const Vector* v) {
  return v && v->Size() == kVariables &&
         beam3d_graph::Finite(v->Data(), kVariables);
}
bool ModeValid(BeamTimeDerivative mode) {
  return mode == BeamTimeDerivative::kNone ||
         mode == BeamTimeDerivative::kDisplacement ||
         mode == BeamTimeDerivative::kVelocity ||
         mode == BeamTimeDerivative::kAcceleration;
}
}  // namespace

struct BeamTimeKinematicWorkspace::Storage {
  TimeArena* arena = nullptr;
  TimePair* fields[kFields] = {};
  bool valid = false;
  bool differentiated = false;
  Storage() = default;
  Storage(const Storage&) = delete;
  Storage& operator=(const Storage&) = delete;
  Storage(Storage&&) = delete;
  Storage& operator=(Storage&&) = delete;
  ~Storage() { delete arena; }
};

BeamTimeKinematicWorkspace::BeamTimeKinematicWorkspace() = default;
BeamTimeKinematicWorkspace::~BeamTimeKinematicWorkspace() { delete storage_; }

StatusCode BeamTimeKinematicWorkspace::Create(BeamTimeKinematicWorkspace** out,
                                              Diagnostic* d) {
  if (!out) return beam3d_graph::Invalid(d);
  *out = nullptr;
  std::unique_ptr<BeamTimeKinematicWorkspace> work(
      new (std::nothrow) BeamTimeKinematicWorkspace);
  if (!work) return Allocation(d);
  work->storage_ = new (std::nothrow) Storage;
  if (!work->storage_) return Allocation(d);
  work->storage_->arena = new (std::nothrow) TimeArena;
  if (!work->storage_->arena) return Allocation(d);
  const StatusCode status = work->storage_->arena->Initialize(d);
  if (status != StatusCode::kOk) return status;
  *out = work.release();
  return StatusCode::kOk;
}

StatusCode BeamTimeKinematicWorkspace::Evaluate(
    const BeamKinematicInput* pose, const Vector* velocity14,
    const Vector* acceleration14, BeamTimeDerivative derivative, double xi,
    TrialEvaluationFailure* failure, Diagnostic* d) {
  storage_->valid = false;
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (!failure || !beam3d_graph::Valid(pose, xi) || !MotionValid(velocity14) ||
      !MotionValid(acceleration14) || !ModeValid(derivative))
    return beam3d_graph::Invalid(d);
  const double length = pose->length;
  const double la = length + pose->internal_position->Data()[0];
  const double lb = length + pose->internal_position->Data()[1];
  if (!(la > 0) || !(lb > 0) || !std::isfinite(la) || !std::isfinite(lb))
    return beam3d_graph::Geometry(TrialEvaluationFailure::kGeometry,
                                  "BeamKinematics nonpositive endpoint stretch",
                                  failure, d);
  TimePairOps ops(storage_->arena, derivative, d);
  ops.Motion(velocity14, acceleration14);
  const TimePair* xa = ops.Seed(pose->positions, 0, false);
  const TimePair* xb = ops.Seed(pose->positions, 1, false);
  const TimePair* relative_xa = pose->relative_positions
      ? ops.Seed(pose->relative_positions, 0, false) : xa;
  const TimePair* relative_xb = pose->relative_positions
      ? ops.Seed(pose->relative_positions, 1, false) : xb;
  const TimePair* ta = ops.Seed(pose->tangents, 0, true);
  const TimePair* tb = ops.Seed(pose->tangents, 1, true);
  const TimePair* na = ops.Seed(pose->normals, 0, true);
  const TimePair* nb = ops.Seed(pose->normals, 1, true);
  const TimePair* length_a = ops.LengthSeed(la, 0);
  const TimePair* length_b = ops.LengthSeed(lb, 1);
  const StatusCode status =
      beam3d_graph::Evaluate(&ops, xa, xb, relative_xa, relative_xb,
                             ta, tb, na, nb, length_a, length_b,
                             length, xi, storage_->fields, failure, d);
  if (status != StatusCode::kOk) return status;
  storage_->valid = true;
  storage_->differentiated = derivative != BeamTimeDerivative::kNone;
  return StatusCode::kOk;
}

const Matrix* BeamTimeKinematicWorkspace::Coefficient(BeamKinematicField field,
                                                      int order) const {
  const int index = static_cast<int>(field);
  return storage_->valid && index >= 0 && index < kFields && order >= 0 &&
                 order < kOrders
             ? storage_->fields[index]->value[order]
             : nullptr;
}
const Matrix* BeamTimeKinematicWorkspace::Derivative(BeamKinematicField field,
                                                     int order) const {
  const int index = static_cast<int>(field);
  return storage_->valid && storage_->differentiated && index >= 0 &&
                 index < kFields && order >= 0 && order < kOrders
             ? storage_->fields[index]->derivative[order]
             : nullptr;
}
}  // namespace fesa
