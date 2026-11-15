#include "analysis/static_assembler.h"

#include <algorithm>
#include <limits>
#include <memory>
#include <new>

#include "adapters/mkl/sparse_ops.h"
#include "adapters/openmp/thread_runtime.h"
#include "analysis/numerical_support.h"
#include "model/model_storage.h"
namespace fesa {
namespace {
using numerical_support::Buffers;
StatusCode Flush(Buffers* buffers, int matrix_count, int vector_count,
                 SparseMatrix** accumulated, Vector* rhs, Diagnostic* d) {
  if (vector_count > 0) {
    const StatusCode status = rhs->ScatterAdd(
        buffers->vector_values, buffers->vector_indices, vector_count, d);
    if (status != StatusCode::kOk) return status;
  }
  if (matrix_count == 0) return StatusCode::kOk;
  SparseMatrix* chunk = nullptr;
  StatusCode status =
      SparseMatrix::Create(rhs->Size(), matrix_count, buffers->rows,
                           buffers->columns, buffers->matrix_values, &chunk, d);
  std::unique_ptr<SparseMatrix> chunk_owner(chunk);
  if (status != StatusCode::kOk) return status;
  SparseMatrix* merged = nullptr;
  status = SparseMatrix::Add(*accumulated, chunk, &merged, d);
  if (status != StatusCode::kOk) return status;
  delete *accumulated;
  *accumulated = merged;
  return StatusCode::kOk;
}
struct AssemblyTaskContext {
  const HhtSolidAssemblyInput* dynamic = nullptr;
  const DofManager* dofs = nullptr;
  const Vector* prescribed = nullptr;
  Buffers* buffers = nullptr;
  const BeamLineLoadState* line_loads = nullptr;
};
StatusCode PrepareElement(Element* element, const DofManager* dofs, Buffers* b,
                          int item, Diagnostic* d) {
  b->elements[item] = element;
  int* global = b->global_indices + b->element_dofs * item;
  auto status = numerical_support::ElementIndices(element, dofs, global, d);
  if (status != StatusCode::kOk) return status;
  int matrix = b->matrix_offsets[item], vector = b->vector_offsets[item];
  for (int i = 0; i < element->dof_count(); ++i) {
    const int row = dofs->Equation(global[i]);
    if (row < 0) continue;
    b->vector_indices[vector++] = row;
    for (int j = 0; j < element->dof_count(); ++j) {
      const int column = dofs->Equation(global[j]);
      if (column < row) continue;
      b->rows[matrix] = row;
      b->columns[matrix++] = column;
    }
  }
  b->matrix_offsets[item + 1] = matrix;
  b->vector_offsets[item + 1] = vector;
  return StatusCode::kOk;
}
StatusCode ComputeAssembly(int item, int worker, void* context, Diagnostic* d) {
  auto* task = static_cast<AssemblyTaskContext*>(context);
  Buffers* b = task->buffers;
  auto* local = b->workers[worker];
  const Element* element = b->elements[item];
  d->element_id = element->Id();
  const int size = element->dof_count();
  std::unique_ptr<numerical_support::WorkerBuffers> owned;
  if (local->u->Size() != size) {
    numerical_support::WorkerBuffers* raw = nullptr;
    const auto status =
        numerical_support::WorkerBuffers::Create(local->element, &raw, d, size);
    owned.reset(raw);
    if (status != StatusCode::kOk) return status;
    local = raw;
  }
  StatusCode status = StatusCode::kOk;
  if (!task->dynamic) {
    status = element->ComputeStiffness(local->element, local->stiffness, d);
    if (status != StatusCode::kOk) return status;
  }
  const int* global = b->global_indices + b->element_dofs * item;
  int entries[Buffers::kMaxUpperEntries] = {};
  int local_free[Buffers::kMaxElementDofs] = {};
  int count = 0, free_count = 0;
  for (int i = 0; i < element->dof_count(); ++i) {
    const int row = task->dofs->Equation(global[i]);
    if (row < 0) continue;
    local_free[free_count++] = i;
    for (int j = 0; j < element->dof_count(); ++j)
      if (task->dofs->Equation(global[j]) >= row)
        entries[count++] = size * i + j;
  }
  if (task->dynamic) {
    Vector *raw_committed = nullptr, *raw_predictor = nullptr;
    status = Vector::Create(size, &raw_committed, d);
    std::unique_ptr<Vector> committed(raw_committed);
    if (status != StatusCode::kOk) return status;
    status = Vector::Create(size, &raw_predictor, d);
    std::unique_ptr<Vector> predictor(raw_predictor);
    if (status != StatusCode::kOk) return status;
    status =
        committed->GatherFrom(task->dynamic->committed, global, size, 0, d);
    if (status != StatusCode::kOk) return status;
    status =
        predictor->GatherFrom(task->dynamic->predictor, global, size, 0, d);
    if (status != StatusCode::kOk) return status;
    status = local->u->GatherFrom(task->prescribed, global, size, 0, d);
    if (status != StatusCode::kOk) return status;
    std::unique_ptr<Vector> velocities[2];
    if (HasMaterialRayleighDamping(element)) {
      const Vector* sources[] = {task->dynamic->previous_velocity,
                                 task->dynamic->predicted_velocity};
      for (int i = 0; i < 2; ++i) {
        Vector* raw_v = nullptr;
        status = Vector::Create(size, &raw_v, d);
        velocities[i].reset(raw_v);
        if (status == StatusCode::kOk)
          status = raw_v->GatherFrom(sources[i], global, size, 0, d);
        if (status != StatusCode::kOk) return status;
      }
    }
    status = EvaluateHhtSolidContribution(
        element, local->element, task->dynamic->mode, task->dynamic->alpha,
        task->dynamic->c0, committed.get(), predictor.get(), local->u,
        local->stiffness, local->force, d, velocities[0].get(),
        velocities[1].get(), task->dynamic->c1);
    if (status != StatusCode::kOk) return status;
    status = local->stiffness->GatherTo(entries, count, b->matrix_values,
                                        b->matrix_offsets[item], d);
    if (status != StatusCode::kOk || free_count == 0) return status;
    return b->vector_values->GatherFrom(local->force, local_free, free_count,
                                        b->vector_offsets[item], d);
  }
  status = local->stiffness->GatherTo(entries, count, b->matrix_values,
                                      b->matrix_offsets[item], d);
  if (status != StatusCode::kOk || free_count == 0) return status;
  status = local->u->GatherFrom(task->prescribed, global, size, 0, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyVector(local->stiffness, false, local->u,
                                  local->force, d);
  if (status != StatusCode::kOk) return status;
  status = local->force->Scale(-1, d);  // E1: -Kfc*uc.
  if (status != StatusCode::kOk) return status;
  if (task->line_loads) {
    status = numerical_support::EvaluateLineLoad(task->line_loads, element,
                                                local, d);
    if (status != StatusCode::kOk) return status;
    status = local->force->AddScaled(1, local->line_global, d);
    if (status != StatusCode::kOk) return status;
  }
  return b->vector_values->GatherFrom(local->force, local_free, free_count,
                                      b->vector_offsets[item], d);
}
}  // namespace
StatusCode AccumulateAssemblyCount(int free_count, int element_free_count,
                                   std::size_t* count, Diagnostic* d) {
  const auto maximum =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  if (!count || free_count < 0 || element_free_count < 0 ||
      element_free_count > Buffers::kMaxElementDofs ||
      (free_count == 0 && (*count != 0 || element_free_count != 0)))
    return model_storage::Invalid("Assembly conservative count", d);
  if (free_count == 0) return StatusCode::kOk;
  const auto extra = static_cast<std::size_t>(element_free_count *
                                              (element_free_count + 1) / 2);
  if (*count > maximum || extra > maximum - *count ||
      !sparse_ops::ValidCsrStorageSize(free_count, *count + extra))
    return model_storage::Invalid("Assembly conservative count", d);
  *count += extra;
  return StatusCode::kOk;
}
StatusCode ValidateAssemblyCounts(const AnalysisModel* model,
                                  const DofManager* dofs, Diagnostic* d,
                                  int* maximum_element_dofs, bool* beam_only) {
  if (!model || !dofs) return model_storage::Invalid("Assembly count input", d);
  std::size_t count = static_cast<std::size_t>(dofs->FreeCount());
  auto status = AccumulateAssemblyCount(dofs->FreeCount(), 0, &count, d);
  if (status != StatusCode::kOk) return status;
  int maximum = 0;
  bool beams = true;
  for (const auto& entry : *model->Elements()) {
    int indices[Buffers::kMaxElementDofs] = {};
    status = numerical_support::ElementIndices(entry.second, dofs, indices, d);
    if (status != StatusCode::kOk) return status;
    int free = 0;
    maximum = std::max(maximum, entry.second->dof_count());
    beams = beams && entry.second->family() == ElementFamily::kBeam;
    for (int i = 0; i < entry.second->dof_count(); ++i)
      if (dofs->Equation(indices[i]) >= 0) ++free;
    status = AccumulateAssemblyCount(dofs->FreeCount(), free, &count, d);
    if (status != StatusCode::kOk) return status;
  }
  if (maximum_element_dofs) *maximum_element_dofs = maximum ? maximum : 12;
  if (beam_only) *beam_only = beams;
  return StatusCode::kOk;
}
StatusCode AssembleReducedSystem(const AnalysisModel* model,
                                 const DofManager* dofs,
                                 const EffectiveConditions* conditions,
                                 StaticWorkspace* workspace,
                                 EquationSystem** out, Diagnostic* d,
                                 ElementWorkspace** workers, int worker_count,
                                 const HhtSolidAssemblyInput* dynamic) {
  if (!out) return model_storage::Invalid("AssembleReducedSystem", d);
  *out = nullptr;
  StatusCode status =
      numerical_support::Validate(model, dofs, conditions, workspace, d);
  if (status != StatusCode::kOk) return status;
  int maximum_element_dofs = 12;
  bool beam_only = true;
  status =
      ValidateAssemblyCounts(model, dofs, d, &maximum_element_dofs, &beam_only);
  if (status != StatusCode::kOk) return status;
  if (dynamic) {
    status = ValidateHhtSolidAssembly(model, dofs, conditions, workspace,
                                      dynamic, workers, worker_count, d);
    if (status != StatusCode::kOk) return status;
  }
  if (dofs->FreeCount() == 0) return StatusCode::kOk;
  if (!workspace->ReducedRhs() ||
      workspace->ReducedRhs()->Size() != dofs->FreeCount())
    return model_storage::Invalid("Static assembly reduced RHS", d);
  status = dynamic ? workspace->ExternalForce()->CopyFrom(dynamic->external, d)
                   : numerical_support::ExternalForce(
                         dofs, conditions, workspace->ExternalForce(), d);
  if (status != StatusCode::kOk) return status;
  status = workspace->Increment()->Scale(0, d);
  if (status != StatusCode::kOk) return status;
  for (const auto* boundary : *conditions->Boundaries()) {
    const int global =
        dofs->GlobalDof(boundary->node_id, boundary->definition->Dof());
    status =
        workspace->Increment()->Set(global, boundary->definition->Value(), d);
    if (status != StatusCode::kOk) return status;
  }
  const int size = dofs->FreeCount();
  std::unique_ptr<int[]> indices(new (std::nothrow) int[size]);
  if (!indices) return model_storage::AllocationFailure("Reduced indices", d);
  for (int i = 0; i < size; ++i) indices[i] = dofs->GlobalDofForEquation(i);
  status = workspace->ReducedRhs()->GatherFrom(workspace->ExternalForce(),
                                               indices.get(), size, 0, d);
  if (status != StatusCode::kOk) return status;
  SparseMatrix* matrix = nullptr;
  status = SparseMatrix::Create(size, 0, nullptr, nullptr, nullptr, &matrix, d);
  std::unique_ptr<SparseMatrix> owner(matrix);
  if (status != StatusCode::kOk) return status;
  Buffers* raw = nullptr;
  const int capacity = static_cast<int>(std::min(
      model->Elements()->size(), static_cast<std::size_t>(Buffers::kElements)));
  status = Buffers::Create(true, std::max(1, capacity), &raw, d, workers,
                           worker_count, maximum_element_dofs, beam_only);
  std::unique_ptr<Buffers> buffers(raw);
  if (status != StatusCode::kOk) return status;
  status = numerical_support::PrepareLineLoads(conditions->LineLoadState(),
                                                raw, d);
  if (status != StatusCode::kOk) return status;
  AssemblyTaskContext task;
  task.dynamic = dynamic;
  task.dofs = dofs;
  task.prescribed = workspace->Increment();
  task.buffers = raw;
  task.line_loads = conditions->LineLoadState();
  auto next = model->Elements()->begin();
  while (next != model->Elements()->end()) {
    int count = 0;
    raw->matrix_offsets[0] = raw->vector_offsets[0] = 0;
    while (count < capacity && next != model->Elements()->end()) {
      status = PrepareElement(next->second, dofs, raw, count, d);
      if (status != StatusCode::kOk) {
        if (d) d->element_id = next->first;
        return status;
      }
      ++count;
      ++next;
    }
    status = ExecuteElementTasks(count, worker_count, ComputeAssembly, &task,
                                 raw->diagnostics, nullptr, d);
    if (status != StatusCode::kOk) return status;
    matrix = owner.release();
    status = Flush(raw, raw->matrix_offsets[count], raw->vector_offsets[count],
                   &matrix, workspace->ReducedRhs(), d);
    owner.reset(matrix);
    if (status != StatusCode::kOk) return status;
  }
  Vector* rhs = nullptr;
  status = Vector::Create(size, &rhs, d);
  std::unique_ptr<Vector> rhs_owner(rhs);
  if (status != StatusCode::kOk) return status;
  status = rhs->CopyFrom(workspace->ReducedRhs(), d);
  if (status != StatusCode::kOk) return status;
  status = EquationSystem::Create(matrix, rhs, MatrixType::kSymmetricPositiveDefinite,
                                  out, d);
  if (status == StatusCode::kOk) {
    owner.release();
    rhs_owner.release();
  }
  return status;
}
}  // namespace fesa
