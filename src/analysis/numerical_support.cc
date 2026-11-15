#include "analysis/numerical_support.h"

#include <memory>
#include <new>

#include "model/model_storage.h"
namespace fesa {
namespace numerical_support {
WorkerBuffers::~WorkerBuffers() {
  delete line_workspace;
  delete line_density;
  delete line_global;
  delete line_local;
  delete stiffness;
  delete u;
  delete force;
}
StatusCode WorkerBuffers::Create(ElementWorkspace* element, WorkerBuffers** out,
                                 Diagnostic* d, int element_dofs) {
  if (!out) return model_storage::Invalid("Worker scratch output", d);
  *out = nullptr;
  if (element_dofs <= 0 || element_dofs > Buffers::kMaxElementDofs)
    return model_storage::Invalid("Worker scratch layout", d);
  std::unique_ptr<WorkerBuffers> result(new (std::nothrow) WorkerBuffers);
  if (!result) return model_storage::AllocationFailure("Worker scratch", d);
  result->element = element;
  auto status =
      Matrix::Create(element_dofs, element_dofs, &result->stiffness, d);
  if (status != StatusCode::kOk) return status;
  for (auto** vector : {&result->u, &result->force}) {
    status = Vector::Create(element_dofs, vector, d);
    if (status != StatusCode::kOk) return status;
  }
  *out = result.release();
  return StatusCode::kOk;
}
Buffers::~Buffers() {
  for (int i = 0; i < worker_count; ++i) {
    if (workers) delete workers[i];
    if (owned_elements_) delete owned_elements_[i];
  }
  delete[] workers;
  delete[] owned_elements_;
  delete[] diagnostics;
  delete[] elements;
  delete[] global_indices;
  delete[] matrix_offsets;
  delete[] vector_offsets;
  delete vector_values;
  delete matrix_values;
  delete[] vector_indices;
  delete[] rows;
  delete[] columns;
}
StatusCode Buffers::Create(bool matrix_storage, int capacity, Buffers** out,
                           Diagnostic* d, ElementWorkspace** elements,
                           int worker_count, int element_dofs,
                           bool beam_workspace) {
  *out = nullptr;
  if (capacity <= 0 || capacity > kElements || worker_count <= 0 ||
      element_dofs <= 0 || element_dofs > kMaxElementDofs)
    return model_storage::Invalid("Assembly capacity", d);
  std::unique_ptr<Buffers> result(new (std::nothrow) Buffers);
  if (!result) return model_storage::AllocationFailure("Assembly buffers", d);
  result->worker_count = worker_count;
  result->element_dofs = element_dofs;
  const int upper = element_dofs * (element_dofs + 1) / 2;
  StatusCode status = StatusCode::kOk;
  if (!elements && beam_workspace) {
    status =
        PrepareWorkerWorkspaces(worker_count, 12, &result->owned_elements_, d);
    if (status != StatusCode::kOk) return status;
    elements = result->owned_elements_;
  }
  result->workers = new (std::nothrow) WorkerBuffers*[worker_count]();
  result->diagnostics = new (std::nothrow) Diagnostic[worker_count];
  result->elements = new (std::nothrow) Element*[capacity];
  result->global_indices = new (std::nothrow) int[element_dofs * capacity];
  result->matrix_offsets = new (std::nothrow) int[capacity + 1];
  result->vector_offsets = new (std::nothrow) int[capacity + 1];
  result->vector_indices = new (std::nothrow) int[element_dofs * capacity];
  if (!result->workers || !result->diagnostics || !result->elements ||
      !result->global_indices || !result->matrix_offsets ||
      !result->vector_offsets || !result->vector_indices)
    return model_storage::AllocationFailure("Assembly indices/workers", d);
  for (int i = 0; i < worker_count; ++i) {
    if (beam_workspace && (!elements || !elements[i]))
      return model_storage::Invalid("Missing worker", d);
    status = WorkerBuffers::Create(elements ? elements[i] : nullptr,
                                   &result->workers[i], d, element_dofs);
    if (status != StatusCode::kOk) return status;
  }
  status = Vector::Create(element_dofs * capacity, &result->vector_values, d);
  if (status != StatusCode::kOk) return status;
  if (matrix_storage) {
    status = Vector::Create(upper * capacity, &result->matrix_values, d);
    if (status != StatusCode::kOk) return status;
    result->rows = new (std::nothrow) int[upper * capacity];
    result->columns = new (std::nothrow) int[upper * capacity];
    if (!result->rows || !result->columns)
      return model_storage::AllocationFailure("COO indices", d);
  }
  *out = result.release();
  return StatusCode::kOk;
}
StatusCode Validate(const AnalysisModel* model, const DofManager* dofs,
                    const EffectiveConditions* conditions,
                    StaticWorkspace* workspace, Diagnostic* d) {
  if (!model || !dofs || !conditions || !workspace ||
      model->Conditions() != conditions || dofs->NodeCount() <= 0 ||
      model->Nodes()->size() != static_cast<std::size_t>(dofs->NodeCount()) ||
      workspace->ExternalForce()->Size() != dofs->GlobalCount())
    return model_storage::Invalid("Static numerical inputs", d);
  return StatusCode::kOk;
}
StatusCode PrepareLineLoads(const BeamLineLoadState* state, Buffers* buffers,
                            Diagnostic* d) {
  if (!state) return StatusCode::kOk;
  for (int i = 0; i < buffers->worker_count; ++i) {
    WorkerBuffers* worker = buffers->workers[i];
    StatusCode status =
        BeamLineLoadWorkspace::Create(&worker->line_workspace, d);
    if (status != StatusCode::kOk) return status;
    status = Vector::Create(3, &worker->line_density, d);
    if (status != StatusCode::kOk) return status;
    status = Vector::Create(12, &worker->line_global, d);
    if (status != StatusCode::kOk) return status;
    status = Vector::Create(12, &worker->line_local, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode EvaluateLineLoad(const BeamLineLoadState* state,
                            const Element* element, WorkerBuffers* worker,
                            Diagnostic* d) {
  const StatusCode status =
      state->CopyDensity(element->Id(), worker->line_density, d);
  return status == StatusCode::kOk
             ? ComputeBeamLineLoad(element, worker->line_density,
                                    worker->line_workspace, worker->line_global,
                                    worker->line_local, d)
             : status;
}
StatusCode ElementIndices(const Element* element, const DofManager* dofs,
                          int* indices, Diagnostic* d) {
  if (!element || !dofs || !indices || element->node_count() <= 0 ||
      (element->dofs_per_node() != 3 && element->dofs_per_node() != 6) ||
      element->node_count() >
          Buffers::kMaxElementDofs / element->dofs_per_node() ||
      element->dof_count() != element->node_count() * element->dofs_per_node())
    return model_storage::Invalid("Element DOF layout", d);
  const int components = element->dofs_per_node();
  for (int local = 0; local < element->node_count(); ++local) {
    const Node* node = element->GetNode(local);
    if (!node || dofs->NodeDofCount(node->Id()) != components)
      return model_storage::Invalid("Element node layout", d);
    for (int component = 0; component < components; ++component) {
      const int index = local * components + component;
      indices[index] = dofs->GlobalDof(node->Id(), component + 1);
      if (indices[index] < 0 || dofs->Equation(indices[index]) == -2)
        return model_storage::Invalid("Element DOF mapping", d);
    }
  }
  return StatusCode::kOk;
}
StatusCode ExternalForce(const DofManager* dofs,
                         const EffectiveConditions* conditions, Vector* out,
                         Diagnostic* d) {
  StatusCode status = out->Scale(0, d);
  if (status != StatusCode::kOk) return status;
  for (const auto* load : *conditions->LoadTotals()) {
    status =
        out->Set(dofs->GlobalDof(load->node_id, load->dof), load->value, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
}  // namespace numerical_support
}  // namespace fesa
