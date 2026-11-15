#include "adapters/abaqus/initial_velocity_reader.h"

#include <map>
#include <memory>
#include <new>
#include <string>
#include <vector>

#include "adapters/abaqus/input_syntax.h"
#include "fesa/model/index_set.h"

namespace fesa {
namespace abaqus {
StatusCode ResolveInitialVelocity(const InputDeck* deck, const Domain* domain,
                                  const DofManager* dofs, Vector** out,
                                  Diagnostic* d) {
  if (out) *out = nullptr;
  if (!out || !deck || !domain || !dofs)
    return Error(nullptr, 0, 0, "Initial velocity arguments", d);
  if (domain->Nodes()->size() != static_cast<size_t>(dofs->NodeCount()))
    return Error(deck->SourceFile(), 0, 0, "Initial velocity node mapping", d);
  for (const auto& entry : *domain->Nodes()) {
    const int count = dofs->NodeDofCount(entry.first);
    if (count != 3 && count != 6)
      return Error(deck->SourceFile(), 0, 0, "Initial velocity requires 3N/6N",
                   d, StatusCode::kUnsupported);
  }
  Vector* raw = nullptr;
  auto status = Vector::Create(dofs->GlobalCount(), &raw, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> velocity(raw);
  try {
    std::map<int, int> first_lines;
    for (const auto* block : *deck->Blocks()) {
      if (*block->keyword != "INITIALCONDITIONS") continue;
      for (const auto* row : *block->rows) {
        std::vector<std::string> fields(0);
        Fields(row->text->c_str(), &fields);
        int component = 0;
        double value = 0;
        if (fields.size() != 3 ||
            !Integer(fields[1].c_str(), 1, 6, &component) ||
            !Real(fields[2].c_str(), &value))
          return Error(deck->SourceFile(), row->line, 0,
                       "Invalid initial velocity row", d);
        if (d) {
          d->node_id = 0;
          d->condition_dof = component;
        }
        std::string scope = "ASSEMBLY.";
        scope.append(domain->AssemblyName());
        std::string target = fields[0];
        const auto dot = target.find('.');
        if (dot != std::string::npos) {
          if (!Equal(target.substr(0, dot).c_str(), domain->InstanceName()))
            return Error(deck->SourceFile(), row->line, 0,
                         "Undefined initial velocity Instance", d);
          scope.assign("PART.");
          scope.append(domain->PartName());
          target = target.substr(dot + 1);
        }
        int node_id = 0;
        const int* ids = &node_id;
        int count = 1;
        if (!Integer(target.c_str(), 1, 999999999, &node_id)) {
          const IndexSet* set = nullptr;
          status = domain->FindIndexSetByName(scope.c_str(), target.c_str(),
                                              false, &set, d);
          if (status != StatusCode::kOk)
            return Error(deck->SourceFile(), row->line, 0,
                         "Undefined initial velocity node set", d,
                         status == StatusCode::kAllocationFailure
                             ? status
                             : StatusCode::kInputError);
          ids = set->Ids();
          count = set->Size();
        }
        for (int i = 0; i < count; ++i) {
          if (dofs->NodeDofCount(ids[i]) > 0 &&
              component > dofs->NodeDofCount(ids[i])) {
            status = Error(deck->SourceFile(), row->line, 0,
                           "Initial velocity component exceeds node DOFs", d,
                           StatusCode::kUnsupported);
            if (d) {
              d->node_id = ids[i];
              d->condition_dof = component;
            }
            return status;
          }
          const int global = dofs->GlobalDof(ids[i], component);
          const auto inserted = first_lines.emplace(global, row->line);
          const char* failure = nullptr;
          if (global < 0)
            failure = "Undefined initial velocity node";
          else if (!inserted.second)
            failure = "Duplicate initial velocity DOF";
          else if (dofs->Equation(global) == -1 && value != 0)
            failure = "Nonzero initial velocity on fixed DOF";
          if (failure) {
            status = Error(deck->SourceFile(), row->line, 0, failure, d,
                           StatusCode::kInputError,
                           inserted.second ? 0 : inserted.first->second);
            if (d) {
              d->node_id = ids[i];
              d->condition_dof = component;
            }
            return status;
          }
          status = velocity->Set(global, value, d);
          if (status != StatusCode::kOk) return status;
        }
      }
    }
    *out = velocity.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return Error(deck->SourceFile(), 0, 0, "Initial velocity allocation", d,
                 StatusCode::kAllocationFailure);
  }
}
}  // namespace abaqus
}  // namespace fesa
