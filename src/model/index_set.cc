#include "fesa/model/index_set.h"

#include <algorithm>
#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {

IndexSet::IndexSet(int id, char* scope, char* name, bool element_set,
                   int count, int* ids)
    : id_(id), scope_(scope), name_(name), element_set_(element_set),
      count_(count), ids_(ids) {}
IndexSet::~IndexSet() {
  delete[] scope_;
  delete[] name_;
  delete[] ids_;
}
StatusCode IndexSet::Create(int id, const char* scope, const char* name,
                            bool element_set, int count, const int* ids,
                            IndexSet** out, Diagnostic* d) {
  if (out == nullptr) return model_storage::Invalid("IndexSet::Create", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || scope == nullptr || name == nullptr
      || name[0] == '\0' || count < 0 || (count > 0 && ids == nullptr)) {
    return model_storage::Invalid("IndexSet::Create", d);
  }
  for (int i = 0; i < count; ++i) {
    if (!model_storage::ValidId(ids[i])) {
      return model_storage::Invalid("IndexSet::Create member ID", d);
    }
  }
  std::unique_ptr<char[]> scope_copy(model_storage::CopyText(scope));
  std::unique_ptr<char[]> name_copy(model_storage::CopyText(name));
  std::unique_ptr<int[]> members(count > 0 ? new (std::nothrow) int[count]
                                         : nullptr);
  if (!scope_copy || !name_copy || (count > 0 && !members)) {
    return model_storage::AllocationFailure("IndexSet::Create", d);
  }
  int unique_count = 0;
  if (count > 0) {
    std::copy_n(ids, count, members.get());
    std::sort(members.get(), members.get() + count);
    unique_count = static_cast<int>(std::unique(members.get(),
                    members.get() + count) - members.get());
  }
  *out = new (std::nothrow) IndexSet(id, scope_copy.get(), name_copy.get(),
                                    element_set, unique_count, members.get());
  if (*out == nullptr) {
    return model_storage::AllocationFailure("IndexSet::Create", d);
  }
  scope_copy.release();
  name_copy.release();
  members.release();
  return StatusCode::kOk;
}
int IndexSet::Id() const { return id_; }
const char* IndexSet::Scope() const { return scope_; }
const char* IndexSet::Name() const { return name_; }
bool IndexSet::IsElementSet() const { return element_set_; }
int IndexSet::Size() const { return count_; }
const int* IndexSet::Ids() const { return ids_; }

}  // namespace fesa
