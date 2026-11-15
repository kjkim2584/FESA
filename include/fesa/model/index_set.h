/// @file
/// @brief Linear-static model ownership and immutable input definitions.
#ifndef FESA_MODEL_INDEX_SET_H_
#define FESA_MODEL_INDEX_SET_H_

#include "fesa/base/diagnostic.h"

namespace fesa {

/// Immutable sorted set of original IDs with owned scope/name strings.
class IndexSet {
 public:
  /// Copies strings/IDs and sorts/deduplicates membership. Empty sets allowed.
  /// Caller owns *out on success; failure returns null. IDs must be positive.
  static StatusCode Create(int id, const char* scope, const char* name,
                           bool element_set, int count, const int* ids,
                           IndexSet** out, Diagnostic* d);
  ~IndexSet();
  IndexSet(const IndexSet&) = delete;
  IndexSet& operator=(const IndexSet&) = delete;
  IndexSet(IndexSet&&) = delete;
  IndexSet& operator=(IndexSet&&) = delete;
  /// Internal identifier, distinct from name and member IDs.
  int Id() const;
  /// Borrowed interpretation scope key: PART.<name> or ASSEMBLY.<name>.
  /// The suffix retains original spelling; valid until destruction.
  const char* Scope() const;
  /// Original set name, retaining input spelling.
  const char* Name() const;
  /// True for an element set, false for a node set.
  bool IsElementSet() const;
  /// Number of sorted unique identifiers.
  int Size() const;
  /// Borrowed read-only member IDs, valid until destruction; null if empty.
  const int* Ids() const;

 private:
  IndexSet(int id, char* scope, char* name, bool element_set,
           int count, int* ids);
  int id_;
  char* scope_;  ///< Owned.
  char* name_;  ///< Owned original spelling.
  bool element_set_;
  int count_;
  int* ids_;  ///< Owned index array, not numerical data.
};

}  // namespace fesa

#endif  // FESA_MODEL_INDEX_SET_H_
