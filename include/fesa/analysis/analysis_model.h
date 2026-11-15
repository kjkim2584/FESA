/// @file
/// @brief Linear-static model ownership and immutable input definitions.
#ifndef FESA_ANALYSIS_ANALYSIS_MODEL_H_
#define FESA_ANALYSIS_ANALYSIS_MODEL_H_

#include <map>

#include "fesa/analysis/effective_conditions.h"
#include "fesa/model/domain.h"

namespace fesa {

/// Non-owning snapshot of Domain nodes/elements. Clear before Domain mutation,
/// replacement or destruction. Bind is atomic: failure preserves previous view.
class AnalysisModel {
 public:
  /// Caller owns *out on success; failure nulls it.
  static StatusCode Create(AnalysisModel** out, Diagnostic* d);
  ~AnalysisModel();
  AnalysisModel(const AnalysisModel&) = delete;
  AnalysisModel& operator=(const AnalysisModel&) = delete;
  AnalysisModel(AnalysisModel&&) = delete;
  AnalysisModel& operator=(AnalysisModel&&) = delete;
  /// Borrows domain and its objects until Clear/destruction or next Bind.
  /// Caller must keep Domain alive. Failure leaves this view unchanged.
  StatusCode Bind(const Domain* domain, Diagnostic* d);
  /// Clears old view/condition references FIRST, then binds the full model and
  /// borrows conditions. Both inputs must outlive this view. Failure leaves an
  /// empty/unbound view. No mapping or condition ownership is transferred.
  StatusCode BindStep(const Domain* domain,
                      const EffectiveConditions* conditions, Diagnostic* d);
  /// Borrowed conditions; null after Clear/successful Bind or failed BindStep.
  const EffectiveConditions* Conditions() const;
  /// Drops references only, never deletes Domain objects. Idempotent.
  void Clear();
  /// Borrowed view maps, valid while bound and Domain remains unchanged.
  const std::map<int, Node*>* Nodes() const;
  /// Borrowed view map; this object never owns the elements.
  const std::map<int, Element*>* Elements() const;

 private:
  AnalysisModel() = default;
  const EffectiveConditions* conditions_ = nullptr;  ///< Borrowed.
  const Domain* domain_ = nullptr;  ///< Borrowed until Clear.
  std::map<int, Node*> nodes_;  ///< Non-owning entries.
  std::map<int, Element*> elements_;  ///< Non-owning entries.
};

}  // namespace fesa

#endif  // FESA_ANALYSIS_ANALYSIS_MODEL_H_
