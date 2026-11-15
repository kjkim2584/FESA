/// @file
/// @brief Isolated incremental-elastic solid definitions and official results.
#ifndef FESA_ELEMENTS_INCREMENTAL_SOLID_ELEMENT_H_
#define FESA_ELEMENTS_INCREMENTAL_SOLID_ELEMENT_H_

#include "fesa/elements/solid_reference_data.h"
#include "fesa/model/element.h"
#include "fesa/model/solid_section.h"

namespace fesa {
/// Domain-owned element; borrows immutable nodes/section/material through
/// destruction. Owns official results only, never committed material history.
class IncrementalSolidElement : public Element {
 public:
  ~IncrementalSolidElement() override;
  IncrementalSolidElement(const IncrementalSolidElement&) = delete;
  IncrementalSolidElement& operator=(const IncrementalSolidElement&) = delete;
  IncrementalSolidElement(IncrementalSolidElement&&) = delete;
  IncrementalSolidElement& operator=(IncrementalSolidElement&&) = delete;
  /// Immutable physical family and original input TYPE.
  ElementFamily family() const override;
  const char* type_name() const override;
  /// Immutable full formulation; independent of constitutive identity.
  SolidFormulation formulation() const;
  /// Connectivity size and three global translations per node.
  int node_count() const override;
  int dofs_per_node() const override;
  int dof_count() const override;
  /// Borrowed node; invalid local index returns null.
  const Node* GetNode(int local) const override;
  /// Initial ordered-corner convergence length, L. Failure preserves out.
  StatusCode CharacteristicLength(double* out, Diagnostic* d) const override;
  /// Physical constitutive IP count.
  int integration_point_count() const;
  /// Borrowed global node-major3n resisting force, F. State controls validity.
  const Vector* internal_force() const;
  /// Borrowed IP ng×6 global tensor S/LE, [11,22,33,12,13,23].
  /// Stress has F/L² units; logarithmic strain is dimensionless.
  const Matrix* ip_stress() const;
  const Matrix* ip_strain() const;
  /// Borrowed n×6 element-node tensors, never shared-node averaged.
  const Matrix* nodal_stress() const;
  const Matrix* nodal_strain() const;
  /// Borrowed immutable Domain definitions, valid through destruction.
  const SolidSection* section() const;
  const ElasticMaterial* material() const;
  /// Incremental evaluation uses the immutable initial gradient and the call's
  /// committed baseline. Writes only requested outputs and candidate stress;
  /// preflight preserves outputs, subsequent failures invalidate them.
  bool SupportsTrialEvaluation() const override;
  StatusCode EvaluateTrial(const ElementTrialState* state,
                           const ElementWorkspace* workspace,
                           ElementEvaluationRequest request, Vector* force,
                           Matrix* tangent, TrialEvaluationFailure* failure,
                           Diagnostic* d) const override;
  /// Checks current normalized sampled Jacobians without constitutive recovery.
  bool SupportsTrialGeometryCheck() const override;
  StatusCode CheckTrialGeometry(const ElementTrialState* state,
                                const ElementWorkspace* workspace,
                                TrialEvaluationFailure* failure,
                                Diagnostic* d) const override;
  /// Recomputes force and S/LE at the candidate. Failure leaves official
  /// payload unconfirmed; Analysis owns validation and history commitment.
  StatusCode UpdateTrialResults(const ElementTrialState* state,
                                const ElementWorkspace* workspace,
                                TrialEvaluationFailure* failure,
                                Diagnostic* d) override;

 protected:
  /// Full formulation supplied by a validating factory; owns no definitions.
  IncrementalSolidElement(int id, const SolidSection* section,
                         SolidFormulation formulation);
  /// Validates and borrows n distinct3DOF nodes, allocates exact result storage.
  /// Failure leaves partial owned storage for destruction by the factory.
  StatusCode Initialize(const Node* const* nodes, Diagnostic* d);
  /// Factory-only initial H/q initialization and initial-form geometry guard.
  /// Failure invalidates definitions; the factory destroys partial storage.
  StatusCode ValidateGeometry(Diagnostic* d) const;

 private:
  const Node** nodes_ = nullptr;    ///< Owned array, borrowed Domain nodes.
  const SolidSection* section_;     ///< Borrowed immutable Domain section.
  const SolidFormulation formulation_;  ///< Immutable original full TYPE.
  const int node_count_, point_count_;  ///< Immutable, checked against tables.
  Matrix** gradients_ = nullptr;   ///< Owned array and ng initial3×n H, 1/L.
  Vector* weight_ = nullptr;       ///< Owned ng initial volume weights, L³.
  Vector* force_ = nullptr;         ///< Owned official3n force, initially zero.
  Matrix* stress_ = nullptr;        ///< Owned ng×6 official Cauchy stress.
  Matrix* strain_ = nullptr;        ///< Owned ng×6 official logarithmic strain.
  Matrix* nodal_stress_ = nullptr;  ///< Owned n×6 recovered stress.
  Matrix* nodal_strain_ = nullptr;  ///< Owned n×6 recovered strain.
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_INCREMENTAL_SOLID_ELEMENT_H_
