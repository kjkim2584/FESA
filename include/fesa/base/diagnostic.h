/// @file
/// @brief Caller-owned reusable failure metadata with no dynamic storage.
#ifndef FESA_BASE_DIAGNOSTIC_H_
#define FESA_BASE_DIAGNOSTIC_H_

#include "fesa/base/status_code.h"

namespace fesa {

/// Stage where a checked operation failed.
enum class DiagnosticStage { kNone, kValidation, kAllocation, kComputation };

/// Dimensionless geometry quantity; never inferred from the numeric value.
enum class GeometryQuantity {
  kNone,
  kDeformationJacobian,       ///< det(F), current/initial volume ratio.
  kNormalizedCurrentJacobian, ///< det(Jcurrent/hcurrent), natural mapping.
};

/// Interpretation of borrowed source_file; I/O consumers compose derived paths.
enum class DiagnosticPathKind {
  kLiteralFile,
  kDefaultConfigDirectory,
  kAnalysisInputFile
};

/// One worker owns this record. Pointers are borrowed, never released here.
/// The operation name is static; source_file must outlive diagnostic reporting.
/// Source/context fields are filled by the caller and preserved by Record.
struct Diagnostic {
  Diagnostic() = default;
  Diagnostic(const Diagnostic&) = delete;
  Diagnostic& operator=(const Diagnostic&) = delete;
  Diagnostic(Diagnostic&&) = delete;
  Diagnostic& operator=(Diagnostic&&) = delete;

  /// Changes failure fields and clears active-flux/geometry context;
  /// other source/context fields remain unchanged.
  void Record(StatusCode value, DiagnosticStage failed_stage,
              const char* failed_operation, int native_code = 0);

  StatusCode code = StatusCode::kOk;  ///< Most recently reported failure.
  DiagnosticStage stage = DiagnosticStage::kNone;  ///< Failure stage.
  const char* operation = nullptr;  ///< Borrowed static operation name.
  const char* source_file =
      nullptr;  ///< Borrowed UTF-8 path/base; see path_kind.
  DiagnosticPathKind path_kind = DiagnosticPathKind::kLiteralFile;
  ///< source_file is literal, exe directory for solver.cfg, or input file for
  ///< <base>_analysis.log. Base storage must outlive diagnostic reporting.
  int source_line = 0;  ///< One-based input line; zero means unavailable.
  int first_source_line =
      0;               ///< First definition for a duplicate; otherwise 0.
  int step = 0;        ///< One-based input step; zero means unavailable.
  int node_id = 0;     ///< Original node ID; zero means unavailable.
  int element_id = 0;  ///< Original element ID; zero means unavailable.
  /// Global load/BC component 1..6; zero unavailable. Source-context owners
  /// clear this on unrelated failures; Record preserves supplied context.
  int condition_dof = 0;
  bool duplicate_boundary = false;  ///< True only for expanded BC collision.
  int boundary_dof = 0;  ///< Global component 1..6 when duplicate_boundary.
  double first_boundary_value = 0;  ///< First prescribed value, model units.
  double boundary_value = 0;  ///< Conflicting prescribed value, model units.
  bool active_flux_failure = false;  ///< Undefined selected-active mean.
  int flux_field = 0;                ///< 0 force, 1 moment.
  int flux_count = 0, active_flux_count = 0;  ///< Contribution counts.
  double flux_epsilon_l = 0;  ///< Dimensionless active threshold multiplier.
  double flux_history_max = 0;  ///< Previous committed maximum, field units.
  double flux_threshold = 0;  ///< Active cutoff, field units.
  double flux_current_max = 0;  ///< Same field units as contributions.
  bool geometry_point_available = false;  ///< Failure has a natural location.
  double natural_xi = 0, natural_eta = 0, natural_zeta = 0;
  ///< Dimensionless element coordinates; valid only with the flag above.
  bool geometry_value_available = false;  ///< True only for a computed value.
  GeometryQuantity geometry_quantity = GeometryQuantity::kNone;
  double geometry_value = 0;  ///< Dimensionless, valid only with the flag.
  double geometry_limit = 0;  ///< Strict lower guard for the named quantity.
  int vendor_code = 0;          ///< Native backend error when available.
};

/// Records failure if d is non-null; borrows d and operation during the call.
StatusCode ReportFailure(StatusCode code, DiagnosticStage stage,
                         const char* operation, Diagnostic* d);

}  // namespace fesa

#endif  // FESA_BASE_DIAGNOSTIC_H_
