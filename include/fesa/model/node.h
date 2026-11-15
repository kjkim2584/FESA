/// @file
/// @brief Linear-static model ownership and immutable input definitions.
#ifndef FESA_MODEL_NODE_H_
#define FESA_MODEL_NODE_H_

#include "fesa/math/vector.h"
#include "fesa/model/nodal_temperature_result.h"

namespace fesa {

/// Immutable original ID/global coordinates; owns mechanical or thermal results.
/// Coordinate units are the input model's consistent length unit.
class Node {
 public:
  /// Caller owns *out on success. On failure *out is null. d is borrowed.
  static StatusCode Create(int id, double x, double y, double z,
                           Node** out, Diagnostic* d);
  /// Creates an immutable 3- or 6-component node; same ownership/failure
  /// contract as the six-component convenience overload above.
  static StatusCode Create(int id, double x, double y, double z, int dof_count,
                           Node** out, Diagnostic* d);
  /// Creates a thermal-only node with 1..20 temperature components (11..30).
  /// Caller owns success; failure nulls *out. No mechanical slots exist.
  /// Invalid component count is kInputError; d is optional and borrowed.
  static StatusCode CreateThermal(int id, double x, double y, double z,
                                  int temperature_dof_count, Node** out,
                                  Diagnostic* d);
  /// Number of mechanical components; zero for a thermal-only node.
  int dof_count() const;
  /// Number of thermal components; zero for a mechanical node.
  int temperature_dof_count() const;
  ~Node();
  Node(const Node&) = delete;
  Node& operator=(const Node&) = delete;
  Node(Node&&) = delete;
  Node& operator=(Node&&) = delete;
  /// Original positive input identifier.
  int Id() const;
  /// Copies the selected global coordinate; invalid axis preserves *out.
  StatusCode Coordinate(int axis, double* out, Diagnostic* d) const;
  /// Overwrites owned global [U,UR,RF,RM] (12) or [U,RF] (6) from values by
  /// explicit CopyFrom. Units: length, radians, force, force*length.
  /// Caller marks common AnalysisState unconfirmed immediately before update.
  /// Failure may invalidate stored results; no restoration or commit occurs.
  /// Thermal-only nodes are rejected without modifying their thermal results.
  StatusCode UpdateResults(const Vector* values, Diagnostic* d);
  /// Borrowed mechanical results, initially zero, valid until Node destruction.
  /// Null for thermal-only nodes.
  /// Six-DOF order: U[3], UR[3], RF[3], RM[3]; rotations in radians.
  /// Three-DOF order: U[3], RF[3], no rotational storage.
  const Vector* Results() const;
  /// Borrowed thermal results until destruction; null for mechanical nodes.
  const NodalTemperatureResult* TemperatureResults() const;
  /// Explicit NT/RFL copy in label order. Mechanical node or invalid input
  /// (including cross-slot NT/RFL alias) is rejected before mutation.
  /// Backend failure may leave partial unconfirmed results.
  /// Caller marks ThermalState unconfirmed before update; no commit/rollback.
  StatusCode UpdateTemperatureResults(const Vector* nt, const Vector* rfl,
                                      Diagnostic* d);

 private:
  Node(int id, double x, double y, double z, Vector* results);
  int id_;
  double x_;
  double y_;
  double z_;
  Vector* results_;  ///< Owned; caller reads only.
  NodalTemperatureResult* temperature_results_ = nullptr;  ///< Owned or null.
};

}  // namespace fesa

#endif  // FESA_MODEL_NODE_H_
