/// @file
/// @brief Source-preserving approved Abaqus input, owned independently of
/// Domain.
#ifndef FESA_IO_INPUT_DECK_H_
#define FESA_IO_INPUT_DECK_H_

#include <map>
#include <string>
#include <vector>

#include "fesa/base/diagnostic.h"

namespace fesa {

/// Lexical scope category, distinct even when Part and Assembly names coincide.
enum class InputScope { kModel, kPart, kAssembly, kStep };

/// Original data row. Text is owned; line is one-based in SourceFile().
struct InputLine {
  InputLine() = default;
  ~InputLine();
  InputLine(const InputLine&) = delete;
  InputLine& operator=(const InputLine&) = delete;
  InputLine(InputLine&&) = delete;
  InputLine& operator=(InputLine&&) = delete;
  std::string* text = nullptr;  ///< Owned original row, excluding CR/LF.
  int line = 0;                 ///< Original physical line.
};

/// Validated logical ELEMENT record; original physical InputLine rows remain.
/// InputBlock owns these records. All IDs and source lines are immutable after
/// successful read. Ten nodes is the maximum of the approved element subset.
struct InputElement {
  InputElement() = default;
  ~InputElement() = default;
  InputElement(const InputElement&) = delete;
  InputElement& operator=(const InputElement&) = delete;
  InputElement(InputElement&&) = delete;
  InputElement& operator=(InputElement&&) = delete;
  int id = 0;               ///< Original element ID.
  int node_count = 0;       ///< Validated physical connectivity size.
  int node_ids[10] = {};    ///< Original node order; only node_count are used.
  int node_lines[10] = {};  ///< Physical source line of each node token.
  int line = 0;             ///< First physical data line (element ID).
  int last_line = 0;        ///< Last physical continuation line.
};

/// Parsed keyword block; owns all pointer members and row entries.
/// Keyword/option keys are canonical ASCII; values retain original spelling.
/// Scope is the original Part/Assembly name, empty for model/step data.
/// Read-only after publication by ReadInputDeck; no concurrent mutation.
struct InputBlock {
  InputBlock() = default;
  ~InputBlock();
  InputBlock(const InputBlock&) = delete;
  InputBlock& operator=(const InputBlock&) = delete;
  InputBlock(InputBlock&&) = delete;
  InputBlock& operator=(InputBlock&&) = delete;
  std::string* keyword = nullptr;   ///< Owned canonical full keyword.
  std::string* original = nullptr;  ///< Owned keyword text with continuations.
  std::string* scope = nullptr;     ///< Owned original scope spelling.
  std::map<std::string, std::string>* options = nullptr;  ///< Owned options.
  std::vector<InputLine*>* rows = nullptr;         ///< Owned list and entries.
  std::vector<InputElement*>* elements = nullptr;  ///< Owned ELEMENT records.
  int line = 0;  ///< First physical keyword line.
  int step = 0;  ///< Zero for model, otherwise one-based step index.
  InputScope scope_kind = InputScope::kModel;  ///< Original lexical scope.
  /// Accepted metadata/output request not applied as written. The log text
  /// is keyword-specific: fixed output, inactive surface or stored constants.
  bool compatibility_warning = false;
  /// Borrowed original option value, or empty string if missing.
  const char* Option(const char* key) const;
  /// Tests presence, including a flag with an empty value.
  bool HasOption(const char* key) const;
};

/// Owns ordered lexical blocks. All borrowed accessors expire at destruction.
class InputDeck {
 public:
  ~InputDeck();
  InputDeck(const InputDeck&) = delete;
  InputDeck& operator=(const InputDeck&) = delete;
  InputDeck(InputDeck&&) = delete;
  InputDeck& operator=(InputDeck&&) = delete;
  /// True if any successfully validated step contains CONTROLS, including
  /// RESET.
  bool HasControls() const;
  /// Original UTF-8 path copied on successful read.
  const char* SourceFile() const;
  /// Borrowed ordered blocks, preserving empty NEW and original definitions.
  const std::vector<InputBlock*>* Blocks() const;

 private:
  friend StatusCode ReadInputDeck(const char*, InputDeck**, Diagnostic*);
  InputDeck(std::string* source, std::vector<InputBlock*>* blocks);
  std::string* source_;               ///< Owned input path.
  std::vector<InputBlock*>* blocks_;  ///< Owned list and entries.
};

/// Reads the approved B33/C3D4/S4 linear subset with optional density.
/// HEAT TRANSFER selects the separate thermal whitelist and single-instance
/// DC3D4/DC3D8 model preparation subset. Thermal steps, IC, conditions and
/// inactive metadata retain source context; model success does not imply a
/// connected HeatTransferAnalysis execution path. See docs/HeatTransfer/README.md.
/// S4 requires homogeneous SHELL SECTION, SECTION INTEGRATION=GAUSS and
/// explicit thickness,2 data; offsets, layers and nodal thickness are rejected.
/// Forward references resolve in BuildDomain; S4 maps to FESA MITC4. Parsed
/// NLGEOM S4 parses the approved static increment/time-control subset; model
/// support does not imply an implemented nonlinear Analysis execution path.
/// Approved HHT DYNAMIC syntax is preserved without resolving increment
/// defaults. BuildDomain connects linear DIRECT definitions; FESASolver still
/// gates unsupported execution before result creation.
/// Global INITIAL CONDITIONS, TYPE=VELOCITY preserves target/DOF/value rows
/// and source context; expanded targets and A26 constraints are not resolved.
/// BuildDomain rejects initial velocity in static-only decks; dynamic velocity
/// resolution is separate. Parsing does not imply initialization capability.
/// See docs/HHT/io_contract.md for the staged dynamic input boundary.
/// @see docs/MITC4/nonlinear_input_contract.md.
/// Solid subset also accepts ten-node continuation records and the types
/// documented in docs/SolidElement/specification.md (input contract).
/// Caller owns *out on success; failure sets null and releases partial state.
/// out is creation-only. Failure d borrows input_path, which must outlive d.
/// Successful calls preserve d; file bytes are never modified.
StatusCode ReadInputDeck(const char* input_path, InputDeck** out,
                         Diagnostic* d);

}  // namespace fesa
#endif  // FESA_IO_INPUT_DECK_H_
