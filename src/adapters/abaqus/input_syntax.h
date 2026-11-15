/// @file
/// @brief Shared lexical primitives for the approved single-file input adapter.
#ifndef FESA_ADAPTERS_ABAQUS_INPUT_SYNTAX_H_
#define FESA_ADAPTERS_ABAQUS_INPUT_SYNTAX_H_
#include "fesa/io/input_deck.h"
#include "fesa/model/nonlinear_static_parameters.h"
namespace fesa {
namespace abaqus {
/// Writes ASCII uppercase with whitespace removed; inputs borrowed.
void Canonical(const char* text, std::string* out);
/// Writes text with only outer ASCII whitespace removed.
void Trim(const char* text, std::string* out);
/// Splits CSV, preserving internal blanks; permits one empty trailing field.
void Fields(const char* text, std::vector<std::string>* out);
/// Full signed decimal integer conversion within inclusive limits.
bool Integer(const char* text, int minimum, int maximum, int* out);
/// Full finite decimal conversion; E/e/D/d exponents only, no hex/NaN/Inf.
bool Real(const char* text, double* out);
/// Nonempty ASCII unquoted name; period is allowed only for qualified targets.
bool Name(const char* text, bool qualified = false);
/// Canonical ASCII name equality, ignoring whitespace.
bool Equal(const char* left, const char* right);
/// Records input context; source borrowed. Zero first_line means no duplicate.
StatusCode Error(const char* source, int line, int step, const char* operation,
                 Diagnostic* d, StatusCode code = StatusCode::kInputError,
                 int first_line = 0);
/// Shared NODE/ELEMENT/NSET/ELSET row validation, called only after each
/// procedure whitelist and scope validator. Owns parsed element records in block.
StatusCode ValidateInputMeshRows(const char* source, InputBlock* block,
                                  Diagnostic* d);
/// Reads approved STATIC data, returning effective period; input is validated.
bool StaticPeriod(const InputBlock* block, double* period);
/// Parses four STATIC scalars with documented defaults, preserving out on
/// invalid data. Caller fills STEP INC/extrapolation separately. No ownership.
bool NonlinearStaticData(const InputBlock* block,
                         NonlinearStaticParameters* out);
}  // namespace abaqus
}  // namespace fesa
#endif  // FESA_ADAPTERS_ABAQUS_INPUT_SYNTAX_H_
