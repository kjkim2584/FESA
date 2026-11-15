#include "fesa/elements/beam3d.h"
#include "fesa/analysis/beam_internal_state.h"
#include "fesa/elements/mitc4_results.h"
#include "fesa/elements/mitc4_shell.h"
#include "elements/solid_element_read_view.h"
#include "fesa/io/result_writer.h"
#include "fesa/analysis/rotation_state.h"

#define NOMINMAX
#include <Windows.h>

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <memory>
#include <new>
#include <vector>

#include "adapters/hdf5/hdf5_result_writer.h"
#include "adapters/hdf5/hdf5_write_support.h"
#include "adapters/hdf5/hdf5_thermal_schema.h"
#include "fesa/analysis/thermal_state.h"
#include "hdf5.h"

namespace fesa {
namespace hdf5 {
namespace {
thread_local FaultPoint next_fault = FaultPoint::kNone;
thread_local int error_depth = 0;
}  // namespace
namespace internal {
bool Inject(FaultPoint point) {
  if (next_fault != point) return false;
  next_fault = FaultPoint::kNone;
  return true;
}
StatusCode Failure(const char* operation, Diagnostic* d,
                   StatusCode code,
                   bool native) {
  if (native) {
    // Capture before cleanup replaces the stack. All vendor storage stays here.
    const hid_t stack = H5Eget_current_stack();
    const ssize_t depth = stack < 0 ? 0 : H5Eget_num(stack);
    error_depth = depth > INT_MAX ? INT_MAX : static_cast<int>(depth);
    if (stack >= 0) H5Eclose_stack(stack);
  }
  if (d)
    d->Record(code,
              code == StatusCode::kAllocationFailure
                  ? DiagnosticStage::kAllocation
                  : DiagnosticStage::kValidation,
              operation, native ? -1 : 0);
  return code;
}
}  // namespace internal
namespace {
using internal::Handle;
using internal::Failure;
using internal::Inject;

bool Check(hid_t result, const char* operation, Diagnostic* d) {
  if (result >= 0) return true;
  Failure(operation, d, StatusCode::kIoFailure, true);
  return false;
}
bool Text(hid_t object, const char* name, const char* const* values, int count,
          bool scalar, Diagnostic* d) {
  Handle type(H5Tcopy(H5T_C_S1), H5Tclose);
  if (!Check(type.get(), "H5Tcopy string", d) ||
      !Check(H5Tset_size(type.get(), H5T_VARIABLE), "H5Tset_size", d) ||
      !Check(H5Tset_cset(type.get(), H5T_CSET_UTF8), "H5Tset_cset", d))
    return false;
  const hsize_t size = static_cast<hsize_t>(count);
  Handle space(
      scalar ? H5Screate(H5S_SCALAR) : H5Screate_simple(1, &size, nullptr),
      H5Sclose);
  if (!Check(space.get(), "H5Screate attribute", d)) return false;
  Handle attr(H5Acreate2(object, name, type.get(), space.get(), H5P_DEFAULT,
                         H5P_DEFAULT),
              H5Aclose);
  if (!Check(attr.get(), "H5Acreate string", d) ||
      !Check(H5Awrite(attr.get(), type.get(), values), "H5Awrite string", d))
    return false;
  return attr.Close(d) && space.Close(d) && type.Close(d);
}
bool ScalarText(hid_t object, const char* name, const char* value,
                Diagnostic* d) {
  return Text(object, name, &value, 1, true, d);
}
bool Number(hid_t object, const char* name, hid_t file_type, hid_t memory_type,
            const void* value, Diagnostic* d) {
  Handle space(H5Screate(H5S_SCALAR), H5Sclose);
  if (!Check(space.get(), "H5Screate scalar", d)) return false;
  const htri_t exists = H5Aexists(object, name);
  if (!Check(exists, "H5Aexists", d)) return false;
  Handle attr(exists ? H5Aopen(object, name, H5P_DEFAULT)
                     : H5Acreate2(object, name, file_type, space.get(),
                                  H5P_DEFAULT, H5P_DEFAULT),
              H5Aclose);
  if (!Check(attr.get(), "H5Aopen/create number", d) ||
      !Check(H5Awrite(attr.get(), memory_type, value), "H5Awrite number", d))
    return false;
  return attr.Close(d) && space.Close(d);
}
bool Integer(hid_t object, const char* name, std::int64_t value,
             Diagnostic* d) {
  return Number(object, name, H5T_STD_I64LE, H5T_NATIVE_INT64, &value, d);
}
bool State(hid_t object, const char* name, unsigned char value, Diagnostic* d) {
  return Number(object, name, H5T_STD_U8LE, H5T_NATIVE_UCHAR, &value, d);
}
bool Real(hid_t object, const char* name, double value, Diagnostic* d) {
  return Number(object, name, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, &value, d);
}
const char* const kXyz[] = {"X", "Y", "Z"};
const char* const kAxes[] = {"t", "n1", "n2"};
const char* const kForces[] = {"N", "V1", "V2", "Ttor", "M1", "M2"};
const char* const kForceUnits[] = {"F", "F", "F", "FL", "FL", "FL"};
const char* const kStress[] = {"sigma_ax"};
const char* const kStrain[] = {"epsilon_ax"};
// Solid nonlinear_contracts.md: identifiers describe recovered kinematics.
const char* NonlinearSolidRule(SolidFormulation formulation) {
  switch (formulation) {
    case SolidFormulation::kTetra4Full:
      return "TET1";
    case SolidFormulation::kTetra10Full:
      return "TET4";
    case SolidFormulation::kHexa8Full:
      return "HEX8";
    case SolidFormulation::kWedge6:
      return "WEDGE2";
    case SolidFormulation::kPyramid5Full:
      return "PYR5_FELIPPA";
    default:
      return nullptr;
  }
}
bool Metadata(hid_t ds, const char* name, Diagnostic* d, bool shell = false,
              bool nonlinear = false, bool initial = false) {
  const char* unit = nullptr;
  const char* evaluation = "increment_end_node";
  const char* coordinates = "global_cartesian";
  const char* const* components = kXyz;
  int count = 3;
  if (std::strcmp(name, "coordinates") == 0 ||
      std::strcmp(name, "node_coordinates") == 0) {
    unit = "L";
    evaluation = "initial_node";
  } else if (std::strcmp(name, "orientation_quaternion") == 0) {
    unit = "1";
    static const char* const quaternion[] = {"w", "x", "y", "z"};
    components = quaternion;
    count = 4;
    if (!ScalarText(ds, "representation", "unit_wxyz_initial_to_current", d))
      return false;
  } else if (std::strcmp(name, "local_axes") == 0) {
    unit = "1";
    evaluation = "initial_element";
    components = kAxes;
  } else if (std::strcmp(name, "U") == 0)
    unit = "L";
  else if (std::strcmp(name, "UR") == 0)
    unit = "rad";
  else if (std::strcmp(name, "V") == 0)
    unit = "L/T";
  else if (std::strcmp(name, "A") == 0)
    unit = "L/T²";
  else if (std::strcmp(name, "angular_velocity") == 0)
    unit = "rad/T";
  else if (std::strcmp(name, "angular_acceleration") == 0)
    unit = "rad/T²";
  else if (std::strcmp(name, "RF") == 0)
    unit = "F";
  else if (std::strcmp(name, "RM") == 0)
    unit = "FL";
  else if (std::strcmp(name, "end_forces") == 0) {
    evaluation = "element_end_a_b";
    coordinates = nonlinear ? "current_material_end_t_n1_n2"
                            : "element_local_t_n1_n2";
    components = kForces;
    count = 6;
    if (!Text(ds, "unit_dimension", kForceUnits, 6, false, d)) return false;
  } else if (std::strcmp(name, "axial_stress") == 0 ||
             std::strcmp(name, "axial_strain") == 0) {
    const bool stress = std::strcmp(name, "axial_stress") == 0;
    unit = stress ? "F/L²" : "1";
    evaluation = "section_average_at_end";
    coordinates = "element_local_t_n1_n2";
    components = stress ? kStress : kStrain;
    count = 1;
    if (nonlinear &&
        (!ScalarText(ds, "measure", stress ? "energetic_material_E_times_Green"
                                          : "Green_centroid", d) ||
         !ScalarText(ds, "section_location", "centroid", d))) return false;
  } else if (std::strcmp(name, "initial_output_axes") == 0 ||
             std::strcmp(name, "current_output_axes") == 0) {
    unit = "1";
    evaluation = std::strcmp(name, "initial_output_axes") == 0
                     ? "initial_element"
                     : "increment_end_element";
    static const char* const output_axes[] = {"a", "b", "n"};
    components = output_axes;
  } else if (std::strcmp(name, "ip_natural_coordinates") == 0) {
    unit = "1";
    coordinates = "natural";
    evaluation = "reference_integration_point";
    static const char* const natural[] = {"xi", "eta", "zeta"};
    components = natural;
  } else if (std::strcmp(name, "element_internal_forces") == 0) {
    unit = "F";
    evaluation = "increment_end_element_node";
    if (shell) {
      unit = nullptr;
      static const char* const forces[] = {"F1", "F2", "F3", "M1", "M2", "M3"};
      components = forces;
      count = 6;
      if (!Text(ds, "unit_dimension", kForceUnits, 6, false, d)) return false;
    }
  } else if (std::strcmp(name, "ip_stress") == 0 ||
             std::strcmp(name, "ip_strain") == 0 ||
             std::strcmp(name, "element_nodal_stress") == 0 ||
             std::strcmp(name, "element_nodal_strain") == 0 ||
             std::strcmp(name, "surface_stress") == 0 ||
             std::strcmp(name, "surface_strain") == 0) {
    const bool stress = std::strstr(name, "stress") != nullptr;
    unit = stress ? "F/L²" : "1";
    evaluation = std::strncmp(name, "ip_", 3) == 0 ? "integration_point"
                 : std::strncmp(name, "surface_", 8) == 0
                     ? "surface_direct"
                     : "element_nodal_extrapolated";
    if (shell) coordinates = "current_element_frame";
    static const char* const tensor[] = {"11", "22", "33", "12", "13", "23"};
    components = tensor;
    count = 6;
    if (!ScalarText(ds, "measure",
                    (shell || nonlinear)
                        ? (stress ? "cauchy" : "logarithmic")
                        : (stress ? "cauchy_small_strain" : "infinitesimal"),
                    d) ||
        (!stress && !ScalarText(ds, "shear_convention", "tensor", d)))
      return false;
  } else
    return true;
  if (initial && std::strcmp(evaluation, "increment_end_node") == 0)
    evaluation = "initial_node";
  if (initial && std::strcmp(evaluation, "increment_end_element_node") == 0)
    evaluation = "initial_element_node";
  if ((shell || nonlinear) && std::strcmp(name, "UR") == 0 &&
      !ScalarText(ds, "representation", "axis_angle", d))
    return false;
  if ((std::strcmp(name, "A") == 0 ||
       std::strcmp(name, "angular_acceleration") == 0) &&
      !ScalarText(ds, "acceleration_representation", "hht_algorithmic", d))
    return false;
  return (!unit || ScalarText(ds, "unit_dimension", unit, d)) &&
         ScalarText(ds, "coordinate_system", coordinates, d) &&
         ScalarText(ds, "evaluation", evaluation, d) &&
         Text(ds, "components", components, count, false, d);
}
bool Dataset(hid_t group, const char* name, hid_t file_type, hid_t memory_type,
             int rank, const hsize_t* dims, const void* data, Diagnostic* d,
             bool nonlinear = false, bool initial = false,
             const internal::DatasetMetadata* metadata = nullptr) {
  Handle space(H5Screate_simple(rank, dims, nullptr), H5Sclose);
  if (!Check(space.get(), "H5Screate dataset", d)) return false;
  Handle ds(H5Dcreate2(group, name, file_type, space.get(), H5P_DEFAULT,
                       H5P_DEFAULT, H5P_DEFAULT),
            H5Dclose);
  if (!Check(ds.get(), "H5Dcreate2", d)) return false;
  if (Inject(FaultPoint::kDataset)) {
    Failure("H5Dwrite payload", d);
    return false;
  }
  if ((dims[0] && !Check(H5Dwrite(ds.get(), memory_type, H5S_ALL, H5S_ALL,
                                  H5P_DEFAULT, data),
                         "H5Dwrite payload", d)) ||
      !(metadata
            ? (ScalarText(ds.get(), "unit_dimension", metadata->unit_dimension,
                          d) &&
               ScalarText(ds.get(), "coordinate_system",
                          metadata->coordinate_system, d) &&
               ScalarText(ds.get(), "evaluation", metadata->evaluation, d) &&
               Text(ds.get(), "components", metadata->components,
                    metadata->component_count, false, d))
            : Metadata(ds.get(), name, d, false, nonlinear, initial)))
    return false;
  return ds.Close(d) && space.Close(d);
}
// Shell payloads write borrowed producer arrays directly. No solver math
// buffers or recovery are allocated here; engineering shear conversion belongs
// to the HDF5 transfer, which leaves the producer arrays unchanged.
hid_t ShellDataset(hid_t group, const char* name, int rank, const hsize_t* dims,
                   Diagnostic* d, bool initial = false) {
  Handle space(H5Screate_simple(rank, dims, nullptr), H5Sclose);
  if (!Check(space.get(), "Shell H5Screate dataset", d)) return -1;
  const hid_t ds = H5Dcreate2(group, name, H5T_IEEE_F64LE, space.get(),
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (!Check(ds, "Shell H5Dcreate", d)) return -1;
  if (Inject(FaultPoint::kDataset)) {
    Failure("Shell H5Dwrite payload", d);
    H5Dclose(ds);
    return -1;
  }
  if (!Metadata(ds, name, d, true, false, initial) || !space.Close(d)) {
    H5Dclose(ds);
    return -1;
  }
  return ds;
}
bool ShellSelection(hid_t dataset, int memory_rank, const hsize_t* dimensions,
                    const hsize_t* memory_start, const hsize_t* memory_count,
                    const hsize_t* file_start, const hsize_t* file_count,
                    const double* data, bool half_shear, Diagnostic* d) {
  Handle memory(H5Screate_simple(memory_rank, dimensions, nullptr), H5Sclose);
  if (!Check(memory.get(), "Shell memory dataspace", d)) return false;
  Handle file(H5Dget_space(dataset), H5Sclose);
  if (!Check(file.get(), "Shell file dataspace", d) ||
      !Check(H5Sselect_hyperslab(memory.get(), H5S_SELECT_SET, memory_start,
                                 nullptr, memory_count, nullptr),
             "Shell memory selection", d) ||
      !Check(H5Sselect_hyperslab(file.get(), H5S_SELECT_SET, file_start,
                                 nullptr, file_count, nullptr),
             "Shell file selection", d))
    return false;
  Handle transfer(H5Pcreate(H5P_DATASET_XFER), H5Pclose);
  if (!Check(transfer.get(), "Shell transfer properties", d) ||
      (half_shear && !Check(H5Pset_data_transform(transfer.get(), "0.5*x"),
                            "Shell tensor shear transform", d)) ||
      !Check(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memory.get(), file.get(),
                      transfer.get(), data),
             "Shell H5Dwrite selection", d))
    return false;
  return transfer.Close(d) && file.Close(d) && memory.Close(d);
}
bool OutputAxesReady(const Matrix* axes) {
  if (!axes || axes->Rows() != 3 || axes->Columns() != 3) return false;
  const double* data = axes->Data();
  for (int i = 0; i < 9; ++i)
    if (!std::isfinite(data[i])) return false;
  for (int row = 0; row < 3; ++row)
    if (data[3 * row] == 0 && data[3 * row + 1] == 0 && data[3 * row + 2] == 0)
      return false;
  return true;
}
bool StringDataset(hid_t group, const char* name, int count,
                   const char* const* values, Diagnostic* d) {
  Handle type(H5Tcopy(H5T_C_S1), H5Tclose);
  if (!Check(type.get(), "H5Tcopy dataset string", d) ||
      !Check(H5Tset_size(type.get(), H5T_VARIABLE), "H5Tset_size dataset", d) ||
      !Check(H5Tset_cset(type.get(), H5T_CSET_UTF8), "H5Tset_cset dataset", d))
    return false;
  const hsize_t dims[] = {static_cast<hsize_t>(count)};
  return Dataset(group, name, type.get(), type.get(), 1, dims, values, d) &&
         type.Close(d);
}
bool Flush(hid_t file, FaultPoint point, const char* operation, Diagnostic* d) {
  if (Inject(point)) {
    Failure(operation, d);
    return false;
  }
  return Check(H5Fflush(file, H5F_SCOPE_GLOBAL), operation, d);
}
bool Utf8(const char* text) {
  return text && MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, text, -1,
                                     nullptr, 0) > 0;
}
bool AbsoluteInput(const char* input) {
  // MSVC Debug u8path uses a noexcept default wstring construction internally;
  // its debug proxy allocation terminates instead of returning bad_alloc.
  // The native pointer-path constructor uses a throwing string copy instead.
  const int count =
      MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, input, -1, nullptr, 0);
  std::unique_ptr<wchar_t[]> wide(new wchar_t[count]);
  if (!MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, input, -1, wide.get(),
                           count))
    return false;
  const std::filesystem::path path(wide.get());
  return path.is_absolute();
}
bool NonAscii(const char* text) {
  for (const unsigned char* p = reinterpret_cast<const unsigned char*>(text);
       *p; ++p)
    if (*p >= 128) return true;
  return false;
}
bool CodePageOverride() {
  char* value = nullptr;
  std::size_t length = 0;
  const errno_t status =
      _dupenv_s(&value, &length, "HDF5_PREFER_WINDOWS_CODE_PAGE");
  std::unique_ptr<char, decltype(&std::free)> owner(value, &std::free);
  if (status == ENOMEM) throw std::bad_alloc();
  return value &&
         (std::strcmp(value, "1") == 0 || _stricmp(value, "true") == 0);
}
StatusCode Context(const char* path, const char* operation, Diagnostic* d,
                   StatusCode code = StatusCode::kInvalidArgument) {
  if (d) {
    d->source_file = path;
    d->duplicate_boundary = false;
    d->condition_dof = 0;
    d->path_kind = DiagnosticPathKind::kLiteralFile;
    d->source_line = d->first_source_line = 0;
  }
  return Failure(operation, d, code);
}
// Opens an attribute-only identity guard without delete sharing. The temporary
// handle excludes path replacement while HDF5 verifies/updates the same file.
HANDLE IdentityGuard(const char* path, BY_HANDLE_FILE_INFORMATION* identity,
                     Diagnostic* d) {
  const int count =
      MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, nullptr, 0);
  std::unique_ptr<wchar_t[]> wide(new wchar_t[count]);
  if (!MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, wide.get(),
                           count)) {
    Failure("Result identity UTF-8 conversion", d);
    return INVALID_HANDLE_VALUE;
  }
  HANDLE handle = CreateFileW(wide.get(), FILE_READ_ATTRIBUTES,
                              FILE_SHARE_READ | FILE_SHARE_WRITE, nullptr,
                              OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
  if (handle == INVALID_HANDLE_VALUE ||
      !GetFileInformationByHandle(handle, identity)) {
    const auto error = GetLastError();
    if (handle != INVALID_HANDLE_VALUE) CloseHandle(handle);
    Failure("Result file identity", d);
    if (d) d->vendor_code = static_cast<int>(error);
    return INVALID_HANDLE_VALUE;
  }
  return handle;
}
bool CloseIdentity(HANDLE handle, Diagnostic* d) {
  if (CloseHandle(handle)) return true;
  const auto error = GetLastError();
  Failure("Close result identity handle", d);
  if (d) d->vendor_code = static_cast<int>(error);
  return false;
}
template <typename T>
bool ScalarEquals(hid_t file, const char* name, hid_t disk_type,
                  hid_t memory_type, T expected, Diagnostic* d) {
  Handle attr(H5Aopen(file, name, H5P_DEFAULT), H5Aclose);
  if (!Check(attr.get(), "Late failure H5Aopen", d)) return false;
  Handle type(H5Aget_type(attr.get()), H5Tclose);
  Handle space(H5Aget_space(attr.get()), H5Sclose);
  if (!Check(type.get(), "Late failure attribute type", d) ||
      !Check(space.get(), "Late failure attribute space", d))
    return false;
  if (H5Tequal(type.get(), disk_type) <= 0 ||
      H5Sget_simple_extent_type(space.get()) != H5S_SCALAR) {
    Failure("Late failure attribute schema mismatch", d);
    return false;
  }
  T actual = 0;
  if (!Check(H5Aread(attr.get(), memory_type, &actual), "Late failure H5Aread",
             d))
    return false;
  if (actual != expected) {
    Failure("Late failure metadata mismatch", d);
    return false;
  }
  return space.Close(d) && type.Close(d) && attr.Close(d);
}
bool TextEquals(hid_t file, const char* name, const char* expected,
                Diagnostic* d) {
  Handle attr(H5Aopen(file, name, H5P_DEFAULT), H5Aclose);
  if (!Check(attr.get(), "Late failure input attribute", d)) return false;
  Handle type(H5Aget_type(attr.get()), H5Tclose);
  Handle space(H5Aget_space(attr.get()), H5Sclose);
  if (!Check(type.get(), "Late failure input type", d) ||
      !Check(space.get(), "Late failure input space", d))
    return false;
  if (H5Tget_class(type.get()) != H5T_STRING ||
      H5Tis_variable_str(type.get()) <= 0 ||
      H5Tget_cset(type.get()) != H5T_CSET_UTF8 ||
      H5Sget_simple_extent_type(space.get()) != H5S_SCALAR) {
    Failure("Late failure input schema mismatch", d);
    return false;
  }
  char* actual = nullptr;
  const bool read = Check(H5Aread(attr.get(), type.get(), &actual),
                          "Late failure read input", d);
  const bool matches = read && actual && std::strcmp(actual, expected) == 0;
  const bool reclaimed =
      !actual ||
      Check(H5Treclaim(type.get(), space.get(), H5P_DEFAULT, &actual),
            "Late failure reclaim input", read ? d : nullptr);
  if (!read || !reclaimed) return false;
  if (!matches) {
    Failure("Late failure input mismatch", d);
    return false;
  }
  return space.Close(d) && type.Close(d) && attr.Close(d);
}
class Hdf5ResultWriter final : public ResultWriter {
 public:
  Hdf5ResultWriter(const char* output, const char* input)
      : output_(output), input_(input) {}
  ~Hdf5ResultWriter() override {
    delete thermal_;
    if (file_ >= 0) {
      Repair();
      H5Fclose(file_);
    }
  }
  Hdf5ResultWriter(const Hdf5ResultWriter&) = delete;
  Hdf5ResultWriter& operator=(const Hdf5ResultWriter&) = delete;
  Hdf5ResultWriter(Hdf5ResultWriter&&) = delete;
  Hdf5ResultWriter& operator=(Hdf5ResultWriter&&) = delete;
  StatusCode Open(Diagnostic* d) {
    file_ = H5Fcreate(output_, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (!Check(file_, "H5Fcreate", d)) return StatusCode::kIoFailure;
    if (!ScalarText(file_, "format", "FESA_RESULT", d) ||
        !Integer(file_, "last_completed_frame", 0, d) ||
        !State(file_, "run_state", 0, d) ||
        !Integer(file_, "last_completed_step", 0, d) ||
        !ScalarText(file_, "input_path", input_, d) ||
        !ScalarText(file_, "unit_convention",
                    "consistent_user_units; rotation=rad; no_conversion", d))
      return Abort(d);
    BY_HANDLE_FILE_INFORMATION identity = {};
    HANDLE guard = IdentityGuard(output_, &identity, d);
    if (guard == INVALID_HANDLE_VALUE) return StatusCode::kIoFailure;
    volume_ = identity.dwVolumeSerialNumber;
    index_high_ = identity.nFileIndexHigh;
    index_low_ = identity.nFileIndexLow;
    return CloseIdentity(guard, d) ? StatusCode::kOk : StatusCode::kIoFailure;
  }
  StatusCode WriteModel(const Domain* domain, Diagnostic* d) override;
  StatusCode WriteFrame(const FrameInfo* frame, const Domain* domain,
                        const AnalysisState* state, Diagnostic* d) override;
  StatusCode WriteThermalFrame(const FrameInfo* frame, const Domain* domain,
                                const ThermalState* state, Diagnostic* d) override;
  StatusCode Finish(bool success, int last_completed_step,
                    Diagnostic* d) override;

 private:
  void SetContext(Diagnostic* d) const {
    if (d) {
      d->source_file = output_;
      d->duplicate_boundary = false;
      d->condition_dof = 0;
      d->path_kind = DiagnosticPathKind::kLiteralFile;
      d->source_line = d->first_source_line = 0;
      d->node_id = d->element_id = 0;
    }
  }
  StatusCode Invalid(const char* operation, Diagnostic* d) {
    Failure(operation, d, StatusCode::kInvalidArgument);
    return Abort(d, StatusCode::kInvalidArgument);
  }
  StatusCode Abort(Diagnostic* d, StatusCode code = StatusCode::kIoFailure) {
    failed_ = true;
    Repair();
    return d ? d->code : code;
  }
  // Best effort only. Close may still persist dirty completion data; this does
  // not overwrite original diagnostics or promise disk rollback/recoverability.
  void Repair() {
    if (file_ < 0) return;
    if (*active_frame_) {
      const hid_t group = H5Gopen2(file_, active_frame_, H5P_DEFAULT);
      if (group >= 0) {
        BestState(group, "complete", 0);
        H5Gclose(group);
      }
    }
    BestState(file_, "run_state", 2);
    const hid_t attr = H5Aopen(file_, "last_completed_step", H5P_DEFAULT);
    if (attr >= 0) {
      H5Awrite(attr, H5T_NATIVE_INT64, &last_step_);
      H5Aclose(attr);
    }
    {
      const hid_t frame_attr =
          H5Aopen(file_, "last_completed_frame", H5P_DEFAULT);
      if (frame_attr >= 0) {
        H5Awrite(frame_attr, H5T_NATIVE_INT64, &frames_);
        H5Aclose(frame_attr);
      }
    }
    H5Fflush(file_, H5F_SCOPE_GLOBAL);
  }
  static void BestState(hid_t object, const char* name, unsigned char value) {
    const hid_t attr = H5Aopen(object, name, H5P_DEFAULT);
    if (attr >= 0) {
      H5Awrite(attr, H5T_NATIVE_UCHAR, &value);
      H5Aclose(attr);
    }
  }
  /// Creates the incomplete frame and metadata, transferring its native ID.
  hid_t BeginFrame(const FrameInfo* frame, const char* procedure, Diagnostic* d);
  /// Shared payload flush/completion/highwater flush; promotes writer history
  /// only after both flushes and checked child close succeed.
  StatusCode CompleteFrame(const FrameInfo* frame, Handle* group, Diagnostic* d);
  StatusCode MarkClosedFailure(int last_completed_step, Diagnostic* d);
  StatusCode ModelPayload(hid_t model, const Domain* domain, Diagnostic* d);
  StatusCode ShellModelPayload(hid_t model, const Domain* domain,
                               Diagnostic* d);
  StatusCode ShellFramePayload(const FrameInfo* frame, hid_t group,
                               hid_t elements, const AnalysisState* state,
                               Diagnostic* d);
  StatusCode SolidModelPayload(hid_t model, const Domain* domain,
                               Diagnostic* d);
  StatusCode SolidFramePayload(hid_t group, bool nonlinear, bool initial,
                               Diagnostic* d);
  StatusCode FramePayload(const FrameInfo* frame, hid_t group, hid_t elements,
                          const AnalysisState* state, Diagnostic* d);
  const char* output_;              // Borrowed through diagnostic reporting.
  const char* input_;               // Borrowed through diagnostic reporting.
  ThermalSchema* thermal_ = nullptr;  // Owned thermal selection/scratch.
  const Domain* domain_ = nullptr;  // Borrowed immutable definitions.
  hid_t file_ = -1;
  int n_ = 0;
  int e_ = 0;
  int node_dofs_ = 6;
  ElementFamily family_ = ElementFamily::kBeam;
  bool solid_nonlinear_ = false;  // Immutable model-wide Solid geometry.
  int connection_count_ = 0;
  int ip_count_ = 0;
  FrameProcedure frame_procedure_ = FrameProcedure::kStatic;
  std::int64_t frame_step_ = 0;
  std::int64_t frame_increment_ = 0;
  double step_start_time_ = 0;
  double step_time_ = 0;
  std::int64_t frames_ = 0;
  std::int64_t last_step_ = 0;
  double total_time_ = 0;
  bool failed_ = false;
  bool closed_success_ = false;  // One optional downgrade after known success.
  DWORD volume_ = 0;
  DWORD index_high_ = 0;
  DWORD index_low_ = 0;
  char active_frame_[64] = {};
};
StatusCode Hdf5ResultWriter::ModelPayload(hid_t model, const Domain* domain,
                                          Diagnostic* d) {
  if (family_ == ElementFamily::kShell)
    return ShellModelPayload(model, domain, d);
  std::vector<int> ids(static_cast<std::size_t>(n_));
  Vector* raw = nullptr;
  auto status = Vector::Create(3 * n_, &raw, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> coordinates(raw);
  int row = 0;
  for (const auto& entry : *domain->Nodes()) {
    ids[row] = entry.first;
    for (int axis = 0; axis < 3; ++axis) {
      double value = 0;
      status = entry.second->Coordinate(axis, &value, d);
      if (status == StatusCode::kOk)
        status = coordinates->Set(3 * row + axis, value, d);
      if (status != StatusCode::kOk) return status;
    }
    ++row;
  }
  const hsize_t nodes[] = {static_cast<hsize_t>(n_)};
  const hsize_t xyz[] = {nodes[0], 3};
  if (!Dataset(model, "node_ids", H5T_STD_I32LE, H5T_NATIVE_INT, 1, nodes,
               ids.data(), d) ||
      !Dataset(model, "node_coordinates", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, 2,
               xyz, coordinates->Data(), d))
    return StatusCode::kIoFailure;
  if (family_ == ElementFamily::kSolid)
    return SolidModelPayload(model, domain, d);
  std::vector<int> element_ids(static_cast<std::size_t>(e_));
  std::vector<int> connectivity(static_cast<std::size_t>(2 * e_));
  raw = nullptr;
  if (e_ > 0) {
    status = Vector::Create(9 * e_, &raw, d);
    if (status != StatusCode::kOk) return status;
  }
  std::unique_ptr<Vector> axes(raw);
  const int index[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  row = 0;
  for (const auto& entry : *domain->Elements()) {
    element_ids[row] = entry.first;
    const auto* beam = dynamic_cast<const Beam3D*>(entry.second);
    const Matrix* frame = beam ? beam->LocalAxes() : nullptr;
    if (!frame || frame->Rows() != 3 || frame->Columns() != 3)
      return Failure("Unsupported element output axes", d,
                     StatusCode::kUnsupported);
    status = entry.second->GetDofNodes(connectivity.data() + 2 * row, 2, d);
    if (status == StatusCode::kOk)
      status = frame->GatherTo(index, 9, axes.get(), 9 * row, d);
    if (status != StatusCode::kOk) return status;
    ++row;
  }
  const hsize_t elements[] = {static_cast<hsize_t>(e_)};
  const hsize_t connected[] = {2 * elements[0]},
                offset_shape[] = {elements[0] + 1},
                local[] = {elements[0], 3, 3};
  std::vector<std::int64_t> offsets(e_ + 1);
  for (int i = 0; i <= e_; ++i) offsets[i] = 2 * i;
  // Preserve the historical linear schema label; finite Beam is actual B33.
  const bool finite_beam = !domain->Steps()->empty() &&
                          domain->Steps()->begin()->second->HasFiniteGeometry();
  std::vector<const char*> types(e_, finite_beam ? "B33" : "B31");
  Handle beam(H5Gcreate2(model, "beam", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
              H5Gclose);
  if (!Check(beam.get(), "H5Gcreate beam model", d))
    return StatusCode::kIoFailure;
  return Dataset(model, "element_ids", H5T_STD_I32LE, H5T_NATIVE_INT, 1,
                 elements, element_ids.data(), d) &&
                 StringDataset(model, "element_types", e_, types.data(), d) &&
                 Dataset(model, "connectivity_offsets", H5T_STD_I64LE,
                         H5T_NATIVE_INT64, 1, offset_shape, offsets.data(),
                         d) &&
                 Dataset(model, "connectivity_node_ids", H5T_STD_I32LE,
                         H5T_NATIVE_INT, 1, connected, connectivity.data(),
                         d) &&
                 Dataset(beam.get(), "local_axes", H5T_IEEE_F64LE,
                         H5T_NATIVE_DOUBLE, 3, local,
                         axes ? axes->Data() : nullptr, d) &&
                 beam.Close(d)
             ? StatusCode::kOk
             : StatusCode::kIoFailure;
}
StatusCode Hdf5ResultWriter::ShellModelPayload(hid_t model,
                                               const Domain* domain,
                                               Diagnostic* d) {
  Handle shell(
      H5Gcreate2(model, "shell", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
      H5Gclose);
  if (!Check(shell.get(), "H5Gcreate shell model", d))
    return StatusCode::kIoFailure;
  std::vector<int> node_ids(n_), element_ids(e_), connections(4 * e_),
      nodal_ids(16 * e_);
  std::vector<std::int64_t> connection_offsets(e_ + 1, 0),
      ip_offsets(e_ + 1, 0), nodal_offsets(e_ + 1, 0),
      surface_offsets(e_ + 1, 0);
  std::vector<const char*> types(e_, "S4"), rules(e_, "MITC4_2X2X2"),
      nodal_locations(16 * e_), surface_locations(8 * e_);
  const hsize_t node_shape[] = {static_cast<hsize_t>(n_), 3},
                axis_shape[] = {static_cast<hsize_t>(e_), 3, 3},
                element_shape[] = {static_cast<hsize_t>(e_)},
                ip_shape[] = {static_cast<hsize_t>(8 * e_), 3},
                nodal_shape[] = {static_cast<hsize_t>(16 * e_)};
  Handle coordinates(ShellDataset(model, "node_coordinates", 2, node_shape, d),
                     H5Dclose);
  if (coordinates.get() < 0) return StatusCode::kIoFailure;
  Handle axes(
      ShellDataset(shell.get(), "initial_output_axes", 3, axis_shape, d),
      H5Dclose);
  if (axes.get() < 0) return StatusCode::kIoFailure;
  Handle thickness(ShellDataset(shell.get(), "thickness", 1, element_shape, d),
                   H5Dclose);
  if (thickness.get() < 0) return StatusCode::kIoFailure;
  Handle ip(ShellDataset(model, "ip_natural_coordinates", 2, ip_shape, d),
            H5Dclose);
  if (ip.get() < 0) return StatusCode::kIoFailure;
  Handle surface(ShellDataset(shell.get(), "surface_point_natural_coordinates",
                              2, ip_shape, d),
                 H5Dclose);
  if (surface.get() < 0) return StatusCode::kIoFailure;
  Handle nodal_z(
      ShellDataset(shell.get(), "element_nodal_thickness_coordinates", 1,
                   nodal_shape, d),
      H5Dclose);
  if (nodal_z.get() < 0) return StatusCode::kIoFailure;
  const hsize_t three[] = {3}, zero[] = {0, 0, 0}, one[] = {1},
                row_count[] = {1, 3}, axis_memory[] = {3, 3},
                axis_count[] = {1, 3, 3}, point_memory[] = {8, 3},
                sixteen[] = {16};
  int row = 0;
  for (const auto& entry : *domain->Nodes()) {
    node_ids[row] = entry.first;
    double xyz[3] = {};
    for (int axis = 0; axis < 3; ++axis) {
      const auto status = entry.second->Coordinate(axis, &xyz[axis], d);
      if (status != StatusCode::kOk) return status;
    }
    const hsize_t start[] = {static_cast<hsize_t>(row), 0};
    if (!ShellSelection(coordinates.get(), 1, three, zero, three, start,
                        row_count, xyz, false, d))
      return StatusCode::kIoFailure;
    ++row;
  }
  const double a = 1 / std::sqrt(3.);
  const double plane[4][2] = {{-a, -a}, {a, -a}, {a, a}, {-a, a}};
  double physical[8][3] = {}, surfaces[8][3] = {}, z[16] = {};
  for (int point = 0; point < 4; ++point)
    for (int layer = 0; layer < 2; ++layer) {
      physical[2 * point + layer][0] = plane[point][0];
      physical[2 * point + layer][1] = plane[point][1];
      physical[2 * point + layer][2] = layer == 0 ? -a : a;
      surfaces[4 * layer + point][0] = plane[point][0];
      surfaces[4 * layer + point][1] = plane[point][1];
      surfaces[4 * layer + point][2] = layer == 0 ? -1 : 1;
    }
  const double layers[] = {-a, a, -1, 1};
  for (int layer = 0; layer < 4; ++layer)
    for (int node = 0; node < 4; ++node) z[4 * layer + node] = layers[layer];
  row = 0;
  for (const auto& entry : *domain->Elements()) {
    const auto* shell_element = static_cast<const Mitc4Shell*>(entry.second);
    element_ids[row] = entry.first;
    const auto status =
        shell_element->GetDofNodes(connections.data() + 4 * row, 4, d);
    if (status != StatusCode::kOk) return status;
    for (int layer = 0; layer < 4; ++layer)
      for (int node = 0; node < 4; ++node) {
        nodal_ids[16 * row + 4 * layer + node] = connections[4 * row + node];
        nodal_locations[16 * row + 4 * layer + node] = layer < 2
                                                           ? "thickness_gauss"
                                                       : layer == 2 ? "SNEG"
                                                                    : "SPOS";
      }
    for (int point = 0; point < 8; ++point)
      surface_locations[8 * row + point] = point < 4 ? "SNEG" : "SPOS";
    connection_offsets[row + 1] = 4 * (row + 1);
    ip_offsets[row + 1] = 8 * (row + 1);
    surface_offsets[row + 1] = 8 * (row + 1);
    nodal_offsets[row + 1] = 16 * (row + 1);
    const hsize_t element_start[] = {static_cast<hsize_t>(row), 0, 0},
                  point_start[] = {static_cast<hsize_t>(8 * row), 0},
                  nodal_start[] = {static_cast<hsize_t>(16 * row)};
    const double h = shell_element->section()->thickness();
    if (!ShellSelection(axes.get(), 2, axis_memory, zero, axis_memory,
                        element_start, axis_count,
                        shell_element->initial_output_axes()->Data(), false, d) ||
        !ShellSelection(thickness.get(), 1, one, zero, one, element_start, one,
                        &h, false, d) ||
        !ShellSelection(ip.get(), 2, point_memory, zero, point_memory,
                        point_start, point_memory, &physical[0][0], false, d) ||
        !ShellSelection(surface.get(), 2, point_memory, zero, point_memory,
                        point_start, point_memory, &surfaces[0][0], false, d) ||
        !ShellSelection(nodal_z.get(), 1, sixteen, zero, sixteen, nodal_start,
                        sixteen, z, false, d))
      return StatusCode::kIoFailure;
    ++row;
  }
  const hsize_t nodes[] = {static_cast<hsize_t>(n_)},
                elements[] = {static_cast<hsize_t>(e_)},
                offsets[] = {static_cast<hsize_t>(e_ + 1)},
                connected[] = {static_cast<hsize_t>(4 * e_)},
                nodal[] = {static_cast<hsize_t>(16 * e_)};
  if (!Dataset(model, "node_ids", H5T_STD_I32LE, H5T_NATIVE_INT, 1, nodes,
               node_ids.data(), d) ||
      !Dataset(model, "element_ids", H5T_STD_I32LE, H5T_NATIVE_INT, 1, elements,
               element_ids.data(), d) ||
      !Dataset(model, "connectivity_offsets", H5T_STD_I64LE, H5T_NATIVE_INT64,
               1, offsets, connection_offsets.data(), d) ||
      !Dataset(model, "connectivity_node_ids", H5T_STD_I32LE, H5T_NATIVE_INT, 1,
               connected, connections.data(), d) ||
      !Dataset(model, "ip_offsets", H5T_STD_I64LE, H5T_NATIVE_INT64, 1, offsets,
               ip_offsets.data(), d) ||
      !Dataset(shell.get(), "element_nodal_result_offsets", H5T_STD_I64LE,
               H5T_NATIVE_INT64, 1, offsets, nodal_offsets.data(), d) ||
      !Dataset(shell.get(), "element_nodal_node_ids", H5T_STD_I32LE,
               H5T_NATIVE_INT, 1, nodal, nodal_ids.data(), d) ||
      !Dataset(shell.get(), "surface_point_offsets", H5T_STD_I64LE,
               H5T_NATIVE_INT64, 1, offsets, surface_offsets.data(), d) ||
      !StringDataset(model, "element_types", e_, types.data(), d) ||
      !StringDataset(model, "stiffness_rule_ids", e_, rules.data(), d) ||
      !StringDataset(shell.get(), "element_nodal_locations", 16 * e_,
                     nodal_locations.data(), d) ||
      !StringDataset(shell.get(), "surface_point_locations", 8 * e_,
                     surface_locations.data(), d))
    return StatusCode::kIoFailure;
  return coordinates.Close(d) && axes.Close(d) && thickness.Close(d) &&
                 ip.Close(d) && surface.Close(d) && nodal_z.Close(d) &&
                 shell.Close(d)
             ? StatusCode::kOk
             : StatusCode::kIoFailure;
}
StatusCode Hdf5ResultWriter::ShellFramePayload(const FrameInfo* frame,
                                               hid_t group, hid_t elements,
                                               const AnalysisState* state,
                                               Diagnostic* d) {
  for (const auto& entry : *domain_->Elements()) {
    const auto* shell_element = static_cast<const Mitc4Shell*>(entry.second);
    if (!OutputAxesReady(shell_element->results()->output_axes())) {
      if (d) d->element_id = entry.first;
      return Failure("Shell frame axes require successful recovery", d,
                     StatusCode::kInvalidArgument);
    }
  }
  const bool dynamic = frame->procedure == FrameProcedure::kDynamic;
  const bool initial = frame->kind == FrameKind::kInitial;
  const char* const names[] = {"U",  "UR",
                               "RF", "RM",
                               "V",  "angular_velocity",
                               "A",  "angular_acceleration"};
  const hsize_t xyz[] = {static_cast<hsize_t>(n_), 3},
                memory_shape[] = {static_cast<hsize_t>(n_), 6},
                zero[] = {0, 0, 0};
  for (int field = 0; field < (dynamic ? 8 : 4); ++field) {
    Handle dataset(ShellDataset(group, names[field], 2, xyz, d, initial),
                   H5Dclose);
    if (dataset.get() < 0) return StatusCode::kIoFailure;
    const hsize_t start[] = {0, static_cast<hsize_t>(3 * (field % 2))};
    const auto* values = field < 2   ? state->CandidateDisplacement()
                         : field < 4 ? state->CandidateReaction()
                         : field < 6
                             ? state->Dynamics()->CandidateVelocity()
                             : state->Dynamics()->CandidateAcceleration();
    if (!ShellSelection(dataset.get(), 2, memory_shape, start, xyz, zero, xyz,
                        values->Data(), false, d) ||
        !dataset.Close(d))
      return StatusCode::kIoFailure;
  }
  Handle payload(
      H5Gcreate2(elements, "shell", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
      H5Gclose);
  if (!Check(payload.get(), "Shell H5Gcreate payload", d))
    return StatusCode::kIoFailure;
  const hsize_t axis_shape[] = {static_cast<hsize_t>(e_), 3, 3},
                axis_memory[] = {3, 3}, axis_count[] = {1, 3, 3},
                force_shape[] = {static_cast<hsize_t>(4 * e_), 6},
                force_memory[] = {4, 6};
  Handle axes(
      ShellDataset(payload.get(), "current_output_axes", 3, axis_shape, d),
      H5Dclose);
  if (axes.get() < 0) return StatusCode::kIoFailure;
  Handle forces(ShellDataset(payload.get(), "element_internal_forces", 2,
                             force_shape, d, initial),
                H5Dclose);
  if (forces.get() < 0) return StatusCode::kIoFailure;
  int row = 0;
  for (const auto& entry : *domain_->Elements()) {
    const auto* results =
        static_cast<const Mitc4Shell*>(entry.second)->results();
    const hsize_t axis_start[] = {static_cast<hsize_t>(row), 0, 0},
                  force_start[] = {static_cast<hsize_t>(4 * row), 0};
    if (!ShellSelection(axes.get(), 2, axis_memory, zero, axis_memory,
                        axis_start, axis_count, results->output_axes()->Data(),
                        false, d) ||
        !ShellSelection(forces.get(), 2, force_memory, zero, force_memory,
                        force_start, force_memory,
                        results->internal_force()->Data(), false, d))
      return StatusCode::kIoFailure;
    ++row;
  }
  const char* const fields[] = {"ip_stress",
                                "ip_strain",
                                "surface_stress",
                                "surface_strain",
                                "element_nodal_stress",
                                "element_nodal_strain"};
  for (int field = 0; field < 6; ++field) {
    const int rows = field < 4 ? 8 : 16;
    const hsize_t dims[] = {static_cast<hsize_t>(rows * e_), 6},
                  memory[] = {static_cast<hsize_t>(rows), 6};
    Handle dataset(ShellDataset(payload.get(), fields[field], 2, dims, d),
                   H5Dclose);
    if (dataset.get() < 0) return StatusCode::kIoFailure;
    row = 0;
    for (const auto& entry : *domain_->Elements()) {
      const auto* results =
          static_cast<const Mitc4Shell*>(entry.second)->results();
      const Matrix* source[] = {
          results->ip_stress(),      results->ip_strain(),
          results->surface_stress(), results->surface_strain(),
          results->nodal_stress(),   results->nodal_strain()};
      for (int part = 0; part < (field % 2 ? 2 : 1); ++part) {
        const hsize_t start[] = {0, static_cast<hsize_t>(3 * part)},
                      count[] = {static_cast<hsize_t>(rows),
                                 static_cast<hsize_t>(field % 2 ? 3 : 6)},
                      target[] = {static_cast<hsize_t>(row * rows),
                                  static_cast<hsize_t>(3 * part)};
        if (!ShellSelection(dataset.get(), 2, memory, start, count, target,
                            count, source[field]->Data(), part == 1, d))
          return StatusCode::kIoFailure;
      }
      ++row;
    }
    if (!dataset.Close(d)) return StatusCode::kIoFailure;
  }
  return axes.Close(d) && forces.Close(d) && payload.Close(d)
             ? StatusCode::kOk
             : StatusCode::kIoFailure;
}
StatusCode Hdf5ResultWriter::SolidModelPayload(hid_t model,
                                               const Domain* domain,
                                               Diagnostic* d) {
  std::vector<int> ids(static_cast<std::size_t>(e_)),
      nodes(static_cast<std::size_t>(connection_count_));
  std::vector<std::int64_t> connections(static_cast<std::size_t>(e_) + 1, 0),
      points(static_cast<std::size_t>(e_) + 1, 0);
  std::vector<const char*> types(static_cast<std::size_t>(e_)),
      rules(static_cast<std::size_t>(e_));
  Vector* raw = nullptr;
  auto status = Vector::Create(3 * ip_count_, &raw, d);
  std::unique_ptr<Vector> coordinates(raw);
  if (status != StatusCode::kOk) return status;
  // Own only the active formulation tables for this metadata write.
  std::unique_ptr<SolidReferenceData> references[8];
  int row = 0;
  for (const auto& entry : *domain->Elements()) {
    SolidElementReadView view;
    status = SolidElementReadView::Bind(entry.second, &view, d);
    if (status != StatusCode::kOk) return status;
    const auto* solid = &view;
    const int formulation = static_cast<int>(solid->formulation());
    if (!references[formulation]) {
      SolidReferenceData* raw_reference = nullptr;
      status = SolidReferenceData::Create(solid->formulation(),
                                          IntegrationPurpose::kStiffness,
                                          &raw_reference, d);
      references[formulation].reset(raw_reference);
      if (status != StatusCode::kOk) return status;
    }
    const auto* reference = references[formulation].get();
    ids[row] = entry.first;
    types[row] = solid->element()->type_name();
    rules[row] = solid_nonlinear_ ? NonlinearSolidRule(solid->formulation())
                                  : reference->rule_id();
    if (!rules[row])
      return Invalid("Unsupported nonlinear Solid output rule", d);
    status = solid->element()->GetDofNodes(nodes.data() + connections[row],
                                solid->node_count(), d);
    if (status != StatusCode::kOk) return status;
    for (int point = 0; point < reference->point_count(); ++point) {
      const int xyz[] = {3 * point, 3 * point + 1, 3 * point + 2};
      status = reference->natural_coordinates()->GatherTo(
          xyz, 3, coordinates.get(),
          static_cast<int>(3 * (points[row] + point)), d);
      if (status != StatusCode::kOk) return status;
    }
    connections[row + 1] = connections[row] + solid->node_count();
    points[row + 1] = points[row] + solid->integration_point_count();
    ++row;
  }
  const hsize_t elements[] = {static_cast<hsize_t>(e_)},
                offsets[] = {static_cast<hsize_t>(e_) + 1},
                connectivity[] = {static_cast<hsize_t>(connection_count_)},
                ips[] = {static_cast<hsize_t>(ip_count_), 3};
  return Dataset(model, "element_ids", H5T_STD_I32LE, H5T_NATIVE_INT, 1,
                 elements, ids.data(), d) &&
                 StringDataset(model, "element_types", e_, types.data(), d) &&
                 StringDataset(model, "stiffness_rule_ids", e_, rules.data(),
                               d) &&
                 Dataset(model, "connectivity_offsets", H5T_STD_I64LE,
                         H5T_NATIVE_INT64, 1, offsets, connections.data(), d) &&
                 Dataset(model, "connectivity_node_ids", H5T_STD_I32LE,
                         H5T_NATIVE_INT, 1, connectivity, nodes.data(), d) &&
                 Dataset(model, "ip_offsets", H5T_STD_I64LE, H5T_NATIVE_INT64,
                         1, offsets, points.data(), d) &&
                 Dataset(model, "ip_natural_coordinates", H5T_IEEE_F64LE,
                         H5T_NATIVE_DOUBLE, 2, ips, coordinates->Data(), d)
             ? StatusCode::kOk
             : StatusCode::kIoFailure;
}
StatusCode Hdf5ResultWriter::WriteModel(const Domain* domain, Diagnostic* d) {
  SetContext(d);
  if (failed_ || file_ < 0 || domain_)
    return Invalid("ResultWriter.WriteModel input/order", d);
  const auto preflight = ValidateResultModel(domain, d);
  if (preflight != StatusCode::kOk) return Abort(d, preflight);
  if (!domain->ThermalElements()->empty()) {
    const auto status = ThermalSchema::Create(domain, &thermal_, d);
    if (status != StatusCode::kOk) return Abort(d, status);
  } else {
    n_ = static_cast<int>(domain->Nodes()->size());
    e_ = static_cast<int>(domain->Elements()->size());
    node_dofs_ = domain->Nodes()->begin()->second->dof_count();
    family_ = domain->Elements()->empty()
                  ? ElementFamily::kBeam
                  : domain->Elements()->begin()->second->family();
    if (family_ == ElementFamily::kSolid && !domain->Steps()->empty()) {
      solid_nonlinear_ =
          domain->Steps()->begin()->second->NonlinearParameters() != nullptr;
      for (const auto& step : *domain->Steps())
        if ((step.second->NonlinearParameters() != nullptr) != solid_nonlinear_)
          return Abort(d, Invalid("Mixed Solid output geometry", d));
    }
    for (const auto& entry : *domain->Elements()) {
      connection_count_ += entry.second->node_count();
      if (family_ == ElementFamily::kSolid) {
        SolidElementReadView view;
        const auto view_status =
            SolidElementReadView::Bind(entry.second, &view, d);
        if (view_status != StatusCode::kOk) return Abort(d, view_status);
        ip_count_ += view.integration_point_count();
      } else if (family_ == ElementFamily::kShell)
        ip_count_ += 8;
    }
  }
  try {
    if (!ScalarText(file_, "element_family",
                    thermal_                           ? thermal_->Family()
                    : family_ == ElementFamily::kBeam  ? "beam"
                    : family_ == ElementFamily::kShell ? "shell"
                                                       : "solid",
                    d))
      return Abort(d);
    Handle model(
        H5Gcreate2(file_, "/model", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
        H5Gclose);
    if (!Check(model.get(), "H5Gcreate model", d)) return Abort(d);
    const auto status = thermal_ ? thermal_->WriteModel(model.get(), d)
                                 : ModelPayload(model.get(), domain, d);
    if (status != StatusCode::kOk) return Abort(d, status);
    Handle frames(
        H5Gcreate2(file_, "/frames", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
        H5Gclose);
    if (!Check(frames.get(), "H5Gcreate frames", d) ||
        !ScalarText(file_, "part_name", domain->PartName(), d) ||
        !ScalarText(file_, "instance_name", domain->InstanceName(), d) ||
        !model.Close(d) || !frames.Close(d))
      return Abort(d);
    domain_ = domain;
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    Failure("ResultWriter model allocation", d, StatusCode::kAllocationFailure);
    return Abort(d, StatusCode::kAllocationFailure);
  }
}
StatusCode Hdf5ResultWriter::FramePayload(const FrameInfo* frame, hid_t group,
                                          hid_t elements,
                                          const AnalysisState* state,
                                          Diagnostic* d) {
  if (family_ == ElementFamily::kShell)
    return ShellFramePayload(frame, group, elements, state, d);
  std::vector<int> index(static_cast<std::size_t>(3 * n_));
  Vector* raw = nullptr;
  auto status = Vector::Create(3 * n_, &raw, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> packed(raw);
  const hsize_t xyz[] = {static_cast<hsize_t>(n_), 3};
  const char* const names[] = {"U", "UR", "RF", "RM"};
  for (int field = 0; field < (family_ == ElementFamily::kSolid ? 2 : 4);
       ++field) {
    const int stored_field =
        family_ == ElementFamily::kSolid ? 2 * field : field;
    for (int node = 0; node < n_; ++node)
      for (int c = 0; c < 3; ++c)
        index[3 * node + c] = node_dofs_ * node + 3 * (stored_field % 2) + c;
    status =
        packed->GatherFrom(stored_field < 2 ? state->CandidateDisplacement()
                                            : state->CandidateReaction(),
                           index.data(), 3 * n_, 0, d);
    if (status != StatusCode::kOk) return status;
    if (!Dataset(group, names[stored_field], H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
                 2, xyz, packed->Data(), d,
                 family_ == ElementFamily::kBeam &&
                     frame->geometry == FrameGeometry::kNonlinear,
                 frame->kind == FrameKind::kInitial))
      return StatusCode::kIoFailure;
  }
  if (frame->procedure == FrameProcedure::kDynamic) {
    const auto* motion = state->Dynamics();
    const Vector* values[] = {motion->CandidateVelocity(),
                              motion->CandidateAcceleration()};
    const bool beam_motion = family_ == ElementFamily::kBeam;
    const char* motion_names[] = {"V", "angular_velocity",
                                  "A", "angular_acceleration"};
    for (int field = 0; field < (beam_motion ? 4 : 2); ++field) {
      const auto* source = values[beam_motion ? field / 2 : field];
      const double* data = source->Data();
      if (beam_motion) {
        for (int node = 0; node < n_; ++node)
          for (int c = 0; c < 3; ++c)
            index[3 * node + c] = 6 * node + 3 * (field % 2) + c;
        status = packed->GatherFrom(source, index.data(), 3 * n_, 0, d);
        if (status != StatusCode::kOk) return status;
        data = packed->Data();
      }
      if (!Dataset(group, motion_names[beam_motion ? field : 2 * field],
                   H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, 2, xyz, data, d, false,
                   frame->kind == FrameKind::kInitial))
        return StatusCode::kIoFailure;
    }
  }
  if (family_ == ElementFamily::kSolid)
    return SolidFramePayload(elements, state->FieldHistory() != nullptr,
                             frame->kind == FrameKind::kInitial, d);
  Handle beam(
      H5Gcreate2(elements, "beam", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
      H5Gclose);
  if (!Check(beam.get(), "H5Gcreate beam results", d))
    return StatusCode::kIoFailure;
  const bool finite_beam = frame->geometry == FrameGeometry::kNonlinear;
  if (finite_beam) {
    const auto* rotations = state->Rotations();
    const Matrix* q = rotations ? rotations->CandidateOrientations() : nullptr;
    if (!q || q->Rows() != n_ || q->Columns() != 4)
      return Invalid("Nonlinear Beam output orientations", d);
    const hsize_t shape[] = {static_cast<hsize_t>(n_), 4};
    if (!Dataset(group, "orientation_quaternion", H5T_IEEE_F64LE,
                 H5T_NATIVE_DOUBLE, 2, shape, q->Data(), d, true,
                 frame->kind == FrameKind::kInitial))
      return StatusCode::kIoFailure;
  }
  const char* const names_e[] = {"end_forces", "axial_stress", "axial_strain"};
  for (int field = 0; field < 3; ++field) {
    const int count = field == 0 ? 12 : 2;
    raw = nullptr;
    if (e_ > 0) {
      status = Vector::Create(count * e_, &raw, d);
      if (status != StatusCode::kOk) return status;
    }
    std::unique_ptr<Vector> values(raw);
    const int start = field == 0 ? 0 : 10 + 2 * field;
    int indices[12] = {};
    for (int i = 0; i < count; ++i) indices[i] = start + i;
    int row = 0;
    for (const auto& entry : *domain_->Elements()) {
      const auto* beam_element = dynamic_cast<const Beam3D*>(entry.second);
      const Vector* results =
          beam_element ? beam_element->EndResults() : nullptr;
      if (!results || results->Size() != 16)
        return Failure("Unsupported element end results", d,
                       StatusCode::kUnsupported);
      status = values->GatherFrom(results, indices, count, count * row, d);
      if (status != StatusCode::kOk) return status;
      ++row;
    }
    const hsize_t dims[] = {static_cast<hsize_t>(e_), 2, 6};
    if (!Dataset(beam.get(), names_e[field], H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
                 field == 0 ? 3 : 2, dims, values ? values->Data() : nullptr,
                 d, finite_beam))
      return StatusCode::kIoFailure;
  }
  return beam.Close(d) ? StatusCode::kOk : StatusCode::kIoFailure;
}
StatusCode Hdf5ResultWriter::SolidFramePayload(hid_t group, bool nonlinear,
                                               bool initial, Diagnostic* d) {
  Handle payload(
      H5Gcreate2(group, "solid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
      H5Gclose);
  if (!Check(payload.get(), "H5Gcreate solid payload", d))
    return StatusCode::kIoFailure;
  const char* const names[] = {"element_internal_forces", "ip_stress",
                               "ip_strain", "element_nodal_stress",
                               "element_nodal_strain"};
  for (int field = 0; field < 5; ++field) {
    const bool ip = field == 1 || field == 2;
    const int columns = field == 0 ? 3 : 6;
    const int rows = ip ? ip_count_ : connection_count_;
    Vector* raw = nullptr;
    auto status = Vector::Create(columns * rows, &raw, d);
    std::unique_ptr<Vector> values(raw);
    if (status != StatusCode::kOk) return status;
    int offset = 0;
    int indices[60] = {};
    for (int i = 0; i < 60; ++i) indices[i] = i;
    for (const auto& entry : *domain_->Elements()) {
      SolidElementReadView view;
      const auto view_status = SolidElementReadView::Bind(entry.second, &view, d);
      if (view_status != StatusCode::kOk) return view_status;
      const auto* solid = &view;
      const int count = columns * (ip ? solid->integration_point_count()
                                      : solid->node_count());
      if (field == 0)
        status = values->GatherFrom(solid->internal_force(), indices, count,
                                    offset, d);
      else {
        const Matrix* source = field == 1   ? solid->ip_stress()
                               : field == 2 ? solid->ip_strain()
                               : field == 3 ? solid->nodal_stress()
                                            : solid->nodal_strain();
        status = source->GatherTo(indices, count, values.get(), offset, d);
      }
      if (status != StatusCode::kOk) return status;
      offset += count;
    }
    const hsize_t dims[] = {static_cast<hsize_t>(rows),
                            static_cast<hsize_t>(columns)};
    if (!Dataset(payload.get(), names[field], H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
                 2, dims, values->Data(), d, nonlinear, initial))
      return StatusCode::kIoFailure;
  }
  return payload.Close(d) ? StatusCode::kOk : StatusCode::kIoFailure;
}
hid_t Hdf5ResultWriter::BeginFrame(const FrameInfo* frame,
                                   const char* procedure, Diagnostic* d) {
  const bool initial = frame->kind == FrameKind::kInitial;
  std::snprintf(active_frame_, sizeof(active_frame_), "/frames/%06lld",
                static_cast<long long>(frame->frame_number));
  Handle group(
      H5Gcreate2(file_, active_frame_, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
      H5Gclose);
  if (!Check(group.get(), "H5Gcreate frame", d)) return -1;
  if (!State(group.get(), "complete", 0, d) ||
      !Integer(group.get(), "frame_number", frame->frame_number, d) ||
      !Integer(group.get(), "step_number", frame->step_number, d) ||
      !Integer(group.get(), "increment_number", frame->increment_number, d) ||
      !ScalarText(group.get(), "step_name", frame->step_name, d) ||
      !Real(group.get(), "step_time", frame->step_time, d) ||
      !Real(group.get(), "total_time", frame->total_time, d) ||
      !State(group.get(), "step_completed", frame->completes_step ? 1 : 0, d) ||
      !ScalarText(group.get(), "analysis_procedure", procedure, d) ||
      !ScalarText(
          group.get(), "geometry_mode",
          frame->geometry == FrameGeometry::kNonlinear ? "nonlinear" : "linear",
          d) ||
      !ScalarText(group.get(), "frame_kind",
                  initial ? "initial" : "increment_end", d))
    return -1;
  return group.Release();
}
StatusCode Hdf5ResultWriter::CompleteFrame(const FrameInfo* frame,
                                           Handle* group, Diagnostic* d) {
  if (!Flush(file_, FaultPoint::kPayloadFlush, "H5Fflush payload", d))
    return Abort(d);
  if (Inject(FaultPoint::kComplete)) {
    Failure("Write frame complete", d);
    return Abort(d);
  }
  if (!State(group->get(), "complete", 1, d)) return Abort(d);
  if (Inject(FaultPoint::kHighWater)) {
    Failure("Write last completed step", d);
    return Abort(d);
  }
  const auto completed_step =
      frame->completes_step ? frame->step_number : last_step_;
  if (!Integer(file_, "last_completed_step", completed_step, d) ||
      !Integer(file_, "last_completed_frame", frame->frame_number, d) ||
      !group->Close(d) ||
      !Flush(file_, FaultPoint::kCompletionFlush, "H5Fflush completion", d))
    return Abort(d);
  frames_ = frame->frame_number;
  last_step_ = completed_step;
  if (frame_step_ != frame->step_number) step_start_time_ = total_time_;
  frame_step_ = frame->step_number;
  frame_increment_ = frame->increment_number;
  frame_procedure_ = frame->procedure;
  step_time_ = frame->step_time;
  total_time_ = frame->total_time;
  *active_frame_ = '\0';
  return StatusCode::kOk;
}
StatusCode Hdf5ResultWriter::WriteThermalFrame(const FrameInfo* frame,
                                               const Domain* domain,
                                               const ThermalState* state,
                                               Diagnostic* d) {
  SetContext(d);
  const bool same = frame && frame_step_ == frame->step_number;
  const bool steady =
      frame && frame->procedure == FrameProcedure::kHeatTransferSteady;
  const bool transient =
      frame && frame->procedure == FrameProcedure::kHeatTransferTransient;
  if (failed_ || file_ < 0 || !thermal_ || domain != domain_ || !frame ||
      !state || frames_ == INT64_MAX || last_step_ >= INT_MAX ||
      frame->frame_number != frames_ + 1 ||
      frame->step_number != last_step_ + 1 ||
      frame->increment_number != (same ? frame_increment_ + 1 : 1) ||
      (!steady && !transient) || frame->geometry != FrameGeometry::kLinear ||
      frame->kind != FrameKind::kIncrementEnd || !Utf8(frame->step_name) ||
      !std::isfinite(frame->step_time) || frame->step_time <= 0 ||
      !std::isfinite(frame->total_time) ||
      frame->total_time !=
          (same ? step_start_time_ : total_time_) + frame->step_time ||
      (same && (frame->procedure != frame_procedure_ ||
                frame->step_time <= step_time_)) ||
      state->Results() != ResultState::kUnconfirmed ||
      !state->HasValidCandidate() ||
      state->ResultStepNumber() != frame->step_number ||
      state->ResultIncrementNumber() != frame->increment_number ||
      !state->CandidateTemperature() || !state->CandidateReactionFlux() ||
      state->CandidateTemperature()->Size() != state->GlobalCount() ||
      state->CandidateReactionFlux()->Size() != state->GlobalCount())
    return Invalid("Thermal frame input/order/state", d);
  const auto layout = state->ValidateDomainLayout(domain, d);
  if (layout != StatusCode::kOk) return Abort(d, layout);
  if (d) d->step = static_cast<int>(frame->step_number);
  Handle group(
      BeginFrame(frame,
                 steady ? "heat_transfer_steady" : "heat_transfer_transient",
                 d),
      H5Gclose);
  if (group.get() < 0) return Abort(d);
  Handle nodal(
      H5Gcreate2(group.get(), "nodal", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
      H5Gclose);
  if (!Check(nodal.get(), "Create thermal nodal group", d)) return Abort(d);
  const auto status = thermal_->WriteFramePayload(nodal.get(), state, d);
  if (status != StatusCode::kOk) return Abort(d, status);
  if (!nodal.Close(d)) return Abort(d);
  return CompleteFrame(frame, &group, d);
}
StatusCode Hdf5ResultWriter::WriteFrame(const FrameInfo* frame,
                                        const Domain* domain,
                                        const AnalysisState* state,
                                        Diagnostic* d) {
  SetContext(d);
  const bool dynamic = frame && frame->procedure == FrameProcedure::kDynamic;
  const bool initial = frame && frame->kind == FrameKind::kInitial;
  const bool same = frame && frame_step_ == frame->step_number;
  if (thermal_ || failed_ || file_ < 0 || !domain_ || domain != domain_ || !frame ||
      !state || frames_ == INT64_MAX || last_step_ >= INT_MAX ||
      frame->frame_number != frames_ + 1 ||
      frame->step_number != last_step_ + 1 ||
      frame->increment_number !=
          (same ? frame_increment_ + 1 : (dynamic ? 0 : 1)) ||
      (same && (initial || frame->procedure != frame_procedure_ ||
                frame->step_time <= step_time_)) ||
      (!same && dynamic && !initial) ||
      (initial &&
       (!dynamic || frame->step_time != 0 || frame->completes_step)) ||
      !Utf8(frame->step_name) || !std::isfinite(frame->step_time) ||
      (!initial && frame->step_time <= 0) ||
      !std::isfinite(frame->total_time) ||
      frame->total_time !=
          (frame_step_ == frame->step_number ? step_start_time_ : total_time_) +
              frame->step_time ||
      state->Results() != ResultState::kUnconfirmed ||
      state->ResultStepNumber() != frame->step_number ||
      state->ResultIncrementNumber() != frame->increment_number ||
      !state->CandidateDisplacement() || !state->CandidateReaction() ||
      state->CandidateDisplacement()->Size() != node_dofs_ * n_ ||
      state->CandidateReaction()->Size() != node_dofs_ * n_)
    return Invalid("ResultWriter.WriteFrame input/order/state", d);
  if (d) d->step = static_cast<int>(frame->step_number);
  if ((frame->procedure != FrameProcedure::kStatic && !dynamic) ||
      (frame->kind != FrameKind::kIncrementEnd && !initial))
    return Invalid("Frame procedure/kind", d);
  if (dynamic &&
      ((family_ != ElementFamily::kSolid && family_ != ElementFamily::kShell &&
        family_ != ElementFamily::kBeam) ||
       (family_ == ElementFamily::kSolid &&
        frame->geometry != FrameGeometry::kLinear) ||
       (family_ == ElementFamily::kBeam &&
        frame->geometry != FrameGeometry::kNonlinear)))
    return Abort(d, Failure("Dynamic frame family/geometry not connected", d,
                            StatusCode::kUnsupported));
  if (dynamic) {
    const auto* motion = state->Dynamics();
    if (!motion || !motion->CandidateVelocity() ||
        !motion->CandidateAcceleration() ||
        motion->CandidateVelocity()->Size() != node_dofs_ * n_ ||
        motion->CandidateAcceleration()->Size() != node_dofs_ * n_ ||
        motion->CandidateStepTime() != frame->step_time ||
        motion->CandidateTotalTime() != frame->total_time)
      return Invalid("Dynamic frame motion/time", d);
  }
  if (dynamic && family_ == ElementFamily::kBeam) {
    const auto* internal = state->BeamInternals();
    const auto* rotations = state->Rotations();
    const auto* q = rotations ? rotations->CandidateOrientations() : nullptr;
    if (!internal || internal->ElementCount() != e_ ||
        !internal->CandidatePosition() || !internal->CandidateVelocity() ||
        !internal->CandidateAcceleration() ||
        internal->CandidatePosition()->Size() != 2 * e_ ||
        internal->CandidateVelocity()->Size() != 2 * e_ ||
        internal->CandidateAcceleration()->Size() != 2 * e_ ||
        !q || q->Rows() != n_ || q->Columns() != 4)
      return Invalid("Dynamic Beam joint internal motion/pose", d);
    int row = 0;
    for (const auto& entry : *domain->Elements())
      if (internal->ElementId(row++) != entry.first)
        return Invalid("Dynamic Beam internal element order", d);
  }
  if ((frame->geometry != FrameGeometry::kLinear &&
       frame->geometry != FrameGeometry::kNonlinear) ||
      (frame->geometry == FrameGeometry::kNonlinear) !=
          (state->FieldHistory() != nullptr) ||
      (family_ == ElementFamily::kSolid &&
       (frame->geometry == FrameGeometry::kNonlinear) != solid_nonlinear_))
    return Invalid("Frame geometry/state mismatch", d);
  try {
    Handle group(BeginFrame(frame, dynamic ? "dynamic" : "static", d), H5Gclose);
    if (group.get() < 0) return Abort(d);
    Handle nodal(
        H5Gcreate2(group.get(), "nodal", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
        H5Gclose);
    Handle elements(H5Gcreate2(group.get(), "elements", H5P_DEFAULT,
                               H5P_DEFAULT, H5P_DEFAULT),
                    H5Gclose);
    if (!Check(nodal.get(), "H5Gcreate nodal", d) ||
        !Check(elements.get(), "H5Gcreate elements", d))
      return Abort(d);
    const auto status =
        FramePayload(frame, nodal.get(), elements.get(), state, d);
    if (status != StatusCode::kOk) return Abort(d, status);
    if (!nodal.Close(d) || !elements.Close(d)) return Abort(d);
    return CompleteFrame(frame, &group, d);
  } catch (const std::bad_alloc&) {
    Failure("ResultWriter frame allocation", d, StatusCode::kAllocationFailure);
    return Abort(d, StatusCode::kAllocationFailure);
  }
}
StatusCode Hdf5ResultWriter::MarkClosedFailure(int last_completed_step,
                                               Diagnostic* d) {
  closed_success_ =
      false;  // Consume allowance even if validation/repair fails.
  failed_ = true;
  if (last_completed_step != last_step_)
    return Failure("Late failure history mismatch", d,
                   StatusCode::kInvalidArgument);
  try {
    if (!Utf8(output_) || !Utf8(input_) ||
        (NonAscii(output_) && CodePageOverride()))
      return Failure("Late failure UTF-8/code-page policy", d);
    if (Inject(FaultPoint::kLateReopen))
      return Failure("Late failure reopen", d);
    BY_HANDLE_FILE_INFORMATION identity = {};
    HANDLE raw = IdentityGuard(output_, &identity, d);
    if (raw == INVALID_HANDLE_VALUE) return StatusCode::kIoFailure;
    std::unique_ptr<void, decltype(&CloseHandle)> guard(raw, CloseHandle);
    if (identity.dwVolumeSerialNumber != volume_ ||
        identity.nFileIndexHigh != index_high_ ||
        identity.nFileIndexLow != index_low_)
      return Failure("Late failure file identity mismatch", d);
    Handle file(H5Fopen(output_, H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
    if (!Check(file.get(), "Late failure H5Fopen", d))
      return StatusCode::kIoFailure;
    if (!TextEquals(file.get(), "format", "FESA_RESULT", d) ||
        !TextEquals(file.get(), "element_family",
                    thermal_ ? thermal_->Family()
                    : family_ == ElementFamily::kBeam    ? "beam"
                    : family_ == ElementFamily::kShell ? "shell"
                                                       : "solid",
                    d) ||
        !ScalarEquals<std::int64_t>(file.get(), "last_completed_frame",
                                    H5T_STD_I64LE, H5T_NATIVE_INT64, frames_,
                                    d) ||
        !ScalarEquals<std::int64_t>(file.get(), "last_completed_step",
                                    H5T_STD_I64LE, H5T_NATIVE_INT64, last_step_,
                                    d) ||
        !ScalarEquals<unsigned char>(file.get(), "run_state", H5T_STD_U8LE,
                                     H5T_NATIVE_UCHAR, 1, d) ||
        !TextEquals(file.get(), "input_path", input_, d))
      return StatusCode::kIoFailure;
    if (Inject(FaultPoint::kLateState))
      return Failure("Late failure run state", d);
    if (!State(file.get(), "run_state", 2, d) ||
        !Flush(file.get(), FaultPoint::kFinishFlush, "Late failure flush", d) ||
        !file.Close(d))
      return StatusCode::kIoFailure;
    return CloseIdentity(guard.release(), d) ? StatusCode::kOk
                                             : StatusCode::kIoFailure;
  } catch (const std::bad_alloc&) {
    return Failure("Late failure allocation", d,
                   StatusCode::kAllocationFailure);
  }
}
StatusCode Hdf5ResultWriter::Finish(bool success, int last_completed_step,
                                    Diagnostic* d) {
  SetContext(d);
  if (file_ < 0 && !success && closed_success_)
    return MarkClosedFailure(last_completed_step, d);
  if (file_ < 0)
    return Failure("ResultWriter already closed", d,
                   StatusCode::kInvalidArgument);
  StatusCode status = StatusCode::kOk;
  if (last_completed_step != last_step_ ||
      (success &&
       (!domain_ || frames_ == 0 || failed_ || frame_step_ != last_step_))) {
    status = Failure("ResultWriter Finish failed/inconsistent history", d);
    failed_ = true;
  }
  if (!success || failed_) Repair();
  if (status == StatusCode::kOk) {
    if (Inject(FaultPoint::kFinishState))
      status = Failure("Write run state", d);
    else if (!State(file_, "run_state", success ? 1 : 2, d) ||
             !Integer(file_, "last_completed_step", last_step_, d) ||
             !Flush(file_, FaultPoint::kFinishFlush, "H5Fflush Finish", d))
      status = StatusCode::kIoFailure;
  }
  if (status != StatusCode::kOk) {
    failed_ = true;
    Repair();
  }
  const hid_t file = file_;
  file_ = -1;
  const bool closed = Check(H5Fclose(file), "H5Fclose",
                            status == StatusCode::kOk ? d : nullptr);
  if (!closed && H5Iis_valid(file) > 0) file_ = file;
  const bool injected = Inject(FaultPoint::kFileClose);
  if (status == StatusCode::kOk && (!closed || injected)) {
    failed_ = true;
    if (injected) return Failure("H5Fclose", d);
    return StatusCode::kIoFailure;
  }
  closed_success_ = success && status == StatusCode::kOk;
  return status;
}
 }  // namespace
 namespace internal {
 bool Check(hid_t result, const char* operation, Diagnostic* d) {
   return ::fesa::hdf5::Check(result, operation, d);
 }
 bool ScalarText(hid_t object, const char* name, const char* value,
                 Diagnostic* d) {
   return ::fesa::hdf5::ScalarText(object, name, value, d);
 }
 bool StringDataset(hid_t group, const char* name, int count,
                    const char* const* values, Diagnostic* d) {
   return ::fesa::hdf5::StringDataset(group, name, count, values, d);
 }
 bool Dataset(hid_t group, const char* name, hid_t file_type, hid_t memory_type,
              int rank, const hsize_t* dims, const void* data, Diagnostic* d,
              const DatasetMetadata* metadata) {
   return ::fesa::hdf5::Dataset(group, name, file_type, memory_type, rank, dims,
                                data, d, false, false, metadata);
 }
 }  // namespace internal
void FailNext(FaultPoint point) { next_fault = point; }
int NativeErrorDepth() { return error_depth; }
}  // namespace hdf5
StatusCode ValidateResultModel(const Domain* domain, Diagnostic* d) {
  if (domain && !domain->ThermalElements()->empty()) {
    if (!hdf5::Utf8(domain->PartName()) || !hdf5::Utf8(domain->InstanceName()))
      return hdf5::Failure("Thermal result UTF-8 metadata", d, StatusCode::kInvalidArgument);
    return hdf5::ValidateThermalResultModel(domain, d);
  }
  if (!domain || domain->Nodes()->empty() ||
      domain->Nodes()->size() > INT_MAX / 6 ||
      domain->Elements()->size() > INT_MAX / 16 ||
      domain->Nodes()->size() > SIZE_MAX / (6 * sizeof(double)) ||
      domain->Elements()->size() > SIZE_MAX / (16 * sizeof(double)) ||
      !hdf5::Utf8(domain->PartName()) || !hdf5::Utf8(domain->InstanceName()))
    return hdf5::Failure("Result model counts/UTF-8 metadata", d,
                         StatusCode::kInvalidArgument);
  const int node_dofs = domain->Nodes()->begin()->second->dof_count();
  const auto family = domain->Elements()->empty()
                          ? ElementFamily::kBeam
                          : domain->Elements()->begin()->second->family();
  if (family == ElementFamily::kShell &&
      (domain->Elements()->size() > INT_MAX / 96 ||
       domain->Elements()->size() > SIZE_MAX / (96 * sizeof(double))))
    return hdf5::Failure("Shell result packing range", d,
                         StatusCode::kInvalidArgument);
  if (node_dofs == 3 &&
      (domain->Elements()->size() > INT_MAX / 60 ||
       domain->Elements()->size() > SIZE_MAX / (60 * sizeof(double))))
    return hdf5::Failure("Solid result model packing range", d,
                         StatusCode::kInvalidArgument);
  for (const auto& entry : *domain->Nodes()) {
    if (entry.second->dof_count() == 0)
      return hdf5::Failure("Thermal node in mechanical result model", d,
                           StatusCode::kInvalidArgument);
    if (entry.first < 1 || entry.first > 999999999 ||
        entry.second->dof_count() != node_dofs)
      return hdf5::Failure("Result node ID/layout", d,
                           StatusCode::kInvalidArgument);
  }
  std::size_t connection_count = 0, point_count = 0;
  for (const auto& entry : *domain->Elements()) {
    const Element* element = entry.second;
    const auto* beam = dynamic_cast<const Beam3D*>(element);
    SolidElementReadView solid_view;
    const auto* solid = SolidElementReadView::Bind(
        element, &solid_view, nullptr) == StatusCode::kOk ? &solid_view : nullptr;
    const auto* shell = dynamic_cast<const Mitc4Shell*>(element);
    if (element->family() != family)
      return hdf5::Failure("Mixed result element families", d,
                           StatusCode::kUnsupported);
    bool supported = false;
    if (beam && node_dofs == 6) {
      const Matrix* axes = beam->LocalAxes();
      const Vector* results = beam->EndResults();
      supported = axes && axes->Rows() == 3 && axes->Columns() == 3 &&
                  results && results->Size() == 16;
    } else if (shell && node_dofs == 6) {
      const auto* results = shell->results();
      supported = shell->node_count() == 4 && shell->dof_count() == 24 &&
                  shell->section() && shell->section()->thickness() > 0 &&
                  hdf5::OutputAxesReady(shell->initial_output_axes()) &&
                  results && results->internal_force() &&
                  results->internal_force()->Size() == 24;
      if (results) {
        const Matrix* stores[] = {
            results->ip_stress(),      results->ip_strain(),
            results->surface_stress(), results->surface_strain(),
            results->nodal_stress(),   results->nodal_strain()};
        for (int i = 0; i < 6; ++i)
          supported = supported && stores[i] &&
                      stores[i]->Rows() == (i < 4 ? 8 : 16) &&
                      stores[i]->Columns() == 6;
        supported = supported && results->output_axes() &&
                    results->output_axes()->Rows() == 3 &&
                    results->output_axes()->Columns() == 3;
      }
      point_count += 8;
    } else if (solid && node_dofs == 3) {
      const int nodes = solid->node_count(),
                points = solid->integration_point_count();
      const auto formulation = solid->formulation();
      const bool tetra4 = formulation == SolidFormulation::kTetra4Full &&
                          std::strcmp(solid->element()->type_name(), "C3D4") == 0 &&
                          nodes == 4 && points == 1;
      const bool tetra10 = formulation == SolidFormulation::kTetra10Full &&
                           std::strcmp(solid->element()->type_name(), "C3D10") == 0 &&
                           nodes == 10 && points == 4;
      const bool hexa8 = formulation == SolidFormulation::kHexa8Full &&
                         std::strcmp(solid->element()->type_name(), "C3D8") == 0 &&
                         nodes == 8 && points == 8;
      const bool reduced_hexa8 =
          formulation == SolidFormulation::kHexa8Reduced &&
          std::strcmp(solid->element()->type_name(), "C3D8R") == 0 && nodes == 8 &&
          points == 1;
      const bool wedge6 = formulation == SolidFormulation::kWedge6 &&
                          std::strcmp(solid->element()->type_name(), "C3D6") == 0 &&
                          nodes == 6 && points == 2;
      const bool pyramid5 = formulation == SolidFormulation::kPyramid5Full &&
                            std::strcmp(solid->element()->type_name(), "C3D5") == 0 &&
                            nodes == 5 && points == 5;
      const bool raw_pyramid5 =
          formulation == SolidFormulation::kPyramid5Reduced &&
          std::strcmp(solid->element()->type_name(), "FESA_C3D5R") == 0 && nodes == 5 &&
          points == 1;
      supported = (tetra4 || tetra10 || hexa8 || reduced_hexa8 ||
                   wedge6 || pyramid5 || raw_pyramid5) &&
                  solid->internal_force() &&
                  solid->internal_force()->Size() == 3 * nodes;
      const Matrix* stores[] = {solid->ip_stress(), solid->ip_strain(),
                                solid->nodal_stress(), solid->nodal_strain()};
      for (int i = 0; i < 4; ++i)
        supported = supported && stores[i] &&
                    stores[i]->Rows() == (i < 2 ? points : nodes) &&
                    stores[i]->Columns() == 6;
      point_count += points;
    }
    if (entry.first < 1 || entry.first > 999999999 || !supported)
      return hdf5::Failure("Result element ID/store shape", d,
                           StatusCode::kUnsupported);
    connection_count += element->node_count();
    if (connection_count > INT_MAX / 6 || point_count > INT_MAX / 6)
      return hdf5::Failure("Result packed offset range", d,
                           StatusCode::kInvalidArgument);
    for (int local = 0; local < element->node_count(); ++local) {
      const Node* node = element->GetNode(local);
      const Node* registered = nullptr;
      if (!node ||
          domain->FindNode(node->Id(), &registered, nullptr) !=
              StatusCode::kOk ||
          registered != node)
        return hdf5::Failure("Result connectivity outside Domain", d,
                             StatusCode::kInvalidArgument);
    }
  }
  if (node_dofs == 3 && domain->Elements()->empty())
    return hdf5::Failure("Solid result model requires elements", d,
                         StatusCode::kUnsupported);
  return StatusCode::kOk;
}
StatusCode OpenResultWriter(const char* output_path, const char* input_path,
                            ResultWriter** out, Diagnostic* d) {
  if (out) *out = nullptr;
  if (!out || !output_path || !*output_path || !input_path || !*input_path ||
      !hdf5::Utf8(output_path) || !hdf5::Utf8(input_path))
    return hdf5::Context(output_path, "OpenResultWriter UTF-8 input", d);
  try {
    if (!hdf5::AbsoluteInput(input_path))
      return hdf5::Context(output_path, "Result input_path must be absolute",
                           d);
    if (hdf5::NonAscii(output_path) && hdf5::CodePageOverride())
      return hdf5::Context(
          output_path,
          "HDF5_PREFER_WINDOWS_CODE_PAGE incompatible with UTF-8 output", d,
          StatusCode::kIoFailure);
    std::unique_ptr<hdf5::Hdf5ResultWriter> writer(
        new hdf5::Hdf5ResultWriter(output_path, input_path));
    if (d) {
      d->source_file = output_path;
      d->duplicate_boundary = false;
      d->condition_dof = 0;
      d->path_kind = DiagnosticPathKind::kLiteralFile;
      d->source_line = d->first_source_line = 0;
    }
    const auto status = writer->Open(d);
    if (status != StatusCode::kOk) return status;
    *out = writer.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return hdf5::Context(output_path, "OpenResultWriter allocation", d,
                         StatusCode::kAllocationFailure);
  }
}
}  // namespace fesa
