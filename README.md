# FESA

**Finite Element Solver for Structural Analysis**

FESA is a research finite element solver for Windows PCs. It reads Abaqus-style
`.inp` files, performs the analysis, and saves results in HDF5 format. It runs
from the command line and uses OpenMP for CPU parallelism on a single PC.

This repository provides source code for users who want to build the solver.
Executables and external libraries are not included.

## Supported analyses

| Model | Elements and analyses |
| --- | --- |
| Beam | `B33`: two-node Euler–Bernoulli beam; linear static, geometrically nonlinear static, and HHT dynamic analyses |
| Shell | `S4`: FESA's four-node MITC4 shell; linear and geometrically nonlinear static analyses, and HHT dynamic analyses |
| Solid | `C3D4`, `C3D10`, `C3D6`, `C3D5`, `C3D8`, `C3D8R`, `FESA_C3D5R`: linear static and linear dynamic analyses. The five fully integrated types also support geometrically nonlinear static analysis with an incremental elastic formulation |
| Heat transfer | `DC3D4`, `DC3D8`, planar `DS4`: uncoupled steady-state and transient heat transfer, including conduction, heat flux, convection, and radiation to the ambient environment |

Structural materials are primarily isotropic and linearly elastic. Geometrically
nonlinear analyses target large displacements and rotations with small elastic
strains. Plasticity and contact are not supported. Mixed beam/shell/solid
structural models, nonlinear solid dynamics, and coupled thermal-stress analyses
are outside the supported scope. Beams must be initially straight, have a constant
cross-section with no eccentricity, and satisfy `I12=0`.

The solver is under development. Some MITC4 rotation components in large-rotation
cases have not passed comparisons with reference results, and general HHT dynamic
accuracy validation is incomplete. Dynamic acceleration output represents the
time-integration state; reactions are estimates based on the corresponding
inertia. New models require separate comparisons against analytical solutions or
trusted reference results.

## 1. Prerequisites

The current build target is **Windows x64 + MSVC + C++17**.

| Component | Required configuration |
| --- | --- |
| Visual Studio 2026 or Build Tools | Desktop development with C++, MSVC x64 tools, and Windows SDK |
| CMake | Version 4.2 or later, with support for the `Visual Studio 18 2026` generator |
| Intel oneAPI oneMKL | x64 LP64 DLL libraries |
| Intel OpenMP | `libiomp5md.lib` and runtime DLLs from the oneAPI compiler package |
| HDF5 | x64 C DLL development package containing `include/hdf5.h` and `lib/hdf5.lib` |
| Git | Required to clone the repository; downloading a ZIP is also supported |

The build does not download external SDKs automatically. GoogleTest and Python
are not required to build or run the solver. To read results in Python, you can
install `h5py` separately.

The development environment has used MSVC 19.51, CMake 4.4.0, oneMKL runtime
2026.1, and HDF5 2.1.1. Compatibility with other SDK combinations requires
separate verification.

### Get the source

Run the following commands in PowerShell:

```powershell
git clone https://github.com/kjkim2584/FESA.git
Set-Location FESA
```

If you downloaded a ZIP, start in the extracted repository folder.

### Set the SDK paths

**Replace the paths below with your installation locations** before running
these commands. Environment variables apply only to the current PowerShell
session; set them again when you open a new window.

```powershell
$env:FESA_MKL_ROOT = 'C:\Program Files (x86)\Intel\oneAPI\mkl\latest'
$env:FESA_HDF5_ROOT = 'C:\Program Files\HDF_Group\HDF5\2.1.1'

$intelOpenMpRoot = [IO.Path]::GetFullPath(
    (Join-Path $env:FESA_MKL_ROOT '../../compiler/latest'))
$env:PATH = "$env:FESA_MKL_ROOT\bin;$intelOpenMpRoot\bin;$env:FESA_HDF5_ROOT\bin;$env:PATH"
```

`FESA_MKL_ROOT` and `FESA_HDF5_ROOT` must each point to a folder containing
`include`, `lib`, and `bin`. The current CMake configuration expects
`mkl/<version>` and `compiler/latest` under the same oneAPI installation folder.
Install the compiler package that provides Intel OpenMP as well. If your HDF5
DLL requires additional DLLs, add their folders to `PATH` too.

## 2. Build

Run these commands from the repository root. Use an optimized Release build for
analysis runs.

```powershell
cmake --preset windows-msvc -DBUILD_TESTING=OFF
cmake --build --preset windows-release --parallel
```

Check that each command returns exit code `0`. In PowerShell, inspect
`$LASTEXITCODE` immediately after the command. The executable is generated at:

```text
out/cmake/windows-msvc/Release/FESA.exe
```

For a Debug build, use `cmake --build --preset windows-debug --parallel`.
After changing installation paths or SDK combinations, configure a new build
folder with a different `-B` path, or reconfigure using
`cmake --fresh --preset windows-msvc -DBUILD_TESTING=OFF`.

## 3. Run your first analysis

### Create an input file

Save the following content as `cantilever.inp` in a writable working folder.
This linear static example fixes one end of a 1 m beam and applies −100 N in the
global Z direction at the free end. It uses N–m–s units: Young's modulus is in
Pa, cross-sectional area in m², and second moments of area and the torsion
constant in m⁴.

```text
*Heading
FESA cantilever example: N-m-s
*Part, name=Beam
*Node
1, 0, 0, 0
2, 1, 0, 0
*Element, type=B33
1, 1, 2
*Elset, elset=AllElements
1
*Beam General Section, elset=AllElements, material=Steel, section=GENERAL
0.01, 8.333333333333333e-6, 0, 8.333333333333333e-6, 1.4e-5
0, 0, -1
*End Part
*Assembly, name=Assembly
*Instance, name=BeamInstance, part=Beam
*End Instance
*End Assembly
*Material, name=Steel
*Elastic
2.1e11, 0.3
*Boundary
BeamInstance.1, 1, 6, 0
*Step, name=TipLoad, nlgeom=NO
*Static
1, 1
*Cload
BeamInstance.2, 3, -100
*End Step
```

### Run the solver

From the repository root, specify the actual path to your input file. Put paths
containing spaces in quotes.

```powershell
& ./out/cmake/windows-msvc/Release/FESA.exe 'C:\work\cantilever.inp' --threads 1
$LASTEXITCODE
```

The expected free-end Z displacement is approximately `-1.90476e-5 m`.
The solver does not convert units automatically. Use consistent units for
length, force, time, density, and thermal properties. Structural rotations are
in radians.

## 4. Command-line options and configuration

```text
FESA.exe <input.inp> [--threads N] [--config <path>] [--b33-material-first]
FESA.exe --help
```

| Option | Description |
| --- | --- |
| `<input.inp>` | Input file to analyze. Required |
| `--threads N` | Number of CPU threads. Default: `1`. Must be a positive integer no greater than the number of logical CPUs available to the current process |
| `--config <path>` | Path to the convergence and mass configuration file |
| `--b33-material-first` | Use the material tangent for the first correction in a pure B33 nonlinear static analysis. Disabled by default; it may increase the iteration count for some models |
| `--help` | Print usage information |

For structural analyses, convergence settings in the input file's `*CONTROLS`
take precedence. The cfg file is still read and validated for independent
settings such as the mass formulation. Specify its path with `--config`;
otherwise, the solver uses `solver.cfg` next to the executable. If the selected
cfg file does not exist, the solver creates it with default values, so its
folder must be writable. Heat-transfer convergence controls use the input's
thermal settings. You can start with the default configuration.

Only a subset of Abaqus `.inp` keywords and options is supported; this is not
full Abaqus compatibility. Use a single Part and a single Instance without
transformations, and remove unused output requests and unsupported keywords.
`S4` selects FESA's MITC4 formulation. Matching element names do not guarantee
the same formulation or results as Abaqus.

## 5. Inspect the results

Output files are created **in the same folder as the input file**.

| File | Contents |
| --- | --- |
| `cantilever.h5` | Official analysis results: model information and results by step and frame |
| `cantilever_analysis.log` | Analysis progress, selected settings, and error diagnostics |

Running the same input again may overwrite its results and log. Preserve any
results you need before rerunning. Open HDF5 files with HDFView or an HDF5
library such as Python's `h5py`.

Key paths include `/model/node_ids`, `/model/node_coordinates`, and
`/frames/<number>/nodal/U`. Rows in `U` follow the order of `node_ids`; columns
contain global X/Y/Z displacements. Beam and shell rotations are stored in
`UR`, and heat-transfer temperatures in `NT`. Additional results depend on the
element and analysis type.

The root attribute `run_state` is `0` while running, `1` after successful
completion, and `2` after failure. Also check `last_completed_frame` and each
frame's `complete` attribute. A failed run may leave some completed frames, so
the presence of an output file alone does not indicate success.

| Exit code | Meaning |
| --- | --- |
| `0` | Analysis and result output succeeded, or help was printed successfully |
| `1` | Input, configuration, numerical, convergence, or file I/O error |
| `2` | Invalid command-line usage |

## Troubleshooting

| Symptom | What to check |
| --- | --- |
| `FESA_*_ROOT` or missing SDK file error | Verify that the environment variable points to the actual SDK root containing `include` and `lib` |
| Intel OpenMP import library error | Check that oneAPI `compiler/latest/lib/libiomp5md.lib` is installed |
| Missing DLL when starting the solver | Add the SDK `bin` folders to `PATH` in the same PowerShell session |
| MSVC x64 or generator error | Check the Visual Studio 2026 C++ tools, Windows SDK, and CMake version |
| Input rejected | Check the reported file, line, keyword, element type, and support limitations in the log |
| Singular matrix or failure to converge | Check boundary conditions, connectivity, material and section properties, consistent units, and support for the requested analysis |

When reporting a problem, include a reproducible `.inp` file, the command you
ran, its exit code, the analysis log, and your compiler and SDK versions.
