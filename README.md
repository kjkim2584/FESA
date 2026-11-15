# FESA

**Finite Element Solver for Structural Analysis**

FESA는 Windows PC에서 실행하는 연구용 유한요소 해석 프로그램입니다.
Abaqus 형식의 `.inp` 파일을 읽어 해석하고, 결과를 HDF5 파일에 저장합니다.
명령행에서 실행하며, 단일 PC의 CPU를 OpenMP로 병렬 사용합니다.

이 저장소는 직접 빌드하는 사용자를 위한 소스 배포본입니다.
실행 파일과 외부 라이브러리는 포함하지 않습니다.

## 어떤 해석을 할 수 있나요?

| 대상 | 요소와 해석 |
| --- | --- |
| 보 | `B33`: 2절점 Euler–Bernoulli 보, 선형 정적·기하비선형 정적·HHT 동적 |
| 셸 | `S4`: FESA의 4절점 MITC4 셸, 선형·기하비선형 정적 및 HHT 동적 |
| 솔리드 | `C3D4`, `C3D10`, `C3D6`, `C3D5`, `C3D8`, `C3D8R`, `FESA_C3D5R`: 선형 정적·선형 동적. 이 중 완전적분 5종은 증분 탄성 기하비선형 정적 지원 |
| 열전달 | `DC3D4`, `DC3D8`, 평면 `DS4`: 비연성 정상·비정상 열전달, 전도·열유량·대류·주위 환경과의 복사 |

구조 재료의 기본 범위는 등방 선형탄성입니다. 기하비선형 기능은 큰 변위·회전과
작은 탄성변형률을 대상으로 하며, 소성이나 접촉은 지원하지 않습니다.
보·셸·솔리드를 섞은 구조 모델, 솔리드 비선형 동적, 열응력 연성해석은 지원 범위 밖입니다.
보는 초기 직선·일정 단면·단면 편심 없음·`I12=0` 조건을 사용합니다.

현재 개발 중인 솔버로, MITC4의 일부 대회전 회전 성분은 기준 해석과의 비교를
통과하지 못했으며 일반 HHT 동적 정확도 검증도 완료되지 않았습니다.
동적 출력의 가속도는 시간적분 상태이고 반력은 해당 관성을 사용한 추정값입니다.
새 모델에는 해석해 또는 신뢰할 수 있는 기준 결과와의 별도 비교가 필요합니다.

## 1. 빌드 준비

현재 빌드 대상은 **Windows x64 + MSVC + C++17**입니다.

| 준비 항목 | 필요한 구성 |
| --- | --- |
| Visual Studio 2026 또는 Build Tools | C++ 데스크톱 개발 도구, MSVC x64, Windows SDK |
| CMake | 4.2 이상, `Visual Studio 18 2026` 생성기 지원 |
| Intel oneAPI oneMKL | x64 LP64 DLL 라이브러리 |
| Intel OpenMP | oneAPI compiler 패키지의 `libiomp5md.lib`와 런타임 DLL |
| HDF5 | x64 C 라이브러리의 DLL 개발 패키지, `include/hdf5.h`와 `lib/hdf5.lib` |
| Git | 저장소를 clone할 때 사용. ZIP 다운로드도 가능 |

빌드는 외부 SDK를 자동 다운로드하지 않습니다. GoogleTest와 Python은 솔버 빌드·실행에
필요하지 않습니다. 결과를 Python으로 읽고 싶다면 `h5py`를 별도로 설치할 수 있습니다.

개발 환경에서 사용한 조합은 MSVC 19.51, CMake 4.4.0, oneMKL 런타임 2026.1,
HDF5 2.1.1입니다. 다른 SDK 조합의 호환성은 별도로 확인해야 합니다.

### 소스 받기

PowerShell에서 다음을 실행합니다.

```powershell
git clone https://github.com/kjkim2584/FESA.git
Set-Location FESA
```

ZIP으로 받았다면 압축을 푼 폴더에서 시작합니다.

### SDK 경로 설정

아래 경로를 **자신의 설치 위치로 변경한 뒤** 실행합니다.
환경 변수는 현재 PowerShell 창에만 적용되므로, 새 창을 열면 다시 설정합니다.

```powershell
$env:FESA_MKL_ROOT = 'C:\Program Files (x86)\Intel\oneAPI\mkl\latest'
$env:FESA_HDF5_ROOT = 'C:\Program Files\HDF_Group\HDF5\2.1.1'

$intelOpenMpRoot = [IO.Path]::GetFullPath(
    (Join-Path $env:FESA_MKL_ROOT '../../compiler/latest'))
$env:PATH = "$env:FESA_MKL_ROOT\bin;$intelOpenMpRoot\bin;$env:FESA_HDF5_ROOT\bin;$env:PATH"
```

`FESA_MKL_ROOT`와 `FESA_HDF5_ROOT`는 각각 `include`, `lib`, `bin`을 포함한 폴더입니다.
현재 CMake는 oneAPI의 `mkl/<버전>`과 `compiler/latest`가 같은 oneAPI 폴더 아래에
설치된 구성을 사용합니다. Intel OpenMP를 제공하는 compiler 패키지도 설치해야 합니다.
HDF5 DLL이 추가 DLL에 의존한다면 그 DLL의 폴더 역시 `PATH`에 있어야 합니다.

## 2. 빌드

저장소 루트에서 실행합니다. 실제 해석에는 최적화된 Release 빌드를 권장합니다.

```powershell
cmake --preset windows-msvc -DBUILD_TESTING=OFF
cmake --build --preset windows-release --parallel
```

각 명령의 종료 코드가 `0`인지 확인합니다. PowerShell에서는 실행 직후 `$LASTEXITCODE`로
확인할 수 있습니다. 생성되는 실행 파일은 다음 위치에 있습니다.

```text
out/cmake/windows-msvc/Release/FESA.exe
```

Debug 빌드가 필요하면 `cmake --build --preset windows-debug --parallel`을 사용합니다.
설치 경로·SDK 조합을 변경할 때는 기존 빌드 폴더와 다른 `-B` 경로로 새로 구성하거나,
`cmake --fresh --preset windows-msvc -DBUILD_TESTING=OFF`로 다시 구성합니다.

## 3. 첫 해석 실행

### 입력 파일 만들기

쓰기 가능한 작업 폴더에 아래 내용을 `cantilever.inp`라는 이름으로 저장합니다.
길이 1 m인 보의 한쪽 끝을 고정하고 반대쪽 끝에 전역 Z 방향으로 −100 N을 가하는
선형 정적 예제입니다. 예제의 단위는 N–m–s이며, 탄성계수는 Pa, 단면적은 m²,
단면2차모멘트와 비틀림상수는 m⁴입니다.

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

### 실행하기

저장소 루트에서 입력 파일의 실제 경로를 지정합니다. 공백이 있는 경로는 따옴표로 감쌉니다.

```powershell
& ./out/cmake/windows-msvc/Release/FESA.exe 'C:\work\cantilever.inp' --threads 1
$LASTEXITCODE
```

이 예제에서는 자유단의 Z 변위가 약 `-1.90476e-5 m`입니다.
솔버는 단위를 자동 변환하지 않으므로 길이·힘·시간·밀도·열물성에는 일관된 단위를
사용해야 합니다. 구조 회전은 rad입니다.

## 4. 실행 옵션과 설정

```text
FESA.exe <input.inp> [--threads N] [--config <path>] [--b33-material-first]
FESA.exe --help
```

| 옵션 | 설명 |
| --- | --- |
| `<input.inp>` | 해석할 입력 파일. 필수 |
| `--threads N` | 사용할 CPU 스레드 수. 기본값 `1`. 현재 프로세스가 사용할 수 있는 논리 CPU 수 이내의 양의 정수 |
| `--config <path>` | 수렴·질량 설정 파일 경로 |
| `--b33-material-first` | 순수 B33 비선형 정적해석의 첫 보정에 재료 접선을 사용. 기본값 꺼짐이며 모델에 따라 반복 횟수가 늘 수 있음 |
| `--help` | 사용법 출력 |

구조 해석에서 입력 파일에 `*CONTROLS`가 있으면 입력의 수렴 설정을 우선 사용합니다.
cfg는 질량 방식 등 독립 설정을 위해 이 경우에도 읽고 검증합니다.
cfg 경로는 `--config`로 지정하며, 옵션이 없으면 실행 파일 옆의 `solver.cfg`를
사용합니다. 선택된 cfg가 없으면 기본값으로 생성하므로 해당 폴더에 쓰기 권한이 필요합니다.
열전달의 수렴 제어는 입력의 열 설정을 적용합니다. 처음 사용할 때는 기본 설정으로 시작할 수 있습니다.

입력은 Abaqus `.inp`의 일부 keyword·옵션만 지원하며 Abaqus 전체 호환을 제공하지 않습니다.
단일 Part와 변환 없는 단일 Instance를 사용하고, 사용하지 않는 출력 요청이나
미지원 keyword는 입력에서 제거해야 합니다. `S4`는 FESA MITC4 정식화를 선택합니다.
같은 요소 이름이 Abaqus와 동일한 정식화·결과를 보장하지 않습니다.

## 5. 결과 확인

결과 파일은 **입력 파일이 있는 폴더**에 생성됩니다.

| 파일 | 내용 |
| --- | --- |
| `cantilever.h5` | 공식 해석 결과: 모델 정보와 step·frame별 결과 |
| `cantilever_analysis.log` | 해석 진행, 사용 설정과 오류 진단 |

같은 입력을 다시 실행하면 결과와 로그가 덮어써질 수 있으므로 필요한 결과는 먼저 보관합니다.
HDF5는 HDFView 또는 Python의 `h5py` 같은 HDF5 도구로 열 수 있습니다.

주요 경로는 `/model/node_ids`, `/model/node_coordinates`, `/frames/<번호>/nodal/U`입니다.
`U`의 행은 `node_ids` 순서에 대응하고, 열은 전역 X/Y/Z 변위입니다.
보·셸의 회전은 `UR`, 열전달의 온도는 `NT`에 저장합니다. 제공되는 추가 결과는
요소와 해석 종류에 따라 다릅니다.

루트 attribute `run_state`는 `0` 실행 중, `1` 성공 종료, `2` 실패 종료를 뜻합니다.
`last_completed_frame`과 각 frame의 `complete`도 함께 확인합니다.
실패한 실행에도 일부 완료 frame이 남을 수 있으므로 파일 존재만으로 성공을 판단하지 않습니다.

| 종료 코드 | 의미 |
| --- | --- |
| `0` | 해석·결과 저장 성공, 또는 도움말 출력 성공 |
| `1` | 입력·설정·수치 계산·수렴·파일 입출력 오류 |
| `2` | 명령행 사용법 오류 |

## 자주 발생하는 문제

| 증상 | 확인할 내용 |
| --- | --- |
| `FESA_*_ROOT` 또는 SDK 파일 오류 | 환경 변수가 `include`·`lib`를 포함한 실제 SDK 루트를 가리키는지 확인 |
| Intel OpenMP import library 오류 | oneAPI `compiler/latest/lib/libiomp5md.lib` 설치 여부 확인 |
| 실행 시 DLL을 찾지 못함 | 같은 PowerShell 세션에서 SDK의 `bin` 경로를 `PATH`에 추가했는지 확인 |
| MSVC x64 또는 생성기 오류 | Visual Studio 2026 C++ 도구, Windows SDK, CMake 버전 확인 |
| 입력이 거부됨 | 로그에 표시된 파일·행·keyword·요소 종류와 지원 제한 확인 |
| 특이행렬·미수렴 | 경계조건, 연결, 재료·단면, 일관된 단위와 해당 해석의 지원 범위 확인 |

문제를 보고할 때는 재현 가능한 `.inp`, 실행 명령, 종료 코드, 해석 로그와
컴파일러·SDK 버전을 함께 제공해 주세요.
