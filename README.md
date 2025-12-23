# LCMS Adduct Finder

**Automated Targeted Feature Extraction & Adduct Verification Tool for LC-MS Data.**

This Python tool is designed for the **targeted analysis** of LC-MS data. By providing a list of **Chemical Formulas**,
it automatically performs a comprehensive scan for various adduct forms (e.g., `[M+H]+`, `[M+Na]+`). It extracts
Extracted Ion Chromatograms (EIC) and rigorously evaluates **peak quality** using Gaussian fitting to determine the
reliability of the detected signals.

---

## Key Features

* **Targeted Extraction**: Instantly converts chemical formulas (e.g., `C6H12O6`) into target m/z values, enabling
  precise extraction of specific metabolites or compounds.
* **Multi-Adduct Verification**:
    * Automatically scans for **14+ different adduct types** (Monomers, Dimers, Na/NH4 adducts, etc.) simultaneously.
    * Helps confirm the identity of a substance by checking if multiple adducts elute at the same Retention Time (RT).
* **Peak Quality & Existence Check**:
    * **Gaussian Scoring**: Fits a Gaussian curve to the raw peak data and calculates an R² score.
    * Distinguishes high-quality peaks ("Excellent/Good") from noise or irregular shapes ("Poor/Noise").
* **Precision Mass Calculation**: Uses high-precision logic considering electron mass:
  $$m/z = \frac{(M \times n + \Delta) - (Charge \times m_e)}{|Charge|}$$
* **Visual Inspection (Optional)**: Saves EIC plots as PNG images with Gaussian fit overlay.
* **MS2 Matching**: Links MS1 features to MS2 events for additional confirmation.

---

## Prerequisites

This tool uses `pythonnet` and `fisher-py` to read Thermo `.raw` files, which requires a .NET runtime. Installation varies by operating system:

### Windows

.NET Framework is typically pre-installed on Windows 10/11. If needed, install the [.NET Runtime](https://dotnet.microsoft.com/download) (version 4.7.2 or later recommended).

### Linux (Ubuntu/Debian)

Install Mono runtime:

```bash
sudo apt update
sudo apt install -y mono-complete
```

For other distributions, see the [Mono installation guide](https://www.mono-project.com/download/stable/#download-lin).

### macOS

Install Mono using Homebrew:

```bash
brew install mono
```

Or download the installer from the [Mono project website](https://www.mono-project.com/download/stable/#download-mac).

---

## Installation

### From PyPI (recommended)

```bash
pip install lib_eic
```

With YAML configuration support:

```bash
pip install lib_eic[yaml]
```

### From Source

```bash
git clone https://github.com/SNUFML/lib_eic.git
cd lib_eic
pip install -e .
```

### Dependencies

- pandas, openpyxl (Excel I/O)
- molmass (mass calculations)
- scipy, numpy (numerical processing)
- matplotlib (plotting)
- fisher-py, pythonnet (Thermo .raw file reading)

---

## Supported Adducts

The tool automatically detects the ionization mode (Positive/Negative) and scans for the following adducts:

| Mode             | Adduct Types                                                              |
|:-----------------|:--------------------------------------------------------------------------|
| **Positive (+)** | `[M+H]+`, `[M+Na]+`, `[M+NH4]+`, `[M+ACN+H]+`, `[2M+H]+`, `[M-H2O+H]+`, etc. |
| **Negative (-)** | `[M-H]-`, `[M+FA-H]-`, `[M-H2O-H]-`, `[2M-H]-`, etc.            |

---

## Usage

### Quick Start

1. Place your Thermo `.raw` files in `./raw` folder
   - Nested layouts are supported, e.g. `./raw/{RP,HILIC}/{1st,2nd}/*.raw`
2. Create an input Excel file (`file_list.xlsx`) with your compound list
3. Run the tool:

```bash
lib_eic
# or
python -m lib_eic
```

### Input Excel Format

`lib_eic` supports two input formats (auto-detected by column names).

#### A) Direct m/z format (recommended for EIC plot generation)

The Excel file contains **separate sheets** for chromatography modes:
- `RP` (Reverse Phase)
- `HILIC` (Hydrophilic Interaction Liquid Chromatography)

The Excel file may contain merged cells in row 1; headers/data start from **row 2**.

| num | File name              | mixture | Compound name | Polarity | m/z     |
|:----|:-----------------------|:--------|:--------------|:---------|:--------|
| 1   | `Library_POS_Mix121`   | 121     | Spermine      | POS      | 203.223 |
| 2   | `Library_POS_Mix121`   | 121     | Putrescine    | POS      | 89.107  |
| 3   | `Library_NEG_Mix121`   | 121     | Glucose       | NEG      | 179.056 |

* **num**: Optional ordering number (used to prefix plot filenames for easier sorting)
* **File name**: Partial raw filename prefix used for matching (e.g., matches `File name.raw`, `File name_2nd.raw`, ...)
* **mixture**: Mixture identifier (used in plot filenames)
* **Compound name**: Display name for plots and Excel output
* **Polarity**: `POS` or `NEG`
* **m/z**: Direct target m/z value

EIC plots are saved under:
`EIC_Plots_Export/{LC mode}/{Polarity}/{File name}/[{num}_]{Compound name}_{Polarity}_{mixture}{suffix}.png`

Notes:
- For direct m/z input (separate `RP`/`HILIC` sheets), if `--raw-folder` contains an `{LC mode}` subfolder, the tool searches that first.
- If raw files are further split by run folders (e.g. `1st/`, `2nd/`), the run label is carried into the output (and plot filenames) to avoid overwrites.

#### B) Formula-based format (legacy)

| RawFile         | Mode | Formula        |
|:----------------|:-----|:---------------|
| `sample_01.raw` | POS  | C6H12O6        |
| `sample_01.raw` | POS  | C10H16N5O13P3  |
| `sample_02.raw` | NEG  | C6H12O6        |

* **RawFile**: Filename (extension can be omitted; `.raw` is appended if missing)
* **Mode**: `POS` or `NEG` (uppercase)
* **Formula**: Chemical formula to analyze

---

## Configuration

### Option 1: Command Line Arguments

```bash
lib_eic --raw-folder ./my_raw_files --input compounds.xlsx --output results.xlsx --ppm 10.0 -v
```

Common options:
- `--raw-folder`: Path to folder containing .raw files (default: `./raw`)
- `--input`: Input Excel file path (default: `file_list.xlsx`)
- `--output`: Output Excel file path (default: `Final_Result_With_Plots.xlsx`)
- `--pivots`: Enable per-target pivot table sheets
- `--no-pivots`: Disable per-target pivot table sheets (default; faster for large runs)
- `--ppm`: Mass tolerance in ppm (default: `10.0`)
- `--no-plots`: Disable EIC plot generation
- `--no-fitting`: Disable Gaussian fitting
- `--no-ms2`: Disable MS2 indexing/matching
- `--workers N`: Number of worker processes (default: auto; use `1` for sequential)
- `--sequential`: Force sequential processing (equivalent to `--workers 1`)
- `--no-progress`: Disable the tqdm progress bar
- `-v, --verbose`: Enable verbose output
- `--help`: Show all available options

Performance notes:
- Parallelism is **file-level** (one worker per raw file); if you only have 1 raw file, speedup is limited.
- If CPU usage stays low, the run is likely bottlenecked by disk I/O (`.raw` reads) or output writing (Excel/plots); try fewer workers and/or an SSD, and consider `--no-plots` and leaving pivot sheets disabled.

### Option 2: YAML Configuration File

Generate a default configuration template:

```bash
lib_eic --generate-config config.yaml
```

Then edit and use:

```bash
lib_eic --config config.yaml
```

Example `config.yaml`:

```yaml
raw_data_folder: "./raw"
input_excel: "file_list.xlsx"
input_sheets: ["RP", "HILIC"]
output_excel: "Final_Result_With_Plots.xlsx"
include_pivot_tables: false
show_progress: true
num_workers: 0          # 0 = auto, 1 = sequential, N = N workers
parallel_mode: "auto"   # "auto", "sequential", "file" (file-level multiprocessing)
ppm_tolerance: 10.0
min_peak_intensity: 100000
enable_fitting: true
enable_plotting: true
enable_ms2: true
export_plot_folder: "EIC_Plots_Export"
area_method: "sum"  # or "trapz"
ms2_match_mode: "rt_linked"  # or "global"
```

### Option 3: Python API

```python
from lib_eic.config import Config
from lib_eic.processor import process_all

# Create configuration
config = Config(
    raw_data_folder="./raw",
    input_excel="file_list.xlsx",
    output_excel="results.xlsx",
    ppm_tolerance=10.0,
    enable_plotting=True,
    enable_fitting=True
)

# Run analysis
process_all(config)
```

#### Advanced: Direct Module Access

```python
from lib_eic.chemistry.mass import get_exact_mass, calculate_target_mz
from lib_eic.chemistry.adducts import get_enabled_adducts
from lib_eic.analysis.eic import build_targets, extract_eic
from lib_eic.analysis.fitting import fit_gaussian_and_score
from lib_eic.io.raw_file import RawFileReader

# Calculate exact mass
mass = get_exact_mass("C6H12O6")  # ~180.0634

# Get enabled adducts for positive mode
adducts = get_enabled_adducts("POS")

# Build targets from formulas
targets = build_targets(["C6H12O6"], adducts)

# Read raw file and extract EIC
with RawFileReader("sample.raw") as reader:
    rt, intensity = reader.get_chromatogram(target_mz=181.0707, ppm=10.0)
```

---

## Output

### 1. Excel Report (`Final_Result_With_Plots.xlsx`)

* **All_Features Sheet**: Complete per-target results table (includes targets that were not reported as features)
  - Adds `EICGenerated` (whether chromatogram extraction returned data) and `FilteredOut` (below `--min-intensity`)
  - Formula-based: RawFile, Mode, Formula, Adduct, mz_theoretical, RT_min, Intensity, Area, GaussianScore, PeakQuality, HasMS2, EICGenerated, FilteredOut
  - Direct m/z: num, RawFile, File name, lc_mode, mixture, Compound name, Polarity, mz_target, RT_min, Intensity, Area, GaussianScore, PeakQuality, HasMS2, EICGenerated, FilteredOut

* **Per-Target Sheets**: Pivot tables for each Formula / Compound name
  - **Area Table**: Peak areas across samples and adducts
  - **Retention Time Table**: RT consistency verification across adducts

### 2. EIC Plots (`EIC_Plots_Export/`)

Visual validation of detected peaks (when plotting is enabled):

* Direct m/z: `EIC_Plots_Export/{Polarity}/{File name}/{Compound name}_{Polarity}_{mixture}{suffix}.png`
* **Blue Line**: Raw EIC data
* **Red Marker**: Apex RT indicator
* **Dashed Line**: Gaussian fit curve (when fitting is enabled)

---

## Quality Scoring System

The Gaussian fitting provides an R² score to assess peak quality:

| R² Score | Label         | Interpretation                              |
|----------|---------------|---------------------------------------------|
| > 0.8    | **Excellent** | High-quality Gaussian peak, reliable signal |
| 0.5-0.8  | **Good**      | Adequate fit, signal is likely real         |
| < 0.5    | **Poor**      | Irregular shape, may be artifact or noise   |
| N/A      | **Not Fitted**| Fitting was disabled                        |

---

## Testing

Run the test suite:

```bash
pytest tests/
pytest -v tests/  # verbose output
pytest --cov=lib_eic tests/  # with coverage report
```

---

## Project Structure

```
lib_eic/
├── chemistry/          # Mass calculations and adduct definitions
│   ├── mass.py         # Exact mass and m/z calculations
│   └── adducts.py      # Adduct type definitions
├── analysis/           # Core analysis modules
│   ├── eic.py          # EIC extraction and target building
│   ├── fitting.py      # Gaussian peak fitting
│   └── ms2.py          # MS2 precursor matching
├── io/                 # Input/output operations
│   ├── raw_file.py     # Thermo .raw file reader
│   ├── excel.py        # Excel I/O
│   └── plotting.py     # EIC visualization
├── config.py           # Configuration management
├── processor.py        # Main processing pipeline
├── cli.py              # Command-line interface
└── validation.py       # Input validation
```

---

## Limitations

- **Supported File Format**: Thermo `.raw` files only (via fisher-py)
- **mzML Not Supported**: Use vendor-specific raw files for best results
- **Platform**: Requires .NET runtime (Windows native, or Mono on Linux/macOS)

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Author

**Jihyun Chun** (jihyun5311@snu.ac.kr)

Repository: https://github.com/SNUFML/lib_eic
