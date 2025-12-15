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
    * Automatically scans for **16+ different adduct types** (Monomers, Dimers, Na/NH4 adducts, etc.) simultaneously.
    * Helps confirm the identity of a substance by checking if multiple adducts elute at the same Retention Time (RT).
* **Peak Quality & Existence Check**:
    * **Gaussian Scoring**: Fits a Gaussian curve to the raw peak data and calculates an $R^2$ score.
    * Distinguishes high-quality peaks ("Excellent/Good") from noise or irregular shapes ("Poor/Noise").
* **Precision Mass Calculation**: Uses high-precision logic considering electron mass:
  $$m/z = \frac{(M \times n + \Delta) - (Charge \times m_e)}{|Charge|}$$
* **Visual Inspection (Optional)**: Saves EIC plots as PNG images when enabled (see `ENABLE_PLOTTING` in `main.py`).

---

## Supported Adducts

The tool automatically detects the ionization mode (Positive/Negative) and scans for the following adducts to maximize
detection coverage:

| Mode             | Adduct Types                                                   |
|:-----------------|:---------------------------------------------------------------|
| **Positive (+)** | `[M+H]+`, `[M+Na]+`, `[M+NH4]+`, `[M+ACN+H]+`, `[2M+H]+`, etc. |
| **Negative (-)** | `[M-H]-`, `[M+Cl]-`, `[M+HCOO]-`, `[M+FA-H]-`, `[2M-H]-`, etc. |

---

## Prerequisites

### 1. Input Data

This tool reads **Thermo `.raw`** files directly via `fisher-py`.

### 2. Python Dependencies

Install the required libraries:

```bash
pip install pandas openpyxl molmass scipy matplotlib numpy fisher-py pythonnet
```

---

## Usage

### Step 1. Prepare Data Folder

Place your Thermo `.raw` files in `./raw` (or set `RAW_DATA_FOLDER` to your folder).

### Step 2. Prepare Input Excel

Edit the provided `file_list.xlsx` (or create one) with a sheet named **`Final`** (place it next to `main.py`, or change `INPUT_EXCEL_FILE`).
The script supports two input layouts:

**Long format (one formula per row)**:

| RawFile         | Mode | Formula        |
|:----------------|:-----|:---------------|
| `sample_01.raw` | POS  | C6H12O6        |
| `sample_01.raw` | POS  | C10H16N5O13P3  |
| `sample_02.raw` | NEG  | C6H12O6        |


* **RawFile**: Filename (extension can be omitted; `.raw` is appended if missing).
* **Mode**: `POS` or `NEG` (uppercase).
* **Formula columns**: Each non-empty cell is treated as a formula; formulas are normalized and deduplicated per RAW file (order preserved).

### Step 3. Configure & Run

Open the script (`main.py`) and set your data folder path:

```python
# In main.py
RAW_DATA_FOLDER = r"C:\Path\To\Your\raw_Files"
```

Then run the script:

```bash
python main.py
```

---
## Output Examples

The tool generates reports that help you answer: *"Does my target compound exist in this sample, and is the signal
reliable?"*

### 1. Excel Report (`Final_Result_With_Plots.xlsx`)

* **All_Features**: Raw data for every detected adduct.
* **[Formula_Name] Sheets**:
    * **Area Table**: Pivot table of peak areas for each adduct.
    * **Retention Time Table**: Check if RTs are consistent across different adducts.

### 2. EIC Plots (`EIC_Plots_Export/`)

Visual validation of the detected peaks (only generated when `ENABLE_PLOTTING = True`).

* **Blue Line**: Raw EIC data
* **Red Marker/Line**: Apex RT indicator

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
