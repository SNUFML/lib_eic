# Performance Optimization Plan for lib_eic

**Context**: Processing ~480 raw files with 300 DPI publication-quality plots.
**Goal**: Reduce wall-clock time by optimizing algorithms and eliminating redundancy.

---

## 1. Plotting Configuration (Quality Up, Overhead Down)

### 1A. Increase Resolution to 300 DPI
**File**: `lib_eic/config.py:37`
```python
# Change default from 120 to 300
plot_dpi: int = 300
```

### 1B. Eliminate Redundant Array Calculations in Plotting
**Files**: `lib_eic/io/plotting.py`, `lib_eic/processor.py`

**Current State** (`plotting.py:176-182`):
```python
max_int = np.max(int_arr) if int_arr.size > 0 else 1.0
apex_idx = np.argmax(int_arr) if int_arr.size > 0 else 0
apex_rt = rt_arr[apex_idx] if rt_arr.size > 0 else 0.0
```
These values are already computed in `_summarize_eic_peak` and the calling code.

**Changes**:
1. Modify `save_eic_plot` and `save_eic_plot_direct_mz` signatures to accept:
   - `apex_rt: float` (already computed)
   - `max_intensity: float` (already computed)
   - `apex_idx: int` (computed alongside)

2. Update `_EICPlotter.update()` to accept pre-computed `rel_abundance` (already does) and avoid redundant calculations

3. Modify callers in `process_raw_file` (line ~847) and `process_raw_file_direct_mz` (line ~1051) to pass pre-computed values

---

## 2. MS2 Analysis Optimization (Lazy Initialization)

**Files**: `lib_eic/processor.py`

**Current State** (`processor.py:762`, `processor.py:921`):
```python
ms2_enabled, ms2_index = _build_ms2_context(reader, config)  # Called unconditionally at start
```

**Problem**: MS2 index is built even when all peaks are below threshold or MS2 matching is never used.

**Changes**:
1. Create `_LazyMS2Context` class that defers `build_ms2_index()` until first access:
```python
class _LazyMS2Context:
    def __init__(self, reader: RawFileReader, config: Config):
        self._reader = reader
        self._config = config
        self._index: Optional[MS2Index] = None
        self._built = False

    @property
    def enabled(self) -> bool:
        return self._config.enable_ms2 and self._reader.has_scan_events_api()

    def get_index(self) -> Optional[MS2Index]:
        if not self._built and self.enabled:
            self._index = build_ms2_index(self._reader)
            self._built = True
        return self._index
```

2. Modify `_maybe_match_ms2` to trigger lazy build only when a peak passes threshold

---

## 3. Data Processing & Raw File Interaction

### 3A. Configurable Chromatogram Batch Size
**Files**: `lib_eic/config.py`, `lib_eic/io/raw_file.py`, `lib_eic/processor.py`

1. Add to `Config` class:
```python
chromatogram_batch_size: int = 256  # Can test with 512, 1024
```

2. Update `get_chromatograms_batch` signature to use config value
3. Pass through from `_get_chromatograms_batch` in processor.py

### 3B. RT Array Optimization (Shared Time Axis Detection)
**File**: `lib_eic/io/raw_file.py:243-283`

**Current State**: Each trace processes its RT array independently with normalization and sorting.

**Optimization** (shape-only check for performance):
```python
def get_chromatograms_batch(...) -> List[Tuple[np.ndarray, np.ndarray]]:
    # ... existing batch extraction code ...

    # After extracting batch, check if all RT arrays share same axis
    first_rt = None
    all_same_rt = True

    for rt_arr, int_arr in zip(positions, intensities):
        rt_normalized = self._normalize_rt_to_minutes(np.asarray(rt_arr, dtype=float))
        if first_rt is None:
            first_rt = rt_normalized
            # Sort once if needed
            if first_rt.size > 1 and np.any(np.diff(first_rt) < 0):
                sort_order = np.argsort(first_rt)
                first_rt = first_rt[sort_order]
            else:
                sort_order = None
        else:
            # Shape-only check (faster than value comparison)
            if rt_normalized.shape != first_rt.shape:
                all_same_rt = False
                break

    # If all same RT, reuse first_rt and only apply sort_order to intensities
    if all_same_rt and first_rt is not None:
        for int_arr in intensities:
            intensity = np.asarray(int_arr, dtype=float)[:first_rt.size]
            if sort_order is not None:
                intensity = intensity[sort_order]
            results.append((first_rt, intensity))
    else:
        # Fall back to per-trace processing (current behavior)
        ...
```

### 3C. File Matching Optimization (Prefix Lookup)
**File**: `lib_eic/processor.py:168-223`

**Current State**: Linear scan through all entries for each partial filename.

**Optimization**: Build a prefix index once, use for all lookups:
```python
def _build_raw_entry_index(entries: List[Path]) -> Dict[str, List[Path]]:
    """Build prefix -> paths index for O(1) lookup."""
    index: Dict[str, List[Path]] = {}
    for entry in entries:
        if not entry.name.lower().endswith(".raw"):
            continue
        stem_lower = entry.stem.lower()
        # Index by first N characters for prefix matching
        for prefix_len in range(1, len(stem_lower) + 1):
            prefix = stem_lower[:prefix_len]
            index.setdefault(prefix, []).append(entry)
    return index
```

Modify `_build_formula_work_items` and `_build_direct_mz_work_items` to build index once and reuse.

---

## 4. Curve Fitting & Data Structure

### 4A. Reduce maxfev in curve_fit
**File**: `lib_eic/analysis/fitting.py:97`

```python
# Change from 2000 to 800
popt, _ = curve_fit(
    gaussian_func,
    x_data,
    y_data,
    p0=[a_guess, x0_guess, sigma_guess],
    maxfev=800,  # Reduced from 2000
)
```

### 4B. Remove Duplicate Result Deduplication
**File**: `lib_eic/io/excel.py:230-234`

**Current State**:
```python
df_all = df_results
if not df_status.empty:
    df_all = pd.concat([df_results, df_status], ignore_index=True, sort=False).drop_duplicates()
```

**Problem**: `results` and `status_rows` contain overlapping entries. Results are appended to both lists, then deduplicated.

**Fix in processor.py** (`process_raw_file` ~870-875, `process_raw_file_direct_mz` ~1108-1113):
```python
# Current: appends to both lists
if status_rows is not None:
    status_rows.append(row_out)
if not filtered_out:
    results.append(row_out)

# Change to: status_rows is canonical, results derived at end
if status_rows is not None:
    status_rows.append(row_out)
# Remove: results.append(row_out)
```

Then in `_write_results_if_any`, derive results from status_rows:
```python
results = [r for r in status_rows if not r.get("FilteredOut", False)]
```

This eliminates the concat + drop_duplicates overhead.

---

## Files to Modify (Summary)

| File | Changes |
|------|---------|
| `lib_eic/config.py` | Add `chromatogram_batch_size`, change `plot_dpi` default |
| `lib_eic/io/plotting.py` | Accept pre-computed apex/max values in save functions |
| `lib_eic/io/raw_file.py` | RT array sharing, configurable batch size |
| `lib_eic/analysis/fitting.py` | Reduce `maxfev` to 800 |
| `lib_eic/analysis/ms2.py` | (No changes - lazy wrapper in processor) |
| `lib_eic/processor.py` | Lazy MS2, pass pre-computed plot values, file index, dedup fix |
| `lib_eic/io/excel.py` | Remove drop_duplicates(), derive results from status_rows |

---

## Implementation Order (with testing after each group)

### Group 1: Quick wins
- 4A: Reduce `maxfev` to 800 (1 line)
- 1A: Set `plot_dpi` default to 300 (1 line)
- **Run tests**: `pytest tests/`

### Group 2: Batch size & plot optimization
- 3A: Configurable batch size (3 files, ~10 lines)
- 1B: Pre-computed plot values (2 files, ~30 lines)
- **Run tests**: `pytest tests/`

### Group 3: Data structure cleanup
- 4B: Deduplication fix - use status_rows as single source (2 files, ~20 lines)
- **Run tests**: `pytest tests/`

### Group 4: Lazy initialization & advanced optimizations
- 2: Lazy MS2 initialization (~40 lines)
- 3B: RT array optimization - **shape-only check** (~30 lines)
- 3C: File matching prefix index (~40 lines)
- **Run tests**: `pytest tests/`

---

## Testing Strategy
Run `pytest tests/` after each implementation group to catch regressions early.
