from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from lib_eic.config import Config
from lib_eic.io.excel import write_results_excel
from lib_eic.processor import process_raw_file_direct_mz


class _DummyReader:
    def __init__(self, filename: str, eics: list[tuple[np.ndarray, np.ndarray]]):
        self.filename = filename
        self._eics = eics

    def has_scan_events_api(self) -> bool:
        return False

    def has_multi_chromatogram_api(self) -> bool:
        return True

    def get_chromatograms_batch(self, target_mzs: list[float], ppm: float) -> list[tuple[np.ndarray, np.ndarray]]:
        assert len(target_mzs) == len(self._eics)
        return self._eics


def test_process_raw_file_direct_mz_collects_status_rows() -> None:
    config = Config(
        ppm_tolerance=10.0,
        min_peak_intensity=1000.0,
        enable_fitting=False,
        enable_plotting=False,
        enable_ms2=False,
    )

    rt = np.array([0.0, 1.0, 2.0], dtype=float)
    eic_hi = (rt, np.array([0.0, 2000.0, 0.0], dtype=float))
    eic_lo = (rt, np.array([0.0, 500.0, 0.0], dtype=float))
    reader: Any = _DummyReader("sample.raw", [eic_hi, eic_lo])

    targets_df = pd.DataFrame(
        [
            {
                "num": 1,
                "File name": "Library_POS_Mix121",
                "mixture": 121,
                "Compound name": "A",
                "Polarity": "POS",
                "m/z": 100.0,
            },
            {
                "num": 2,
                "File name": "Library_POS_Mix121",
                "mixture": 121,
                "Compound name": "B",
                "Polarity": "POS",
                "m/z": 200.0,
            },
            {
                "num": 3,
                "File name": "Library_POS_Mix121",
                "mixture": 121,
                "Compound name": "C",
                "Polarity": "POS",
                "m/z": None,
            },
        ]
    )

    status_rows: list[dict] = []
    results = process_raw_file_direct_mz(
        reader=reader,
        targets_df=targets_df,
        lc_mode="RP",
        polarity="POS",
        partial_filename="Library_POS_Mix121",
        file_suffix="",
        config=config,
        raw_file_id="RP/sample.raw",
        status_rows=status_rows,
    )

    assert len(results) == 1
    assert results[0]["Compound name"] == "A"
    assert results[0]["num"] == 1
    assert results[0]["EICGenerated"] is True
    assert results[0]["FilteredOut"] is False

    assert len(status_rows) == 3
    status_by_compound = {r["Compound name"]: r for r in status_rows}

    assert status_by_compound["A"]["EICGenerated"] is True
    assert status_by_compound["A"]["FilteredOut"] is False
    assert status_by_compound["A"]["num"] == 1

    assert status_by_compound["B"]["EICGenerated"] is True
    assert status_by_compound["B"]["FilteredOut"] is True
    assert status_by_compound["B"]["num"] == 2

    assert status_by_compound["C"]["EICGenerated"] is False
    assert status_by_compound["C"]["FilteredOut"] is False
    assert status_by_compound["C"]["num"] == 3


def test_write_results_excel_appends_status_rows_to_all_features(
    tmp_path: Path,
) -> None:
    out_path = tmp_path / "out.xlsx"

    results = [
        {
            "RawFile": "sample.raw",
            "num": 1,
            "Compound name": "A",
            "Polarity": "POS",
            "mz_target": 100.0,
            "Area": 1.0,
            "EICGenerated": True,
            "FilteredOut": False,
        }
    ]
    status_rows = [
        {**results[0]},
        {
            "RawFile": "sample.raw",
            "num": 2,
            "Compound name": "B",
            "Polarity": "POS",
            "mz_target": 200.0,
            "Area": None,
            "EICGenerated": True,
            "FilteredOut": True,
        },
    ]

    write_results_excel(
        results,
        str(out_path),
        include_pivot_tables=False,
        status_rows=status_rows,
    )

    xl = pd.ExcelFile(out_path)
    assert "All_Features" in xl.sheet_names
    assert "Target_Status" not in xl.sheet_names

    df_all = pd.read_excel(out_path, sheet_name="All_Features")
    assert "num" in df_all.columns
    assert df_all["num"].iloc[0] == 1
    assert {"EICGenerated", "FilteredOut"}.issubset(set(df_all.columns))
    assert len(df_all) == 2
