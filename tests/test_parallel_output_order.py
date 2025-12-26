from __future__ import annotations

import time
from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Sequence

import pandas as pd

from lib_eic.config import Config
from lib_eic import processor


@contextmanager
def _thread_pool(
    *,
    max_workers: Optional[int],
    log_level: str = "INFO",
    log_handlers: Optional[Sequence[Any]] = None,
) -> Iterator[ThreadPoolExecutor]:
    del log_level, log_handlers
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        yield executor


def test_process_all_formula_based_preserves_work_item_order(monkeypatch, tmp_path: Path) -> None:
    raw_root = tmp_path / "raw"
    raw_root.mkdir()
    (raw_root / "A.raw").mkdir()
    (raw_root / "B.raw").mkdir()

    input_xlsx = tmp_path / "input.xlsx"
    input_xlsx.write_text("x", encoding="utf-8")

    meta_data = pd.DataFrame(
        {
            "RawFile": ["A", "B"],
            "Mode": ["POS", "POS"],
            "Formula": ["H2O", "H2O"],
        }
    )

    def fake_read_input_excel(*_args: Any, **_kwargs: Any) -> pd.DataFrame:
        return meta_data

    def fake_worker(
        *,
        raw_file_path: str,
        formulas: List[str],
        mode: str,
        config_dict: Dict[str, Any],
        raw_file_id: str,
    ) -> Dict[str, Any]:
        del raw_file_path, mode, config_dict
        if raw_file_id.startswith("A"):
            time.sleep(0.05)
        row = {
            "RawFile": raw_file_id,
            "Formula": formulas[0],
            "Adduct": "M+H",
            "RT_min": 1.0,
            "Intensity": 1.0,
            "Area": 1.0,
        }
        status_row = {**row, "EICGenerated": True, "FilteredOut": False}
        return {
            "results": [row],
            "status_rows": [status_row],
            "raw_file_id": raw_file_id,
        }

    captured: Dict[str, Any] = {}

    def fake_write_results_excel(
        results: List[Dict],
        output_path: str,
        include_pivot_tables: bool = True,
        *,
        status_rows: Optional[List[Dict]] = None,
    ) -> None:
        del output_path, include_pivot_tables
        captured["results"] = list(results)
        captured["status_rows"] = list(status_rows or [])

    monkeypatch.setattr(processor, "read_input_excel", fake_read_input_excel)
    monkeypatch.setattr(processor, "_process_single_file_formula_worker", fake_worker)
    monkeypatch.setattr(processor, "create_process_pool", _thread_pool)
    monkeypatch.setattr(processor, "write_results_excel", fake_write_results_excel)

    config = Config(
        raw_data_folder=raw_root,
        input_excel=str(input_xlsx),
        output_excel=str(tmp_path / "out.xlsx"),
        parallel_mode="file",
        num_workers=2,
        enable_ms2=False,
        enable_plotting=False,
    )

    processor.process_all_formula_based(config)

    assert [r["RawFile"] for r in captured["results"]] == ["A.raw", "B.raw"]
    assert [r["RawFile"] for r in captured["status_rows"]] == ["A.raw", "B.raw"]


def test_process_all_direct_mz_preserves_work_item_order(monkeypatch, tmp_path: Path) -> None:
    raw_root = tmp_path / "raw"
    raw_rp = raw_root / "RP"
    raw_rp.mkdir(parents=True)
    (raw_rp / "Mix1.raw").mkdir()
    (raw_rp / "Mix2.raw").mkdir()

    input_xlsx = tmp_path / "input.xlsx"
    input_xlsx.write_text("x", encoding="utf-8")

    rp_data = pd.DataFrame(
        [
            {
                "File name": "Mix1",
                "mixture": 1,
                "Compound name": "A",
                "Polarity": "POS",
                "m/z": 100.0,
            },
            {
                "File name": "Mix2",
                "mixture": 2,
                "Compound name": "B",
                "Polarity": "POS",
                "m/z": 200.0,
            },
        ]
    )

    def fake_read_all_lc_mode_sheets(*_args: Any, **_kwargs: Any) -> Dict[str, pd.DataFrame]:
        return {"RP": rp_data}

    def fake_worker(
        *,
        raw_file_path: str,
        targets_records: List[Dict[str, Any]],
        config_dict: Dict[str, Any],
        lc_mode: str,
        polarity: str,
        partial_filename: str,
        file_suffix: str,
        raw_file_id: str,
    ) -> Dict[str, Any]:
        del (
            raw_file_path,
            targets_records,
            config_dict,
            lc_mode,
            polarity,
            partial_filename,
            file_suffix,
        )
        if raw_file_id.endswith("Mix1.raw"):
            time.sleep(0.05)
        row = {
            "RawFile": raw_file_id,
            "Compound name": "X",
            "Polarity": "POS",
            "mz_target": 123.0,
            "RT_min": 1.0,
            "Intensity": 1.0,
            "Area": 1.0,
        }
        status_row = {**row, "EICGenerated": True, "FilteredOut": False}
        return {
            "results": [row],
            "status_rows": [status_row],
            "raw_file_id": raw_file_id,
        }

    captured: Dict[str, Any] = {}

    def fake_write_results_excel(
        results: List[Dict],
        output_path: str,
        include_pivot_tables: bool = True,
        *,
        status_rows: Optional[List[Dict]] = None,
    ) -> None:
        del output_path, include_pivot_tables
        captured["results"] = list(results)
        captured["status_rows"] = list(status_rows or [])

    monkeypatch.setattr(processor, "read_all_lc_mode_sheets", fake_read_all_lc_mode_sheets)
    monkeypatch.setattr(processor, "_process_single_file_direct_mz_worker", fake_worker)
    monkeypatch.setattr(processor, "create_process_pool", _thread_pool)
    monkeypatch.setattr(processor, "write_results_excel", fake_write_results_excel)

    config = Config(
        raw_data_folder=raw_root,
        input_excel=str(input_xlsx),
        input_sheets=["RP"],
        output_excel=str(tmp_path / "out.xlsx"),
        parallel_mode="file",
        num_workers=2,
        enable_ms2=False,
        enable_plotting=False,
    )

    processor.process_all_direct_mz(config)

    assert [r["RawFile"] for r in captured["results"]] == [
        "RP/Mix1.raw",
        "RP/Mix2.raw",
    ]
    assert [r["RawFile"] for r in captured["status_rows"]] == [
        "RP/Mix1.raw",
        "RP/Mix2.raw",
    ]
