from pathlib import Path

from lib_eic.processor import (
    _build_raw_file_id,
    _collect_raw_entries_recursive,
    _infer_run_label,
    _resolve_lc_mode_raw_folder,
    find_matching_raw_files,
)


def test_collect_raw_entries_recursive_prunes_raw_dirs(tmp_path: Path) -> None:
    raw_root = tmp_path / "raw"

    (raw_root / "HILIC" / "1st").mkdir(parents=True)
    raw_file = raw_root / "HILIC" / "1st" / "Library_HILIC_NEG_Mixture_1.raw"
    raw_file.write_text("x", encoding="utf-8")

    (raw_root / "RP" / "1st").mkdir(parents=True)
    raw_dir = raw_root / "RP" / "1st" / "Bundle.raw"
    raw_dir.mkdir()
    inner = raw_dir / "inner.raw"
    inner.write_text("x", encoding="utf-8")

    entries = _collect_raw_entries_recursive(raw_root)

    assert raw_file in entries
    assert raw_dir in entries
    assert inner not in entries


def test_resolve_lc_mode_raw_folder_case_insensitive(tmp_path: Path) -> None:
    raw_root = tmp_path / "Raw"
    (raw_root / "HiLiC").mkdir(parents=True)

    assert _resolve_lc_mode_raw_folder(raw_root, "HILIC") == raw_root / "HiLiC"


def test_infer_run_label_from_folder(tmp_path: Path) -> None:
    search_root = tmp_path / "HILIC"
    (search_root / "1st").mkdir(parents=True)
    raw_path = search_root / "1st" / "sample.raw"
    raw_path.write_text("x", encoding="utf-8")

    assert _infer_run_label(raw_path, search_root) == "1st"


def test_build_raw_file_id_uses_relative_posix_path(tmp_path: Path) -> None:
    raw_root = tmp_path / "Raw"
    raw_path = raw_root / "HILIC" / "1st" / "sample.raw"
    raw_path.parent.mkdir(parents=True)
    raw_path.write_text("x", encoding="utf-8")

    assert _build_raw_file_id(raw_path, raw_root) == "HILIC/1st/sample.raw"


def test_find_matching_raw_files_avoids_numeric_prefix_overmatch(
    tmp_path: Path,
) -> None:
    root = tmp_path / "HILIC"
    root.mkdir()

    f1 = root / "Library_HILIC_NEG_Mixture_1.raw"
    f2 = root / "Library_HILIC_NEG_Mixture_1_2nd.raw"
    f10 = root / "Library_HILIC_NEG_Mixture_10.raw"
    for path in (f1, f2, f10):
        path.write_text("x", encoding="utf-8")

    matches = find_matching_raw_files(
        "Library_HILIC_NEG_Mixture_1",
        root,
        raw_entries=[f1, f2, f10],
    )

    assert matches == [f1, f2]
