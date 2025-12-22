from pathlib import Path

from openpyxl import Workbook

from lib_eic.io.excel import read_input_excel_direct_mz


def test_read_input_excel_direct_mz_keeps_optional_num_column(tmp_path: Path) -> None:
    path = tmp_path / "direct_mz.xlsx"

    wb = Workbook()
    ws = wb.active
    ws.title = "RP"

    # The reader uses skiprows=1 because some templates have merged cells in row 1.
    ws.append(["merged header row"])
    ws.append(["Num", "File name", "mixture", "Compound name", "Polarity", "m/z", "extra"])
    ws.append([1, "Library_POS_Mix121", 121, "Spermine", "POS", 203.223, "x"])
    ws.append([2, "Library_POS_Mix121", 121, "Putrescine", "POS", 89.107, "y"])

    wb.save(path)

    df = read_input_excel_direct_mz(str(path), sheet_name="RP")

    assert list(df.columns) == [
        "num",
        "File name",
        "mixture",
        "Compound name",
        "Polarity",
        "m/z",
    ]
    assert df["num"].tolist() == [1, 2]
