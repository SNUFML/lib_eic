from lib_eic.io.plotting import build_direct_mz_plot_filename


def test_build_direct_mz_plot_filename_no_num_prefix() -> None:
    name = build_direct_mz_plot_filename(
        compound_name="Spermine",
        polarity="POS",
        mixture="121",
    )

    assert name == "Spermine_POS_121.png"


def test_build_direct_mz_plot_filename_with_num_prefix_and_suffix() -> None:
    name = build_direct_mz_plot_filename(
        compound_name="1-Methylhistamine (dihydrochloride)",
        polarity="NEG",
        mixture="1",
        file_suffix="_2nd",
        num_prefix="003",
    )

    assert name == "003_1-Methylhistamine_dihydrochloride_NEG_1_2nd.png"


def test_build_direct_mz_plot_filename_ignores_nan_num_prefix() -> None:
    name = build_direct_mz_plot_filename(
        compound_name="Spermine",
        polarity="POS",
        mixture="121",
        num_prefix="NaN",
    )

    assert name == "Spermine_POS_121.png"
