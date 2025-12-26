from lib_eic.cli import build_config_from_args, parse_args


def test_cli_no_ms2_sets_config_flag() -> None:
    args = parse_args(["--no-ms2"])
    config = build_config_from_args(args)
    assert config.enable_ms2 is False
