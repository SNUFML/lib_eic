"""Command-line interface for LCMS Adduct Finder."""

import argparse
import sys
from pathlib import Path
from typing import List, Optional

from . import __version__


def create_parser() -> argparse.ArgumentParser:
    """Create argument parser for CLI.

    Returns:
        Configured ArgumentParser.
    """
    parser = argparse.ArgumentParser(
        prog="lib_eic",
        description=(
            "Automated Targeted Feature Extraction & Adduct Verification Tool "
            "for LC-MS Data."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default settings
  lib_eic

  # Run with custom input/output files
  lib_eic --input compounds.xlsx --output results.xlsx

  # Run with custom raw data folder
  lib_eic --raw-folder /path/to/raw/files

  # Run with YAML config file
  lib_eic --config settings.yaml

  # Run with verbose output
  lib_eic -v

  # Generate default config file
  lib_eic --generate-config config.yaml
""",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    # Input/output options
    io_group = parser.add_argument_group("Input/Output")
    io_group.add_argument(
        "-i",
        "--input",
        dest="input_excel",
        metavar="FILE",
        help="Input Excel file with compound list (default: file_list.xlsx)",
    )
    io_group.add_argument(
        "-o",
        "--output",
        dest="output_excel",
        metavar="FILE",
        help="Output Excel file for results (default: Final_Result_With_Plots.xlsx)",
    )
    io_group.add_argument(
        "-r",
        "--raw-folder",
        dest="raw_data_folder",
        metavar="DIR",
        help="Folder containing .raw files (default: ./raw)",
    )
    io_group.add_argument(
        "-s",
        "--sheet",
        dest="input_sheet",
        metavar="NAME",
        help=(
            "Single sheet name in input Excel file (legacy; also sets --sheets). "
            "Default: Final"
        ),
    )
    io_group.add_argument(
        "--sheets",
        dest="input_sheets",
        metavar="NAME",
        nargs="+",
        help="Sheet names to read for direct m/z input (default: RP HILIC)",
    )

    # Config options
    config_group = parser.add_argument_group("Configuration")
    config_group.add_argument(
        "-c",
        "--config",
        dest="config_file",
        metavar="FILE",
        help="YAML configuration file",
    )
    config_group.add_argument(
        "--generate-config",
        dest="generate_config",
        metavar="FILE",
        help="Generate default config file and exit",
    )

    # Processing options
    proc_group = parser.add_argument_group("Processing")
    proc_group.add_argument(
        "--ppm",
        dest="ppm_tolerance",
        type=float,
        metavar="N",
        help="Mass tolerance in ppm (default: 10.0)",
    )
    proc_group.add_argument(
        "--min-intensity",
        dest="min_peak_intensity",
        type=float,
        metavar="N",
        help="Minimum peak intensity threshold (default: 100000)",
    )
    proc_group.add_argument(
        "--no-fitting",
        dest="disable_fitting",
        action="store_true",
        help="Disable Gaussian peak fitting",
    )
    proc_group.add_argument(
        "--no-plots",
        dest="disable_plotting",
        action="store_true",
        help="Disable EIC plot generation",
    )
    proc_group.add_argument(
        "--no-ms2",
        dest="disable_ms2",
        action="store_true",
        help="Disable MS2 indexing/matching",
    )
    proc_group.add_argument(
        "--area-method",
        dest="area_method",
        choices=["sum", "trapz"],
        help="Peak area calculation method (default: sum)",
    )

    # Logging options
    log_group = parser.add_argument_group("Logging")
    log_group.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose (DEBUG) logging",
    )
    log_group.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress non-error output",
    )
    log_group.add_argument(
        "--log-file",
        dest="log_file",
        metavar="FILE",
        help="Write log to file",
    )

    return parser


def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        args: List of arguments (defaults to sys.argv).

    Returns:
        Parsed arguments namespace.
    """
    parser = create_parser()
    return parser.parse_args(args)


def build_config_from_args(args: argparse.Namespace):
    """Build Config object from parsed arguments.

    Args:
        args: Parsed arguments namespace.

    Returns:
        Config instance with CLI overrides applied.
    """
    from .config import Config, load_config

    # Start with config file or defaults
    if args.config_file:
        config = load_config(args.config_file)
    else:
        config = Config()

    # Apply CLI overrides
    if args.input_excel:
        config.input_excel = args.input_excel
    if args.output_excel:
        config.output_excel = args.output_excel
    if args.raw_data_folder:
        config.raw_data_folder = Path(args.raw_data_folder)
    if args.input_sheet:
        config.input_sheet = args.input_sheet
        config.input_sheets = [args.input_sheet]
    if getattr(args, "input_sheets", None):
        config.input_sheets = list(args.input_sheets)
    if args.ppm_tolerance is not None:
        config.ppm_tolerance = args.ppm_tolerance
    if args.min_peak_intensity is not None:
        config.min_peak_intensity = args.min_peak_intensity
    if args.disable_fitting:
        config.enable_fitting = False
    if args.disable_plotting:
        config.enable_plotting = False
    if getattr(args, "disable_ms2", False):
        config.enable_ms2 = False
    if args.area_method:
        config.area_method = args.area_method
    if args.log_file:
        config.log_file = args.log_file

    # Set log level based on verbosity
    if args.verbose:
        config.log_level = "DEBUG"
    elif args.quiet:
        config.log_level = "ERROR"

    return config


def generate_default_config(output_path: str) -> None:
    """Generate default configuration file.

    Args:
        output_path: Path to write config file.
    """
    from .config import Config, save_config

    config = Config()
    save_config(config, output_path)
    print(f"Generated default config: {output_path}")


def main(args: Optional[List[str]] = None) -> int:
    """Main entry point for CLI.

    Args:
        args: Command-line arguments (defaults to sys.argv).

    Returns:
        Exit code (0 for success, non-zero for errors).
    """
    parsed_args = parse_args(args)

    # Handle --generate-config
    if parsed_args.generate_config:
        try:
            generate_default_config(parsed_args.generate_config)
            return 0
        except Exception as e:
            print(f"Error generating config: {e}", file=sys.stderr)
            return 1

    # Build config from args
    try:
        config = build_config_from_args(parsed_args)
    except Exception as e:
        print(f"Configuration error: {e}", file=sys.stderr)
        return 1

    # Setup logging
    from .logging_setup import setup_logging

    setup_logging(level=config.log_level, log_file=config.log_file)

    # Run processing
    from .processor import process_all

    try:
        process_all(config)
        return 0
    except KeyboardInterrupt:
        print("\nInterrupted by user", file=sys.stderr)
        return 130
    except Exception as e:
        import logging

        logger = logging.getLogger(__name__)
        logger.exception("Processing failed: %s", e)
        return 1


if __name__ == "__main__":
    sys.exit(main())
