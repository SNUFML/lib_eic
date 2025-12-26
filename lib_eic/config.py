"""Configuration management for LCMS Adduct Finder."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Any


@dataclass
class Config:
    """Configuration settings for LCMS Adduct Finder.

    All settings that were previously global constants are now encapsulated here.
    This allows for different configurations per run and eliminates global state.
    """

    # File paths
    raw_data_folder: Path = field(default_factory=lambda: Path("./raw"))
    input_excel: str = "file_list.xlsx"
    input_sheets: List[str] = field(default_factory=lambda: ["RP", "HILIC"])
    input_sheet: str = "Final"
    output_excel: str = "Final_Result_With_Plots.xlsx"
    export_plot_folder: str = "EIC_Plots_Export"
    include_pivot_tables: bool = False
    show_progress: bool = True

    # Mass tolerance
    ppm_tolerance: float = 10.0

    # Chromatogram extraction
    chromatogram_batch_size: int = 256

    # Peak filtering
    min_peak_intensity: float = 100000.0

    # Gaussian fitting controls
    enable_fitting: bool = True
    fit_rt_window_min: float = 0.30

    # Plotting controls
    enable_plotting: bool = True
    plot_dpi: int = 300

    # MS2 matching settings
    enable_ms2: bool = True
    ms2_match_mode: str = "rt_linked"  # "rt_linked" or "global"
    ms2_rt_window_min: float = 0.30
    store_ms2_match_details: bool = True

    # Area calculation method
    area_method: str = "sum"  # "sum" or "trapz"

    # Extension features (disabled by default)
    export_ms2_mgf: bool = False

    # Adduct names to include (if empty, use all from ADDUCT_DEFINITIONS)
    enabled_adducts: List[str] = field(default_factory=list)

    # Logging
    log_level: str = "INFO"
    log_file: Optional[str] = None

    # Parallelization
    # - num_workers <= 0: auto
    # - num_workers == 1: sequential
    # - num_workers > 1: that many worker processes
    num_workers: int = 0
    # "auto" (default), "sequential", "file", "task"
    parallel_mode: str = "auto"

    def __post_init__(self) -> None:
        """Convert string paths to Path objects if needed."""
        if isinstance(self.raw_data_folder, str):
            self.raw_data_folder = Path(self.raw_data_folder)

        # Normalize input_sheets
        input_sheets_value: Any = self.input_sheets
        if input_sheets_value is None:
            self.input_sheets = []
        elif isinstance(input_sheets_value, str):
            self.input_sheets = [s.strip() for s in str(input_sheets_value).split(",") if s.strip()]
        else:
            self.input_sheets = [str(s).strip() for s in input_sheets_value if str(s).strip()]

        # Validate ms2_match_mode
        if self.ms2_match_mode not in ("rt_linked", "global"):
            raise ValueError(f"Invalid ms2_match_mode: {self.ms2_match_mode!r}. " "Expected 'rt_linked' or 'global'.")

        # Validate area_method
        if self.area_method not in ("sum", "trapz"):
            raise ValueError(f"Invalid area_method: {self.area_method!r}. " "Expected 'sum' or 'trapz'.")

        try:
            chromatogram_batch_size = int(self.chromatogram_batch_size)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid chromatogram_batch_size: {self.chromatogram_batch_size!r}. Expected a positive integer."
            ) from exc
        if chromatogram_batch_size <= 0:
            raise ValueError(
                f"Invalid chromatogram_batch_size: {self.chromatogram_batch_size!r}. Expected a positive integer."
            )
        self.chromatogram_batch_size = chromatogram_batch_size

        # Normalize and validate parallel_mode
        parallel_mode_value: Any = self.parallel_mode
        if parallel_mode_value is None:
            self.parallel_mode = "auto"
        self.parallel_mode = str(self.parallel_mode).strip().lower()
        if self.parallel_mode not in ("auto", "sequential", "file", "task"):
            raise ValueError(
                f"Invalid parallel_mode: {self.parallel_mode!r}. " "Expected 'auto', 'sequential', 'file', or 'task'."
            )

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Config":
        """Create a Config instance from a dictionary."""
        # Filter out unknown keys
        valid_fields = {f.name for f in cls.__dataclass_fields__.values()}
        filtered_data = {k: v for k, v in data.items() if k in valid_fields}
        return cls(**filtered_data)

    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary."""
        return {
            "raw_data_folder": str(self.raw_data_folder),
            "input_excel": self.input_excel,
            "input_sheets": self.input_sheets,
            "input_sheet": self.input_sheet,
            "output_excel": self.output_excel,
            "export_plot_folder": self.export_plot_folder,
            "include_pivot_tables": bool(self.include_pivot_tables),
            "show_progress": bool(self.show_progress),
            "ppm_tolerance": self.ppm_tolerance,
            "chromatogram_batch_size": self.chromatogram_batch_size,
            "min_peak_intensity": self.min_peak_intensity,
            "enable_fitting": self.enable_fitting,
            "fit_rt_window_min": self.fit_rt_window_min,
            "enable_plotting": self.enable_plotting,
            "plot_dpi": self.plot_dpi,
            "enable_ms2": self.enable_ms2,
            "ms2_match_mode": self.ms2_match_mode,
            "ms2_rt_window_min": self.ms2_rt_window_min,
            "store_ms2_match_details": self.store_ms2_match_details,
            "area_method": self.area_method,
            "export_ms2_mgf": self.export_ms2_mgf,
            "enabled_adducts": self.enabled_adducts,
            "log_level": self.log_level,
            "log_file": self.log_file,
            "num_workers": self.num_workers,
            "parallel_mode": self.parallel_mode,
        }


def load_config(config_path: Optional[str] = None) -> Config:
    """Load configuration from a YAML file.

    Args:
        config_path: Path to YAML config file. If None, returns default config.

    Returns:
        Config instance with loaded settings.

    Raises:
        FileNotFoundError: If config file doesn't exist.
        ValueError: If config file has invalid format.
    """
    if config_path is None:
        return Config()

    config_file = Path(config_path)
    if not config_file.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    try:
        import yaml
    except ImportError:
        raise ImportError("PyYAML is required for config file support. " "Install with: pip install pyyaml")

    with open(config_file, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    if data is None:
        return Config()

    if not isinstance(data, dict):
        raise ValueError(f"Config file must contain a YAML dictionary, got: {type(data)}")

    return Config.from_dict(data)


def save_config(config: Config, config_path: str) -> None:
    """Save configuration to a YAML file.

    Args:
        config: Config instance to save.
        config_path: Path to save YAML config file.
    """
    try:
        import yaml
    except ImportError:
        raise ImportError("PyYAML is required for config file support. " "Install with: pip install pyyaml")

    config_file = Path(config_path)
    config_file.parent.mkdir(parents=True, exist_ok=True)

    with open(config_file, "w", encoding="utf-8") as f:
        yaml.dump(config.to_dict(), f, default_flow_style=False, sort_keys=False)
