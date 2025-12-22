"""Logging configuration for LCMS Adduct Finder."""

import logging
import sys
from typing import Optional


class TqdmLoggingHandler(logging.Handler):
    """Logging handler that plays nicely with tqdm progress bars."""

    def __init__(self, stream=None) -> None:
        super().__init__()
        self.stream = stream if stream is not None else sys.stderr

    def emit(self, record: logging.LogRecord) -> None:
        try:
            msg = self.format(record)
            try:
                from tqdm.auto import tqdm
            except Exception:
                stream = self.stream
                stream.write(msg + "\n")
                stream.flush()
                return

            tqdm.write(msg, file=self.stream)
        except Exception:
            self.handleError(record)


def setup_logging(
    level: str = "INFO",
    log_file: Optional[str] = None,
    name: str = "root",
    *,
    use_tqdm: bool = False,
) -> logging.Logger:
    """Configure logging for the application.

    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
        log_file: Optional path to log file. If None, logs only to console.
        name: Ignored (kept for backward compatibility).
        use_tqdm: Use a tqdm-aware console handler (recommended when showing a
            progress bar).

    Returns:
        Configured logger instance.
    """
    # Convert string level to logging constant
    numeric_level = getattr(logging, level.upper(), logging.INFO)

    # Configure root logger so module loggers (lib_eic.*) inherit handlers.
    logger = logging.getLogger()
    logger.setLevel(numeric_level)

    # Remove existing handlers to avoid duplicates
    logger.handlers.clear()

    # Create formatter
    formatter = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console handler
    console_handler: logging.Handler
    if use_tqdm:
        console_handler = TqdmLoggingHandler(stream=sys.stderr)
    else:
        console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(numeric_level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # File handler (optional)
    if log_file:
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(numeric_level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def get_logger(name: str = "lcms_adduct_finder") -> logging.Logger:
    """Get a logger instance.

    Args:
        name: Logger name. Use __name__ for module-specific loggers.

    Returns:
        Logger instance.
    """
    return logging.getLogger(name)
