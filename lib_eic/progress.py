"""Progress bar helpers.

This module keeps ``tqdm`` usage optional and centralized.
"""

from __future__ import annotations

import sys
from contextlib import contextmanager
from typing import Any, Iterator


class _NullProgressBar:
    def update(self, _n: int = 1) -> None:
        return None

    def set_postfix_str(self, _s: str = "", refresh: bool = True) -> None:
        del refresh
        return None


def should_show_progress(enabled: bool) -> bool:
    """Return True if a progress bar should be shown for this run."""
    if not enabled:
        return False
    try:
        return sys.stderr.isatty()
    except Exception:
        # If isatty() is not available, default to enabled.
        return True


@contextmanager
def progress_bar(
    *,
    total: int,
    desc: str,
    enabled: bool,
    unit: str = "it",
) -> Iterator[Any]:
    """Context manager that yields a tqdm-like progress bar (or a no-op)."""
    if not enabled:
        yield _NullProgressBar()
        return

    try:
        from tqdm.auto import tqdm as _tqdm
    except Exception:
        yield _NullProgressBar()
        return

    tqdm: Any = _tqdm
    with tqdm(
        total=int(total),
        desc=str(desc),
        unit=str(unit),
        file=sys.stderr,
        dynamic_ncols=True,
        leave=True,
    ) as pbar:
        yield pbar
