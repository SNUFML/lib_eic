"""Multiprocessing helpers for lib_eic.

This module centralizes cross-platform (Windows/Linux) multiprocessing behavior.
We always use the ``spawn`` start method for consistency and safety.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from logging.handlers import QueueHandler, QueueListener
from typing import Iterable, Iterator, Optional, Sequence


def _shutdown_executor(executor: ProcessPoolExecutor, *, wait: bool, cancel_futures: bool) -> None:
    """Shutdown an executor with best-effort Python version compatibility."""
    try:
        executor.shutdown(wait=wait, cancel_futures=cancel_futures)
    except TypeError:
        # Python < 3.9 does not support cancel_futures.
        executor.shutdown(wait=wait)


def _iter_executor_processes(executor: ProcessPoolExecutor):
    processes = getattr(executor, "_processes", None)
    if isinstance(processes, dict):
        return [p for p in processes.values() if p is not None]
    if isinstance(processes, (list, tuple, set)):
        return [p for p in processes if p is not None]
    return []


def _terminate_processes(processes) -> None:
    """Best-effort terminate worker processes to avoid Ctrl+C hangs.

    Notes:
        ``ProcessPoolExecutor`` may block on exit waiting for workers to finish a
        long-running C/.NET call. On interrupts, we forcefully terminate workers
        so the parent can exit promptly.

        We prefer to let the executor manager thread do the actual joins to
        avoid double-joining processes from multiple threads.
    """
    proc_list = [p for p in (processes or []) if p is not None]

    for proc in proc_list:
        try:
            if getattr(proc, "is_alive", None) and proc.is_alive():
                proc.terminate()
        except Exception:
            pass

    for proc in proc_list:
        try:
            if getattr(proc, "is_alive", None) and proc.is_alive():
                kill = getattr(proc, "kill", None)
                if callable(kill):
                    kill()
        except Exception:
            pass


def resolve_max_workers(num_workers: int, *, num_items: int) -> Optional[int]:
    """Resolve an executor ``max_workers`` value from config.

    Rules:
    - ``num_workers <= 0``: auto (cap by ``num_items``)
    - ``num_workers == 1``: sequential (caller should bypass executor)
    - ``num_workers > 1``: use that value (cap by ``num_items``)
    """
    if num_items <= 0:
        return 1

    if num_workers <= 0:
        return max(1, min(mp.cpu_count(), num_items))

    return max(1, min(int(num_workers), num_items))


def should_use_process_pool(*, parallel_mode: str, num_workers: int, num_items: int) -> bool:
    """Return True if we should use file-level multiprocessing."""
    mode = str(parallel_mode or "").strip().lower()
    if mode in {"sequential", "off", "none"}:
        return False
    if mode == "task":
        logging.getLogger(__name__).warning(
            "parallel_mode='task' is not implemented; using file-level multiprocessing."
        )
        mode = "file"
    if num_workers == 1:
        return False
    return num_items > 1 and mode in {"auto", "file"}


def _init_worker(log_queue, log_level: str) -> None:
    """Initialize worker process: logging and expensive imports.

    This function runs once per worker when the pool is created.
    Pre-importing fisher_py triggers pythonnet/.NET CLR initialization,
    avoiding this overhead on the first file processed by each worker.
    """
    # Let the parent process handle Ctrl+C. This avoids workers raising
    # KeyboardInterrupt while inside long-running C/.NET calls.
    try:
        import signal

        signal.signal(signal.SIGINT, signal.SIG_IGN)
    except Exception:
        pass

    # Configure logging first
    if log_queue is not None:
        numeric_level = getattr(logging, str(log_level).upper(), logging.INFO)
        root = logging.getLogger()
        root.handlers.clear()
        root.setLevel(numeric_level)
        root.addHandler(QueueHandler(log_queue))

    # Pre-import fisher_py modules to initialize pythonnet/.NET CLR once per worker.
    # Without this, CLR initialization happens on first RawFileReader use,
    # adding 1-3 seconds overhead per file instead of per worker.
    #
    # We import all commonly used modules to ensure the CLR loads the required
    # .NET assemblies upfront, and access a class attribute to force full
    # initialization (pythonnet can be lazy).
    try:
        from fisher_py.raw_file_reader import RawFileReaderAdapter
        from fisher_py.data import Device  # noqa: F401
        from fisher_py.data.business import ChromatogramTraceSettings  # noqa: F401
        from fisher_py.data.business import TraceType  # noqa: F401
        from fisher_py.data.business.mass_options import MassOptions  # noqa: F401

        # Force CLR initialization by accessing class - pythonnet may defer
        # loading .NET types until first access
        _ = RawFileReaderAdapter.__name__
    except ImportError:
        pass  # fisher_py not installed; workers will fail later with clear error


def _noop() -> None:
    """No-op function used to pre-warm workers."""
    pass


def _prewarm_workers(executor: ProcessPoolExecutor, max_workers: int) -> None:
    """Force all workers to start immediately by submitting dummy tasks.

    ProcessPoolExecutor creates workers lazily. Without pre-warming, the first
    N real tasks would each wait for worker creation + CLR initialization.
    By submitting dummy tasks first, we parallelize worker creation.
    """
    if max_workers <= 0:
        return

    # Submit dummy tasks to force worker creation
    futures = [executor.submit(_noop) for _ in range(max_workers)]

    # Wait for all workers to be ready
    for f in futures:
        f.result()


@contextmanager
def create_process_pool(
    *,
    max_workers: Optional[int],
    log_level: str = "INFO",
    log_handlers: Optional[Sequence[logging.Handler]] = None,
) -> Iterator[ProcessPoolExecutor]:
    """Create a ``ProcessPoolExecutor`` using the ``spawn`` start method.

    If ``log_handlers`` are provided, worker logs are forwarded to the parent
    process via a queue listener using those handlers.

    Workers are pre-warmed immediately after pool creation to ensure all
    CLR initialization happens in parallel before any real work begins.

    On ``KeyboardInterrupt`` (Ctrl+C), workers are terminated to avoid hangs
    while waiting for long-running C/.NET calls to return.
    """
    ctx = mp.get_context("spawn")

    listener = None
    log_queue = None

    handlers = list(log_handlers or [])
    if handlers:
        log_queue = ctx.Queue()
        listener = QueueListener(log_queue, *handlers, respect_handler_level=True)
        listener.start()

    # Always use _init_worker to pre-import fisher_py and initialize CLR.
    # log_queue may be None if no handlers; _init_worker handles this.
    initializer = _init_worker
    initargs = (log_queue, log_level)

    executor: Optional[ProcessPoolExecutor] = None
    aborted = False
    try:
        executor = ProcessPoolExecutor(
            max_workers=max_workers,
            mp_context=ctx,
            initializer=initializer,
            initargs=initargs,
        )

        # Pre-warm workers: force all workers to start and run their initializer
        # NOW, in parallel, before any real work is submitted.
        # Without this, workers are created lazily (one at a time as tasks are
        # submitted), serializing the CLR initialization overhead.
        prewarm_n = max_workers
        if prewarm_n is None:
            prewarm_n = getattr(executor, "_max_workers", 0)
        _prewarm_workers(executor, int(prewarm_n or 0))

        yield executor
    except BaseException:
        aborted = True
        if executor is not None:
            processes = _iter_executor_processes(executor)
            _terminate_processes(processes)
            try:
                _shutdown_executor(executor, wait=True, cancel_futures=True)
            except Exception:
                pass
        raise
    finally:
        if executor is not None and not aborted:
            _shutdown_executor(executor, wait=True, cancel_futures=False)
        if listener is not None:
            listener.stop()
        if log_queue is not None:
            try:
                log_queue.close()
            except Exception:
                pass
            try:
                log_queue.join_thread()
            except Exception:
                pass


def iter_log_handlers() -> Iterable[logging.Handler]:
    """Return the currently configured root logger handlers."""
    return list(logging.getLogger().handlers)
