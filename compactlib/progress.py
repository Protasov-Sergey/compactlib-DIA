from __future__ import annotations

import time
from contextlib import contextmanager
from typing import Iterator

try:  # optional nice terminal output
    from rich.console import Console
    from rich.status import Status
except Exception:  # pragma: no cover - fallback for minimal environments
    Console = None
    Status = None


class ProgressLogger:
    """Small logging/progress helper for CLI commands.

    The core library operations are mostly vectorized pandas transformations, so a
    row-level progress bar would be misleading and can slow processing. Instead,
    compactlib reports major stages with elapsed times and optionally displays a
    terminal spinner via rich.
    """

    def __init__(self, verbose: bool = True, progress: bool = True) -> None:
        self.verbose = verbose
        self.progress = progress
        self.console = Console(stderr=True) if Console is not None else None
        self.t0 = time.perf_counter()

    def log(self, message: str) -> None:
        if not self.verbose:
            return
        if self.console is not None:
            self.console.print(message)
        else:
            print(message)

    @contextmanager
    def step(self, message: str) -> Iterator[None]:
        """Log a timed processing step."""
        if not self.verbose:
            yield
            return

        t0 = time.perf_counter()
        status_obj = None
        if self.progress and self.console is not None and Status is not None:
            status_obj = self.console.status(f"[bold]{message}[/]", spinner="dots")
            status_obj.start()
        else:
            self.log(f"[compactlib] {message}...")

        try:
            yield
        finally:
            if status_obj is not None:
                status_obj.stop()
            dt = time.perf_counter() - t0
            self.log(f"[compactlib] done: {message} ({dt:.2f} s)")

    def total_elapsed(self) -> float:
        return time.perf_counter() - self.t0


def format_shape(df) -> str:
    try:
        return f"{len(df):,} rows × {len(df.columns):,} columns"
    except Exception:
        return "unknown shape"
