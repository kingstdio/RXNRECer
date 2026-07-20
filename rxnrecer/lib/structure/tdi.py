"""Utilities for working with 3Di token files."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class Tdi:
    """Container for a single 3Di record."""

    name: str = ""
    seq: str = ""
    token_3di: str = ""
    matrix_3di: np.ndarray = field(default_factory=lambda: np.array([]))

    def read_3di_file(self, path: str | Path) -> None:
        """Load a tab-separated 3Di record from disk."""
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"3Di file does not exist: {path}")

        content = path.read_text(encoding="utf-8").strip().split("\t")
        if len(content) < 4:
            raise ValueError(f"Expected 4 tab-separated fields, got {len(content)}: {content}")

        self.name = content[0]
        self.seq = content[1]
        self.token_3di = content[2]
        matrix = np.array(content[3].split(","), dtype=float)
        self.matrix_3di = matrix.reshape(-1, 10)

    def show(self) -> None:
        """Print the record for quick debugging."""
        print(f"name: {self.name}")
        print(f"seq: {self.seq}")
        print(f"token_3di: {self.token_3di}")
        print(f"matrix_3di:\n{self.matrix_3di}")
