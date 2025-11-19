"""Test fixtures that stub external dependencies."""

import sys
from types import SimpleNamespace


class _DummyFoldCompound:
    def __init__(self, sequence: str) -> None:
        self.sequence = sequence

    def pf(self):
        return None, -5.0

    def bpp(self):
        # Use 1-based indices that match the boundary positions used in the tests.
        return [
            (8, 50, 0.3),
            (9, 51, 0.2),
            (1, 10, 0.1),
        ]


def _fold_compound(sequence: str) -> _DummyFoldCompound:
    return _DummyFoldCompound(sequence)


sys.modules["RNA"] = SimpleNamespace(fold_compound=_fold_compound)
