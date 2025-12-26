from __future__ import annotations

from typing import Iterable, Iterator, List

from lib_eic.analysis.eic import build_targets


class _OneShotIterable(Iterable[str]):
    def __init__(self, items: List[str]) -> None:
        self._items = items
        self._iterated = False

    def __iter__(self) -> Iterator[str]:
        if self._iterated:
            raise RuntimeError("Iterable was iterated more than once")
        self._iterated = True
        return iter(self._items)


def test_build_targets_does_not_iterate_formulas_twice() -> None:
    formulas = _OneShotIterable(["H2O", "C2H6O"])

    targets = build_targets(formulas, mode="POS")

    assert targets
