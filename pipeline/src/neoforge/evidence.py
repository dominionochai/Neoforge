"""Evidence-kind discipline: the single source of truth for what kind of
quantity any value in NEOFORGE is."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any


class EvidenceKind(str, Enum):
    MEASURED = "measured"
    PREDICTED = "predicted"
    DERIVED = "derived"
    MISSING = "missing"


@dataclass(frozen=True)
class FeatureValue:
    kind: EvidenceKind
    value: Any = None
    source: str | None = None
    note: str | None = None

    def as_dict(self) -> dict:
        return {
            "kind": self.kind.value,
            "value": self.value,
            "source": self.source,
            "note": self.note,
        }


def prediction(value: Any, source: str) -> FeatureValue:
    return FeatureValue(EvidenceKind.PREDICTED, value=value, source=source)


def measurement(value: Any, source: str) -> FeatureValue:
    return FeatureValue(EvidenceKind.MEASURED, value=value, source=source)


def derivation(value: Any, source: str) -> FeatureValue:
    return FeatureValue(EvidenceKind.DERIVED, value=value, source=source)


def missing(note: str) -> FeatureValue:
    return FeatureValue(EvidenceKind.MISSING, note=note)
