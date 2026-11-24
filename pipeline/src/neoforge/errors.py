"""Typed errors so failures surface with their reason, never as silent drops."""
from __future__ import annotations


class NeoforgeError(Exception):
    """Root of all deliberately raised pipeline errors."""


class InputIntegrityError(NeoforgeError):
    """A reference/input file failed its pinned SHA-256 check or was unreadable."""


class CaseBundleError(NeoforgeError):
    """A case bundle is malformed, incomplete, or internally inconsistent."""


class ArtifactChainError(NeoforgeError):
    """A stage was asked to run on artifacts a prerequisite stage never wrote."""
