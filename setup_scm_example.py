"""Compatibility placeholder kept for public release layout stability.

The current package version is managed through ``pyproject.toml`` and
``rxnrecer.__version__``. This file remains in the release tree so existing
users of older release layouts do not see an unnecessary structural change.
"""

from rxnrecer import __version__


def get_release_version() -> str:
    return f"RXNRECer v{__version__}"
