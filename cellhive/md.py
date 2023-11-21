"""Helper module for easy import of metadata related modules."""

import importlib

__version__ = importlib.metadata.version("cellhive")

# not importing this in __init__.py to not slow down cli functions if not required.

from .metadata_tools import (  # pylint: disable=wrong-import-position,unused-import
    check,
    check_2,
    layers,
    md,
    obs,
    obsm,
)
