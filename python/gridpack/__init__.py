# -------------------------------------------------------------
# file: gridpack/__init__.py
# -------------------------------------------------------------
# GridPACK Python package.
#
# The compiled pybind11 bindings live in gridpack._gridpack (private).
# This module re-exports them so that:
#
#   import gridpack
#   gridpack.Environment()          # works
#   gridpack.hadrec.Module()        # works
#   import gridpack.hadrec          # works (bare-name import statement)
#
# High-level API classes (Session, PowerFlow, DynamicSim, ...) are added
# alongside as pure-Python wrappers.
# -------------------------------------------------------------

from ._gridpack import *          # noqa: F401,F403  re-export top-level symbols
from . import _gridpack as _ext

import sys as _sys

# Republish each pybind11 def_submodule as both an attribute on this
# package (so `gridpack.hadrec.Module()` works without an explicit
# import) and a sys.modules entry (so `import gridpack.hadrec` as a
# statement also works). Do not silently skip a missing submodule --
# the compiled extension is expected to publish all five.
for _name in (
    "hadrec",
    "dynamic_simulation",
    "powerflow",
    "state_estimation",
    "emt",
):
    _sub = getattr(_ext, _name)
    globals()[_name] = _sub
    _sys.modules[__name__ + "." + _name] = _sub

del _name, _sub, _sys

# High-level Python API.  Low-level pybind11 access continues to work
# through the re-exports above.
from .session import Session  # noqa: E402

__all__ = ["Session"]
