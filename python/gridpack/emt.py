# -------------------------------------------------------------
# file: gridpack/emt.py
# -------------------------------------------------------------
# Re-export the compiled pybind11 emt submodule.
# High-level EMT wrapper is added here in a later commit.
# -------------------------------------------------------------

from ._gridpack.emt import *  # noqa: F401,F403
from ._gridpack import emt as _mod  # noqa: F401


def __getattr__(name):
    return getattr(_mod, name)
