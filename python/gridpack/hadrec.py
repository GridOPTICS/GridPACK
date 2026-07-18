# -------------------------------------------------------------
# file: gridpack/hadrec.py
# -------------------------------------------------------------
# Re-export the compiled pybind11 hadrec submodule.  High-level HADREC
# wrapper classes are added here (see class Hadrec).
# -------------------------------------------------------------

from ._gridpack.hadrec import *  # noqa: F401,F403
from ._gridpack import hadrec as _mod  # noqa: F401


def __getattr__(name):
    # Any pybind11 attribute not explicitly re-exported is still
    # accessible: `gridpack.hadrec.<anything>` falls back to the
    # compiled submodule.
    return getattr(_mod, name)
