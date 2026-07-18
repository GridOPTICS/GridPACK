# -------------------------------------------------------------
# file: gridpack/state_estimation.py
# -------------------------------------------------------------
# Re-export the compiled pybind11 state_estimation submodule.
# High-level StateEstimation wrapper is added here in a later commit.
# -------------------------------------------------------------

from ._gridpack.state_estimation import *  # noqa: F401,F403
from ._gridpack import state_estimation as _mod  # noqa: F401


def __getattr__(name):
    return getattr(_mod, name)
