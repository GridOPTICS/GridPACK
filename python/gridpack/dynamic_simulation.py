# -------------------------------------------------------------
# file: gridpack/dynamic_simulation.py
# -------------------------------------------------------------
# Re-export the compiled pybind11 dynamic_simulation submodule.
# High-level DynamicSim wrapper lives in gridpack/dynamic_sim.py to keep
# the compat namespace clean.
# -------------------------------------------------------------

from ._gridpack.dynamic_simulation import *  # noqa: F401,F403
from ._gridpack import dynamic_simulation as _mod  # noqa: F401


def __getattr__(name):
    return getattr(_mod, name)
