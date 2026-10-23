# -------------------------------------------------------------
# file: gridpack/dynamic_simulation.py
# -------------------------------------------------------------
# Re-export the compiled pybind11 dynamic_simulation submodule.
# High-level DynamicSim wrapper lives in gridpack/dynamic_sim.py to keep
# the compat namespace clean.
# -------------------------------------------------------------

from ._gridpack.dynamic_simulation import *  # noqa: F401,F403
from ._gridpack import dynamic_simulation as _mod  # noqa: F401

# High-level Python classes live in gridpack.dynamic_sim; re-export them
# here so `from gridpack.dynamic_simulation import DynamicSim` works too.
from .dynamic_sim import (  # noqa: F401,E402
    DynamicSim,
    DynamicSimStepper,
    DSFResult,
)


def __getattr__(name):
    return getattr(_mod, name)
