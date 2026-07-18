# -------------------------------------------------------------
# file: gridpack/__init__.py
# -------------------------------------------------------------
# GridPACK Python package.
#
# The compiled pybind11 bindings live in gridpack._gridpack (private).
# Backwards-compat surface: existing scripts using
#
#   import gridpack
#   gridpack.Environment()
#   import gridpack.hadrec
#   gridpack.hadrec.Module()
#   gridpack.dynamic_simulation.Event()          # bare attribute access
#
# all continue to work.  Each pybind11 submodule has a matching Python
# module (gridpack/hadrec.py, gridpack/powerflow.py, ...) that re-exports
# the compiled symbols and hosts new high-level classes when we add them.
# The Python modules are imported eagerly so attribute access without an
# explicit `import` statement still resolves.
# -------------------------------------------------------------

from ._gridpack import *          # noqa: F401,F403  re-export top-level symbols
from . import _gridpack           # noqa: F401       keep private extension reachable

# Eager import of each subpackage so that:
#   * `gridpack.hadrec.Module` works after just `import gridpack`
#   * `import gridpack.hadrec` (statement form) also works
# Each of these modules re-exports its pybind11 def_submodule.
from . import hadrec              # noqa: F401,E402
from . import dynamic_simulation  # noqa: F401,E402
from . import powerflow           # noqa: F401,E402
from . import state_estimation    # noqa: F401,E402
from . import emt                 # noqa: F401,E402

# High-level Python API.
from .session import Session                              # noqa: E402
from .powerflow import PowerFlow                          # noqa: E402
from .dynamic_sim import DynamicSim, DynamicSimStepper    # noqa: E402
from .results import PowerFlowResult, DSFResult           # noqa: E402

__all__ = [
    "Session",
    "PowerFlow",
    "PowerFlowResult",
    "DynamicSim",
    "DynamicSimStepper",
    "DSFResult",
]
