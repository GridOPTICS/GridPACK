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
# -------------------------------------------------------------

from . import _gridpack           # noqa: F401  keep private extension reachable

# Bound one by one, NOT via `from ._gridpack import *`: the wildcard also
# bound the pybind11 submodules as attributes of this package, and CPython
# then skipped importing the same-named Python shims below (`from . import X`
# is a no-op when X is already an attribute).  dynamic_simulation and emt
# were the two nothing else imported *from*, so they stayed the raw compiled
# modules and gridpack.dynamic_simulation.DynamicSim raised AttributeError.
from ._gridpack import (          # noqa: F401
    CoarseTimer,
    Communicator,
    Configuration,
    ConfigurationCursor,
    Environment,
    NoPrint,
    TaskCounter,
    TaskManager,
)

# Eager so that `gridpack.hadrec.Module` works after just `import gridpack`.
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
from .state_estimation import (                            # noqa: E402
    StateEstimation,
    StateEstimationResult,
)
from .hadrec import Hadrec                                 # noqa: E402
from .contingency import (                                  # noqa: E402
    Contingency,
    ContingencyAnalysis,
    ContingencyResult,
)
from .exceptions import GridPACKError, PowerFlowDiverged   # noqa: E402

#: Compat shims mirroring the pybind11 submodules.  Public, but not in
#: __all__ -- `from gridpack import *` should not pull in module objects.
COMPAT_MODULES = (
    "hadrec",
    "dynamic_simulation",
    "powerflow",
    "state_estimation",
    "emt",
)

_HIGH_LEVEL = [
    "Session",
    "PowerFlow",
    "PowerFlowResult",
    "DynamicSim",
    "DynamicSimStepper",
    "DSFResult",
    "StateEstimation",
    "StateEstimationResult",
    "Hadrec",
    "Contingency",
    "ContingencyAnalysis",
    "ContingencyResult",
    "GridPACKError",
    "PowerFlowDiverged",
]

# Thin pybind11 wrappers over the C++ objects.  Public because the scripts
# under python/src drive GridPACK entirely through them, but prefer Session
# and the high-level classes: these do no lifetime management of their own.
_LOW_LEVEL = [
    "CoarseTimer",
    "Communicator",
    "Configuration",
    "ConfigurationCursor",
    "Environment",
    "NoPrint",
    "TaskCounter",
    "TaskManager",
]

__all__ = _HIGH_LEVEL + _LOW_LEVEL


def __dir__():
    """The supported surface only; implementation modules stay out of it."""
    return sorted(__all__ + list(COMPAT_MODULES))


def __getattr__(name):
    # Everything supported is bound above, so reaching here is a typo or a
    # private name; point at the curated surface rather than at dir().
    raise AttributeError(
        "module 'gridpack' has no attribute %r; the public API is %s "
        "plus the submodules %s"
        % (name, ", ".join(sorted(__all__)), ", ".join(COMPAT_MODULES))
    )
