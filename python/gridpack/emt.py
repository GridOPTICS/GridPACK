# -------------------------------------------------------------
# file: gridpack/emt.py
# -------------------------------------------------------------
# Re-export the compiled pybind11 emt submodule.
# High-level EMT wrapper is added here in a later commit.
# -------------------------------------------------------------

import warnings

from ._gridpack import emt as _mod  # noqa: F401

__all__ = ["EMT"]

_EXPERIMENTAL = (
    "gridpack.emt is experimental: EMT is the one application with no Session "
    "integration, so it is outside the close()-ordering guarantees the others "
    "have, and its API may change without notice."
)


# Deliberately no `from ._gridpack.emt import *`: EMT has to resolve through
# here for the warning to fire at all.
def __getattr__(name):
    if name.startswith("_"):
        raise AttributeError(name)
    attr = getattr(_mod, name)
    warnings.warn(_EXPERIMENTAL, FutureWarning, stacklevel=2)
    return attr


def __dir__():
    return sorted(n for n in dir(_mod) if not n.startswith("_"))
