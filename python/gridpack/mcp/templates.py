# -------------------------------------------------------------
# file: gridpack/mcp/templates.py
# -------------------------------------------------------------
# Build an input XML around a user's RAW (and DYR) file.
# -------------------------------------------------------------

"""Default Powerflow and Dynamic_simulation blocks for a network file.

The settings are those of the bundled ``input_9b3g.xml`` without its
events; every generator in the RAW is watched.
"""

from ..rawedit import raw_generators, raw_version  # noqa: F401

__all__ = ["raw_version", "raw_generators", "build_xml"]

_LINEAR = """    <LinearSolver>
      <PETScOptions>
        -ksp_type richardson
        -pc_type lu
        -pc_factor_mat_solver_type superlu_dist
        -ksp_max_it 1
      </PETScOptions>
    </LinearSolver>
"""

_POWERFLOW = """  <Powerflow>
    <%(tag)s>%(raw)s</%(tag)s>
    <maxIteration>50</maxIteration>
    <tolerance>1.0e-6</tolerance>
%(linear)s    <UseNonLinear>false</UseNonLinear>
    <UseNewton>false</UseNewton>
    <NewtonRaphsonSolver>
      <SolutionTolerance>1.0E-05</SolutionTolerance>
      <FunctionTolerance>1.0E-05</FunctionTolerance>
      <MaxIterations>50</MaxIterations>
      <LinearSolver>
        <SolutionTolerance>1.0E-08</SolutionTolerance>
        <MaxIterations>50</MaxIterations>
        <PETScOptions>
          -ksp_type bicg
          -pc_type bjacobi
          -sub_pc_type ilu -sub_pc_factor_levels 5 -sub_ksp_type preonly
        </PETScOptions>
      </LinearSolver>
    </NewtonRaphsonSolver>
    <NonlinearSolver>
      <SolutionTolerance>1.0E-05</SolutionTolerance>
      <FunctionTolerance>1.0E-05</FunctionTolerance>
      <MaxIterations>50</MaxIterations>
      <PETScOptions>
        -ksp_type bicg
        -pc_type bjacobi
        -sub_pc_type ilu -sub_pc_factor_levels 5 -sub_ksp_type preonly
      </PETScOptions>
    </NonlinearSolver>
  </Powerflow>
"""

_DYNAMICS = """  <Dynamic_simulation>
    <generatorParameters>%(dyr)s</generatorParameters>
    <simulationTime>10.0</simulationTime>
    <timeStep>0.01</timeStep>
    <Events>
    </Events>
    <generatorWatch>
%(watch)s    </generatorWatch>
    <generatorWatchFrequency>1</generatorWatchFrequency>
    <generatorWatchFileName>gen_watch.csv</generatorWatchFileName>
%(linear)s    <LinearMatrixSolver>
      <Ordering>nd</Ordering>
      <Package>superlu_dist</Package>
      <Iterations>1</Iterations>
      <Fill>5</Fill>
    </LinearMatrixSolver>
  </Dynamic_simulation>
"""

_WATCH = """      <generator>
        <busID>%d</busID>
        <generatorID>%s</generatorID>
      </generator>
"""


def build_xml(raw, version, dyr=None, generators=()):
    """The XML text: a Powerflow block, plus Dynamic_simulation with
    ``dyr`` watching ``generators`` (``[(bus, id)]``)."""
    tag = "networkConfiguration" if version == 23 \
        else "networkConfiguration_v%d" % version
    text = _POWERFLOW % {"tag": tag, "raw": raw, "linear": _LINEAR}
    if dyr:
        watch = "".join(_WATCH % g for g in generators)
        text += _DYNAMICS % {"dyr": dyr, "watch": watch, "linear": _LINEAR}
    return '<?xml version="1.0" encoding="utf-8"?>\n<Configuration>\n%s' \
        '</Configuration>\n' % text
