from grid2op.Action import TopologyAction
from grid2op.Reward import L2RPNReward
from grid2op.Rules import DefaultRules
from grid2op.Chronics import Multifolder
from grid2op.Chronics import GridStateFromFileWithForecasts

import sys
sys.path.append("../../")
from grid2op_backend import GridPACKBackend

config = {
    "backend": GridPACKBackend
}