from grid2op.Action import PlayableAction
from grid2op.Reward import L2RPNReward, RedispReward
from grid2op.Rules import DefaultRules
from grid2op.Chronics import Multifolder, GridStateFromFileWithForecasts
from grid2op.Backend import PandaPowerBackend

config = {}
# config = {
#     "action_class": PlayableAction,
#     "reward_class": L2RPNReward,
#     "gamerules_class": DefaultRules,
#     "chronics_class": None,
#     "grid_value_class": GridStateFromFileWithForecasts,
#     "volagecontroler_class": None,
#     "thermal_limits": None,
#     "names_chronics_to_grid": None,
# }