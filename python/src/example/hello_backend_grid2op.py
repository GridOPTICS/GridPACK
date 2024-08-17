import grid2op
import warnings
import sys
sys.path.append("../")
from grid2op_backend import GridPACKBackend
from grid2op.PlotGrid.PlotMatplot import PlotMatplot

filename = "input_39bus_step005_v33.xml"

# make environment
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    env = grid2op.make("test_grid2op",
                       grid_path=filename,
                       backend=GridPACKBackend()
                       )

print("============ Initialized Environment =============")

# reset environment
obs = env.reset()

# print network analytics
print("Number of buses:  %d" % (obs.n_sub))
print("Number of generators: %d" % (obs.n_gen))
print("Number of loads: %d" % (obs.n_load))
print("Number of lines: %d" % (obs.n_line))
print("Number of storage units: %d" % (obs.n_storage))

## step environment
# obs, reward, done, info = env.step(env.action_space())

# print(env.observation_space)
# plot_helper = PlotMatplot(env.observation_space)
# fig = plot_helper.plot_layout()
