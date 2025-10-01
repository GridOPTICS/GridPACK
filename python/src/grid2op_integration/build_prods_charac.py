import os
import pandas as pd

import gridpack
from gridpack.dynamic_simulation import DSFullApp, Event, EventVector

# inname = "input_9bus.xml"
# inname = "input_9bus_fault.xml"
# inname = "input_39bus_IBR.xml"
# inname = "input_39bus_IBR_fault.xml"
# inname = "input_145.xml"
inname = "input_145_fault.xml"
# inname = "input_300_cmpld.xml"
# inname = "input_3000.xml"

### failed ones
# inname = "input_240.xml"

gridpack.NoPrint().setStatus(False)
env = gridpack.Environment()
comm = gridpack.Communicator()

ds_app = DSFullApp()

ds_app.solvePowerFlowBeforeDynSimu(inname, -1)

conf = gridpack.Configuration()
cursor = conf.getCursor("Configuration.Dynamic_simulation")

ds_app.readGenerators(-1);
ds_app.readSequenceData();
ds_app.initialize();
ds_app.setGeneratorWatch();
ds_app.setObservations(cursor)

busfault = Event()

busfault.start = 1.0 # fault start time
busfault.end = 1.1   # fault end time
busfault.step = 0.005  # fault duration simu time step
busfault.isBus = True
busfault.bus_idx = 22 # bus number of the fault

busfaultlist = EventVector([busfault])

ds_app.solvePreInitialize(busfaultlist[0]);

out = []
nbus = ds_app.totalBuses()
out_dict = {}
i = 0
for bus in range(nbus):
    out_dict[ds_app.getBusInfoInt(bus, "BUS_NUMBER")] = bus
    for g in range(ds_app.numGenerators(bus)):
        gentype = ds_app.getBusInfoString(bus, "GENERATOR_MODEL", g)

        Pmax = ds_app.getBusInfoReal(bus, "GENERATOR_PMAX", g)
        
        # if gentype == "GENROU":
        gentype = "thermal"
        max_ramp_up = Pmax
        max_ramp_down = Pmax
        marginal_cost = 35
        start_cost = 120 * Pmax
        shut_down_cost = 0.10 * start_cost
        min_up_time = 4
        min_down_time = 2
        if gentype == "REGCA1":
            gentype = "wind"
            max_ramp_up = Pmax
            max_ramp_down = Pmax
            marginal_cost = 0
            start_cost = 0
            shut_down_cost = 0
            min_up_time = 0
            min_down_time = 0
        
        out.append({
            "Pmax": Pmax,
            "Pmin": max(0, ds_app.getBusInfoReal(bus, "GENERATOR_PMIN", g)),
            "name": f"gen_{out_dict[ds_app.getBusInfoInt(bus, 'GENERATOR_BUSNUMBER', g)]}_{i}",
            "type": gentype,
            "bus": ds_app.getBusInfoInt(bus, "GENERATOR_BUSNUMBER", g),
            "max_ramp_up": max_ramp_up,
            "max_ramp_down": max_ramp_down,
            "start_cost": start_cost,
            "shut_down_cost": shut_down_cost,
            "marginal_cost": marginal_cost,
            "min_up_time": min_up_time,
            "min_down_time": min_down_time,
            "V": ds_app.getBusInfoReal(bus, "GENERATOR_VS", g) * ds_app.getBusInfoReal(bus, "BUS_BASEKV")
        })
        i += 1

outdir = f"./grid2op_{inname.split('.')[0]}"
os.makedirs(outdir, exist_ok=True)
df_out = pd.DataFrame(out)
print(df_out, df_out["Pmin"].sort_values())
df_out.to_csv(f"{outdir}/prods_charac.csv", index=False)
        
ds_app.close()
del ds_app
