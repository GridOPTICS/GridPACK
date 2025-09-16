import sys
import gridpack

class TestGridPACKReset():
    def __init__(self):
        self.sel_index = 0
        self.full_path = "input_9bus.xml"

    def _init_gridpack(self):
        gridpack.NoPrint().setStatus(False)
        self.env = gridpack.Environment()
        self.comm = gridpack.Communicator()

        # Create hadrec module
        self.dsapp = gridpack.dynamic_simulation.DSFullApp()

        # solve the power flow - to load the grid data
        self.dsapp.solvePowerFlowBeforeDynSimu(self.full_path, -1)  # 0 inidcates that solves the first raw file for power flow, the xml file supports multiple power flow raw files read in

        conf = gridpack.Configuration()
        cursor = conf.getCursor("Configuration.Dynamic_simulation")
        
        # get faults
        faults = self.dsapp.getEvents(cursor)

        self.dsapp.readGenerators(self.sel_index);
        self.dsapp.readSequenceData();
        self.dsapp.initialize();
        self.dsapp.setGeneratorWatch();
        self.dsapp.setObservations(cursor);
        self.dsapp.solvePreInitialize(faults[0])

    def reset(self):
        # breakpoint()
        # del self.dsapp
        self.dsapp = None
        self.env = None
        self.comm = None

    def copy(self):
        res = type(self)()
        res._init_gridpack()
        return res

if __name__ == "__main__":
    test = TestGridPACKReset()
    test._init_gridpack()
    test.reset()
    print("GridPACK reset done, will re-initialize now.")
    test._init_gridpack()
    test.reset()
    res = test.copy()
    res.reset()