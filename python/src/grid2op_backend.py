import os, sys
import numpy as np
import pandas as pd
from typing import Optional, Tuple, Union
from bs4 import BeautifulSoup
import xml.etree.ElementTree as ET

import grid2op
from grid2op.Backend import Backend   # required

import gridpack
from mpi4py import MPI

class Grid:
    def __init__(self):
        pass

class GridPACKBackend(Backend):
    def __init__(self, grid2op_stepsize, log_freq=1) -> None:
        # Run Backend init
        super().__init__(can_be_copied=False)

        # Create GridPACK environment and pass the communicator to it
        self.env = gridpack.Environment()
        comm = gridpack.Communicator()

        np = gridpack.NoPrint()
        sys.stdout.write("%d: NoPrint status: %r\n" % (comm.rank(), np.status()))
        np.setStatus (True)

        # Create hadrec module
        self._dsapp = gridpack.dynamic_simulation.DSFullApp()
        self._grid2op_stepsize = grid2op_stepsize
        print(f"[INFO] Grid2Op Simulation Timestep: {self._grid2op_stepsize} seconds")

        # NOTE: add a timestamp along with the counter - get the time from GridPACK - if possible
        self._log_freq = log_freq

    def _read_bus_gen_load_data(self):
        # load and generator data
        # NOTE: Generator and Load values change only when there are dynamic load/generators. If they are not dynamic, the code currently returns the previous value. 
        bus_data = []
        gen_data = []
        load_data = []
        for bus in range(self.n_sub):
            # Dynamic simulation values are in pu, need to convert them MW and MVar
            CASE_SBASE = self._dsapp.getBusInfoReal(bus, "CASE_SBASE")
            BUS_BASEKV = self._dsapp.getBusInfoReal(bus, "BUS_BASEKV")

            # bus list
            vn_pu = self._dsapp.getBusInfoReal(bus, "BUS_VMAG_CURRENT") # in p.u.
            vn_kv = vn_pu * BUS_BASEKV # in kV

            # BUS_VMAG_CURRENT - latest value
            bus_data.append({
                "name": self._dsapp.getBusInfoString(bus, "BUS_NAME"),
                "id": self._dsapp.getBusInfoInt(bus, "BUS_NUMBER"),
                "tick": self._counter_time,
                "vn_pu": vn_pu,
                "vn_kv": vn_kv,
                "frequency": self._dsapp.getBusInfoReal(bus, "BUS_FREQUENCY"),
                "type": self._dsapp.getBusInfoInt(bus, "BUS_TYPE")
            })
            
            # bus list
            for g in range(self._dsapp.numGenerators(bus)):
                # gen data
                gen_data.append({
                    "tick": self._counter_time,
                    "name": self._dsapp.getBusInfoString(bus, "GENERATOR_ID", g),
                    "bus": bus,
                    "p_mw": self._dsapp.getBusInfoReal(bus, 'GENERATOR_PG_CURRENT', g) * CASE_SBASE,
                    "q_mvar": self._dsapp.getBusInfoReal(bus, 'GENERATOR_QG_CURRENT', g) * CASE_SBASE,
                    "vn_kv": vn_kv,
                    "min_q_mvar": self._dsapp.getBusInfoReal(bus, "GENERATOR_QMIN", g),
                    "max_q_mvar": self._dsapp.getBusInfoReal(bus, "GENERATOR_QMAX", g),
                    "in_service": self._dsapp.getBusInfoInt(bus, "GENERATOR_STAT", g)
                })
            for l in range(self._dsapp.numLoads(bus)):
                # load list
                load_data.append({
                    "tick": self._counter_time,
                    "name": self._dsapp.getBusInfoString(bus, "LOAD_ID", l),
                    "bus": bus,
                    "p_mw_int": self._dsapp.getBusInfoReal(bus, 'LOAD_PL', l),
                    "q_mvar_int": self._dsapp.getBusInfoReal(bus, 'LOAD_QL', l),
                    "p_mw": self._dsapp.getBusInfoReal(bus, 'LOAD_PL_CURRENT', l),
                    "q_mvar": self._dsapp.getBusInfoReal(bus, 'LOAD_QL_CURRENT', l),
                    "scaling": self._dsapp.getBusInfoReal(bus, "LOAD_SCALE", l), 
                    "in_service": self._dsapp.getBusInfoInt(bus, "LOAD_STATUS", l)
                })

        return bus_data, gen_data, load_data

    def _read_line_data(self, branch, elem_num, branch_info):
        data_dict = {}
        
        # record
        data_dict["from_bus"] = branch_info["from_bus"]
        data_dict["to_bus"] = branch_info["to_bus"]
        
        # voltages
        data_dict["vn_from_kv"] = branch_info["vn_from_kv"]
        data_dict["vn_to_kv"] = branch_info["vn_to_kv"]
            
        # p-values
        data_dict["p_from_mw"] = self._dsapp.getBranchInfoReal(branch, 'BRANCH_FROM_P_CURRENT', elem_num) * branch_info["f_case_sbase"]
        data_dict["p_to_mw"] = self._dsapp.getBranchInfoReal(branch, 'BRANCH_TO_P_CURRENT', elem_num) * branch_info["t_case_sbase"]

        # q-values
        data_dict["q_from_mvar"] = self._dsapp.getBranchInfoReal(branch, 'BRANCH_FROM_Q_CURRENT', elem_num) * branch_info["f_case_sbase"]
        data_dict["q_to_mvar"] = self._dsapp.getBranchInfoReal(branch, 'BRANCH_TO_Q_CURRENT', elem_num) * branch_info["t_case_sbase"]

        # current values
        data_dict["i_from_ka"] = self._dsapp.getBranchInfoReal(branch, 'BRANCH_IRFLOW_CURRENT', elem_num) * branch_info["f_case_sbase"]
        data_dict["i_to_ka"] = self._dsapp.getBranchInfoReal(branch, 'BRANCH_IRFLOW_CURRENT', elem_num) * branch_info["t_case_sbase"] # output in amps

        data_dict["length_km"] = self._dsapp.getBranchInfoReal(branch, "BRANCH_LENGTH", elem_num) 
        data_dict["r_ohm_per_km"] = self._dsapp.getBranchInfoReal(branch, "BRANCH_SEQ_RLINZ", elem_num)
        data_dict["x_ohm_per_km"] = self._dsapp.getBranchInfoReal(branch, "BRANCH_SEQ_XLINZ", elem_num)
        # data_dict["max_i_ka"] = 
        
        return data_dict
        
    def _read_transformer_data(self, branch, elem_num, branch_info):
        data_dict = {}
        
        # record
        data_dict["hv_bus"] = branch_info["from_bus"]
        data_dict["lv_bus"] = branch_info["to_bus"]

        # voltages
        data_dict["vn_hv_kv"] = branch_info["vn_from_kv"]
        data_dict["vn_lv_kv"] = branch_info["vn_to_kv"]
            
        # p-values
        data_dict["p_hv_mw"] = self._dsapp.getBranchInfoReal(branch, 'BRANCH_FROM_P_CURRENT', elem_num) * branch_info["f_case_sbase"]
        data_dict["p_lv_mw"] = self._dsapp.getBranchInfoReal(branch, 'BRANCH_TO_P_CURRENT', elem_num) * branch_info["t_case_sbase"]

        # q-values
        data_dict["q_hv_mvar"] = self._dsapp.getBranchInfoReal(branch, 'BRANCH_FROM_Q_CURRENT', elem_num) * branch_info["f_case_sbase"]
        data_dict["q_lv_mvar"] = self._dsapp.getBranchInfoReal(branch, 'BRANCH_TO_Q_CURRENT', elem_num) * branch_info["t_case_sbase"]

        # current values
        data_dict["i_hv_ka"] = self._dsapp.getBranchInfoReal(branch, 'BRANCH_IRFLOW_CURRENT', elem_num) * branch_info["f_case_sbase"]
        data_dict["i_lv_ka"] = self._dsapp.getBranchInfoReal(branch, 'BRANCH_IRFLOW_CURRENT', elem_num) * branch_info["t_case_sbase"] # output in amps
        
        return data_dict
        
    def _read_line_transformer_data(self):
        # objects for data collection
        line_data = []
        transformer_data = []
        
        # line level data - lines and branches could be different. The two end points can have multiple lines but only one branch. 
        nbranch = self._dsapp.totalBranches()
        
        # iterate over each branch
        for branch in range(0, nbranch):
            # branch information
            branch_info = {}
            branch_info["name"] = self._dsapp.getBranchInfoString(branch, "BRANCH_NAME")

            # end points
            (f, t) = self._dsapp.getBranchEndpoints(branch)
            branch_info["from_bus"] = self.BUS_DICT[f]
            branch_info["to_bus"] = self.BUS_DICT[t]

            # from and to bus kv
            branch_info["f_buskv"] = self._dsapp.getBusInfoReal(
                branch_info["from_bus"], "BUS_BASEKV")
            branch_info["t_buskv"] = self._dsapp.getBusInfoReal(
                branch_info["to_bus"], "BUS_BASEKV")
            
            # from and to bus case sbase
            branch_info["f_case_sbase"] = self._dsapp.getBusInfoReal(
                branch_info["from_bus"], "CASE_SBASE")
            branch_info["t_case_sbase"] = self._dsapp.getBusInfoReal(
                branch_info["to_bus"], "CASE_SBASE")
            
            # from and to voltage
            branch_info["vn_from_kv"] = self._dsapp.getBusInfoReal(
                branch_info["from_bus"], 'BUS_VMAG_CURRENT') * branch_info["f_buskv"]
            branch_info["vn_to_kv"] = self._dsapp.getBusInfoReal(
                branch_info["to_bus"], 'BUS_VMAG_CURRENT') * branch_info["t_buskv"]

            # number of lines
            n_elements = self._dsapp.getBranchInfoInt(branch, 'BRANCH_NUM_ELEMENTS')

            # iterate over each line
            for elem_num in range(n_elements):
                # check if line or transformer
                is_transformer = self._dsapp.getBranchInfoReal(branch, "BRANCH_TAP", elem_num)

                if is_transformer == 1: # if transformer
                    data_dict = self._read_transformer_data(branch, elem_num, branch_info)
                else:
                    data_dict = self._read_line_data(branch, elem_num, branch_info)
                
                # common data
                data_dict["tick"] = self._counter_time
                data_dict["branch_num"] = branch
                data_dict["line_num"] = elem_num
                data_dict["name"] = branch_info["name"]
                data_dict["in_service"] = self._dsapp.getBranchInfoInt(branch, "BRANCH_STATUS", elem_num)
                data_dict["thermal_limit_a"] = self._dsapp.getBranchInfoReal(branch, "BRANCH_RATING_A", elem_num) # - each branch - line connects two points and there could be multiple branches

                # append dict to list
                if is_transformer == 1: # if transformer
                    transformer_data.append(data_dict)
                else:
                    line_data.append(data_dict)
                
        return line_data, transformer_data
    
    def _build_grid(self):
        # create grid object and update data
        grid = Grid()
        self._dsapp.updateData()   

        # then fill the "n_sub" and "sub_info"
        self.n_sub = self._dsapp.totalBuses()

        # read bus, load, and gen data
        bus_data, gen_data, load_data = self._read_bus_gen_load_data()
        # convert to dataframes
        # bus and its results
        grid.bus = pd.DataFrame(bus_data)
        grid.res_bus = pd.DataFrame(index=grid.bus.index, columns=grid.bus.columns)
        
        # gen and its results
        grid.gen = pd.DataFrame(gen_data)
        grid.res_gen = pd.DataFrame(index=grid.gen.index, columns=grid.gen.columns)
        
        # load and its results
        grid.load = pd.DataFrame(load_data)
        grid.res_load = pd.DataFrame(index=grid.load.index, columns=grid.load.columns)
        # print("Actual Load")
        # print(grid.load.loc[grid.load["bus"].isin([52, 56, 71, 144, 167])])
        # print(grid.load)
        # print("After Load")
        
        
        # create bus dict for translation
        self.BUS_DICT = grid.bus[["id"]].reset_index().set_index("id").to_dict()['index']

        # read line and transformer data
        line_data, transformer_data = self._read_line_transformer_data()

        # line and its results
        # NOTE: bus dict is needed to translate actual bus number to dataframe index. The from_bus and to_bus columns in grid.line needs to be replaced with this new index to be consistent across the code. 
        grid.line = pd.DataFrame(line_data)
        grid.res_line = pd.DataFrame(index=grid.line.index, columns=grid.line.columns)
        
        # transformers and its results - variabe called BRANCH_TAP in the line - non-zero tap implies transformer
        grid.trafo = pd.DataFrame(transformer_data)
        grid.res_trafo = pd.DataFrame(index=grid.trafo.index, columns=grid.trafo.columns)
        
        # return
        return grid

    def load_grid(self, 
                  path : Union[os.PathLike, str], 
                  filename : Optional[Union[os.PathLike, str]]=None
        ) -> None:
        '''
        # called once   
        This step is called only ONCE, when the grid2op environment is created. In this step, you read a grid file (in the format that you want) and the backend should inform grid2op about the "objects" on this powergrid and their location.
        '''
        # first load the grid from the file
        # self.full_path = path
        # if filename is not None:
        #     self.full_path = os.path.join(self.full_path, filename)
        self.full_path = self.make_complete_path(path, filename)
        
        # read XML to get timestep
        with open(self.full_path, 'r') as f:
            data = f.read()
        bs_data = BeautifulSoup(data, "lxml")
        # timestep
        self._gridpack_stepsize = float(bs_data.find('timestep').text)
        print(f"[INFO] GridPACK simulation step size: {self._gridpack_stepsize} seconds")

        self._counter = 0
        self._counter_time = self._counter * self._gridpack_stepsize
        
        # NOTE: Need to duble check this
        self.cannot_handle_more_than_2_busbar()

        # select index
        sel_index = -1

        # solve the power flow - to load the grid data
        self._dsapp.solvePowerFlowBeforeDynSimu(self.full_path, sel_index)  # 0 inidcates that solves the first raw file for power flow, the xml file supports multiple power flow raw files read in
        self._dsapp.readGenerators(sel_index);
        self._dsapp.readSequenceData();
        self._dsapp.initialize();
        self._dsapp.setGeneratorWatch();
        
        # Building grid object from hadapp
        self._grid = self._build_grid()
        
        # then fill the number and location of loads
        self.n_load = self._grid.load.shape[0]
        self.load_to_subid = np.zeros(self.n_load, dtype=int)
        for load_id in range(self.n_load):
            self.load_to_subid[load_id] = self._grid.load.iloc[load_id]["bus"]
            
        # then fill the number and location of generators
        self.n_gen = self._grid.gen.shape[0]
        self.gen_to_subid = np.zeros(self.n_gen, dtype=int)
        for gen_id in range(self.n_gen):
            self.gen_to_subid[gen_id] = self._grid.gen.iloc[gen_id]["bus"]
            
        # then fill the number and location of storage units
        # self.n_storage = self._grid.storage.shape[0]
        # self.storage_to_subid = np.zeros(self.n_storage, dtype=int)
        # for storage_id in range(self.n_storage):
        #     self.storage_to_subid[storage_id] = self._grid.storage.iloc[storage_id]["bus"]
        
        # WARNING
        # for storage, their description is loaded in a different file (see 
        # the doc of Backend.load_storage_data)
        # to start we recommend you to ignore the storage unit of your grid with:
        self.set_no_storage()
        
        # finally handle powerlines
        # NB: grid2op considers that trafos are powerlines.
        # so we decide here to say: first n "powerlines" of grid2Op
        # will be pandapower powerlines and
        # last k "powerlines" of grid2op will be the trafos of pandapower.
        self.n_line = self._grid.line.shape[0] + self._grid.trafo.shape[0]
        self.line_or_to_subid = np.zeros(self.n_line, dtype=int)
        self.line_ex_to_subid = np.zeros(self.n_line, dtype=int)
        for line_id in range(self._grid.line.shape[0]):
            self.line_or_to_subid[line_id] = self._grid.line.iloc[line_id]["from_bus"]
            self.line_ex_to_subid[line_id] = self._grid.line.iloc[line_id]["to_bus"]

        nb_powerline = self._grid.line.shape[0]
        for trafo_id in range(self._grid.trafo.shape[0]):
            self.line_or_to_subid[trafo_id + nb_powerline] = self._grid.trafo.iloc[trafo_id]["hv_bus"]
            self.line_ex_to_subid[trafo_id + nb_powerline] = self._grid.trafo.iloc[trafo_id]["lv_bus"]
            
        # FIXME: Missing info from gridpack
        # and now the thermal limit
        # self.thermal_limit_a = 1000. * np.concatenate(
        #     (
        #         self._grid.line["max_i_ka"].values,
        #         self._grid.trafo["sn_mva"].values
        #         / (np.sqrt(3) * self._grid.trafo["vn_hv_kv"].values),
        #     )
        # )
        
        # FIXME: Random number
        self.thermal_limit_a = 1000 * np.ones(self.n_line, dtype=int)
            
        self._compute_pos_big_topo()

        # # transfer data from power flow network to dynamic simulation network
        # self._hadapp.transferPFtoDS()
        
        # # define a bus fault
        # busfault = gridpack.dynamic_simulation.Event()
        # busfault.start = 1.0 # fault start time
        # busfault.end = 1.1   # fault end time
        # busfault.step = 0.005  # fault duration simu time step
        # busfault.isBus = True
        # busfault.bus_idx = 22 # bus number of the fault		  

        # busfaultlist = gridpack.dynamic_simulation.EventVector([busfault])

        # # initialize the dynamic simulation
        # self._hadapp.initializeDynSimu(busfaultlist, 0) # 0 inidcates read in the first dyr dynamic parameter file, the xml file supports multiple dyr files read in

        # Remember the input file was read into the Configuration singleton
        conf = gridpack.Configuration()
        cursor = conf.getCursor("Configuration.Dynamic_simulation")

        # get faults
        faults = self._dsapp.getEvents(cursor)

        # solve with faults
        self._dsapp.solvePreInitialize(faults[0])

    def reset(self, 
              path: Union[os.PathLike, str], 
              grid_filename: Optional[Union[os.PathLike, str]]=None
            ) -> None:
        # TODO: Reset GridPACK simulator
        # self._counter = 0
        # self._counter_time = self._counter * self._gridpack_stepsize
        pass

    def apply_action(self, backendAction: Union["grid2op.Action._backendAction._BackendAction", None]) -> None:
        '''
        # called for each "step", thousands of times
        # modify the topology, load, generation etc.
        '''
        "NOTE: Grid2Op registers an event and GridPACK executes that action in the next time step. We need to check can we register an event later in the code."
        # the following few lines are highly recommended
        if backendAction is None:
            return
        
        (
            active_bus,
            (prod_p, prod_v, load_p, load_q, storage),
            _,
            shunts__,
        ) = backendAction()
        print("[GridPACK] Executing action")
        # print(self._grid.res_load)
        print(load_p)

        # change the active values of the loads
        load_dict = {}
        for load_id, new_p in load_p:
            print(load_id, new_p, self._grid.load["p_mw"].iloc[load_id])
            bus = self.BUS_DICT[self._grid.load.loc[load_id, "bus"]]
            # bus = self._grid.load.loc[load_id, "bus"]
            case_sbase = self._dsapp.getBusInfoReal(bus, "CASE_SBASE")
            load_dict[bus] = {"p": new_p / case_sbase}
            # self._grid.load["p_mw"].iloc[load_id] = new_p
            
        # change the reactive values of the loads
        for load_id, new_q in load_q:
            # bus = self.BUS_DICT[]
            bus = self.BUS_DICT[self._grid.load.loc[load_id, "bus"]]
            # bus = self._grid.load.loc[load_id, "bus"]
            case_sbase = self._dsapp.getBusInfoReal(bus, "CASE_SBASE")
            if bus in load_dict:
                load_dict[bus]["q"] = new_q / case_sbase
            # self._grid.load["q_mvar"].iloc[load_id] = new_q

        # update load values
        if len(load_dict) != 0:
            load_frame = pd.DataFrame(load_dict).T
            self._dsapp.scatterInjectionLoadNew_compensateY(
                load_frame.index.values,
                load_frame["p"].values,
                load_frame["q"].values,
            )
        
        # change the active value of generators
        for gen_id, new_p in prod_p:
            self._grid.gen["p_mw"].iloc[gen_id] = new_p
            
        # for the voltage magnitude, pandapower expects pu but grid2op provides kV,
        # so we need a bit of change
        for gen_id, new_v in prod_v:
            self._grid.gen["vn_pu"].iloc[gen_id] = new_v  # but new_v is not pu !

        # disconnected powerlines are indicated because they are
        # connected to bus "-1" in the `get_lines_or_bus()` and
        # `get_lines_ex_bus()`
        # NB : at time of writing, grid2op side a powerline disconnected
        # on a side (*eg* "or" side or "ext" side) is
        # disconnected on both.
        
        # the only difficulty here is that grid2op considers that
        # trafo are also powerline.
        # We already "solved" that by saying that the "k" last "lines"
        # from grid2op point of view will indeed be trafos.
        
        # n_line_pp = self._grid.line.shape[0]
        
        # handle the disconnection on "or" side
        # lines_or_bus = backendAction.get_lines_or_bus()
        # for line_id, new_bus in lines_or_bus:
        #     if line_id < n_line_pp:
        #         # a pandapower powerline has bee disconnected in grid2op
        #         dt = self._grid.line
        #         line_id_db = line_id
        #     else:
        #         # a pandapower trafo has bee disconnected in grid2op
        #         dt = self._grid.trafo
        #         line_id_db = line_id - n_line_pp

        #     if new_bus == -1:
        #         # element was disconnected
        #         dt.iloc[line_id_db]["in_service"] = False
        #     else:
        #         # element was connected
        #         dt.iloc[line_id_db]["in_service"] = True

        # lines_ex_bus = backendAction.get_lines_ex_bus()
        # for line_id, new_bus in lines_ex_bus:
        #     if line_id < n_line_pp:
        #         # a pandapower powerline has bee disconnected in grid2op
        #         dt = self._grid.line
        #         line_id_db = line_id
        #     else:
        #         # a pandapower trafo has bee disconnected in grid2op
        #         dt = self._grid.trafo
        #         line_id_db = line_id - n_line_pp

        #     if new_bus == -1:
        #         # element was disconnected
        #         dt.iloc[line_id_db]["in_service"] = False
        #     else:
        #         # element was connected
        #         dt.iloc[line_id_db]["in_service"] = True

    def _update_bus_gen_load_data(self):
        # read data
        bus_data, gen_data, load_data = self._read_bus_gen_load_data()

        # update the dataframe
        self._grid.res_bus = pd.DataFrame(bus_data)
        self._grid.res_gen = pd.DataFrame(gen_data)
        self._grid.res_load = pd.DataFrame(load_data)
        
    def _update_line_transformer_data(self):
        # read data
        line_data, trafo_data = self._read_line_transformer_data()
        
        # line data
        self._grid.res_line = pd.DataFrame(line_data)
        
        # transformer data
        self._grid.res_trafo = pd.DataFrame(trafo_data)

        # print(self._grid.res_line[self._grid.res_line.isnull().any(axis=1)].tail(5))
        # print(self._grid.res_trafo[self._grid.res_trafo.isnull().any(axis=1)].tail(5))
        
        # sys.exit(1)
        
    def _update_data(self):
        # to get current data from dynamic simulation to data collection object
        self._dsapp.updateData()
        
        # update data for loads and generators
        self._update_bus_gen_load_data()

        # update line and transformer data
        self._update_line_transformer_data()

    def _reset_data_collectors(self):
        # output for plotting
        self.bus_logger = []
        self.gen_logger = []
        self.load_logger = []
        self.line_logger = []
        self.trafo_logger = []

    def _log_data(self):
        # append to out list
        self.bus_logger.append(pd.concat([self._grid.bus[["id"]], self._grid.res_bus], axis=1))
        self.gen_logger.append(pd.concat([self._grid.gen[["name", "bus"]], self._grid.res_gen], axis=1))
        self.load_logger.append(pd.concat([self._grid.load[["name", "bus"]], self._grid.res_load], axis=1))
        self.line_logger.append(self._grid.res_line)
        self.trafo_logger.append(self._grid.res_trafo)
    
    def runpf(self, is_dc : bool=False):
        '''
        # called for each "step", thousands of times
        # run the solver
        '''
        
        # Move GridPACK simulation by n_steps 	
        n_steps = int(self._grid2op_stepsize / self._gridpack_stepsize)
        assert n_steps > 0, f"[ERROR] Grid2Op timestep ({self._grid2op_stepsize} seconds) cannot be smaller that GridPACK timestep ({self._gridpack_stepsize} seconds)."

        # reset data collector object
        self._reset_data_collectors()

        # run GridPACK
        for i in range(n_steps):
            # run GridPACK simulation by one time step
            self._dsapp.executeOneSimuStep()  

            # update data in result dataframes
            self._update_data()

            # log data at log frequency
            if ((self._counter-1) % self._log_freq) == 0:
                self._log_data()

            # counter
            self._counter += 1
            self._counter_time = self._counter * self._gridpack_stepsize

        return True, None    
    
    def _aux_get_topo_vect(self, res, dt, key, pos_topo_vect, add_id=0):
        # we loop through each element of the table
        # (each table representing either the loads, or the generators or the powerlines or the trafos)
        # then we assign the right bus (local - eg 1 or 2) to the right
        # component of the vector "res"  (the component is given by the "pos_topo_vect" - eg self.load_pos_topo_vect
        # when we look at the loads)
        el_id = 0
        for (status, bus_id) in dt[["in_service", key]].values:
            my_pos_topo_vect = pos_topo_vect[el_id + add_id]
            if status:
                local_bus = self.global_bus_to_local_int(bus_id, my_pos_topo_vect)
            else:
                local_bus = -1
            res[my_pos_topo_vect] = local_bus
            el_id += 1

    def get_topo_vect(self):
        '''
        # retrieve the results
        '''
        res = np.full(self.dim_topo, fill_value=-2, dtype=int)
        # read results for load
        self._aux_get_topo_vect(res, self._grid.load, "bus", self.load_pos_topo_vect)
        # then for generators
        self._aux_get_topo_vect(res, self._grid.gen, "bus", self.gen_pos_topo_vect)
        # then each side of powerlines
        self._aux_get_topo_vect(res, self._grid.line, "from_bus", self.line_or_pos_topo_vect)
        self._aux_get_topo_vect(res, self._grid.line, "to_bus", self.line_ex_pos_topo_vect)
        
        # then for the trafos, but remember pandapower trafos are powerlines in grid2Op....
        # so we need to trick it a bit
        # (we can do this trick because we put the trafo "at the end" of the powerline in grid2op
        # in the Step1_loading.py)
        n_line_pp = self._grid.line.shape[0]
        self._aux_get_topo_vect(res, self._grid.trafo, "hv_bus", self.line_or_pos_topo_vect, add_id=n_line_pp)
        self._aux_get_topo_vect(res, self._grid.trafo, "lv_bus", self.line_ex_pos_topo_vect, add_id=n_line_pp)            
        return res

    def loads_info(self)-> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        '''
        # retrieve the results
        '''
        # carefull with copy / deep copy
        load_p = self._grid.res_load["p_mw"].values  # in MW
        load_q = self._grid.res_load["q_mvar"].values  # in MVAr
        load_v = self._grid.res_bus.iloc[self._grid.load["bus"].values]["vn_kv"].values  # in kV

        return load_p, load_q, load_v

    def generators_info(self)-> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        '''
        # retrieve the results
        '''
        prod_p = self._grid.res_gen["p_mw"].values  # in MW
        prod_q = self._grid.res_gen["q_mvar"].values  # in MVAr
        prod_v = self._grid.res_gen["vn_kv"].values  # in kV
        
        return prod_p, prod_q, prod_v
        
    def _aux_get_line_info(self, colname_powerline, colname_trafo):
        """
        concatenate the information of powerlines and trafo using 
        the convention that "powerlines go first" then trafo
        """
        res = np.concatenate(
            (
                self._grid.res_line[colname_powerline].values,
                self._grid.res_trafo[colname_trafo].values,
            )
        )
        # NOTE: Typecasting this variable to address numpy abs issue. Need to comeback to this.
        return np.float32(res)
    
    def lines_or_info(self)-> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Retrieve the results
        
        Main method to retrieve the information at the "origin" side of the powerlines and transformers.

        We simply need to follow the convention we adopted:

        - origin side (grid2op) will be "from" side for pandapower powerline
        - origin side (grid2op) will be "hv" side for pandapower trafo
        - we chose to first have powerlines, then transformers

        (convention chosen in :func:`EducPandaPowerBackend.load_grid`)

        """
        # FIXME: Placeholder values
        # return np.ones(self.n_line), np.ones(self.n_line), np.ones(self.n_line), np.ones(self.n_line)
    
        p_or = self._aux_get_line_info("p_from_mw", "p_hv_mw")
        q_or = self._aux_get_line_info("q_from_mvar", "q_hv_mvar")
        v_or = self._aux_get_line_info("vn_from_kv", "vn_hv_kv")  # in kV
        a_or = self._aux_get_line_info("i_from_ka", "i_hv_ka")  # grid2op expects amps (A) gridpack returns amps (A)
        # * 100 / kV (TAmps)
        
        # there would be a bug in v_or because of the way pandapower
        # internally looks at the extremity of powerlines / trafos.
        # we fix it here:
        status = np.concatenate(
            (
                self._grid.line["in_service"].values,
                self._grid.trafo["in_service"].values,
            )
        )
        
        # NOTE: v_or[~s] doesn't work in this version of python
        v_or[[~s for s in status]] = 0.
        
        return p_or, q_or, v_or, a_or

    def lines_ex_info(self):
        """
        Retrieve the results
        
        Main method to retrieve the information at the "extremity" side of the powerlines and transformers.

        We simply need to follow the convention we adopted:

        - extremity side (grid2op) will be "to" side for pandapower powerline
        - extremity side (grid2op) will be "lv" side for pandapower trafo
        - we chose to first have powerlines, then transformers

        (convention chosen in function :func:`EducPandaPowerBackend.load_grid`)

        """
        # FIXME: Placeholder values
        # return np.ones(self.n_line), np.ones(self.n_line), np.ones(self.n_line), np.ones(self.n_line)
    
        p_ex = self._aux_get_line_info("p_to_mw", "p_lv_mw")
        q_ex = self._aux_get_line_info("q_to_mvar", "q_lv_mvar")
        v_ex = self._aux_get_line_info("vn_to_kv", "vn_lv_kv")  # in pu
        a_ex = self._aux_get_line_info("i_to_ka", "i_lv_ka") # grid2op expects amps (A) gridpack returns amps (A)
        
        # there would be a bug in v_ex because of the way pandapower
        # internally looks at the extremity of powerlines / trafos.
        # we fix it here:
        status = np.concatenate(
            (
                self._grid.line["in_service"].values,
                self._grid.trafo["in_service"].values,
            )
        )

        # NOTE: v_or[~s] doesn't work in this version of python
        v_ex[[~s for s in status]] = 0.
        
        return p_ex, q_ex, v_ex, a_ex
    
    def close(self):
        self._dsapp = None
        self.env = None