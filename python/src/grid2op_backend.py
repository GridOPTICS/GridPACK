import os, sys
import numpy as np
import pandas as pd
from typing import Optional, Tuple, Union

import grid2op
from grid2op.Backend import Backend   # required

import gridpack
from mpi4py import MPI

class Grid:
    def __init__(self):
        pass

class GridPACKBackend(Backend):
    def __init__(self) -> None:
        # Run Backend init
        super().__init__(can_be_copied=False)

        # Create GridPACK environment and pass the communicator to it
        self.env = gridpack.Environment()
        comm = gridpack.Communicator()

        np = gridpack.NoPrint()
        sys.stdout.write("%d: NoPrint status: %r\n" % (comm.rank(), np.status()))
        np.setStatus (True)

        # Create hadrec module
        # FIXME: gridpack.dynamic_simulation.DSFullApp() is unable to retrieve any data
        self._hadapp = gridpack.dynamic_simulation.DSFullApp()

        self._counter = 0

    def build_grid(self):
        grid = Grid()

        # data lists
        bus_list = []
        gen_list = []
        load_list = []
        branch_list = []

        # then fill the "n_sub" and "sub_info"
        self.n_sub = self._hadapp.totalBuses()

        for bus in range(self.n_sub):
            
            # bus list
            bus_voltage = self._hadapp.getBusInfoReal(bus, "BUS_VOLTAGE_MAG")
            bus_list.append({
                "name": self._hadapp.getBusInfoString(bus, "BUS_NAME"),
                "id": self._hadapp.getBusInfoInt(bus, "BUS_NUMBER"),
                "vn_kv": bus_voltage,
                "type": self._hadapp.getBusInfoInt(bus, "BUS_TYPE")
            })
            for g in range(self._hadapp.numGenerators(bus)):
                # generator list
                gen_list.append({
                    "name": self._hadapp.getBusInfoString(bus, "GENERATOR_ID", g),
                    "bus": bus,
                    "p_mw": self._hadapp.getBusInfoReal(bus, "GENERATOR_PG", g),
                    "q_mvar": self._hadapp.getBusInfoReal(bus, "GENERATOR_QG", g),
                    "vm_pu": bus_voltage,
                    "min_q_mvar": self._hadapp.getBusInfoReal(bus, "GENERATOR_QMIN", g),
                    "max_q_mvar": self._hadapp.getBusInfoReal(bus, "GENERATOR_QMAX", g),
                    "in_service": self._hadapp.getBusInfoBool(bus, "GENERATOR_STAT", g)
                })
            for l in range(self._hadapp.numLoads(bus)):
                # load list
                load_list.append({
                    "name": self._hadapp.getBusInfoString(bus, "LOAD_ID", l),
                    "bus": bus,
                    "p_mw": self._hadapp.getBusInfoReal(bus, "LOAD_PL", l) / 100,
                    "q_mvar": self._hadapp.getBusInfoReal(bus, "LOAD_QL", l) / 100,
                    "scaling": self._hadapp.getBusInfoInt(bus, "LOAD_SCALE", l), 
                    "in_service": self._hadapp.getBusInfoBool(bus, "LOAD_STATUS", l)
                })
        
        # branch frame
        nbranch = self._hadapp.totalBranches()
        for branch in range(0, nbranch):
            # branch end points
            (f, t) = self._hadapp.getBranchEndpoints(branch)
            
            # append data to branch list
            branch_list.append({
                    "name": self._hadapp.getBranchInfoString(branch, "BRANCH_NAME"), 
                    "from_bus": f,
                    "to_bus": t,
                    "length_km": self._hadapp.getBranchInfoReal(branch, "BRANCH_LENGTH"), 
                    "r_ohm_per_km": self._hadapp.getBranchInfoReal(branch, "BRANCH_SEQ_RLINZ"),
                    "x_ohm_per_km": self._hadapp.getBranchInfoReal(branch, "BRANCH_SEQ_XLINZ"),
                    "max_i_ka": 0.0, 
                    # FIXME: Need in service value
                    "in_service": True # self._hadapp.getBranchInfoBool(branch, "BRANCH_STATUS")
                })

        # convert to dataframes
        # bus and its results
        grid.bus = pd.DataFrame(bus_list)
        grid.res_bus = pd.DataFrame(index=grid.bus.index, columns=["vm_pu"])
        
        # gen and its results
        grid.gen = pd.DataFrame(gen_list)
        grid.res_gen = pd.DataFrame(index=grid.gen.index, columns=["p_mw", "q_mvar", "vm_pu"])
        
        # load and its results
        grid.load = pd.DataFrame(load_list)
        grid.res_load = pd.DataFrame(index=grid.load.index, columns=["p_mw", "q_mvar"])
        
        # line and its results
        # NOTE: bus dict is needed to translate actual bus number to dataframe index. The from_bus and to_bus columns in grid.line needs to be replaced with this new index to be consistent across the code. 
        bus_dict = grid.bus[["id"]].reset_index().set_index("id").to_dict()['index']
        grid.line = pd.DataFrame(branch_list)
        grid.line.from_bus = [bus_dict[b] for b in grid.line.from_bus]
        grid.line.to_bus = [bus_dict[b] for b in grid.line.to_bus]
        grid.res_line = pd.DataFrame(index=grid.line.index, columns=["p_from_mw", "p_to_mw", "q_from_mvar", "q_to_mvar", "vm_from_pu", "vm_to_pu", "i_from_ka", "i_to_ka"])
        
        # transformers and its results
        grid.trafo = pd.DataFrame(columns=["name", "from_bus", "to_bus", "hv_bus", "lv_bus", "in_service"])
        grid.res_trafo = pd.DataFrame(index=grid.trafo.index, columns=["p_hv_mw", "p_lv_mw", "q_hv_mvar", "q_lv_mvar", "vm_hv_pu", "vm_lv_pu", "i_hv_ka", "i_lv_ka"])
        
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
        full_path = path
        if filename is not None:
            full_path = os.path.join(full_path, filename)
        
        # solve the power flow - to load the grid data
        self._hadapp.solvePowerFlowBeforeDynSimu(full_path, -1)  # 0 inidcates that solves the first raw file for power flow, the xml file supports multiple power flow raw files read in
        self._hadapp.readGenerators();
        self._hadapp.readSequenceData();
        self._hadapp.initialize();
        self._hadapp.setGeneratorWatch();
        
        # Building grid object from hadapp
        self._grid = self.build_grid()
        
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
        self.thermal_limit_a = np.ones(self.n_line, dtype=int)
            
        self._compute_pos_big_topo()

        # Remember the input file was read into the Configuration singleton
        conf = gridpack.Configuration()
        cursor = conf.getCursor("Configuration.Dynamic_simulation")

        # get faults
        faults = self._hadapp.getEvents(cursor)

        # solve with faults
        self._hadapp.solvePreInitialize(faults[0])

    def apply_action(self, backendAction: Union["grid2op.Action._backendAction._BackendAction", None]) -> None:
        '''
        # called for each "step", thousands of times
        # modify the topology, load, generation etc.
        '''
        
        # the following few lines are highly recommended
        if backendAction is None:
            print("===================================")
            return
        
        (
            active_bus,
            (prod_p, prod_v, load_p, load_q, storage),
            _,
            shunts__,
        ) = backendAction()
        
        # change the active values of the loads
        for load_id, new_p in load_p:
            self._grid.load["p_mw"].iloc[load_id] = new_p
        # change the reactive values of the loads
        for load_id, new_q in load_q:
            self._grid.load["q_mvar"].iloc[load_id] = new_q

        # change the active value of generators
        for gen_id, new_p in prod_p:
            self._grid.gen["p_mw"].iloc[gen_id] = new_p
            
        # for the voltage magnitude, pandapower expects pu but grid2op provides kV,
        # so we need a bit of change
        for gen_id, new_v in prod_v:
            self._grid.gen["vm_pu"].iloc[gen_id] = new_v  # but new_v is not pu !
            self._grid.gen["vm_pu"].iloc[gen_id] /= self._grid.bus["vn_kv"][
                self.gen_to_subid[gen_id]
            ]  # now it is :-)

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
        
        n_line_pp = self._grid.line.shape[0]
        
        # handle the disconnection on "or" side
        lines_or_bus = backendAction.get_lines_or_bus()
        for line_id, new_bus in lines_or_bus:
            if line_id < n_line_pp:
                # a pandapower powerline has bee disconnected in grid2op
                dt = self._grid.line
                line_id_db = line_id
            else:
                # a pandapower trafo has bee disconnected in grid2op
                dt = self._grid.trafo
                line_id_db = line_id - n_line_pp

            if new_bus == -1:
                # element was disconnected
                dt.iloc[line_id_db]["in_service"] = False
            else:
                # element was connected
                dt.iloc[line_id_db]["in_service"] = True

        lines_ex_bus = backendAction.get_lines_ex_bus()
        for line_id, new_bus in lines_ex_bus:
            if line_id < n_line_pp:
                # a pandapower powerline has bee disconnected in grid2op
                dt = self._grid.line
                line_id_db = line_id
            else:
                # a pandapower trafo has bee disconnected in grid2op
                dt = self._grid.trafo
                line_id_db = line_id - n_line_pp

            if new_bus == -1:
                # element was disconnected
                dt.iloc[line_id_db]["in_service"] = False
            else:
                # element was connected
                dt.iloc[line_id_db]["in_service"] = True

    def _random_data_generator(self, n_rows, columns):
        if n_rows > 0:
            new_data = np.random.rand(n_rows, len(columns))
            return pd.DataFrame(new_data, columns=columns, dtype="float32")
        else:
            return pd.DataFrame(columns=columns)

    def _update_gen_load_data(self):
        # load and generator data
        # NOTE: Generator and Load values change only when there are dynamic load/generators. If they are not dynamic, the code currently returns the previous value. 
        gen_data = []
        load_data = []
        for bus in range(self.n_sub):
            # Dynamic simulation values are in pu, need to convert them MW and MVar
            CASE_SBASE = 1 # self._hadapp.getBusInfoReal(bus, "CASE_SBASE")

            # bus list
            for g in range(self._hadapp.numGenerators(bus)):
                # gen data
                gen_data.append({
                    "p_mw": self._hadapp.getBusInfoReal(bus, 'GENERATOR_PG_CURRENT', g) * CASE_SBASE,
                    "q_mvar": self._hadapp.getBusInfoReal(bus, 'GENERATOR_QG_CURRENT', g) * CASE_SBASE,
                    "vm_pu": self._hadapp.getBusInfoReal(bus, 'BUS_VOLTAGE_MAG')
                })
            for l in range(self._hadapp.numLoads(bus)):
                # load list
                load_data.append({
                    "p_mw": self._hadapp.getBusInfoReal(bus, 'LOAD_PL_CURRENT', l) * CASE_SBASE / 100,
                    "q_mvar": self._hadapp.getBusInfoReal(bus, 'LOAD_QL_CURRENT', l) * CASE_SBASE / 100
                })
            
        self._grid.res_load = pd.DataFrame(load_data)
        self._grid.res_gen = pd.DataFrame(gen_data) # p and q also gets multiplied by 100 - CASE_SBASE

        print(self._grid.load)
        print(self._grid.res_load)
        sys.exit(1)
    
    def _update_line_transformer_data(self):
        # line level data - lines and branches could be different. The two end points can have multiple lines but only one branch. 
        nbranch = self._hadapp.totalBranches()
        line_data = []
        # 
        for branch in range(0, nbranch):
            n_elements = self._hadapp.getBranchInfoInt(branch, 'BRANCH_NUM_ELEMENTS')
            for elem_num in range(n_elements):
                # empty object
                line_dict = dict()
                line_dict["branch_num"] = branch
                line_dict["line_num"] = elem_num

                # end points
                (f, t) = self._hadapp.getBranchEndpoints(branch)
                line_dict["from_bus"] = f
                line_dict["to_bus"] = t

                # from and to bus kv
                f_buskv = 1 # self._hadapp.getBusInfoReal(bus, "BUS_KV") # "BUS_KV"
                t_buskv = 1 # self._hadapp.getBusInfoReal(bus, "BUS_KV") # "BUS_KV"
                
                # from voltage
                from_bus = self._grid.bus.index[self._grid.bus['id'] == f].values[0]
                line_dict["vm_from_pu"] = self._hadapp.getBusInfoReal(from_bus, 'BUS_VOLTAGE_MAG') * f_buskv
                
                # to voltage
                to_bus = self._grid.bus.index[self._grid.bus['id'] == t].values[0]
                line_dict["vm_to_pu"] = self._hadapp.getBusInfoReal(to_bus, 'BUS_VOLTAGE_MAG') * t_buskv

                # p-values
                line_dict["p_from_mw"] = self._hadapp.getBranchInfoReal(branch, 'BRANCH_FROM_P_CURRENT', elem_num)
                line_dict["p_to_mw"] = self._hadapp.getBranchInfoReal(branch, 'BRANCH_TO_P_CURRENT', elem_num)

                # q-values
                line_dict["q_from_mvar"] = self._hadapp.getBranchInfoReal(branch, 'BRANCH_FROM_Q_CURRENT', elem_num)
                line_dict["q_to_mvar"] = self._hadapp.getBranchInfoReal(branch, 'BRANCH_TO_Q_CURRENT', elem_num)

                # current values
                line_dict["i_from_ka"] = self._hadapp.getBranchInfoReal(branch, 'BRANCH_IRFLOW_CURRENT', elem_num)
                line_dict["i_to_ka"] = self._hadapp.getBranchInfoReal(branch, 'BRANCH_IRFLOW_CURRENT', elem_num) # output in amps

                line_data.append(line_dict)
        
        # res line data using random data generator
        # self._grid.res_line = self._random_data_generator(self._grid.line.shape[0], columns=self._grid.res_line.columns)
        # self._grid.res_line["vm_from_pu"] = vm_from_pu
        # self._grid.res_line["vm_to_pu"] = vm_to_pu # * Bus_basekV
        # print(self._grid.res_line)
        self._grid.res_line = pd.DataFrame(line_data)
        # print(self._grid.line)
        # print(self._grid.res_line)
        # sys.exit(1)
        
        # transformer data
        self._grid.res_trafo = self._random_data_generator(self._grid.trafo.shape[0], columns=self._grid.res_trafo.columns)
    
    def _update_data(self):
        # to get current data from dynamic simulation to data collection object
        self._hadapp.updateData()
        
        # counter
        self._counter += 1
        print(f"Counter value: {self._counter}")

        # update data for loads and generators
        self._update_gen_load_data()

        # update line and transformer data
        self._update_line_transformer_data()

        # print(self._grid.res_load)
        # print(self._grid.res_gen)
        print(self._grid.res_line[["p_from_mw", "p_to_mw", "q_from_mvar", "q_to_mvar"]])
        # sys.exit(1)

    
    def runpf(self, is_dc : bool=False):
        '''
        # called for each "step", thousands of times
        # run the solver
        '''
        # execute one simulation time step	
        self._hadapp.executeOneSimuStep()  

        # update data in result dataframes
        self._update_data()

        # print(
        #     self._grid.res_gen, "\n",
        #     self._grid.res_load, "\n",
        #     self._grid.res_line, "\n",
        #     self._grid.res_trafo
        # )
        
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
        # FIXME: Placeholder values
        # return np.ones(self.n_load), np.ones(self.n_load), np.ones(self.n_load)
    
        # carefull with copy / deep copy
        load_p = self._grid.res_load["p_mw"].values  # in MW
        load_q = self._grid.res_load["q_mvar"].values  # in MVAr
        
        # load_v is the voltage magnitude at the bus at which the load is connected.
        # in pandapower this is not straightforward. We first need to retrieve the
        # voltage in per unit of the bus to which each load is connected.
        # And then we convert the pu to kV. This is what is done below.
        load_v = self._grid.res_bus.iloc[self._grid.load["bus"].values]["vm_pu"].values  # in pu
        load_v *= self._grid.bus.iloc[self._grid.load["bus"].values]["vn_kv"].values # in kv
        return load_p, load_q, load_v

    def generators_info(self)-> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        '''
        # retrieve the results
        '''
        prod_p = self._grid.res_gen["p_mw"].values  # in MW
        prod_q = self._grid.res_gen["q_mvar"].values  # in MVAr
        
        # same as for load, gen_v is not directly accessible in pandapower
        # we first retrieve the per unit voltage, then convert it to kV
        prod_v = self._grid.res_gen["vm_pu"].values
        prod_v *= (
            self._grid.bus["vn_kv"].iloc[self.gen_to_subid].values
        )  # in kV
        return prod_p, prod_q, prod_v
        # return np.ones(self.n_gen), np.ones(self.n_gen), np.ones(self.n_gen)

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
        v_or = self._aux_get_line_info("vm_from_pu", "vm_hv_pu")  # in pu
        a_or = self._aux_get_line_info("i_from_ka", "i_hv_ka") * 1000  # grid2op expects amps (A) pandapower returns kilo-amps (kA)
        # * 100 / kV (TAmps)
        
        # get the voltage in kV (and not in pu)
        bus_id = np.concatenate(
            (
                self._grid.line["from_bus"].values,
                self._grid.trafo["hv_bus"].values,
            )
        )
        v_or *= self._grid.bus.iloc[bus_id]["vn_kv"].values
        
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
        v_ex = self._aux_get_line_info("vm_to_pu", "vm_lv_pu")  # in pu
        a_ex = self._aux_get_line_info("i_to_ka", "i_lv_ka") * 1000  # grid2op expects amps (A) pandapower returns kilo-amps (kA)
        
        # get the voltage in kV (and not in pu)
        bus_id = np.concatenate(
            (
                self._grid.line["to_bus"].values,
                self._grid.trafo["lv_bus"].values,
            )
        )
        v_ex *= self._grid.bus.iloc[bus_id]["vn_kv"].values
        
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
        self._hadapp = None
        self.env = None