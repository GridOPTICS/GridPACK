import os
import numpy as np
from typing import Optional, Tuple, Union

import grid2op
from grid2op.Backend import Backend   # required

import gridpack
from mpi4py import MPI

class GridPACKBackend(Backend):
    def __init__(self) -> None:
        # Run Backend init
        super().__init__()

        # Get World communicator
        comm = MPI.COMM_WORLD
        
        print ('before gridpack ini')
        noprintflag = gridpack.NoPrint()
        noprintflag.setStatus(False) # Setting this to True will disable all printing to stdout

        # Create GridPACK environment and pass the communicator to it
        self.env = gridpack.Environment()
        
        # Create hadrec module
        self._hadapp = gridpack.hadrec.Module()

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
        print('=========== Before Solving Power Flow ===================')
        print(full_path)
        self._hadapp.solvePowerFlowBeforeDynSimu(full_path, 0)  # 0 inidcates that solves the first raw file for power flow, the xml file supports multiple power flow raw files read in
        print ('Finished Solving Power Flow')
        
        # FIXME: Get this grid object from hadapp
        # self._grid = pp.from_json(full_path)
        
        # then fill the "n_sub" and "sub_info"
        self.n_sub = self._hadapp.totalBuses()
        for bus in range(self.n_sub):
            print(bus,
                # self._hadapp.getBusInfoInt(bus, "BUS_NUMBER", -1),
                # self._hadapp.getBusInfoString(bus, "BUS_NAME", -1),
                # self._hadapp.getBusInfoInt(bus, "BUS_TYPE", -1),
                self._hadapp.numGenerators(bus),
                self._hadapp.numLoads(bus))#,
                # self._hadapp.getBusInfoReal(bus, "BUS_VOLTAGE_MAG", -1))
            for g in range(self._hadapp.numGenerators(bus)):
                print(" gen: ", g)#,
                    # self._hadapp.getBusInfoInt(bus, "GENERATOR_NUMBER", g),
                    # self._hadapp.getBusInfoString(bus, "GENERATOR_ID", g),
                    # self._hadapp.getBusInfoReal(bus, "GENERATOR_PG", g),
                    # self._hadapp.getBusInfoReal(bus, "GENERATOR_QG", g))
            for l in range(self._hadapp.numLoads(bus)):
                print("load: ", l)#,
                    # self._hadapp.getBusInfoInt(bus, "LOAD_NUMBER", l),
                    # self._hadapp.getBusInfoString(bus, "LOAD_ID", l),
                    # self._hadapp.getBusInfoReal(bus, "LOAD_PL", l),
                    # self._hadapp.getBusInfoReal(bus, "LOAD_QL", l))
        
        # then fill the number and location of loads
        self.n_load = self._hadapp.numLoads()
        self.load_to_subid = np.zeros(self.n_load, dtype=int)
        for load_id in range(self.n_load):
            # pass
            # FIXME: Missing info from gridpack - check me
            self.load_to_subid[load_id] = self._hadapp.numLoads(load_id)
            
        # then fill the number and location of generators
        # FIXME: Missing info from gridpack
        self.n_gen = self._hadapp.numGenerators()
        self.gen_to_subid = np.zeros(self.n_gen, dtype=int)
        for gen_id in range(self.n_gen):
            # pass
            # FIXME: Missing info from gridpack - check me
            self.gen_to_subid[gen_id] = self._hadapp.numGenerators(gen_id)
            
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
        self.n_line = self._hadapp.numLines() # self._grid.line.shape[0] + self._grid.trafo.shape[0]
        self.line_or_to_subid = np.zeros(self.n_line, dtype=int)
        self.line_ex_to_subid = np.zeros(self.n_line, dtype=int)
        # for line_id in range(self._grid.line.shape[0]):
        #     pass
            # FIXME: Missing info from gridpack
            # self.line_or_to_subid[line_id] = self._grid.line.iloc[line_id]["from_bus"]
            # self.line_ex_to_subid[line_id] = self._grid.line.iloc[line_id]["to_bus"]
        
        # FIXME: Missing info from gridpack
        # nb_powerline = self._grid.line.shape[0]
        # for trafo_id in range(self._grid.trafo.shape[0]):
        #     pass
            # FIXME: Missing info from gridpack
            # self.line_or_to_subid[trafo_id + nb_powerline] = self._grid.trafo.iloc[trafo_id]["hv_bus"]
            # self.line_ex_to_subid[trafo_id + nb_powerline] = self._grid.trafo.iloc[trafo_id]["lv_bus"]
            
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

        # transfer data from power flow network to dynamic simulation network
        self._hadapp.transferPFtoDS()

        # define a bus fault
        busfault = gridpack.dynamic_simulation.Event()
        busfault.start = 1.0 # fault start time
        busfault.end = 1.1   # fault end time
        busfault.step = 0.005  # fault duration simu time step
        busfault.isBus = True
        busfault.bus_idx = 22 # bus number of the fault		  

        busfaultlist = gridpack.dynamic_simulation.EventVector([busfault])

        # initialize the dynamic simulation
        self._hadapp.initializeDynSimu(busfaultlist, 0) # 0 inidcates read in the first dyr dynamic parameter file, the xml file supports multiple dyr files read in

    def apply_action(self, backendAction: Union["grid2op.Action._backendAction._BackendAction", None]) -> None:
        '''
        # called for each "step", thousands of times
        # modify the topology, load, generation etc.
        '''
        return
        # the following few lines are highly recommended
        if backendAction is None:
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
                dt["in_service"].iloc[line_id_db] = False
            else:
                # element was connected
                dt["in_service"].iloc[line_id_db] = True

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
                dt["in_service"].iloc[line_id_db] = False
            else:
                # element was connected
                dt["in_service"].iloc[line_id_db] = True

    def runpf(self, is_dc : bool=False):
        '''
        # called for each "step", thousands of times
        # run the solver
        '''
        # execute one simulation time step	
        self._hadapp.executeDynSimuOneStep()  
        
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
        # # read results for load
        # self._aux_get_topo_vect(res, self._grid.load, "bus", self.load_pos_topo_vect)
        # # then for generators
        # self._aux_get_topo_vect(res, self._grid.gen, "bus", self.gen_pos_topo_vect)
        # # then each side of powerlines
        # self._aux_get_topo_vect(res, self._grid.line, "from_bus", self.line_or_pos_topo_vect)
        # self._aux_get_topo_vect(res, self._grid.line, "to_bus", self.line_ex_pos_topo_vect)
        
        # # then for the trafos, but remember pandapower trafos are powerlines in grid2Op....
        # # so we need to trick it a bit
        # # (we can do this trick because we put the trafo "at the end" of the powerline in grid2op
        # # in the Step1_loading.py)
        # n_line_pp = self._grid.line.shape[0]
        # self._aux_get_topo_vect(res, self._grid.trafo, "hv_bus", self.line_or_pos_topo_vect, add_id=n_line_pp)
        # self._aux_get_topo_vect(res, self._grid.trafo, "lv_bus", self.line_ex_pos_topo_vect, add_id=n_line_pp)            
        return res

    def loads_info(self)-> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        '''
        # retrieve the results
        '''
        # FIXME: Placeholder values
        return np.ones(self.n_load), np.ones(self.n_load), np.ones(self.n_load)
    
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
        # prod_p = self._grid.res_gen["p_mw"].values  # in MW
        # prod_q = self._grid.res_gen["q_mvar"].values  # in MVAr
        
        # # same as for load, gen_v is not directly accessible in pandapower
        # # we first retrieve the per unit voltage, then convert it to kV
        # prod_v = self._grid.res_gen["vm_pu"].values
        # prod_v *= (
        #     self._grid.bus["vn_kv"].iloc[self.gen_to_subid].values
        # )  # in kV
        # return prod_p, prod_q, prod_v
        return np.ones(self.n_gen), np.ones(self.n_gen), np.ones(self.n_gen)

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
        return res
    
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
        return np.ones(self.n_line), np.ones(self.n_line), np.ones(self.n_line), np.ones(self.n_line)
    
        p_or = self._aux_get_line_info("p_from_mw", "p_hv_mw")
        q_or = self._aux_get_line_info("q_from_mvar", "q_hv_mvar")
        v_or = self._aux_get_line_info("vm_from_pu", "vm_hv_pu")  # in pu
        a_or = self._aux_get_line_info("i_from_ka", "i_hv_ka") * 1000  # grid2op expects amps (A) pandapower returns kilo-amps (kA)
        
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
        v_or[~status] = 0.
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
        return np.ones(self.n_line), np.ones(self.n_line), np.ones(self.n_line), np.ones(self.n_line)
    
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
        v_ex[~status] = 0.
        return p_ex, q_ex, v_ex, a_ex
    
    def close(self):
        self._hadapp = None
        self.env = None