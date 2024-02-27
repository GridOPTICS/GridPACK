import os
import numpy as np
from typing import Optional, Tuple, Union

import grid2op
from grid2op.Backend import Backend   # required

import gridpack

class GridPACKBackend(Backend):
    def __init__(self) -> None:
        print ('before gridpack ini')
        noprintflag = gridpack.NoPrint()
        noprintflag.setStatus(False) # Setting this to True will disable all printing to stdout

        # Create GridPACK environment and pass the communicator to it
        env = gridpack.Environment()

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
        self._hadapp.solvePowerFlowBeforeDynSimu(full_path, 0)  # 0 inidcates that solves the first raw file for power flow, the xml file supports multiple power flow raw files read in
        print ('Finished Solving Power Flow')
        
        # FIXME: Get this grid object from hadapp
        # self._grid = pp.from_json(full_path)
        
        # then fill the "n_sub" and "sub_info"
        # FIXME: Missing info from gridpack
        # self.n_sub = self._grid.bus.shape[0]
        
        # then fill the number and location of loads
        # FIXME: Missing info from gridpack
        # self.n_load = self._grid.load.shape[0]
        self.load_to_subid = np.zeros(self.n_load, dtype=int)
        for load_id in range(self.n_load):
            pass
            # FIXME: Missing info from gridpack
            # self.load_to_subid[load_id] = self._grid.load.iloc[load_id]["bus"]
            
        # then fill the number and location of generators
        # FIXME: Missing info from gridpack
        # self.n_gen = self._grid.gen.shape[0]
        self.gen_to_subid = np.zeros(self.n_gen, dtype=int)
        for gen_id in range(self.n_gen):
            pass
            # FIXME: Missing info from gridpack
            # self.gen_to_subid[gen_id] = self._grid.gen.iloc[gen_id]["bus"]
            
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
        # FIXME: Missing info from gridpack
        # self.n_line = self._grid.line.shape[0] + self._grid.trafo.shape[0]
        self.line_or_to_subid = np.zeros(self.n_line, dtype=int)
        self.line_ex_to_subid = np.zeros(self.n_line, dtype=int)
        for line_id in range(self._grid.line.shape[0]):
            pass
            # FIXME: Missing info from gridpack
            # self.line_or_to_subid[line_id] = self._grid.line.iloc[line_id]["from_bus"]
            # self.line_ex_to_subid[line_id] = self._grid.line.iloc[line_id]["to_bus"]
        
        # FIXME: Missing info from gridpack
        # nb_powerline = self._grid.line.shape[0]
        for trafo_id in range(self._grid.trafo.shape[0]):
            pass
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
            
        self._compute_pos_big_topo()

        # transfer data from power flow network to dynamic simulation network
        self._hadapp.transferPFtoDS()

    def apply_action(self, backendAction: Union["grid2op.Action._backendAction._BackendAction", None]) -> None:
        '''
        # called for each "step", thousands of times
        # modify the topology, load, generation etc.
        '''
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

    def runpf(self):
        '''
        # called for each "step", thousands of times
        # run the solver
        '''
        # execute one simulation time step	
        self._hadapp.executeDynSimuOneStep()  

        return True     
    

    def get_topo_vect(self):
        '''
        # retrieve the results
        '''
        pass

    def loads_info(self)-> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        '''
        # retrieve the results
        '''
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