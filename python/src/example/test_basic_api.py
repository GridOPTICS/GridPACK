import sys
import numpy as np
sys.path.append("/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/")
from grid2op_backend import GridPACKBackend

import unittest
import warnings
from grid2op.tests.aaa_test_backend_interface import AAATestBackendAPI
FILE_FORMAT = "xml"  # for this example, put whatever here !

class TestBackendAPI_GridPACKBackend(AAATestBackendAPI, unittest.TestCase):
    def get_path(self):
        return "./"  # or /home/user/Documents/my_path_for_test

    def get_casefile(self):
        # return "input_39bus_step005_v33_itr.xml"
        # return "input_9b3g.xml"
        return "grid.xml"
        #    # or `grid.xml` or any other format

    def make_backend(self, detailed_infos_for_cascading_failures=False):
        # the function that will create your backend
        # do not put "GridPACKBackend" of course, but the class you coded as a backend !
        backend = GridPACKBackend()
        # assert FILE_FORMAT in backend.supported_grid_format, f"your backend does not recognize the '{FILE_FORMAT}' extension, grid2op will not work"
        return backend

    def test_01load_grid(self):
        """Tests the runpf method (AC) without modification
    
        This test supposes that :
        
        - backend.load_grid(...) is implemented        
        - backend.runpf() (DC mode) is implemented
                
        """
        pass

    def test_02modify_load(self):
        """Tests the loads can be modified        

        This test supposes that :
        
        - backend.load_grid(...) is implemented        
        - backend.apply_action(...) for modification of loads is implemented        

        NB: it does not check whether or not the modification is
        consistent with the input. This will be done in a later test"""
        self.skip_if_needed()
        backend = self.aux_make_backend()
        
        np.random.seed(0)
        random_load_p = np.random.uniform(0, 1, size=type(backend).n_load)
        random_load_q = np.random.uniform(0, 1, size=type(backend).n_load)
        
        # try to modify load_p
        action = type(backend)._complete_action_class()
        action.update({"injection": {"load_p": 1.01 * random_load_p}})
        bk_act = type(backend).my_bk_act_class()
        bk_act += action
        backend.apply_action(bk_act)  # modification of load_p only
        
        # try to modify load_q
        action = type(backend)._complete_action_class()
        action.update({"injection": {"load_q": 1.01 * random_load_q}})
        bk_act = type(backend).my_bk_act_class()
        bk_act += action
        backend.apply_action(bk_act)  # modification of load_q only

        print("================== Test Successful!! =====================")

    def test_03modify_gen(self):
        """Tests the generators (including slack !) can be modified 

        This test supposes that :
        
        - backend.load_grid(...) is implemented        
        - backend.apply_action(...) for modification of generators is implemented
                
        NB: it does not check whether or not the modification is
        consistent with the input. This will be done in a later test"""
        self.skip_if_needed()
        print("================== Test-3 Starting!! =====================")
        backend = self.aux_make_backend()
        np.random.seed(0)
        random_gen_p = np.random.uniform(0, 1, size=type(backend).n_gen)
        random_gen_v = np.random.uniform(0, 1, size=type(backend).n_gen)
        
        # try to modify gen_p
        action = type(backend)._complete_action_class()
        action.update({"injection": {"prod_p": 1.01 * random_gen_p}})
        bk_act = type(backend).my_bk_act_class()
        bk_act += action
        backend.apply_action(bk_act)  # modification of prod_p / gen_p only
        
        # try to modify prod_v only
        action = type(backend)._complete_action_class()
        action.update({"injection": {"prod_v": random_gen_v + 0.1}})
        bk_act = type(backend).my_bk_act_class()
        bk_act += action
        backend.apply_action(bk_act)  # modification of prod_v / gen_v only

        backend.close()

        print("================== Test-3 Successful!! =====================")


if __name__ == "__main__":
    unittest.main()