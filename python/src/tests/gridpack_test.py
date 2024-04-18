#!/usr/bin/env python
# -------------------------------------------------------------
# file: hello.py
# -------------------------------------------------------------
# -------------------------------------------------------------
# Copyright (c) 2013 Battelle Memorial Institute
# Licensed under modified BSD License. A copy of this license can be
# found in the LICENSE file in the top level directory of this
# distribution.
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created January 27, 2020 by Perkins
# Last Change: 2024-05-09 10:30:03 d3g096
# -------------------------------------------------------------

import sys, os, time
from unittest import TestCase
import gridpack
import gridpack.hadrec


# -------------------------------------------------------------
# GridPACKTester
# -------------------------------------------------------------
class GridPACKTester(TestCase):
    def hello_test(self):
        c = gridpack.Communicator()
        sys.stdout.write("hello from process %d of %d\n" %
                         (c.rank(), c.size()))
    def hello_comm_test(self):
        c = gridpack.Communicator()
        sys.stdout.write("hello from process %d of %d\n" %
                         (c.rank(), c.size()))
    def task_test(self):
        c = gridpack.Communicator()
        tskmgr = gridpack.TaskManager(c)
        
        task = gridpack.TaskCounter()

        tskmgr.set(100)
        while tskmgr.nextTask(task):
            sys.stdout.write("process %d of %d executing task %d\n" %
                             (c.rank(), c.size(), task.task_id))
        tskmgr = None
    def timer_test(self):
        timer = gridpack.CoarseTimer()
        c1 = timer.createCategory("Category 1")
        c2 = timer.createCategory("Category 2")
        timer.start(c1)
        time.sleep(5)
        timer.start(c2)
        time.sleep(5)
        timer.stop(c2)
        timer.stop(c1)
        timer.dump()

    def configure_test(self):
        comm = gridpack.Communicator()
        conf = gridpack.Configuration()

        print("gridpack.Configuration.KeySep = \"%s\"" % 
              (gridpack.Configuration.KeySep))

        print("conf.KeySep = \"%s\"" % (conf.KeySep))

        self.assertEqual('.', conf.KeySep)

        d = os.path.dirname(os.path.abspath(__file__))
        os.chdir(d)
        
        conf.open("input_tamu500_step005.xml", comm)

        path = ("Configuration%cDynamic_simulation%stimeStep" %
                ( gridpack.Configuration.KeySep,
                  gridpack.Configuration.KeySep))
        dt = conf.get(path)
        self.assertFalse(dt is None)
        self.assertEqual(float(dt), 0.005)
        print("dt = %f" % (float(dt)))

        path = ("Configuration%cDynamic_simulation" %
                ( gridpack.Configuration.KeySep))

        cursor = conf.getCursor(path)
        dt = cursor.get("timeStep")

        self.assertFalse(dt is None)
        self.assertEqual(float(dt), 0.005)
        print("dt = %f" % (float(dt)))

        cursor2 = cursor.getCursor("observations")
        self.assertFalse(cursor2 is None)
        
    # def hadrec_test(self):

        print("Number of buses:  %d" % (hadapp.totalBuses()))
        print("Number of branches: %d" % (hadapp.totalBranches()))
        print("Number of generators: %d" % (hadapp.numGenerators()))
        print("Number of loads: %d" % (hadapp.numLoads()))
        print("Number of storage units: %d" % (hadapp.numStorage()))
        print("Branches connected to bus 1: ", hadapp.getConnectedBranches(1))
        print("Buses connected to branch 1: ", hadapp.getBranchEndpoints(1))

        busfaultlist = gridpack.dynamic_simulation.EventVector()

    #     arg = "input_tamu500_step005.xml"

    #     np = gridpack.NoPrint()
    #     sys.stdout.write("NoPrint status: %r\n" % (np.status()))
    #     np.setStatus(True)
    #     sys.stdout.write("NoPrint status: %r\n" % (np.status()))
        
    #     hadapp = gridpack.hadrec.Module()
    #     hadapp.solvePowerFlowBeforeDynSimu(arg)
    #     hadapp.transferPFtoDS()

    #     busfaultlist = gridpack.dynamic_simulation.EventVector()

    #     hadapp.initializeDynSimu(busfaultlist)

    #     print (hadapp.getBusTotalLoadPower(5))
    #     print (hadapp.getBusTotalLoadPower(7))
    #     print (hadapp.getBusTotalLoadPower(9))
    #     print (hadapp.getBusTotalLoadPower(-1)) # should be None

    #     print (hadapp.getGeneratorPower(1, "1"))
    #     print (hadapp.getGeneratorPower(2, "1"))
    #     print (hadapp.getGeneratorPower(3, "1"))
    #     print (hadapp.getGeneratorPower(3, "14")) # should be None
    #     print (hadapp.getGeneratorPower(-1, "1")) # should be None

    #     # We need some actual tests of get/setState() that work
    #     # print (hadapp.getState(1, "1", "Device", "ANGLE"))
    #     # print (hadapp.setState(1, "1", "Device", "ANGLE", 0.0))

    #     loadshedact = gridpack.hadrec.Action()
    #     loadshedact.actiontype = 0
    #     loadshedact.bus_number = 5
    #     loadshedact.componentID = "1"
    #     loadshedact.percentage = -0.2

    #     loadshedact1 = gridpack.hadrec.Action()
    #     loadshedact1.actiontype = 0
    #     loadshedact1.bus_number = 7
    #     loadshedact1.componentID = "1"
    #     loadshedact1.percentage = -0.2

    #     linetrip = gridpack.hadrec.Action()
    #     linetrip.actiontype = 1;
    #     linetrip.brch_from_bus_number = 6;
    #     linetrip.brch_to_bus_number = 7;
    #     linetrip.branch_ckt = "1 ";

    #     (obs_genBus, obs_genIDs, obs_loadBuses, obs_loadIDs, obs_busIDs) = hadapp.getObservationLists()

    #     print (obs_genBus)
    #     print (obs_genIDs)
    #     print (obs_loadBuses)
    #     print (obs_loadIDs)
    #     print (obs_busIDs)

    #     bApplyAct_LoadShedding = False
    #     bApplyAct_LineTripping = True
    #     isteps = 0

    #     while (not hadapp.isDynSimuDone()):
    #         if (bApplyAct_LoadShedding and
    #             ( isteps == 2500 or isteps == 3000 or
    #               isteps == 3500 or isteps == 4000 or
    #               isteps == 4000 or isteps == 4500 or
    #               isteps == 5000 or isteps == 5500 or
    #               isteps == 5500)):
    #             hadapp.applyAction(loadshedact)
    #             hadapp.applyAction(loadshedact1)

    #         if (bApplyAct_LineTripping and isteps == 400):
    #             hadapp.applyAction(linetrip)
                
    #         hadapp.executeDynSimuOneStep()
    #         ob_vals = hadapp.getObservations()
    #         print (isteps, ob_vals)
    #         isteps = isteps + 1

    #     btest_2dynasimu = True

    #     if (btest_2dynasimu):

    #         busfault = gridpack.dynamic_simulation.Event()
    #         busfault.start = 10.0
    #         busfault.end = 10.2
    #         busfault.step = 0.005
    #         busfault.isBus = True
    #         busfault.bus_idx = 7

    #         busfaultlist = gridpack.dynamic_simulation.EventVector([busfault])

    #         hadapp.transferPFtoDS()
    #         hadapp.initializeDynSimu(busfaultlist)

    #         while (not hadapp.isDynSimuDone()):
    #             if (bApplyAct_LoadShedding and
    #                 ( isteps == 2500 or isteps == 3000 or
    #                   isteps == 3500 or isteps == 4000 )):
    #                 hadapp.applyAction(loadshedact)
    #             hadapp.executeDynSimuOneStep()
    #             ob_vals = hadapp.getObservations()
    #             print (isteps, ob_vals)
    #             isteps = isteps + 1

    #     # See if we could do it again, if we wanted to
    #     hadapp = None
    #     print ('----renke python debug test---before second time ini hadrec module')
    #     hadapp = gridpack.hadrec.Module()
    #     hadapp.solvePowerFlowBeforeDynSimu(arg) 
    #     # It's important to force deallocation order here
    #     hadapp = None
    #     env = None
        
