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
# Last Change: 2020-04-23 14:01:55 d3g096
# -------------------------------------------------------------

import sys, os
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
    def task_test(self):
        c = gridpack.Communicator()
        tskmgr = gridpack.TaskManager(c)
        
        task = gridpack.TaskCounter()

        tskmgr.set(100)
        while tskmgr.nextTask(task):
            sys.stdout.write("process %d of %d executing task %d\n" %
                             (c.rank(), c.size(), task.task_id))
        tskmgr = None
    def hadrec_test(self):

        d = os.path.dirname(os.path.abspath(__file__))
        os.chdir(d)

        arg = "input_tamu500_step005.xml"
        
        hadapp = gridpack.hadrec.Module()
        hadapp.solvePowerFlowBeforeDynSimu(arg)
        hadapp.transferPFtoDS()

        busfaultlist = gridpack.dynamic_simulation.EventVector()

        hadapp.initializeDynSimu(busfaultlist)

        loadshedact = gridpack.hadrec.Action()
        loadshedact.actiontype = 0;
        loadshedact.bus_number = 5;
        loadshedact.componentID = "1";
        loadshedact.percentage = -0.2;

        loadshedact1 = gridpack.hadrec.Action()
        loadshedact1.actiontype = 0;
        loadshedact1.bus_number = 7;
        loadshedact1.componentID = "1";
        loadshedact1.percentage = -0.2;

        (obs_genBus, obs_genIDs, obs_loadBuses, obs_loadIDs, obs_busIDs) = hadapp.getObservationLists()

        print (obs_genBus)
        print (obs_genIDs)
        print (obs_loadBuses)
        print (obs_loadIDs)
        print (obs_busIDs)

        bApplyAct = True
        isteps = 0

        while (not hadapp.isDynSimuDone()):
            if (bApplyAct and
                ( isteps == 2500 or isteps == 3000 or
                  isteps == 3500 or isteps == 4000 )):
                hadapp.applyAction(loadshedact)
                hadapp.applyAction(loadshedact1)
            hadapp.executeDynSimuOneStep()
            ob_vals = hadapp.getObservations();
            print (isteps, ob_vals)
            isteps = isteps + 1

        btest_2dynasimu = True

        if (btest_2dynasimu):

            busfault = gridpack.dynamic_simulation.Event()
            busfault.start = 10.0;
            busfault.end = 10.2;
            busfault.step = 0.005;
            busfault.isBus = True;
            busfault.bus_idx = 7;

            busfaultlist = gridpack.dynamic_simulation.EventVector([busfault])

            hadapp.transferPFtoDS()
            hadapp.initializeDynSimu(busfaultlist)

            while (not hadapp.isDynSimuDone()):
                if (bApplyAct and
                    ( isteps == 2500 or isteps == 3000 or
                      isteps == 3500 or isteps == 4000 )):
                    hadapp.applyAction(loadshedact)
                hadapp.executeDynSimuOneStep()
                ob_vals = hadapp.getObservations();
                print (isteps, ob_vals)
                isteps = isteps + 1

        # See if we could do it again, if we wanted to
        hadapp = None
        print ('----renke python debug test---before second time ini hadrec module')
        hadapp = gridpack.hadrec.Module()
        hadapp.solvePowerFlowBeforeDynSimu(arg)        
        # It's important to force deallocation order here
        hadapp = None
        env = None
        
