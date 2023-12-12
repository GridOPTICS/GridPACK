#!/usr/bin/env python2
# -*- mode: python; py-which-shell: "python2";-*-
# -------------------------------------------------------------
# file: task_manager.py
# -------------------------------------------------------------
# -------------------------------------------------------------
# Battelle Memorial Institute
# Pacific Northwest Laboratory
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created January 27, 2020 by Perkins
# Last Change: 2020-01-27 14:49:13 d3g096
# -------------------------------------------------------------

# RCS ID: $Id$

import sys, os
import gridpack

# -------------------------------------------------------------
# variable initialization
# -------------------------------------------------------------
program = os.path.basename(sys.argv[0])

# -------------------------------------------------------------
# main program
# -------------------------------------------------------------

e = gridpack.Environment()
c = gridpack.Communicator()
tskmgr = gridpack.TaskManager(c)

task = gridpack.TaskCounter()

tskmgr.set(100)
while tskmgr.nextTask(task):
  sys.stdout.write("process %d of %d executing task %d\n" %
                   (c.rank(), c.size(), task.task_id))
tskmgr = None
e = None

