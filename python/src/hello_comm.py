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
# Last Change: 2022-11-03 10:59:25 d3g096
# -------------------------------------------------------------

import sys, os
import gridpack
from mpi4py import MPI

# -------------------------------------------------------------
# variable initialization
# -------------------------------------------------------------
program = os.path.basename(sys.argv[0])

# -------------------------------------------------------------
# main program
# -------------------------------------------------------------

the_comm = MPI.COMM_WORLD

global_me = the_comm.Get_rank()

color = global_me % 2

another_comm = the_comm.Split(color, global_me)

e = gridpack.Environment(another_comm)
c = gridpack.Communicator()

sys.stdout.write("%s: hello from process %d of %d, color %d\n" %
                 (program, c.rank(), c.size(), color))

