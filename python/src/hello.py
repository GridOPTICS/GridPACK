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
# Last Change: 2020-01-27 08:11:21 d3g096
# -------------------------------------------------------------

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

sys.stdout.write("%s: hello from process %d of %d\n" %
                 (program, c.rank(), c.size()))

