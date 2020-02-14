#!/usr/bin/env python
# -------------------------------------------------------------
# file: runme.py
# -------------------------------------------------------------
# -------------------------------------------------------------
# Copyright (c) 2013 Battelle Memorial Institute
# Licensed under modified BSD License. A copy of this license can be found
# in the LICENSE file in the top level directory of this distribution.
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created February 14, 2020 by Perkins
# Last Change: 2020-02-14 09:52:47 d3g096
# -------------------------------------------------------------

import sys, os
from optparse import OptionParser
import Test

# -------------------------------------------------------------
# variable initialization
# -------------------------------------------------------------
program = os.path.basename(sys.argv[0])
usage = "usage: " + program

# -------------------------------------------------------------
# main program
# -------------------------------------------------------------

t = Test.SerialMPI()
s = [ "this", "is", "a", "list", "of", "strings" ]
t.show(s)


