# -------------------------------------------------------------
# file: testit.py
# -------------------------------------------------------------
# -------------------------------------------------------------
# Copyright (c) 2013 Battelle Memorial Institute
# Licensed under modified BSD License. A copy of this license can be found
# in the LICENSE file in the top level directory of this distribution.
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created February 14, 2020 by Perkins
# Last Change: 2020-02-14 07:51:13 d3g096
# -------------------------------------------------------------

from unittest import TestCase

import Test

class TestJoke(TestCase):
    def test_is_string(self):
        s = Test.SerialMPI()
        self.assertTrue(isinstance(s, Test.SerialMPI))
