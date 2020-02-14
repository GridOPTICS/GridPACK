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
# Last Change: 2020-02-14 10:02:36 d3g096
# -------------------------------------------------------------

from unittest import TestCase

import Test


class TestTester(TestCase):
    def test_is_string(self):
        t = Test.SerialMPI()
        self.assertTrue(isinstance(t, Test.SerialMPI))
        s = [ "this", "is", "a", "list", "of", "strings" ]
        t.show(s)
