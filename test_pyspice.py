#!/usr/bin/python
#@+leo-ver=4
#@+node:@file test_pyspice.py
#@@first
#@@language python
#@@tabwidth -4
"""Unit test stuff for pyspice.py"""

__author__ = "Dan White (etihwnad@gmail.com)"
__version__ = "$Revision: 1.3 $"
__date__ = "$Date: 2004/05/05 21:57:20 $"
__copyright__ = "Copyright (c) 2007 Dan White"
__license__ = "GPL"

import pyspice
import unittest
from decimal import Decimal as D

#@+others
#@+node:unit conversions

class unitConversion(unittest.TestCase):
    """Tests conversion between SPICE string and Decimal number."""
    knownValues = ( ('5T',      D('5.0e12')),
                    ('5G',      D('5.0e9')),
                    ('10MEG',   D('10.0e6')),
                    ('342x',    D('342.0e6')),
                    ('15k',     D('15.0e3')),
                    ('1MIL',    D('25.4e-6')),
                    ('435M',    D('435e-3')),
                    ('1U',      D('1.0e-6')),
                    ('67N',     D('67.0e-9')),
                    ('4P',      D('4.0e-12')),
                    ('3F',      D('3.0e-15')) )

    badValues =   ( 'like' )

    def test_unit_knownValues(self):
        """unit() should give known result with known input"""
        for s, num in self.knownValues:
            result = pyspice.unit(s)
            self.assertEqual(num,result)

    def test_unit_badValues(self):
        """unit() should fail with bad input"""
        for s in self.badValues:
            self.assertRaises(pyspice.BadUnitError, pyspice.unit, s)
#@-node:unit conversions
#@+node:netlist parsing

class netlistParsing(unittest.TestCase):
    """Tests the netlist parser"""
    pass
#@nonl
#@-node:netlist parsing
#@-others

if __name__ == "__main__":
    unittest.main()
#@nonl
#@-node:@file test_pyspice.py
#@-leo
