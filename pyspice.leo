<?xml version="1.0" encoding="UTF-8"?>
<?xml-stylesheet ekr_test?>
<leo_file>
<leo_header file_format="2" tnodes="0" max_tnode_index="0" clone_windows="0"/>
<globals body_outline_ratio="0.30390625">
	<global_window_position top="0" left="0" height="687" width="1280"/>
	<global_log_window_position top="0" left="0" height="0" width="0"/>
</globals>
<preferences/>
<find_panel_settings/>
<vnodes>
<v t="etihwnad.20060606092308" a="E"><vh>Notes</vh></v>
<v t="dan.20070114122421" a="E"><vh>pyspice hierarchy</vh>
<v t="dan.20070114122421.1"><vh>readNetlist</vh></v>
<v t="dan.20070114122421.2"><vh>Elements</vh></v>
<v t="dan.20070114122421.3"><vh>writeNetlist</vh></v>
<v t="dan.20070114122421.4" a="E"><vh>helpers</vh>
<v t="dan.20070114122421.5"><vh>read_netlist</vh></v>
</v>
</v>
<v t="etihwnad.20060606092308.1" a="E"
expanded="dan.20061011115934,etihwnad.20060605211347,"><vh>@thin pyspice.py</vh></v>
<v t="dan.20070113143715" tnodeList="dan.20070113143715,dan.20070113145158,dan.20070113154450"><vh>@file pyspicetest.py</vh>
<v t="dan.20070113145158"><vh>converstions</vh></v>
<v t="dan.20070113154450"><vh>netlist parsing</vh></v>
</v>
</vnodes>
<tnodes>
<t tx="dan.20070113143715">@first #!/usr/bin/python
@language python
@tabwidth -4
"""Unit test stuff for pyspice.py"""

__author__ = "Dan White (etihwnad@gmail.com)"
__version__ = "$Revision: 1.3 $"
__date__ = "$Date: 2004/05/05 21:57:20 $"
__copyright__ = "Copyright (c) 2007 Dan White"
__license__ = "GPL"

import pyspice
import unittest
from decimal import Decimal as D

@others

if __name__ == "__main__":
    unittest.main()</t>
<t tx="dan.20070113145158">
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
</t>
<t tx="dan.20070113154450">
class netlistParsing(unittest.TestCase):
    """Tests the netlist parser"""
    pass</t>
<t tx="dan.20070114122421"></t>
<t tx="dan.20070114122421.1"></t>
<t tx="dan.20070114122421.2"></t>
<t tx="dan.20070114122421.3"></t>
<t tx="dan.20070114122421.4"></t>
<t tx="dan.20070114122421.5">
def read_netlist(fname):
    """Reads the entire netlist and returns a list of two-element lists which contain:
        [startNum, spiceLine]
    -whitespace is removed from 'param = val' pairs
    -empty lines are converted to comments
    -lines are lowercased for uniformity since SPICE is case-insensitive"""
    import re
    
    if isinstance(fname,file):
        ifp=fname
    else:
        ifp=open(fname,'rU')

    nline=0
    lines=[]    #raw expanded netlist

    re_param=re.compile(r"(\S*)\s*=\s*(\S*)")
    
    startLine = None
    for line in ifp:
        line=line.strip('\r\n')
        nline+=1
        
        #pass through empty lines
        # done this way in case the line consists of whitespace only and
        # is, in fact, not strictly empty
        if not len(line.split()):
            #convert empty line to comment as a placeholder
            lines.append([nline,'*'])
            continue
        
        #pass through comments, they stay as-is
        elif line[0]=='*':
            lines.append([nline,line])
            continue
        
        #case is unimportant in SPICE, we use names as dictionary keys later
        # so it is important to us
        line=line.lower()
        
        #remove whitespace in parameter assignments
        # to prepare for x.split(' ') that happens later:
        # 'as = 3e-12' =&gt; 'as=3e-12'
        line = re.sub(re_param,r'\1=\2',line)
        
        #beginning of SPICE line
        if line[0] != '+':
            lines.append([nline,line])
        else:            #line continuation
            lines[-1][1] = lines[-1][1]+line[1:]
    return lines</t>
<t tx="etihwnad.20060606092308"></t>
</tnodes>
</leo_file>
