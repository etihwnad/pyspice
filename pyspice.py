#!/usr/bin/env python
#@+leo-ver=4
#@+node:@file pyspice.py
#@@first
#@@language python
#@@tabwidth -4
#@<< head docstring >>
#@+node:<< head docstring>>
"""
#@<< head >>
#@+node:<< head >>
pyspice.py v0.2a

SPICE pre-processor that combines parallel elements (e.g. capacitors, mosfets)
for GREATLY reduced simulation time.  Uses the 'M' parameter of MOSFETS
when combining parallel fets.  Makes simulating fingered transistors easier and
faster.  This module is in a constant state of flux because I am using and
modifying it for my needs.

 Combine parallel capacitors
 Drop combined caps smaller than X fF
 Usage:
  pyspice.py [options] [-i infile] [-o outfile]
 
 Use pyspice.py -h for all the options.
#@-node:<< head >>
#@nl
#@<< copyright >>
#@+node:<< copyright >>
Copyright Dan White 2006

Licensed by the GPL, see http://www.whiteaudio.com/soft/COPYING or the current
GNU GPL license for details.

#@-node:<< copyright >>
#@nl
#@<< release notes >>
#@+node:<< release notes >>
Release Notes, changelog, whatever it turns out as...
-----------------------------------------------------
#@+others
#@+node:v0.2a
pyspice.py v0.2a:
----------------
-added a missing newline before an import statement

#@-node:v0.2a
#@+node:v0.2
pyspice.py v0.2:
----------------
-At least default (pass through) handling of all element types.
NOTE: For combining, this uses a global node name scheme.  In other
    words: subcircuits, libraries, etc. are not in a separate node
    namespace as they should be, beware.
-Changed structure of classes (in LEO), there are base classes that contain
    common attributes and element classes that define the specific behavior.
-This version _should_ work with any netlist and only touch M's and C's, YMMV.
-Work is ongoing on the class structure and most important IMO is getting netlist
    hierarchy implemented. 

#@-node:v0.2
#@+node:v0.1
pyspice v0.1:
-------------
Initial release.
Only worked for netlist containing MOSFETs and Capacitors.
#@-node:v0.1
#@-others
#@nonl
#@-node:<< release notes >>
#@nl
#@+others
#@+node:structures
 Data structures: (outdated)
  lines - list of tuples: (original_line, line.split(), line_type)
  caps  - dictionary of node pairs and capacitance between nodes
#@-node:structures
#@+node:todo
TODO:
xpreserve comments in position
-element class definitions
    xinductor
    xv-source
    xi-source
    -e-source (VCVS)
-preserve node namespaces (within subckts, libraries, etc.)
    -find illegal SPICE node name character to prepend namespace
     e.g. @sub1.node1
-handle control statements specially
    -.model
    -.dc, .ac, .tran
    -.lib/.endl handling
    -.meas
    -.print
    -.probe
    -.param
    -.option
    -.global
    -.ic
    -.subckt/.ends
    -.prot/.unprot
    -.alter blocks
    -.end
    -others?
-option to find/evaluate .param statements?
    -HSPICE takes params, ngspice doesn't
    -logic to find when an explicit number is specified (e.g. 10k) or
     when it is a parameter (e.g. 'value')
    -create dict of .params and use as substitution keys in:
        k=v model parameters
        node names
        element values (Rxx n1 n2 'value'; Cxx n1 n2 'fc/2')
        others?
-option to make a flat output (inline everything):
    .lib statements
    .include statements
    .model statements
    others?
-store node dictionary to keep track of elements connected to that node
    -allows tracking of R-C nodes to drop C or R based on time constant
    -allows more sophisticated dropping/combining with series elements
        specifically R+R+R+R chains -> equivalent R+- (roughly) for faster sims
                      C C C                        C
#@-node:todo
#@-others
"""
#@@nocolor
#@-node:<< head docstring>>
#@nl
#@<< global imports >>
#@+node:<< global imports >>
import sys, getopt, textwrap, warnings
from optparse import OptionParser
import re
from decimal import Decimal

try:
    from cStringIO import StringIO as StringIO
except:
    from StringIO import StringIO as StringIO
#@-node:<< global imports >>
#@nl
#@<< global vars >>
#@+node:<< global vars >>
##
# Global Shorthands
##
stderr=sys.stderr
#set line continuation style and max width
wrapper=textwrap.TextWrapper(subsequent_indent='+ ',width=75)

#not used?
#element lists
_defined_types=['comment',
                'control',
                'capacitor',
                'mosfet',
                'resistor']
capacitors=[]
resistors=[]
mosfets=[]
#end not used?

#global variables
_counts=dict() #keyed by spice element letter
_ncombine_capacitors=0
_ncombine_inductors=0
_ncombine_mosfets=0
_ncombine_res=0

#namespace tracking
# -not implemented yet
_current_scope=''
#@-node:<< global vars >>
#@nl

#@+at
# Debug options: string of debugging elements to print
# * - comments
# follow SPICE element names (c-capacitor, l-inductor)
#@-at
#@@c
dbg=''

#@+others
#@+node:Option processing
#@+node:options
##
# Fancy option processing
##
def options(args=sys.argv):
    """Define options and parse argument list.
    """
    global ifp, ofp
    usage="""%prog [options]"""
    desc="""This script reads a SPICE input files and processes it according
    to the options given.  It is especially useful for processing netlists
    from the output of layout extractors by combining parallel C's and
    FET's.  More features added upon request."""
    parser = OptionParser(usage=usage,description=desc)
    parser.add_option('-i','--infile',dest='infile',default='stdin',
                      help='Input SPICE file to be processed'
                           ', (default: %default)')
    parser.add_option('-o','--outifle',dest='outfile',default='stdout',
                      help='Output file for changes'
                           ', (default: %default)')
    parser.add_option('-d','--dropcap',dest='dropcap',
                      default='10', metavar='X',
                      help='Drop all capacitors smaller than X fF'
                           ', (default: %default)')
    parser.add_option('--no-combine-c',dest='combine_c',action='store_false',
                      default=True,help='Do not combine parallel capacitors')
    parser.add_option('--no-combine-m',dest='combine_m',action='store_false',
                      default=True,help='Do not combine parallel MOSFETs')
    parser.add_option('-v','--verbose',dest='v',action='store_true',
                      default=True,
                      help='Show info and debugging messages (default)')
    parser.add_option('-q','--quiet',dest='v',action='store_false',
                      help='Suppress all messages on stderr')
    parser.add_option('-w','--linewidth',dest='linewidth', default=75,
                      help='Max. line width for netlist (default: %default)')

    (opt, args) = parser.parse_args()
    #infile
    if opt.infile=='stdin':
        ifp=sys.stdin
        if opt.v: info('Input: stdin')
    else:
        try:
            ifp=open(opt.infile,'rU') #python's universal line ending mode
            if opt.v: info('Input: '+opt.infile)
        except IOError, (errno,strerror):
            print>>stderr, "IOError(%s): %s '%s'" % (errno,strerror,opt.infile)
    #outfile
    if opt.outfile=='stdout':
        ofp=sys.stdout
        if opt.v: info('Output: stdout')
    else:
        try:
            ofp=open(opt.outfile,'w')
            if opt.v: info('Output: '+opt.outfile)
        except IOError, (errno,strerror):
            print>>stderr, "IOError(%s): %s '%s'" % (errno,strerror,opt.infile)
    #dropcap
    opt.dropcap=float(opt.dropcap)
    opt.dropcap=opt.dropcap*1e-15
    if opt.v: info('Dropping caps < '+str(opt.dropcap)+' F')
    info("Combining c's:"+str(opt.combine_c))
    info("Combining m's:"+str(opt.combine_m))
    #return the option object
    return opt
#@-node:options
#@-node:Option processing
#@+node:classes
#@+node:exceptions

class PyspiceError(Exception):
    pass
#@+node:class UnitError
class BadUnitError(PyspiceError):
    pass
#@-node:class UnitError
#@+node:class ElementError

class ElementError(LookupError):
    #@	@+others
    #@+node:__init__
    def __init__(self,elm='???'):
        self.elm=elm
    #@-node:__init__
    #@+node:__str__
    def __str__(self):
        return str('No class defined for this element: '+self.elm)
    #@-node:__str__
    #@-others
#@-node:class ElementError
#@-node:exceptions
#@+node:class Netlist

class Netlist:
    """Base class that holds a netlist.
    Notes:
        -this will eventually hold the entire shebang
        -providing __init__ with a file name will:
            -read in netlist
            -classify the lines
            -take care of hierarchy
            -source other files
    """
    #@    @+others
    #@+node:__init__
    def __init__(self,fname=None):
        """Optionally reads a netlist from a file"""
        self.lines = []
        self.deck = []
        if fname:
            self.readfile(fname)
    #@-node:__init__
    #@+node:add_element
    def add_element(self,line):
        """Takes a line (usually split on ' '), classifies it and adds the new
        object to the netlist"""
        self.deck.append( self.classify(line) )
    #@nonl
    #@-node:add_element
    #@+node:add_line
    def add_line(self,line):
        """TODO add_line docstring"""
        self.lines.append(line)
        
    #@nonl
    #@-node:add_line
    #@+node:readfile
    def readfile(self,fname):
        """Read a SPICE netlist from the open file pointer
    
        readfile(filename) -> array lines
        
        Returns a list of expanded lines (without continuation '+')
        Keeps case of comments, all other lines are lowercased
        
        return:
            netlist (list) of SPICE netlist lines
        
        Notes:
            -we need to read at least a full line with continuations before we can
             add the line to the netlist
        """
        import re
        
        if isinstance(fname,file):
            ifp=fname
        else:
            ifp=open(fname,'rU')
        n=0
        #finds a "name = value" pair for shrinking
        re_param=re.compile(r"(\S*)\s*=\s*(\S*)")
        #read in whole file to a list of lines
        
        #
        #
        # NOTE:
        #   This doesn't work yet, want to modify list in-place so we don't
        #   need a copy to iterate over.  The purpose is to be able to
        #   un-continue lines and just delete the continuation.
        #
        #
        lines = ifp.readlines()
        n = -1
        for line in lines:
            n += 1
            line = line.strip('\r\n') #handles any type of line ending
            #pass through empty lines
            if not len(line.split()):
                lines[n] = '*' #convert empty line to comment as a placeholder
                continue
            #pass through comments, they stay as-is
            elif line[0]=='*':
                lines[n] = '*'
                continue
            #case is unimportant in SPICE
            line = line.lower()
            #remove whitespace in parameter assignments
            # to prepare for x.split(' ') that happens next:
            #  'as = 3e-12' => 'as=3e-12'
            line = re.sub(re_param,r'\1=\2',line)
            #handle line continuations
            if line[0]=='+':
                lines[-1].extend(line[1:].split())
                del lines[n]
            else:
                lines[n] = line.split()
        #read in file, now classify
        self.deck.append(map(self.add_element,self.lines))
    #@-node:readfile
    #@+node:classify
    def classify(self,line):
        """Takes a line and creates an appropriate SpiceElement"""
        return elements.handler[line[0][0]](line)
    #@nonl
    #@-node:classify
    #@-others
#@nonl
#@-node:class Netlist
#@+node:class ElementHandler
class ElementHandler:
    #@    @+others
    #@+node:__init__
    def __init__(self):
        # the valid element types
        valid_types = '*.abcdefghijklmnopqrstuvwxyz'
        #setup default element handler
        self.handler = dict()
        for t in valid_types:
            self.handler[t] = SpiceElement
    #@nonl
    #@-node:__init__
    #@+node:add_handler
    def add_handler(self,type,handler):
        """Replaces the existing element object definition with the
        given one"""
        self.handler[type] = handler
    #@nonl
    #@-node:add_handler
    #@-others
#@nonl
#@-node:class ElementHandler
#@+node:base classes
#@+at
# These classes are used to break down the spectrum of SPICE elements
# into 'classes' of elements.  I.e. 2 node passives, 2-node sources, 4-node 
# sources,
# and so on.  The elements found in a real netlist are based on these types.
#@-at
#@@c
#@+node:class SpiceElement

class SpiceElement:
    """Base class for SPICE elements.
    Methods:
        __init__(self,line,num) -> SpiceElement
        __str__(self) -> string spice line
        drop() -> False
    """
    #@	@+others
    #@+node:__init__
    def __init__(self,line,num=None):
        """SpiceElement constructor
        line - netlist expanded line
        type - first character
        num  - input netlist line number (for keeping roughly the same
               order when printing modified netlist)
        """
        #accept lists of 'words' also; BE CAREFUL with this, though
        self.line = line
        self.type = line[0]
        self.num = num
    #@-node:__init__
    #@+node:__str__
    def __str__(self):
        return wrapper.fill(self.line)
    #@-node:__str__
    #@+node:drop
    def drop(self,val=0,mode='<'):
        """Template for dropping elements that defaults to NO if
        not overidden in the element class"""
        return False
    #@-node:drop
    #@-others
#@-node:class SpiceElement
#@+node:class Passive2NodeElement

class Passive2NodeElement(SpiceElement):
    """Base class for 2-node elements.
    Assumes SPICE element line:
    xXXX n1 n2 value p1=val p2=val ...
    Inherits:
        None
    Redefines:
        drop(self,val,mode) -> bool
    """
    #@	@+others
    #@+node:__init__
    def __init__(self,line,num):
        SpiceElement.__init__(self,line,num)
        self.type='passive2'
        arr = line.split()
        self.name = arr[0]
        self.n1 = _current_scope + arr[1]
        self.n2 = _current_scope + arr[2]
        self.value = unit(arr[3])
        self.param = dict() #store x=y as dictionary
        for p in arr[4:]:
            k,v=p.split('=')
            self.param[k]=unit(v)
    #@-node:__init__
    #@+node:__str__
    def __str__(self):
        """Returns the netlist-file representation of this element"""
        s=StringIO()
        
        print>>s, self.name,self.n1,self.n2,self.value,
        
        for k,v in self.param.iteritems():
            print>>s, k+'='+str(v),
            
        return wrapper.fill(s.getvalue())
    #@-node:__str__
    #@+node:drop
    def drop(self,val=0.0,mode='<'):
        """Indicate whether to drop the element from the list.
        Occurs iff (val 'mode' self.value)
        
        Can this be converted to specifying an arbitrary binary function?
          This may allow a more elegant comparison.
          e.g. mode=< instead of mode='<' or mode=cap_smaller(x,y)
        """
        if mode=='<':
            if self.value<val: return True
        elif mode=='<=':
            if self.value<=val: return True
        elif mode=='>':
            if self.value>val: return True
        elif mode=='>=':
            if self.value>=val: return True
        else:
            return False
        return False #shouldn't get here, but...
    #@-node:drop
    #@-others
#@-node:class Passive2NodeElement
#@+node:class Active2NodeElement

class Active2NodeElement(SpiceElement):
    """Base class for active 2-node elements.
    Assumes SPICE element line:
    xXXX n1 n2 value p1=val p2=val ...
    Inherits:
        None
    Redefines:
        None
    """
    #@	@+others
    #@+node:__init__
    def __init__(self,line,num):
        SpiceElement.__init__(self,line,num)
        self.type = 'active2'
        arr = line.split()
        self.name = arr[0]
        self.n1 = _current_scope + arr[1]
        self.n2 = _current_scope + arr[2]
        self.value = unit(arr[3])
        self.param = dict() #store x=y as dictionary
        for p in arr[4:]:
            k,v=p.split('=')
            self.param[k]=unit(v)
    #@-node:__init__
    #@+node:__str__
    def __str__(self):
        """Returns the netlist-file representation of this element"""
        s=StringIO()
        
        print>>s, self.name,self.n1,self.n2,self.value,
        
        for k,v in self.param.iteritems():
            print>>s, k+'='+str(v),
            
        return wrapper.fill(s.getvalue())
    #@-node:__str__
    #@-others
#@-node:class Active2NodeElement
#@+node:class Active4NodeElement

class Active4NodeElement(SpiceElement):
    """Base class for active 4-node elements (xCyS).
    Assumes SPICE element line:
    xXXX n1 n2 value p1=val p2=val ...
    Inherits:
        None
    Redefines:
        None
    """
    #@	@+others
    #@+node:__init__
    def __init__(self,line,num):
        SpiceElement.__init__(self,line,num)
        self.type='active4'
        arr = line.split()
        self.name = arr[0]
        self.n1 = _current_scope + arr[1]
        self.n2 = _current_scope + arr[2]
        self.n3 = _current_scope + arr[3]
        self.n4 = _current_scope + arr[4]
        self.value=unit(arr[5])
        self.param=dict() #store x=y as dictionary
        for p in arr[6:]:
            k,v=p.split('=')
            self.param[k]=unit(v)
    #@-node:__init__
    #@+node:__str__
    def __str__(self):
        s=StringIO()
        print>>s, self.name,self.n1,self.n2,self.n3,self.n4,self.value,
        for k,v in self.param.iteritems():
            #are there instances when 0 is (in)significant?
            #if v==0: continue
            print>>s, k+'='+str(v),
        return wrapper.fill(s.getvalue())
    #@-node:__str__
    #@-others
#@-node:class Active4NodeElement
#@-node:base classes
#@+node:element classes
#@+at
# This is a(n incomplete) definition of the various SPICE elements.
# 
# NOTE: When adding a new element type definition, be sure to add a handler
#   for the new class after defining the class using:
#       elements.add_handler('x',Xdevice)
#@-at
#@@c
#make a repository for element handlers
elements=ElementHandler()
#@+node:class CommentLine

class CommentLine(SpiceElement):
    """SPICE Comment line (/^\*.*/)
    """
    #@	@+others
    #@+node:__init__
    def __init__(self,line,num):
        SpiceElement.__init__(self,line,num)
        self.type='comment'
    #@-node:__init__
    #@-others
elements.add_handler('*',CommentLine)
#@-node:class CommentLine
#@+node:class ControlElement

class ControlElement(SpiceElement):
    """Control statement object, no processing for now.
    Note: this will eventially be a base class for the real control elements
    
    Note: currently has no knowledge of blocks (.lib/.endl, .subckt/.ends)
    has only ONE node namespace, make sure subckt's have unique node names!
    """
    #@	@+others
    #@+node:__init__
    def __init__(self,line,num):
        SpiceElement.__init__(self,line,num)
        self.type='control'
    #@-node:__init__
    #@-others
elements.add_handler('.',ControlElement)
#@-node:class ControlElement
#@+node:class Capacitor

class Capacitor(Passive2NodeElement):
    """Assumes SPICE element line:

    cXXX n1 n2 value p1=val p2=val ...
    Provides:
        isparallel(other)
        combine(other)
    """
    #@	@+others
    #@+node:__init__
    def __init__(self,line,num):
        Passive2NodeElement.__init__(self,line,num)
        self.type='capacitor'
    #@-node:__init__
    #@+node:isparallel
    def isparallel(self,other):
        """Returns True if instance is parallel with other instance
        """
        if self.n1==other.n1 and self.n2==other.n2:
            return True
        elif self.n1==other.n2 and self.n2==other.n1:
            return True
        else:
            return False
    #@-node:isparallel
    #@+node:combine
    def combine(self,other):
        """Adds values if capacitors are in parallel, returns True if
        it combined them.
    
        NOTE: Does not currently touch param dictionary when combining,
          just the values.  How should this be done?  Maybe combine iff
          params are identical to avoid problems?
        """
        global _ncombine_capacitors
        if self.isparallel(other):
            self.value+=other.value
            _ncombine_capacitors+=1
            return True
        else:
            return False
    #@-node:combine
    #@-others
elements.add_handler('c',Capacitor)
#@-node:class Capacitor
#@+node:class Inductor

class Inductor(Passive2NodeElement):
    """Assumes SPICE element line:

    cXXX n1 n2 value p1=val p2=val ...
    Provides:
        isparallel(other)
        combine(other)
    """
    #@	@+others
    #@+node:__init__
    def __init__(self,line,num):
        Passive2NodeElement.__init__(self,line,num)
        self.type='inductor'
    #@-node:__init__
    #@+node:isparallel
    def isparallel(self,other):
        """Returns True if instance is parallel with other instance
        """
        if self.n1==other.n1 and self.n2==other.n2:
            return True
        elif self.n1==other.n2 and self.n2==other.n1:
            return True
        else:
            return False
    #@-node:isparallel
    #@+node:combine
    def combine(self,other):
        """Combines values if inductors are in parallel, returns True if
        it combined them.
    
        NOTE:
         -Does not currently touch param dictionary when combining,
          just the values.  How should this be done?
        """
        global _ncombine_inductors
        if self.isparallel(other):
            self.value=(self.value*other.value)/(self.value+other.value)
            _ncombine_inductors+=1
            return True
        else:
            return False
    #@-node:combine
    #@-others
elements.add_handler('l',Inductor)
#@-node:class Inductor
#@+node:class Mosfet

class Mosfet(SpiceElement):
    """Mosfet constructor takes an array derived from the
    netlist line
    """
    #@	@+others
    #@+node:__init__
    def __init__(self,line,num):
        if isinstance(line,str):
            line=line.split()
        self.line=line
        self.type='mosfet'
        self.num=num
        self.d=line[1]
        self.g=line[2]
        self.s=line[3]
        self.b=line[4]
        self.model=line[5]
        self.param=dict()
        for p in line[6:]:
            k,v=p.split('=')
            self.param[k]=unit(v)
        self.w=self.param['w']
        self.l=self.param['l']
    #@-node:__init__
    #@+node:__str__
    def __str__(self):
        s=StringIO()
        print>>s, self.line[0], self.d, self.g, self.s, self.b, self.model,
        for k,v in self.param.iteritems():
            if v==0: continue
            print>>s, k+'='+str(v),
        return wrapper.fill(s.getvalue())
    #@-node:__str__
    #@+node:isparallel
    def isparallel(self,other):
        """Returns True if transistors are parallel
        """
        # check gate, substrate, and model first
        if self.g==other.g and self.b==other.b and self.model==other.model:
            #source and drain can be reversed
            if self.d==other.d and self.s==other.s:
                return True
            elif self.d==other.s and self.s==other.d:
                return True
            else:
                return False
        else:
            return False
    
    #@-node:isparallel
    #@+node:combine
    def combine(self,other):
        """Combines adds other to self iff the transistors are identical,
        will NOT combine if W/L is different.  Parameter 'M' is incremented on
        self, other is left alone.
    
        Returns True if it combined the transistors.
        
        Increments global _ncombine_mosfets for information.
        
        NOTE: This currently merely adds the parameters (except w, l, and m)
        without regard to their meaning.  Here is the place to specially handle
        certain FET parameters.  Average certain parameters?
        """
        global _ncombine_mosfets
        if self.isparallel(other):
            #combine iff W/L (for original FET) is same also
            if self.w==other.w and self.l==other.l:
                for k,v in other.param.iteritems():
                    if k=='w' or k=='l' or k=='m': continue
                    self.param[k]+=v
                if ('m' in self.param.keys()) or ('m' in other.param.keys()):
                    #add other's M parameter or increment
                    self.param['m']+=other.param.get('m',1)
                else:
                    self.param['m']=2
                
                _ncombine_mosfets+=1
                return True
        else:
            return False
    #@-node:combine
    #@-others
elements.add_handler('m',Mosfet)
#@-node:class Mosfet
#@+node:class Resistor

class Resistor(Passive2NodeElement):
    """Assumes SPICE element line:

    rXXX n1 n2 value p1=val p2=val ...
    """
    #@	@+others
    #@+node:__init__
    def __init__(self,line,num):
        Passive2NodeElement.__init__(self,line,num)
        self.type='resistor'
    #@-node:__init__
    #@-others
elements.add_handler('r',Resistor)
#@-node:class Resistor
#@+node:class Vsource

class Vsource(Active2NodeElement):
    """Assumes SPICE element line:

    vXXX n1 n2 value p1=val p2=val ...
    """
    #@	@+others
    #@+node:__init__
    def __init__(self,line,num):
        Active2NodeElement.__init__(self,line,num)
        self.type='vsource'
    #@-node:__init__
    #@-others
elements.add_handler('v',Vsource)
#@-node:class Vsource
#@+node:class Isource

class Isource(Active2NodeElement):
    """Assumes SPICE element line:

    iXXX n1 n2 value p1=val p2=val ...
    """
    #@	@+others
    #@+node:__init__
    def __init__(self,line,num):
        Active2NodeElement.__init__(self,line,num)
        self.type='isource'
    #@-node:__init__
    #@-others
elements.add_handler('i',Isource)
#@-node:class Isource
#@-node:element classes
#@-node:classes
#@+node:helpers
#@+node:unit
def unit(s):
    """Takes a string and returns the equivalent float.
    '3.0u' -> 3.0e-6"""
    mult={'t'  :Decimal('1.0e12'),
          'g'  :Decimal('1.0e9'),
          'meg':Decimal('1.0e6'),
          'x'  :Decimal('1.0e6'),
          'k'  :Decimal('1.0e3'),
          'mil':Decimal('25.4e-6'),
          'm'  :Decimal('1.0e-3'),
          'u'  :Decimal('1.0e-6'),
          'n'  :Decimal('1.0e-9'),
          'p'  :Decimal('1.0e-12'),
          'f'  :Decimal('1.0e-15')}
    m=re.search('^([0-9e\+\-\.]+)(t|g|meg|x|k|mil|m|u|n|p|f)?',s.lower())
    
    try:
        if m.group(2):
            return Decimal(Decimal(m.group(1)))*mult[m.group(2)]
        else:
            return Decimal(m.group(1))
    except:
        raise BadUnitError
#@nonl
#@-node:unit
#@+node:debug
def debug(message):
    """Print debugging info to stderr."""
    print>>stderr,'Debug:',message
#@-node:debug
#@+node:info
def info(message):
    """Print information to stderr."""
    print>>stderr,'Info:',message
#@-node:info
#@+node:warning
def warning(message,elm=None,num=None):
    """Print warning to stderr.  If elm and num defined,
    print different message"""
    if elm and num:
        message=_opt.infile+":"+str(num)+" '"+elm+\
                "' type not defined yet, passing through..."
    print>>stderr,'Warning:',message
#@-node:warning
#@-node:helpers
#@+node:drop_2node
def drop_2node(elm,val,mode='<',type=None, verbose=True):
    """Drop 2-node elements in elm list according to val.

    mode = '<' | '>'
    type - print info on what and how many it dropped
    Note: only works correctly for capacitors for now
    """
    new_elm=[]
    val=float(val)
    for i in elm:
        if mode=='<' and i.value>=val:
            new_elm.append(i)
        elif mode=='>' and i.value<=val:
            new_elm.append(i)
        else:
            continue
    #Print info about how many it dropped
    if verbose:
        infostr='Dropped '+str((len(elm)-len(new_elm)))+' '
        if type and type[-1]=='s':
            info(infostr+type)
        elif type:
            info(infostr+type+'s')
        else:
            info(infostr+'elements')
    return new_elm
#@-node:drop_2node
#@+node:write_2node
##
# write 2-node SPICE elements to file
# node pair is specified by key in dict:
#  "node0,node1"
##
#NOTE:
# this function is dying a slow, painful death.  Each element is
# getting its own custom __str__() method for writing to netlists.
#
def write_2node(elm, type=None, ofp=sys.stdout, comment=None):
    """Write 2-node SPICE elements to file

    elm - dictionary of elements
    type - SPICE element name
    ofp  - output file pointer
    comment - comment string at head of elm list

    Note: "type" must begin with a SPICE element letter, the rest is printed
    as identifying information, e.g. linductor, capacitor, mosfet.
    """
    if not type:
        raise SyntaxError('Must define a SPICE element name')
    i=1
    if comment:
        print>>ofp,'\n**\n* '+comment+'\n**'
    for k,v in elm.iteritems():
        node=k.split(',')
        print>>ofp, type[0]+str(i).rjust(3,'0'),node[0],node[1],v
        i+=1
    infostr='Wrote '+str(i)+' '+type
    if type[-1]=='s':
        info(infostr)
    else:
        info(infostr+'s')
#@-node:write_2node
#@+node:read_netlist
##
# Read netlist from open file object
#  -make sure if reading from stdin to read only once and
#   use this function to make sure you read the entire netlist
##
def read_netlist(fname):
    #@    << docstring >>
    #@+node:<< docstring >>
    """Read a SPICE netlist from the open file pointer
    
    read_netlist(filename) -> array lines
    
    Returns a list of expanded lines (without continuation '+')
    Keeps case of comments, all other lines are lowercased
    
    return:
        netlist (list) of SPICE netlist lines
    """
    #@-node:<< docstring >>
    #@nl
    #@    << imports >>
    #@+node:<< imports >>
    import re
    #@nonl
    #@-node:<< imports >>
    #@nl

    if isinstance(fname,file):
        ifp=fname
    else:
        ifp=open(fname,'rU')
#    netlist=ifp.readlines()
    nline=0
    lines=[]    #raw expanded netlist
    #
    re_param=re.compile(r"(\S*)\s*=\s*(\S*)") 
    for line in ifp:
        line=line.strip('\r\n')
        #pass through empty lines
        if not len(line.split()):
            #convert empty line to comment as a placeholder
            lines.append('*')
            nline+=1
            continue #next please...
        #pass through comments, they stay asis
        elif line[0]=='*':
            lines.append(line)
            nline+=1
            continue #next please...
        #case is unimportant in SPICE
        line=line.lower()
        #remove whitespace in parameter assignments
        # to prepare for x.split(' ') that happens later:
        # 'as = 3e-12' => 'as=3e-12'
        line=re.sub(re_param,r'\1=\2',line)
        if line[0]!='+': #beginning of SPICE line
            lines.append(line)
            nline+=1
        else:            #line continuation
            line=line[1:]
            lines[-1]=lines[-1]+line
    return lines
#@-node:read_netlist
#@+node:classify
def classify(net):
    """Reads expanded netlist and classifies each line,
    calls appropriate function to combine nodes.
    
    net - list of unwrapped SPICE netlist lines
    """
    # Is there a better (faster) way of doing this?  Maybe generating a dict
    # of handler functions and calling based on the first character.  At any
    # rate, it would make the classification a constant-time operation.
    import string
    elements=dict()
    global _counts
    global _opt
    opt=_opt
    #initialize to all element types
    for ch in '*.'+string.lowercase:
        elements[ch]=[]
        _counts[ch]=0
    num=-1 #line counter
    for line in net:
        num+=1
        arr=line.split()
        x=arr[0][0]

        # Comment
        if x=='*':
            elements[x].append(CommentLine(line,num))
            _counts[x]+=1

        # Control line
        elif x=='.':
            elements[x].append(ControlElement(line,num))
            _counts[x]+=1

        elif x=='a':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        elif x=='b':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        # Capacitor
        elif x=='c':
            elm=Capacitor(line,num)
            #don't combine if it's the first encountered
            if _counts[x]==0:
                elements[x].append(elm)
            else:
                #combine if parallel or add if unique
                if opt.combine_c and elements[x][-1].combine(elm):
                    pass
                else:
                    elements[x].append(elm)
            _counts[x]+=1

        elif x=='d':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        elif x=='e':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        elif x=='f':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        elif x=='g':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        elif x=='h':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        # Current Source
        elif x=='i':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        elif x=='j':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        elif x=='k':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        # Inductor
        elif x=='l':
            elm=Inductor(line,num)
            #don't combine if it's the first encountered
            if _counts[x]==0:
                elements[x].append(elm)
            else:
                #combine if parallel or add if unique
                if elements[x][-1].combine(elm):
                    pass
                else:
                    elements[x].append(elm)
            _counts[x]+=1

        # Mosfet
        elif x=='m':
            elm=Mosfet(line,num)
            if _counts[x]==0:
                #don't combine if it's the first encountered
                elements[x].append(elm)
                mosfets.append(elm)
            elif not opt.combine_m:
                #do not combine elements
                elements[x].append(elm)
            else:
                #search list backwards to find a parallel one
                #in case they aren't adjacent, faster?
                i=len(elements[x])-1
                while i>=0:
                    if elements[x][i].combine(elm):
                        break
                    else:
                        i-=1
                #found unique MOSFET
                if i==-1:
                    elements[x].append(elm)
            _counts['m']+=1

        elif x=='n':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        elif x=='o':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        elif x=='p':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        elif x=='q':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        # Resistor
        elif x=='r':
            elm=Resistor(line,num)
            elements[x].append(elm)
            _counts['r']+=1

        elif x=='s':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        elif x=='t':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        elif x=='u':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        # Voltage Source
        elif x=='v':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        elif x=='w':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        # Subcircuit
        elif x=='x':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        elif x=='y':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        elif x=='z':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
            _counts[x]+=1

        elif x=='+':
            raise Error('No line continuations (+) allowed.')

        else:
            raise Error('Encountered unknown element: '+line)

    #add 'special' elements to netlist
    #debugging/information header
    print>>ofp,'* Input Element counts'
    for k,v in _counts.iteritems():
        if v!=0: print>>ofp, '*',k,'-',v
    #return classified netlist
    return elements

#@-node:classify
#@+node:main
####
# 
# test code
# 
####
def main():
    global _opt
    opt = options()
    _opt=opt
    
    import textwrap
    #set line continuation style and max width
    wrapper=textwrap.TextWrapper(subsequent_indent='+ ',width=opt.linewidth)

    print>>ofp,"* pyspice.py: by Dan White <etihwnad at gmail dot com>"
    print>>ofp,"* mail me bug reports, fixes, and comments if you find this useful"
    print>>ofp,"* ----------------------------------------------------------------"

    netlist = read_netlist(ifp)
    nlist_new=classify(netlist)

    info('Combined %i parallel capacitors' % _ncombine_capacitors)
    info('Combined %i parallel mosfets' % _ncombine_mosfets)

    nlist_new['c']=drop_2node(nlist_new['c'], opt.dropcap, type='capacitor')


    #create a new list of netlist objects
    #in order of first appearance
    all=[]
    print '*\n* Output Element counts'
    for type in nlist_new.keys():
        num=len(nlist_new[type])
        if num:
            print '*',type,'-',num
        for elm in nlist_new[type]:
            all.append(elm)

    #python2.4 specific operation
    #sort netlist by order of appearance in orig. netlist
    all.sort(cmp=lambda x,y: cmp(x.num,y.num))
    for x in all:
        print>>ofp,x
#@-node:main
#@-others

#magic script-maker
if __name__ == '__main__':
    main()
#@-node:@file pyspice.py
#@-leo
