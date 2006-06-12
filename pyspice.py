#!/usr/bin/env python
#@+leo-ver=4-thin
#@+node:etihwnad.20060606092308.1:@thin pyspice.py
#@@first
#@@language python
#@@tabwidth -4
#@<< head docstring >>
#@+node:etihwnad.20060605200356.2:<< head docstring>>
"""
#@<< head >>
#@+node:etihwnad.20060605202504:<< head >>
pyspice.py v0.1

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
#@-node:etihwnad.20060605202504:<< head >>
#@nl
#@<< copyright >>
#@+node:etihwnad.20060605201903:<< copyright >>
Copyright Dan White 2006

Licensed by the GPL, see http://www.whiteaudio.com/soft/COPYING for details.

Some non-binding requests, however:
    -Send an email to the author <etihwnad AT gmail DOT com> if you
     find this useful.
    -Send improvements back to maybe turn this into a full python module
Basically, I just want to see if this is useful to anyone else besides myself
and conrtibute back to the python (and SPICE) community.  I found John
Sheahan's spicepp Perl script a great help and wrote this when its
functionality wasn't what I needed.
#@nonl
#@-node:etihwnad.20060605201903:<< copyright >>
#@nl
#@+others
#@+node:etihwnad.20060605202632:structures
 Data structures: (inaccurate)
  lines - list of tuples: (original_line, line.split(), line_type)
  caps  - dictionary of node pairs and capacitance between nodes
#@-node:etihwnad.20060605202632:structures
#@+node:etihwnad.20060605202632.1:todo
TODO:
xpreserve comments in position
-element class definitions
    -inductor
    -v-source
    -i-source
    -e-source (VCVS)
-preserve node namespaces (within subckts, libraries, etc.)
    -find illegal SPICE node name character to prepend namespace
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
#@nonl
#@-node:etihwnad.20060605202632.1:todo
#@-others
"""
#@@nocolor
#@-node:etihwnad.20060605200356.2:<< head docstring>>
#@nl
#@<< global imports >>
#@+node:etihwnad.20060605210612:<< global imports >>
import sys, getopt, textwrap, warnings
from optparse import OptionParser
import re

try:
    from cStringIO import StringIO as StringIO
except:
    from StringIO import StringIO as StringIO
#@-node:etihwnad.20060605210612:<< global imports >>
#@nl
#@<< global vars >>
#@+node:etihwnad.20060605205852:<< global vars >>
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
_ncombine_mosfets=0
_ncombine_res=0

#namespace tracking
# -not implemented yet
_current_scope=''
#@-node:etihwnad.20060605205852:<< global vars >>
#@nl

#@+at
# Debug options: string of debugging elements to print
# * - comments
# follow SPICE element names (c-capacitor, l-inductor)
#@-at
#@@c
dbg=''

#@+others
#@+node:etihwnad.20060609195838:Option processing
#@+node:etihwnad.20060605200356.36:options
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
    parser.add_option('-i',dest='infile',default='stdin',
                      help='Input SPICE file to be processed'
                           ', (default: %default)')
    parser.add_option('-o',dest='outfile',default='stdout',
                      help='Output file for changes'
                           ', (default: %default)')
    parser.add_option('-d',dest='dropcap',
                      default='10', metavar='X',
                      help='Drop all capacitors smaller than X fF'
                           ', (default: %default)')
    parser.add_option('-v',dest='v',action='store_true',
                      default=True,
                      help='Show info and debugging messages (default)')
    parser.add_option('-q',dest='v',action='store_false',
                      help='Suppress all messages on stderr')
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
    #return the option object
    return opt
#@-node:etihwnad.20060605200356.36:options
#@-node:etihwnad.20060609195838:Option processing
#@+node:etihwnad.20060605211347:classes
#@+node:etihwnad.20060605200356.3:class SpiceElement
class SpiceElement:
    """Base class for SPICE elements.
    Methods:
        __init__(self,line,num) -> SpiceElement
        __str__(self) -> string spice line
        drop() -> False
    """
    #@	@+others
    #@+node:etihwnad.20060605200356.4:__init__
    def __init__(self,line,num):
        """SpiceElement constructor
        line - netlist expanded line
        type - first 'word'
        num  - netlist line number (for keeping roughly the same
               order when printing modified netlist)
        """
        #accept lists of 'words' also; BE CAREFUL with this, though
        if isinstance(line,list):
            line=' '.join(line)
        self.line=line
        self.type=line[0]
        self.num=num
    #@-node:etihwnad.20060605200356.4:__init__
    #@+node:etihwnad.20060605200356.5:__str__
    def __str__(self):
        return wrapper.fill(self.line)
    #@-node:etihwnad.20060605200356.5:__str__
    #@+node:etihwnad.20060605200356.6:drop
    def drop(self,val=0,mode='<'):
        return False
    #@-node:etihwnad.20060605200356.6:drop
    #@-others
#@-node:etihwnad.20060605200356.3:class SpiceElement
#@+node:etihwnad.20060605200356.7:class CommentLine

class CommentLine(SpiceElement):
    """SPICE Comment line (/^\*.*/)
    """
    #@	@+others
    #@+node:etihwnad.20060605200356.8:__init__
    def __init__(self,line,num):
        SpiceElement.__init__(self,line,num)
        self.type='comment'
    #@-node:etihwnad.20060605200356.8:__init__
    #@-others
#@-node:etihwnad.20060605200356.7:class CommentLine
#@+node:etihwnad.20060605200356.9:class ControlElement

class ControlElement(SpiceElement):
    """Control statement object, no processing for now.
    
    Note: currently has no knowledge of blocks (.lib/.endl, .subckt/.ends)
    has only ONE node namespace, make sure subckt's have unique node names!
    """
    #@	@+others
    #@+node:etihwnad.20060605200356.10:__init__
    def __init__(self,line,num):
        SpiceElement.__init__(self,line,num)
        self.type='control'
    #@-node:etihwnad.20060605200356.10:__init__
    #@-others
#@-node:etihwnad.20060605200356.9:class ControlElement
#@+node:etihwnad.20060605200356.11:class Passive2NodeElement

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
    #@+node:etihwnad.20060605200356.12:__init__
    def __init__(self,line,num):
        if isinstance(line,str):
            line=line.split()
        self.line=line
        self.type='passive2'
        self.num=num
        self.n1=_current_scope+line[1]
        self.n2=_current_scope+line[2]
        self.value=unit(line[3])
        self.param=dict() #store x=y as dictionary
        for p in line[4:]:
            k,v=p.split('=')
            self.param[k]=unit(v)
    #@-node:etihwnad.20060605200356.12:__init__
    #@+node:etihwnad.20060605200356.13:__str__
    def __str__(self):
        s=StringIO()
        print>>s, self.line[0],self.n1,self.n2,self.value,
        for k,v in self.param.iteritems():
            if v==0: continue
            print>>s, k+'='+str(v),
        return wrapper.fill(s.getvalue())
    #@-node:etihwnad.20060605200356.13:__str__
    #@+node:etihwnad.20060605200356.14:drop
    def drop(self,val=0.0,mode='<'):
        """Indicate whether to drop the element from the list.
        Occurs iff (val 'mode' self.value)
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
    #@-node:etihwnad.20060605200356.14:drop
    #@-others
#@-node:etihwnad.20060605200356.11:class Passive2NodeElement
#@+node:etihwnad.20060605200356.15:class Capacitor

class Capacitor(Passive2NodeElement):
    """Assumes SPICE element line:

    cXXX n1 n2 value p1=val p2=val ...
    Provides:
        isparallel(other)
        combine(other)
    """
    #@	@+others
    #@+node:etihwnad.20060605200356.16:__init__
    def __init__(self,line,num):
        Passive2NodeElement.__init__(self,line,num)
        self.type='capacitor'
    #@-node:etihwnad.20060605200356.16:__init__
    #@+node:etihwnad.20060605200356.17:isparallel
    def isparallel(self,other):
        """Returns True if instance is parallel with other instance
        """
        if self.n1==other.n1 and self.n2==other.n2:
            return True
        elif self.n1==other.n2 and self.n2==other.n1:
            return True
        else:
            return False
    #@-node:etihwnad.20060605200356.17:isparallel
    #@+node:etihwnad.20060605200356.18:combine
    def combine(self,other):
        """Adds values if capacitors are in parallel, returns True if
        it combined them.
    
        NOTE:
         -Does not currently touch param dictionary when combining,
          just the values.  How should this be done?
        """
        global _ncombine_capacitors
        if self.isparallel(other):
            self.value+=other.value
            _ncombine_capacitors+=1
            return True
        else:
            return False
    #@-node:etihwnad.20060605200356.18:combine
    #@-others
#@-node:etihwnad.20060605200356.15:class Capacitor
#@+node:etihwnad.20060605200356.19:class Mosfet

class Mosfet(SpiceElement):
    """Mosfet constructor takes an array derived from the
    netlist line
    """
    #@	@+others
    #@+node:etihwnad.20060605200356.20:__init__
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
    #@-node:etihwnad.20060605200356.20:__init__
    #@+node:etihwnad.20060605200356.21:__str__
    def __str__(self):
        s=StringIO()
        print>>s, self.line[0], self.d, self.g, self.s, self.b, self.model,
        for k,v in self.param.iteritems():
            if v==0: continue
            print>>s, k+'='+str(v),
        return wrapper.fill(s.getvalue())
    #@-node:etihwnad.20060605200356.21:__str__
    #@+node:etihwnad.20060605200356.22:isparallel
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
    
    #@-node:etihwnad.20060605200356.22:isparallel
    #@+node:etihwnad.20060605200356.23:combine
    def combine(self,other):
        """Combines adds other to self iff the transistors are identical,
        will NOT combine if W/L is different.
    
        Returns True if it combined the transistors
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
    #@-node:etihwnad.20060605200356.23:combine
    #@-others
#@-node:etihwnad.20060605200356.19:class Mosfet
#@+node:etihwnad.20060605200356.24:class Resistor

class Resistor(Passive2NodeElement):
    """Assumes SPICE element line:

    rXXX n1 n2 value p1=val p2=val ...
    """
    #@	@+others
    #@+node:etihwnad.20060605200356.25:__init__
    def __init__(self,line,num):
        Passive2NodeElement.__init__(self,line,num)
        self.type='resistor'
    #@-node:etihwnad.20060605200356.25:__init__
    #@-others
#@-node:etihwnad.20060605200356.24:class Resistor
#@-node:etihwnad.20060605211347:classes
#@+node:etihwnad.20060609195838.1:helpers
#@+node:etihwnad.20060612075426:unit
def unit(s):
    """Takes a string and returns the equivalent float.
    '3.0u' -> 3.0e-6"""
    mult={'t':1.0e12,
          'g':1.0e9,
          'meg':1.0e6,
          'k':1.0e3,
          'mil':25.4e-6,
          'm':1.0e-3,
          'u':1.0e-6,
          'n':1.0e-9,
          'p':1.0e-12,
          'f':1.0e-15}
    m=re.search('^([0-9e\+\-\.]+)(t|g|meg|k|mil|m|u|n|p|f)?',s.lower())
    if m.group(2):
        return float(m.group(1))*mult[m.group(2)]
    else:
        return float(m.group(1))
#@-node:etihwnad.20060612075426:unit
#@+node:etihwnad.20060605200356.27:debug
def debug(message):
    """Print debugging info to stderr."""
    print>>stderr,'Debug:',message
#@-node:etihwnad.20060605200356.27:debug
#@+node:etihwnad.20060605200356.26:info
def info(message):
    """Print information to stderr."""
    print>>stderr,'Info:',message
#@-node:etihwnad.20060605200356.26:info
#@+node:etihwnad.20060605200356.28:warning
def warning(message,elm=None,num=None):
    """Print warning to stderr.  If elm and num defined,
    print different message"""
    if elm and num:
        message=_opt.infile+":"+str(num)+" '"+elm+\
                "' type not defined yet, passing through..."
    print>>stderr,'Warning:',message
#@-node:etihwnad.20060605200356.28:warning
#@-node:etihwnad.20060609195838.1:helpers
#@+node:etihwnad.20060605200356.29:drop_2node
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
#@-node:etihwnad.20060605200356.29:drop_2node
#@+node:etihwnad.20060605200356.30:write_2node
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
#@-node:etihwnad.20060605200356.30:write_2node
#@+node:etihwnad.20060605200356.31:read_netlist
##
# Read netlist from open file object
#  -make sure if reading from stdin to read only once and
#   use this function to make sure you read the entire netlist
##
def read_netlist(fname):
    #@    << docstring >>
    #@+node:etihwnad.20060609200142:<< docstring >>
    """Read a SPICE netlist from the open file pointer
    
    read_netlist(filename) -> array lines
    
    Returns a list of expanded lines (without continuation '+')
    Keeps case of comments, all other lines are lowercased
    
    return:
        netlist (list) of SPICE netlist lines
    """
    #@nonl
    #@-node:etihwnad.20060609200142:<< docstring >>
    #@nl
    #@    << imports >>
    #@+node:etihwnad.20060609200142.1:<< imports >>
    import re
    #@nonl
    #@-node:etihwnad.20060609200142.1:<< imports >>
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
         #pass through empty lines and comments
        if not len(line.split()):
            #convert empty line to comment
            lines.append('*')
            nline+=1
            continue
        elif line[0]=='*':
            lines.append(line)
            nline+=1
            continue
        line=line.lower()
        # 'as = 3e-12' => 'as=3e-12'
        line=re.sub(re_param,r'\1=\2',line)
        if line[0]!='+': #beginning of SPICE line
            lines.append(line)
            nline+=1
        else:            #line continuation
            line=line[1:]
            lines[-1]=lines[-1]+line
    return lines
#@-node:etihwnad.20060605200356.31:read_netlist
#@+node:etihwnad.20060605200356.32:class ElementError

class ElementError(LookupError):
    #@	@+others
    #@+node:etihwnad.20060605200356.33:__init__
    def __init__(self,elm='???'):
        self.elm=elm
    #@-node:etihwnad.20060605200356.33:__init__
    #@+node:etihwnad.20060605200356.34:__str__
    def __str__(self):
        return str('No class defined for this element: '+self.elm)
    #@-node:etihwnad.20060605200356.34:__str__
    #@-others
#@-node:etihwnad.20060605200356.32:class ElementError
#@+node:etihwnad.20060605200356.35:classify
def classify(net):
    """Reads expanded netlist and classifies each line,
    calls appropriate function to combine nodes.
    
    net - list of unwrapped SPICE netlist lines
    """
    import string
    elements=dict()
    global _counts
    #initialize to all element types
    for ch in '*.'+string.lowercase:
        elements[ch]=[]
        _counts[ch]=0
    num=-1 #line counter
    for line in net:
        num+=1
        arr=line.split()
        x=arr[0][0]
        ##
        # Comment
        ##
        if x=='*':
            elements[x].append(CommentLine(line,num))
            _counts[x]+=1
        ##
        # Control line
        ##
        elif x=='.':
            elements[x].append(ControlElement(line,num))
            _counts[x]+=1
        elif x=='a':
            elements[x].append(SpiceElement(line,num))
            warning('',elm=arr[0],num=num)
#            raise ElementError(x)
            _counts[x]+=1
        elif x=='b':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        ##
        # Capacitor
        ##
        elif x=='c':
            elm=Capacitor(line,num)
            #don't combine if it's the first encountered
            if _counts[x]==0:
                elements[x].append(elm)
#                capacitors.append(elm)
            else:
                #combine if parallel or add if unique
                if elements[x][-1].combine(elm):
                    pass
                else:
                    elements[x].append(elm)
            _counts[x]+=1
        elif x=='d':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        elif x=='e':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        elif x=='f':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        elif x=='g':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        elif x=='h':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        ##
        # Current Source
        ##
        elif x=='i':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        elif x=='j':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        elif x=='k':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        ##
        # Inductor
        ##
        elif x=='l':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        ##
        # Mosfet
        ##
        elif x=='m':
            elm=Mosfet(line,num)
            if _counts[x]==0:
                #don't combine if it's the first encountered
                elements[x].append(elm)
                mosfets.append(elm)
            else:
                #search list backwards to find a parallel one
                #in case they aren't adjacent
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
            raise ElementError(x)
            _counts[x]+=1
        elif x=='o':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        elif x=='p':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        elif x=='q':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        ##
        # Resistor
        ##
        elif x=='r':
            elm=Resistor(line,num)
            elements[x].append(elm)
            _counts['r']+=1
        elif x=='s':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        elif x=='t':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        elif x=='u':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        ##
        # Voltage Source
        ##
        elif x=='v':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        elif x=='w':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        ##
        # Subcircuit
        ##
        elif x=='x':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        elif x=='y':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
            _counts[x]+=1
        elif x=='z':
            elements[x].append(SpiceElement(line,num))
            raise ElementError(x)
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
#@-node:etihwnad.20060605200356.35:classify
#@+node:etihwnad.20060605200356.37:main
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
    all.sort(cmp=lambda x,y: cmp(x.num,y.num))
    for x in all:
        print>>ofp,x
#@-node:etihwnad.20060605200356.37:main
#@-others

#magic script-maker
if __name__ == '__main__':
    main()
#@-node:etihwnad.20060606092308.1:@thin pyspice.py
#@-leo
