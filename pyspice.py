#!/usr/bin/env python

__version__ = '0.3'

"""
pyspice.py v0.3

SPICE pre-processor that combines parallel elements (e.g. capacitors, mosfets)
for GREATLY reduced simulation time.  Uses the 'M' parameter of MOSFETS
when combining parallel FETs.  Makes simulating fingered transistors easier and
faster.

 Combine parallel capacitors
 Drop combined caps smaller than X fF
 Usage:
  pyspice.py [options] [-i infile] [-o outfile]

 Use pyspice.py -h for all the options.
Copyright Dan White 2006-8

Licensed by the GPL, see http://www.whiteaudio.com/soft/COPYING or the current
GNU GPL license for details.

Release Notes, changelog
-----------------------------------------------------
pyspice.py v0.3:
----------------
-entire netlist is held in the top level Netlist object
-temporarily does not drop small capacitors
-code cleanup

pyspice.py v0.2a:
----------------
-added a missing newline before an import statement

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

pyspice v0.1:
-------------
Initial release.
Only worked for netlist containing MOSFETs and Capacitors.
"""

##
# Global imports
##
import getopt
import re
import sys
import textwrap
import warnings

from decimal import Decimal
from optparse import OptionParser

try:
    from cStringIO import StringIO as StringIO
except:
    from StringIO import StringIO as StringIO
##
# Global Shorthands
##
stderr = sys.stderr

#set line continuation style and max width
_wrapper = textwrap.TextWrapper(subsequent_indent = '+ ', width = 75)


#global variables
_counts = dict() #keyed by spice element letter
_ncombine_capacitors = 0
_ncombine_inductors = 0
_ncombine_mosfets = 0
_ncombine_res = 0

#namespace tracking
# -not implemented yet
_current_scope = ''

# Debug options: string of debugging elements to print
# * - comments
# follow SPICE element names (c-capacitor, l-inductor)
dbg = ''

##
# Fancy option processing
##
def options(args=sys.argv):
    """Define options and parse argument list.
    """
    global ifp, ofp
    usage = """%prog [options]"""
    desc = """This script reads a SPICE input files and processes it according
    to the options given.  It is especially useful for processing netlists
    from the output of layout extractors by combining parallel C's and
    FET's.  More features added upon request."""

    parser = OptionParser(usage=usage, description=desc)

    parser.add_option('-i', '--infile', dest='infile', default='stdin',
                      help='Input SPICE file to be processed'
                           ', (default: %default)')

    parser.add_option('-o', '--outifle', dest='outfile', default='stdout',
                      help='Output file for changes'
                           ', (default: %default)')

    parser.add_option('-d', '--dropcap', dest='dropcap',
                      default='10', metavar='X',
                      help='Drop all capacitors smaller than X fF'
                           ', (default: %default)')

    parser.add_option('--no-combine-c', dest='combine_c', action='store_false',
                      default=True, help='Do not combine parallel capacitors')

    parser.add_option('--no-combine-m', dest='combine_m', action='store_false',
                      default=True, help='Do not combine parallel MOSFETs')

    parser.add_option('-v', '--verbose', dest='v', action='store_true',
                      default=True, help='Show info and debugging messages (default)')

    parser.add_option('-q', '--quiet', dest='v', action='store_false',
                      help='Suppress all messages on stderr')

    parser.add_option('-w', '--linewidth', dest='linewidth', default=75,
                      help='Max. line width for netlist (default: %default)')

    parser.add_option('--case', dest='case', action='store',
                      choices=('keep', 'lower', 'upper'), default='keep',
                      help='Modify capitalization of lines [keep, lower, upper] (default: %default)')

    (opt, args) = parser.parse_args()

    #infile
    if opt.infile == 'stdin':
        ifp = sys.stdin
        if opt.v: info('Input: stdin')
    else:
        try:
            ifp = open(opt.infile, 'rU') #python's universal line ending mode
            if opt.v: info('Input: '+opt.infile)
        except IOError, (errno, strerror):
            print>>stderr, "IOError(%s): %s '%s'" % (errno, strerror, opt.infile)      

    #outfile
    if opt.outfile == 'stdout':
        ofp = sys.stdout
        if opt.v: info('Output: stdout')
    else:
        try:
            ofp = open(opt.outfile, 'w')
            if opt.v: info('Output: ' + opt.outfile)
        except IOError, (errno, strerror):
            print>>stderr, "IOError(%s): %s '%s'" % (errno, strerror, opt.outfile)

    #dropcap
    opt.dropcap = float(opt.dropcap)
    opt.dropcap = opt.dropcap * 1e-15
    if opt.v:
        info('Dropping caps < ' + str(opt.dropcap) + ' F')
        info("Combining c's: " + str(opt.combine_c))
        info("Combining m's: " + str(opt.combine_m))

    #linewidth
    _wrapper.line_width = opt.linewidth

    return opt


class PyspiceError(Exception):
    '''Base exception for pyspice'''
    pass

class BadUnitError(PyspiceError):
    pass


class ElementError(LookupError):
    def __init__(self, elm='???'):
        self.elm = elm

    def __str__(self):
        return str('No class defined for this element: '+self.elm)



class Netlist:
    """Base class that holds an entire netlist.

    Notes:
        -this will eventually hold the entire shebang
        -providing __init__ with a file name will:
            -read in netlist
            -classify the lines
            -take care of hierarchy
            -source other files
    """
    def __init__(self, fname=None, title=None):
        """Optionally reads a netlist from a file"""
        self.deck = []
        self.lines = []
        self.title = title

        self.elements = dict()
        for e in _elementHandler.validTypes:
            self.elements[e] = []

        if fname:
            self.readfile(fname)
    def _addMassagedLine(self, line):
        """Add the given non-empty line to netlist.  Assumes the line has already
        been massaged."""
        self.lines.append(line)
    def addElement(self, element):
        self.deck.append(element)
        self.elements[element.type].append(element)
    def addLine(self, line):
        """Add the given non-empty line to netlist after massaging"""
        if line:
            #fail on line continuations
            if line[0] == '+':
                raise PyspiceError('addLine does not handle line continuations')
            else:
                line = self.massageLine(line)
                self._addMassagedLine(line)
    def classify(self, line, num=None):
        """Takes a line and creates an appropriate SpiceElement"""
        return _elementHandler.handler[line[0][0].lower()](line, num=num)
    #finds a "name = value" pair for shrinking
    RE_PARAM = re.compile(r"(\S*)\s*=\s*(\S*)")

    #finds only whitespace
    RE_WHITESPACE_EMPTY = re.compile(r'^\s*$')

    def massageLine(self, line):
        #remove trailing newline
        line = line.strip('\r\n')

        #pass through empty lines and convert to comments
        if self.RE_WHITESPACE_EMPTY.search(line):
            return '*'
        # and pass through comments, they stay as-is
        elif line[0] == '*':
            return line

        if _opt == 'keep':
            pass
        elif _opt == 'lower':
            #case is unimportant in SPICE
            #lowercase all non-comment lines
            line = line.lower()
        elif _opt == 'upper':
            line = line.upper()

        #remove whitespace in parameter assignments
        # to prepare for x.split(' ') that happens next:
        #  'as = 3e-12' => 'as=3e-12'
        line = self.RE_PARAM.sub(r'\1=\2', line)

        return line
    def readfile(self, fname):
        """Read a SPICE netlist from the open file pointer into the netlist.

        Reads the file as a netlist into the deck.  Appends the line's text
        and adds a classified SpiceElement to the deck.

        return:
            None

        Note:
            -we need to read at least a full line with continuations before we can
             add the line to the netlist
        """

        if isinstance(fname, file):
            ifp = fname
        else:
            ifp = open(fname, 'rU')

        #first line of any file is ignored
        # typically a title line
        # set title iff title is not set
        line = ifp.readline()
        if not self.title:
            self.title = line

        currentCard = ifp.readline()
        n = 2 #1-indexed line numbers
        for line in ifp:
            n += 1
            #handle line continuations here
            if line[0] == '+':
                currentCard += line[1:]
                continue

            #a new card is started, the previous card is
            #unambiguously finished
            mLine = self.massageLine(currentCard)
            self._addMassagedLine(mLine)
            self.addElement(self.classify(mLine, num=n))
            currentCard = line
    def removeElement(self, element):
        self.deck.remove(element)
        self.elements[element.type].remove(element)
class ElementHandler:
    def __init__(self):
        self.validTypes = '*.abcdefghijklmnopqrstuvwxyz'
        self.handler = dict()
        for t in self.validTypes:
            self.handler[t] = SpiceElement
    def add_handler(self, type, handler):
        """Replaces the existing element object definition with the
        given one"""
        self.handler[type] = handler
# These classes are used to break down the spectrum of SPICE elements
# into 'classes' of elements.  I.e. 2 node passives, 2-node sources, 4-node 
# sources,
# and so on.  The elements found in a real netlist are based on these types.

class SpiceElement:
    """Base class for SPICE elements.
    Methods:
        __init__(self, line, num) -> SpiceElement
        __str__(self) -> string spice line
        drop() -> False
    """
    def __init__(self, line, num=None):
        """SpiceElement constructor
        line - netlist expanded line
        type - first character
        num  - input netlist line number (for keeping roughly the same
               order when printing modified netlist)
        """
        #accept lists of 'words' also; BE CAREFUL with this, though
        self.line = line
        self.type = 'spice'
        self.typeName = 'SpiceElement'
        self.num = num
    def __str__(self):
        return _wrapper.fill(self.line)
    def drop(self, val=0, mode='<'):
        """Template for dropping elements that defaults to NO if
        not overidden in the element class"""
        return False

class Passive2NodeElement(SpiceElement):
    """Base class for 2-node elements.
    Assumes SPICE element line:
    xXXX n1 n2 value p1=val p2=val ...
    Inherits:
        None
    Redefines:
        drop(self, val, mode) -> bool
    """
    def __init__(self, line, num):
        SpiceElement.__init__(self, line, num)
        self.type = 'passive2'
        self.typeName = 'Passive2NodeElement'
        arr = line.split()
        self.name = arr[0]
        self.n1 = _current_scope + arr[1]
        self.n2 = _current_scope + arr[2]
        self.value = unit(arr[3])
        self.param = dict() #store x = y as dictionary
        for p in arr[4:]:
            k, v = p.split('=')
            self.param[k] = unit(v)
    def __str__(self):
        """Returns the netlist-file representation of this element"""
        s = StringIO()

        print>>s, self.name, self.n1, self.n2, self.value,

        for k, v in self.param.iteritems():
            print>>s, k+'='+str(v),

        return _wrapper.fill(s.getvalue())
    def drop(self, val=0.0, mode='<'):
        """Indicate whether to drop the element from the list.
        Occurs iff (val 'mode' self.value)

        Can this be converted to specifying an arbitrary binary function?
          This may allow a more elegant comparison.
          e.g. mode = < instead of mode = '<' or mode = cap_smaller(x, y)
        """
        if mode == '<':
            if self.value < val: return True
        elif mode == '<=':
            if self.value <= val: return True
        elif mode == '>':
            if self.value > val: return True
        elif mode == '>=':
            if self.value >= val: return True
        else:
            return False
        return False #shouldn't get here, but...

class Active2NodeElement(SpiceElement):
    """Base class for active 2-node elements.
    Assumes SPICE element line:
    xXXX n1 n2 value p1=val p2=val ...
    Inherits:
        None
    Redefines:
        None
    """
    def __init__(self, line, num):
        SpiceElement.__init__(self, line, num)
        self.type = 'active2'
        self.typeName = 'Active2NodeElement'
        arr = line.split()
        self.name = arr[0]
        self.n1 = _current_scope + arr[1]
        self.n2 = _current_scope + arr[2]
        self.value = unit(arr[3])
        self.param = dict() #store x = y as dictionary
        for p in arr[4:]:
            k, v = p.split('=')
            self.param[k] = unit(v)
    def __str__(self):
        """Returns the netlist-file representation of this element"""
        s = StringIO()

        print>>s, self.name, self.n1, self.n2, self.value,

        for k, v in self.param.iteritems():
            print>>s, k + '=' + str(v),

        return _wrapper.fill(s.getvalue())

class Active4NodeElement(SpiceElement):
    """Base class for active 4-node elements (xCyS).
    Assumes SPICE element line:
    xXXX n1 n2 value p1=val p2=val ...
    Inherits:
        None
    Redefines:
        None
    """
    def __init__(self, line, num):
        SpiceElement.__init__(self, line, num)
        self.type = 'active4'
        self.typeName = 'Active4NodeElement'
        arr = line.split()
        self.name = arr[0]
        self.n1 = _current_scope + arr[1]
        self.n2 = _current_scope + arr[2]
        self.n3 = _current_scope + arr[3]
        self.n4 = _current_scope + arr[4]
        self.value = unit(arr[5])
        self.param = dict() #store x = y as dictionary
        for p in arr[6:]:
            k, v = p.split('=')
            self.param[k] = unit(v)
    def __str__(self):
        s = StringIO()
        print>>s, self.name, self.n1, self.n2, self.n3, self.n4, self.value,
        for k, v in self.param.iteritems():
            #are there instances when 0 is (in)significant?
            #if v == 0: continue
            print>>s, k + '=' + str(v),
        return _wrapper.fill(s.getvalue())
# This is a(n incomplete) definition of the various SPICE elements.
# 
# NOTE: When adding a new element type definition, be sure to add a handler
#   for the new class after defining the class using:
#       elements.add_handler('x', Xdevice)
#make a repository for element handlers
_elementHandler = ElementHandler()

class CommentLine(SpiceElement):
    """SPICE Comment line (/^\*.*/)
    """
    def __init__(self, line, num):
        SpiceElement.__init__(self, line, num)
        self.type = '*'
        self.typeName = ''
_elementHandler.add_handler('*', CommentLine)

class ControlElement(SpiceElement):
    """Control statement object, no processing for now.
    Note: this will eventially be a base class for the real control elements

    Note: currently has no knowledge of blocks (.lib/.endl, .subckt/.ends)
    has only ONE node namespace, make sure subckt's have unique node names!
    """
    def __init__(self, line, num):
        SpiceElement.__init__(self, line, num)
        self.type = '.'
        self.typeName = 'ControlElement'
_elementHandler.add_handler('.', ControlElement)

class Capacitor(Passive2NodeElement):
    """Assumes SPICE element line:

    cXXX n1 n2 value p1=val p2=val ...
    Provides:
        isparallel(other)
        combine(other)
    """
    def __init__(self, line, num):
        Passive2NodeElement.__init__(self, line, num)
        self.type = 'c'
        self.typeName = 'Capacitor'
    def isparallel(self, other):
        """Returns True if instance is parallel with other instance
        """
        if self.n1 == other.n1 and self.n2 == other.n2:
            return True
        elif self.n1 == other.n2 and self.n2 == other.n1:
            return True
        else:
            return False
    def combine(self, other):
        """Adds values if capacitors are in parallel, returns True if
        it combined them.

        NOTE: Does not currently touch param dictionary when combining,
          just the values.  How should this be done?  Maybe combine iff
          params are identical to avoid problems?
        """
        global _ncombine_capacitors
        if self.isparallel(other):
            self.value += other.value
            _ncombine_capacitors += 1
            return True
        else:
            return False
_elementHandler.add_handler('c', Capacitor)

class Inductor(Passive2NodeElement):
    """Assumes SPICE element line:

    cXXX n1 n2 value p1=val p2=val ...
    Provides:
        isparallel(other)
        combine(other)
    """
    def __init__(self, line, num):
        Passive2NodeElement.__init__(self, line, num)
        self.type = 'l'
        self.typeName = 'Inductor'
    def isparallel(self, other):
        """Returns True if instance is parallel with other instance
        """
        if self.n1 == other.n1 and self.n2 == other.n2:
            return True
        elif self.n1 == other.n2 and self.n2 == other.n1:
            return True
        else:
            return False
    def combine(self, other):
        """Combines values if inductors are in parallel, returns True if
        it combined them.

        NOTE:
         -Does not currently touch param dictionary when combining,
          just the values.  How should this be done?
        """
        global _ncombine_inductors
        if self.isparallel(other):
            self.value = (self.value*other.value)/(self.value+other.value)
            _ncombine_inductors += 1
            return True
        else:
            return False
_elementHandler.add_handler('l', Inductor)

class Mosfet(SpiceElement):
    """Mosfet constructor takes an array derived from the
    netlist line
    """
    def __init__(self, line, num):
        if isinstance(line, str):
            line = line.split()
        self.line = line
        self.type = 'm'
        self.typeName = 'Mosfet'
        self.num = num
        self.name = line[0]
        self.d = line[1]
        self.g = line[2]
        self.s = line[3]
        self.b = line[4]
        self.model = line[5]
        self.param = dict()
        for p in line[6:]:
            k, v = p.split('=')
            self.param[k.lower()] = unit(v)
        self.w = self.param['w']
        self.l = self.param['l']
    def __str__(self):
        s = StringIO()
        print>>s, self.name, self.d, self.g, self.s, self.b, self.model,
        for k, v in self.param.iteritems():
            #TODO really kosher to ignore 0-values, careful of non-zero defaults
            if v == 0: continue
            print>>s, k + '=' + str(v),
        return _wrapper.fill(s.getvalue())
    def isparallel(self, other):
        """Returns True if transistors are parallel
        """
        # check gate, substrate, and model first
        if self.g == other.g and self.b == other.b and self.model == other.model:
            #source and drain can be reversed
            if self.d == other.d and self.s == other.s:
                return True
            elif self.d == other.s and self.s == other.d:
                return True
            else:
                return False
        else:
            return False

    def combine(self, other):
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
            if self.w == other.w and self.l == other.l:
                for k, v in other.param.iteritems():
                    if k == 'w' or k == 'l' or k == 'm': continue
                    self.param[k] += v
                if ('m' in self.param.keys()) or ('m' in other.param.keys()):
                    #add other's M parameter or increment
                    self.param['m'] += other.param.get('m', 1)
                else:
                    self.param['m'] = 2

                _ncombine_mosfets += 1
                return True
        else:
            return False
_elementHandler.add_handler('m', Mosfet)

class Resistor(Passive2NodeElement):
    """Assumes SPICE element line:

    rXXX n1 n2 value p1=val p2=val ...
    """
    def __init__(self, line, num):
        Passive2NodeElement.__init__(self, line, num)
        self.type = 'r'
        self.typeName = 'Resistor'
_elementHandler.add_handler('r', Resistor)

class Vsource(Active2NodeElement):
    """Assumes SPICE element line:

    vXXX n1 n2 value p1=val p2=val ...
    """
    def __init__(self, line, num):
        Active2NodeElement.__init__(self, line, num)
        self.type = 'v'
        self.typeName = 'Vsource'
_elementHandler.add_handler('v', Vsource)

class Isource(Active2NodeElement):
    """Assumes SPICE element line:

    iXXX n1 n2 value p1=val p2=val ...
    """
    def __init__(self, line, num):
        Active2NodeElement.__init__(self, line, num)
        self.type = 'i'
        self.typeName = 'Isource'
_elementHandler.add_handler('i', Isource)
def combineCapacitorsInplace(nlist):
    '''Finds all parallel capacitors and replaces each with a single element
    of equivalent value.  The capacitor is named by the first-occuring name.
    Returns the number of combined capacitors.'''

    caps = [c for c in nlist.deck if c.type == 'c']

    #this modifies the list being iterated over in place
    #usually this is BAD, here it is our way of only checking capacitor
    #combinations for parallel-isity once, the netlist is modified in parallel.
    #
    #This works because different instances of the same class never compare
    #equal
    n = 0
    for c in caps:
        for x in caps[caps.index(c)+1:]:
            if c.combine(x):
                n += 1
                caps.remove(x)
                nlist.removeElement(x)

    return n
def combineMosfetsInplace(nlist):
    '''TODO'''

    fets = [m for m in nlist.deck if m.type == 'm']

    n = 0
    for m in fets:
        for x in fets[fets.index(m)+1:]:
            if m.combine(x):
                n += 1
                fets.remove(x)
                nlist.removeElement(x)

    return n
RE_UNIT = re.compile(r'^([0-9e\+\-\.]+)(t|g|meg|x|k|mil|m|u|n|p|f|a)?')
def unit(s):
    """Takes a string and returns the equivalent float.
    '3.0u' -> 3.0e-6"""
    mult = {'t'  :Decimal('1.0e12'),
            'g'  :Decimal('1.0e9'),
            'meg':Decimal('1.0e6'),
            'x'  :Decimal('1.0e6'),
            'k'  :Decimal('1.0e3'),
            'mil':Decimal('25.4e-6'),
            'm'  :Decimal('1.0e-3'),
            'u'  :Decimal('1.0e-6'),
            'n'  :Decimal('1.0e-9'),
            'p'  :Decimal('1.0e-12'),
            'f'  :Decimal('1.0e-15'),
            'a'  :Decimal('1.0e-18')}

    m = RE_UNIT.search(s.lower())
    try:
        if m.group(2):
            return Decimal(Decimal(m.group(1)))*mult[m.group(2)]
        else:
            return Decimal(m.group(1))
    except:
        raise BadUnitError
def debug(message):
    """Print debugging info to stderr."""
    for m in message.split('\n'):
        if m:
            print>>stderr, 'Debug:', m
def info(message):
    """Print information to stderr."""
    for m in message.split('\n'):
        if m:
            print>>stderr, 'Info:', m
def warning(message, elm=None, num=None):
    """Print warning to stderr.  If elm and num defined,
    print different message"""
    if elm and num:
        message = _opt.infile+":"+str(num)+" '"+elm+\
                "' type not defined yet, passing through..."
    print>>stderr, 'Warning:', message
def main():
    global _opt
    opt = options()
    _opt = opt

    # output file header
    print>>ofp, "* pyspice.py %s: by Dan White <dan@whiteaudio.com>" % __version__
    print>>ofp, "* mail me bug reports, fixes, and comments if you find this useful"
    print>>ofp, "* ----------------------------------------------------------------"

    # Read and parse given input file (as top-level)
    netlist = Netlist(ifp)

    # Show input statistics
    if opt.v:
        info('Read in %i elements' % len(netlist.deck))
        s = StringIO()
        print>>s, 'Input Element counts:'
        for t, v in netlist.elements.iteritems():
            if len(v):
                print>>s, '%s: %i' % (t, len(v))
        info(s.getvalue())


    # Combine elements if requested
    nCombined = dict()
    if opt.combine_c:
        nCombined['c'] = combineCapacitorsInplace(netlist)
        if opt.v:
            info('Combined %i capacitors' % nCombined['c'])

    if opt.combine_m:
        nCombined['m'] = combineMosfetsInplace(netlist)
        if opt.v:
            info('Combined %i mosfets' % nCombined['m'])

    # Show output statistics
    if opt.v:
        s = StringIO()
        print>>s, 'Output Element counts:'
        for t, v in netlist.elements.iteritems():
            if len(v):
                print>>s, '%s: %i' % (t, len(v))
        info(s.getvalue())

    for card in netlist.deck:
        print>>ofp, card


#magic script-maker
if __name__ == '__main__':
    main()
