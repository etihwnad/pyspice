#!/usr/bin/env python

import os
import numpy as np
import subprocess as sp
import struct
from pylab import *


NaN = float('NaN')

def zc(x):
    """Return the flat index array of zero-crossings of x.  Index
    is floor(true-zc-idx)."""
    # entries identical to zero
    ez = list(np.flatnonzero(x == 0))

    # sign change, pass thru zero
    s = np.sign(x)
    cz = list(np.flatnonzero(s[1:] != s[:-1]))

    # there is an entry in both lists for any ez elements, remove from cz
    for e in ez:
        if e-1 in cz:
            cz.remove(e-1)
            cz.remove(e)
    y = np.concatenate((ez, cz))
    y.sort()
    return y


class SimulationData:
    """Base class for holding simulation data.  Data held in a numpy.array.
    Access data via attributes or d.data array.
    """
    def __init__(self, infile=None):
        #defaults only, all set by subclasses
        self.infile = infile
        self.cols = None
        self.colset = None
        self.siglist = []
        self._sig2idx = {}
        self.sweep = {}
        self.sweepvals = []
        self.sweepvar = None

        if infile:
            self.loadData(infile)

    def __getattr__(self, attr):
        """Returns the named signal, sliced with the current self.xrange."""
        if attr in self.siglist:
            idx = self._sig2idx[attr]
            return self.data[self._slice, idx]
        elif attr + '_0' in self.siglist:
            # auto return a complex vector if 'attr' is really attr_0 attr_1
            idx0 = self._sig2idx[attr+'_0']
            idx1 = self._sig2idx[attr+'_1']
            return (self.data[self._slice, idx0] + 
                    1j*self.data[self._slice, idx1])
        else:
            raise AttributeError


    #def __getstate__(self):
        #print '*** in getstate'
        #odict = self.__dict__.copy()
        #print odict['siglist']

        #npzfile = '.' + self.infile + '.npz'
        #odict['npzfile'] = npzfile
        #savestr = 'np.savez(npzfile'
        #savenames = []
        #for k,v in self.__dict__.iteritems():
            #print k, type(v)
            #if isinstance(v, np.ndarray):
                #print '* saving', k
                #savenames.append(k)
                #savestr += ', %s=self.%s' % (k,k)
                #del odict[k]
        #savestr += ')'
        #print savestr
        #exec(savestr)
        #return odict

    #def __setstate__(self, dict):
        #print '*** in setstate'
        #npz = np.load(dict['npzfile'])
        #print 'saved arrays:', npz.files
        #print dict.keys()
        #for v in npz.files:
            #setattr(self, v, npz[v])
        #self.__dict__.update(dict)

    #def __getinitargs__(self):
        #return (self.infile, )

    def getSweep(self, name):
        """Return another SimulationData object with a custom view of the
        rows of self.data and corresponding self._ivar.  Does not
        copy data, just references a slice of self.data."""
        # return cached value
        #sweep index
        if name in self.sweep:
            return self.sweep[name]
        #str(sweep_value)
        elif str(name) in self.sweep:
            return self.sweep[str(name)]

        # cache and return new sweep instance
        if isinstance(name, int) and name >= 0 and name < len(self.sweepvals):
            idx = name
        elif str(name) in map(str, self.sweepvals):
            idx = self.sweepvals.index(name)
            name = str(name)
        else:
            raise ValueError('OOPS: %s is not a valid sweep name.' % repr(name))

        c = SimulationData()
        c.cols = self.cols
        c.colset = self.colset
        c.siglist = self.siglist
        c._sig2idx = self._sig2idx
        c.sweepvar = self.sweepvar
        c.sweepval = self.sweepvals[idx]
        rowidx = (self.data[:,0] == c.sweepval)
        c.data = self.data[rowidx,:]
        c._ivar = self._ivar[rowidx]
        #reset to new range
        #c.xrange()
        c.x = c._ivar
        self.sweep[name] = c
        return c

    def xrange(self, xr=None):
        """Set the current independent axis range as tuple (xmin, xmax).  If
        xr=None, reset the data slice to the full range."""
        if not xr:
            self._xlims = [0, len(self._ivar)+1]
            self._slice = slice(self._xlims[0], self._xlims[1])
        else:
            xmin, xmax = xr
            self._xlims = [xmin, xmax]
            zmin =  zc(self._ivar - xmin)
            zmax = zc(self._ivar - xmax)
            imin = zmin[0] if zmin.size else 0
            imax = zmax[0] if zmax.size else len(self._ivar+1)
            self._slice = slice(imin, imax)
        self.x = self._ivar[self._slice]




class GnucapData(SimulationData):
    """Collect and present a Gnucap simulation output file as a numpy.array with
    access to column labels.

    This is for the ASCII default format, the SBSO v0.1 format is easy
    but unimplemented.
    """
    def loadData(self, infile):
        fin = open(infile, 'rb')

        header = fin.readline()
        cols = header.split()
        self.cols = cols

        #load cache if exists
        npy = '.' + infile + '.npy'
        mtime = os.path.getmtime
        if os.path.exists(npy) and mtime(npy) >= mtime(infile):
            print 'GnucapData: loading cached', npy
            data = np.load(npy, 'r')
        else:
            print 'GnucapData: reading', infile, ' caching to', npy
            data = np.genfromtxt(fin)
            np.save(npy, data)
        
        # handle duplicate columns
        # only use the last-defined name position
        colset = []
        for i,c in enumerate(cols):
            if c not in cols[i+1:]:
                colset.append([i, c])

        self.colset = colset

        #make signal vectors numpy arrays
        self.siglist = []
        self._sig2idx = {}
        for i,name in colset:
            # TODO: np.float only, handle complex also
            #data[name] = np.array(data[name], dtype=np.float)

            #set as an attributs also
            n = name.replace('(', '')
            n = n.replace(')', '')
            n = n.replace('.', '_')
            self.siglist.append(n)
            self._sig2idx[n] = i

            # special independent variable column
            if name.startswith('#'):
                n = '#Sweep' if name == '#' else n
                setattr(self, n[1:], data[:,i])
                setattr(self, '_ivar', data[:,i])
            else:
                pass
                #setattr(self, n, data[:,i])

        self.data = data

        #default to full display range, init relevant attributes
        self.xrange()

        return self





class HspiceData(SimulationData):
    """Convert HSPICE sim data of .option post_version=9601 via sp2sp and
    import into numpy.array with column labels.
    """
    def loadData(self, infile):
        mtime = os.path.getmtime
        exists = os.path.exists

        #intialize vars
        colset = []
        self.siglist = []
        self._sig2idx = {}
        self.sweepvar = None

        # use gwave's converter to help out
        # make cache file with custom sp2sp from gwave-svn
        npy = infile + '.npy'
        if not exists(npy) or mtime(infile) > mtime(npy):
            print 'HspiceData: caching data to', npy
            sp2sp = sp.call(['sp2sp', '-c', 'numpy', '-o', npy, infile])

        #load cache file footer
        print 'HspiceData: loading cached', npy
        fnpy = open(npy, 'rb')
        fnpy.seek(-2, os.SEEK_END)
        dictlen = struct.unpack('<H', fnpy.read(2))[0]
        fnpy.seek(-dictlen, os.SEEK_END)
        s = fnpy.readline().lstrip()
        npyinfo = safe_eval(s)
        self.npy = npyinfo
        fnpy.close()

        cols = npyinfo['cols']
        self.cols = cols

        # load data as a numpy memmap
        data = np.load(npy, 'r')

        # check for sweep
        #TODO: only handles one sweep var, no nested
        nsweepvars = len(npyinfo['sweepvars'])
        if npyinfo['sweepvars']:
            #TODO: sweeprows has lots of (0,0)'s for a MC sweep
            self.sweepvar = npyinfo['sweepvars'][0]
            self.sweepvals = [data[first,0] for first,last in npyinfo['sweeprows'] if last != 0]
        #NOTE: the following is unreachable??, if sweep: 'sweepvars' exists from cache
        elif '.' in cols[0]:
            #parameter sweep
            self.sweepvar = 'param'
            self.sweepvals = sorted(set(data[:,0]))
        elif cols[0].lower() == cols[0]:
            self.sweepvar = cols[0]
            self.sweepvals = sorted(set(data[:,0]))
        elif 'MONTE_CARLO' == cols[0]:
            self.sweepvar = 'MC'
            self.sweepvals = sorted(set(data[:,0]))

        if self.sweepvar:
            print 'Contained sweeps:', self.sweepvar, map(str, self.sweepvals)

        # handle duplicate columns
        # only use the last-defined name position
        # XXX: does this happen with HSPICE data?
        for i,c in enumerate(cols):
            if c not in cols[i+1:]:
                #account for prepended sweep values
                colset.append([i+nsweepvars, c])

        self.colset = colset

        #make signal vectors numpy arrays
        for i,name in colset:
            # TODO: np.float only, handle complex also
            #data[name] = np.array(data[name], dtype=np.float)

            #set as an attribute also
            n = name.replace('(', '')
            #n = n.replace(')', '')
            n = n.replace('.', '_')
            self.siglist.append(n)
            self._sig2idx[n] = i

            if (self.sweepvar and i == 1) or (not self.sweepvar and i == 0):
                setattr(self, n, data[:,i])
                setattr(self, '_ivar', data[:,i])

        self.data = data

        #default to full display range, init relevant attributes
        #self.xrange()
        self.x = self._ivar

        # separate out sweeps into self.sweep[0] and self.sweep['0,0']
        if 0: #use getSweep() instead
            if self.sweepvar:
                for i,val in enumerate(self.sweepvals):
                    d = self.getSweep(i)
                    self.sweep[i] = d
                    self.sweep[str(self.sweepvals[i])] = self.sweep[i]

        return self


def loadSimData(dfile):
    if dfile.endswith('.dat'):
        return GnucapData(dfile)
    else:
        return HspiceData(dfile)


class SignalPlotter():
    def __init__(self, gcdata=None):
        self.gcdata = gcdata

    def __call__(self, ys, x=None, *args, **kwargs):
        """Plot the signal named by string ys, label the curve by the true
        signal name UOS.  Additional args are passed to plot().
        """
        if 'label' not in kwargs:
            idx = self.gcdata._sig2idx[ys]
            kwargs['label'] = self.gcdata.cols[idx]
        y = getattr(self.gcdata, ys)
        plot(self.gcdata.x, y, *args, **kwargs)
        legend(loc='best')

def plotsweep(d, exp, vals=None, ivar=None, globals=None, labelprefix='',
              plotter=None):
    interact = isinteractive()
    if interact:
        interactive(False)

    if not vals:
        vals = d.sweepvals

    for v in vals:
        s = d.getSweep(v)
        
        if isinstance(exp, str):
            if globals:
                y = eval(exp, globals, locals())
            else:
                y = eval(exp)
        else:
            y = exp(s)

        if plotter:
            p = plotter
        else:
            p = plot

        if ivar:
            if isinstance(ivar, str):
                if globals:
                    x = eval(exp, globals, locals())
                else:
                    x = eval(exp)
            else:
                x = exp(s)
        else:
            x = s.x

        p(x, y, label='%s%s=%g' % (labelprefix, s.sweepvar, s.sweepval))

    interactive(interact)
    legend(loc='best')


if __name__ == "__main__":
    import optparse
    import sys

    from pylab import *

    usage = 'usage: %prog [options] simdata.dat'
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-I', '--noninteractive', dest='interactive',
                      help='Do not drop into IPython shell after loading',
                      action='store_false', default=True)

    opts, args = parser.parse_args()

    fname = args[0]
    d = loadSimData(fname)

    #if fname.endswith('.dat'):
        #d = GnucapData(fname)
    #else:
        #d = HspiceData(fname)

    p = SignalPlotter(d)

    interactive(True)
    rcParams['axes.grid'] = True


    if opts.interactive:
        print d.siglist

        from IPython.Shell import IPShellEmbed
        ipshell = IPShellEmbed(banner='*** Dropping into IPython. ***')
        ipshell(header='', global_ns=globals())

