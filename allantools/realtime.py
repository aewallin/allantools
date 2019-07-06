"""
    Real-time statistics 
    ====================
    
    **Author:** Anders Wallin (anders.e.e.wallin "at" gmail.com)
    
    Statistics are computed 'on-the-fly' from a stream of phase/frequency samples
    Initial version 2019 July, Anders E. E. Wallin.
    
    This file is part of allantools, see https://github.com/aewallin/allantools
    
    
    License
    -------
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy

class oadev_realtime(object):
    """
    Overlapping Allan deviation in real-time from a stream of phase samples.
    
    Dobrogowski & Kasznia
    https://doi.org/10.1109/FREQ.2007.4319204
    """
    def __init__(self, afs=[1], tau0=1.0, auto_afs=False, pts_per_decade=4):
        self.x = []                         # phase
        self.afs = afs                      # averaging factor
        self.auto_afs = auto_afs
        self.af_taus = numpy.logspace(0,1,pts_per_decade+1)[:-1] # will fail at >=6 (?), need to remove duplicates?
        self.af_idx = 0
        self.af_decade = 0
        if auto_afs:
            self.afs = numpy.array([1])
        self.S = numpy.zeros(len(afs))      # sum-of-squares
        self.dev = numpy.zeros(len(afs))    # resulting xDEV
        self.tau = tau0                     # time-interval between points

    def update_af(self):
        """ used in auto-AF mode, 
            - check if we can add another AF
            - if yes, add it.
        """
        next_idx = self.af_idx+1
        next_decade = self.af_decade
        if next_idx == len(self.af_taus):
            next_idx = 0
            next_decade = self.af_decade + 1
        
        next_af = int( numpy.round( pow(10.0, next_decade) * self.af_taus[ next_idx ] ) ) # next possible AF
        if len(self.x) >= (2*next_af+1): # can compute next AF
            #print "new AF ", next_af
            self.afs = numpy.append( self.afs, next_af ) # new AF
            self.S = numpy.append( self.S, 0 ) # new S
            self.dev = numpy.append( self.dev, 0 ) # new tau
            self.af_idx = next_idx #self.af_idx+1
            self.af_decade = next_decade
        else:
            pass
            #print "no new AF "

    def add_phase(self, x):
        """ add new phase point, in units of seconds """
        self.x.append(x)
        #print self.x
        for idx, af in enumerate(self.afs):
            if len(self.x) >= (2*af+1):
                self.update_S(idx)
        if self.auto_afs:
            self.update_af()
            
    def update_S(self, idx):
        """ update S, sum-of-squares """
        af = self.afs[idx]
        
        i = len(self.x)-1 # last pt
        S_new = pow( self.x[i] - 2*self.x[i-af] + self.x[i-2*af], 2)
        #print "updateS af=",af," i=",i," Snew=",S_new
        
        self.S[idx] = self.S[idx] + S_new
        self.dev[idx] = math.sqrt( (1.0/(2*pow(af*self.tau,2)*(i+1-2*af))) * self.S[idx] )
        
    def taus(self):
        """ return taus, in unit of seconds """
        return self.tau*self.afs
    
    def devs(self):
        """ return deviation """
        return self.dev
        
    def __str__(self):
        msg = "n_pts: %d " % len(self.x)
        for idx, af in enumerate(self.afs):
            msg += "ADEV(%d)=%f " %(af, self.dev[idx])
        return msg

# end of file realtime.py
