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

class dev_realtime(object):
    """ Base-class for real-time statistics """
    def __init__(self, afs=[1], tau0=1.0, auto_afs=False, pts_per_decade=4):
        self.x = []                         # phase time-series
        self.afs = afs                      # averaging factor, tau = af*tau0
        self.auto_afs = auto_afs
        self.af_taus = numpy.logspace(0, 1, pts_per_decade+1)[:-1] # will fail at >=6 (?), need to remove duplicates?
        self.af_idx = 0       # logspace index, to keep track of afs in auto-af mode
        self.af_decade = 0    # logspace decade
        if auto_afs:
            self.afs = numpy.array([1])
        self.dev = numpy.zeros(len(afs))    # resulting xDEV
        self.tau0 = tau0                     # time-interval between points

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

        next_af = int(numpy.round(pow(10.0, next_decade) * self.af_taus[next_idx])) # next possible AF
        if len(self.x) >= (2*next_af+1): # can compute next AF
            self.afs = numpy.append(self.afs, next_af) # new AF
            self.S = numpy.append(self.S, 0) # new S
            self.dev = numpy.append(self.dev, 0) # new dev
            self.af_idx = next_idx
            self.af_decade = next_decade
        else:
            pass
            #print "no new AF "

    def add_frequency(self, f):
        """ add new frequency point, in units of Hz """
        if not self.x: # empty sequence
            self.add_phase(0) # initialize
        self.add_phase(self.x[-1] + f) # integration

    def taus(self):
        """ return taus, in unit of seconds """
        return self.tau0*self.afs

    def devs(self):
        """ return deviation """
        return self.dev

class oadev_realtime(dev_realtime):
    """ Overlapping Allan deviation in real-time from a stream of phase/frequency samples.

        Dobrogowski & Kasznia
        https://doi.org/10.1109/FREQ.2007.4319204
    """
    def __init__(self, afs=[1], tau0=1.0, auto_afs=False, pts_per_decade=4):
        super(oadev_realtime, self).__init__(afs=afs, tau0=tau0, auto_afs=auto_afs, pts_per_decade=pts_per_decade)
        self.S = numpy.zeros(len(afs))      # sum-of-squares

    def add_phase(self, xnew):
        """ add new phase point, in units of seconds """
        self.x.append(xnew)
        for idx, af in enumerate(self.afs):
            if len(self.x) >= (2*af+1):
                self.update_S(idx)
        if self.auto_afs:
            self.update_af()

    def update_S(self, idx):
        """ update S, sum-of-squares """
        af = self.afs[idx]
        i = len(self.x)-1 # last pt
        S_new = pow(self.x[i] - 2*self.x[i-af] + self.x[i-2*af], 2)
        self.S[idx] = self.S[idx] + S_new
        self.dev[idx] = numpy.sqrt((1.0/(2*pow(af*self.tau0, 2)*(i+1-2*af))) * self.S[idx])

    def __str__(self):
        msg = "n_pts: %d " % len(self.x)
        for idx, af in enumerate(self.afs):
            msg += "ADEV(%d)=%f " %(af, self.dev[idx])
        return msg

class ohdev_realtime(dev_realtime):
    """ Overlapping Hadamard deviation in real-time from a stream of phase/frequency samples.

        Dobrogowski & Kasznia
        https://doi.org/10.1109/FREQ.2007.4319204
    """
    def __init__(self, afs=[1], tau0=1.0, auto_afs=False, pts_per_decade=4):
        super(ohdev_realtime, self).__init__(afs=afs, tau0=tau0, auto_afs=auto_afs, pts_per_decade=pts_per_decade)
        self.S = numpy.zeros(len(afs))      # sum-of-squares

    def add_phase(self, xnew):
        """ add new phase point """
        self.x.append(xnew)
        for idx, af in enumerate(self.afs):
            if len(self.x) > 3*af:
                self.update_S(idx)

    def update_S(self, idx):
        """ update S, sum-of-squares """
        af = self.afs[idx]
        i = len(self.x)-1 # last pt
        #print i,self.x
        S_new = pow(self.x[i] - 3*self.x[i-af] + 3*self.x[i-2*af] - self.x[i-3*af], 2)
        self.S[idx] = self.S[idx] + S_new
        self.dev[idx] = numpy.sqrt((1.0/(6.0*pow(af*self.tau0, 2)*(i+1.0-3*af))) * self.S[idx])

    def __str__(self):
        msg = "n_pts: %d " % len(self.x)
        for idx, af in enumerate(self.afs):
            msg += "OHDEV(%d)=%f " %(af, self.dev[idx])
        return msg


class tdev_realtime(dev_realtime):
    """ Time deviation and Modified Allan deviation in real-time from a stream of phase/frequency samples.

        Dobrogowski & Kasznia
        https://doi.org/10.1109/FREQ.2007.4319204
    """
    def __init__(self, afs=[1], tau0=1.0, auto_afs=False, pts_per_decade=4):
        super(tdev_realtime, self).__init__(afs=afs, tau0=tau0, auto_afs=auto_afs, pts_per_decade=pts_per_decade)
        self.S = numpy.zeros(len(afs))      # sum-of-squares
        self.So = numpy.zeros(len(afs))      # overall sum-of-squares

    def add_phase(self, xnew):
        """ add new phase point """
        self.x.append(xnew)
        for idx, af in enumerate(self.afs):
            if len(self.x) >= 3*af+1:  # 3n+1 samples measured
                self.update_S(idx)
            elif len(self.x) >= 2*af+1: # 2n+1 samples measured
                self.update_S3n(idx)

    def update_S3n(self, idx):
        """ eqn (13) of paper """
        af = self.afs[idx]
        j = len(self.x)-1 # last pt
        self.S[idx] = self.S[idx] + self.x[j] - 2*self.x[j-af] + self.x[j-2*af]
        if len(self.x) == 3*af:
            # last call to this fctn
            self.So[idx] = pow(self.S[idx], 2)
            self.update_dev(idx)

    def update_dev(self, idx):
        # Eqn (14)
        num_pts = len(self.x)
        af = self.afs[idx]
        self.dev[idx] = numpy.sqrt((1.0/6.0)*(1.0/(num_pts-3*af+1.0))*(1.0/pow(af, 2))*(self.So[idx]))

    def update_S(self, idx):
        """ update S, sum-of-squares """
        af = self.afs[idx]
        assert(len(self.x) >= 3*af+1)
        i = len(self.x)-1 # last pt
        # Eqn (12)
        S_new = -1*self.x[i-3*af] + 3*self.x[i-2*af] - 3*self.x[i-af] + self.x[i]

        self.S[idx] = self.S[idx] + S_new
        # Eqn (11)
        self.So[idx] = self.So[idx] + pow(self.S[idx], 2) #??? S_(i-1) in paper for TDEV-sqrt?
        self.update_dev(idx)

    def mdev(self):
        """ scale tdev to output mdev """
        mdev = self.dev.copy()
        for idx, af in enumerate(self.afs):
            mdev[idx] = mdev[idx]*numpy.sqrt(3)/(af*self.tau0)
        return mdev

    def __str__(self):
        msg = "n_pts: %d " % len(self.x)
        for idx, af in enumerate(self.afs):
            msg += "TDEV(%d)=%f " %(af, self.dev[idx])
        return msg

if __name__ == "__main__":
    nbs14_phase = [0.00000, 103.11111, 123.22222, 157.33333, 166.44444, 48.55555, -96.33333, -2.22222, 111.88889, 0.00000]
    nbs14_f = [892.0, 809.0, 823.0, 798.0, 671.0, 644.0, 883.0, 903.0, 677.0]

    nbs14_devs = [(91.22945, 115.8082),  # 0, ADEV(tau=1,tau=2)
                  (91.22945, 85.95287), # 1, OADEV
                  (91.22945, 74.78849),  # 2, MDEV
                  #(91.22945,98.31100),  # TOTDEV, http://www.ieee-uffc.org/frequency-control/learning-riley.asp
                  (91.22945, 93.90379), # 3, TOTDEV, http://tf.nist.gov/general/pdf/2220.pdf page 107
                  (70.80608, 116.7980),  # 4, HDEV
                  (52.67135, 86.35831),  # 5, TDEV
                  (70.80607, 85.61487), # 6, OHDEV
                  #(75.50203, 75.83606),  # 7, MTOTDEV (published table, WFM bias correction)
                  (6.4509e+01, 6.4794e+01), # MTOTDEV Stable32 v1.60, no bias-correction
                  #(43.59112, 87.56794 ), # 8, TTOTDEV (published table, WFM bias correction)
                  (3.7244e+01, 7.4818e+01), # TTOTDEV Stable32 v1.60, no bias-correction
                  (70.80607, 91.16396), # 9 HTOTDEV (published table)
                  (45.704, 81.470)] # 10 HTOTDEV Stable32 (not sure what these are!??)
                  # (100.9770, 102.6039)  # Standard Deviation (sample, not population)
    oadev_rt = oadev_realtime(afs=[1, 2], tau0=1.0)
    for x in nbs14_phase:
        oadev_rt.add_phase(x)
    print("OADEV rt ", oadev_rt.dev[0], " == OADEV ", nbs14_devs[1][0])
    print("OADEV rt ", oadev_rt.dev[1], " == OADEV ", nbs14_devs[1][1])
    assert(numpy.isclose(oadev_rt.dev[0], nbs14_devs[1][0]))
    assert(numpy.isclose(oadev_rt.dev[1], nbs14_devs[1][1]))

    ohdev_rt = ohdev_realtime(afs=[1, 2], tau0=1.0)
    for x in nbs14_phase:
        ohdev_rt.add_phase(x)
    print("OHDEV rt ", ohdev_rt.dev[0], " == OADEV ", nbs14_devs[6][0])
    print("OHDEV rt ", ohdev_rt.dev[1], " == OADEV ", nbs14_devs[6][1])
    assert(numpy.isclose(ohdev_rt.dev[0], nbs14_devs[6][0]))
    assert(numpy.isclose(ohdev_rt.dev[1], nbs14_devs[6][1]))

    tdev_rt = tdev_realtime(afs=[1, 2], tau0=1.0)
    for x in nbs14_phase:
        tdev_rt.add_phase(x)
    print("mDEV rt ", tdev_rt.mdev()[0], " == mDEV ", nbs14_devs[2][0])
    print("mDEV rt ", tdev_rt.mdev()[1], " == mDEV ", nbs14_devs[2][1])
    assert(numpy.isclose(tdev_rt.mdev()[0], nbs14_devs[2][0]))
    assert(numpy.isclose(tdev_rt.mdev()[1], nbs14_devs[2][1]))

# end of file realtime.py
