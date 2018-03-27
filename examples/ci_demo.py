"""
 Allan deviation tools
 Anders Wallin (anders.e.e.wallin "at" gmail.com)
 2016-11-18

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
import matplotlib.pyplot as plt  # only for plotting, not required for calculations

import allantools as at

# this demonstrates how to calculate confidence intervals for ADEV
# using the algorithms from Greenhall2004
data_file = '../tests/phasedat/PHASE.DAT'

def read_datafile(filename):
    p = []
    with open(filename) as f:
        for line in f:
            if not line.startswith("#"):  # skip comments
                p.append(float(line))
    return p

# read input data from file
phase = read_datafile(data_file)
# normal ADEV computation, giving naive 1/sqrt(N) errors
(taus,devs,errs,ns) = at.adev(phase, taus='octave')

# Confidence-intervals for each (tau,adev) pair separately.
cis=[]
for (t,dev) in zip(taus,devs):
    # Greenhalls EDF (Equivalent Degrees of Freedom)
    # alpha     +2,...,-4   noise type, either estimated or known
    # d         1 first-difference variance, 2 allan variance, 3 hadamard variance
    #           we require: alpha+2*d >1     (is this ever false?)
    # m         tau/tau0 averaging factor
    # N         number of phase observations
    edf = at.edf_greenhall( alpha=0, d=2, m=t, N=len(phase), overlapping = False, modified=False )
    # with the known EDF we get CIs 
    (lo,hi) = at.confidence_interval( dev=dev,  edf=edf )
    cis.append( (lo,hi) )

# now we are ready to print and plot the results
print "Tau\tmin Dev\t\tDev\t\tMax Dev"
for (tau,dev,ci) in zip(taus,devs,cis):
    print "%d\t%f\t%f\t%f" % (tau, ci[0], dev, ci[1] )
""" output is
Tau	min Dev		Dev		Max Dev
1	0.285114	0.292232	0.299910
2	0.197831	0.205102	0.213237
4	0.141970	0.149427	0.158198
8	0.102541	0.110135	0.119711
16	0.056510	0.062381	0.070569
32	0.049153	0.056233	0.067632
64	0.027109	0.032550	0.043536
128	0.026481	0.033855	0.055737
256	0.007838	0.010799	0.031075
"""
plt.figure(figsize=(12,8))
#fig, ax = plt.subplots(111)
plt.gca().set_yscale('log')
plt.gca().set_xscale('log')
err_lo = [ d-ci[0] for (d,ci) in zip(devs,cis)]
err_hi = [ ci[1]-d for (d,ci) in zip(devs,cis)]

plt.errorbar(taus, devs, yerr=[ err_lo, err_hi ] ,fmt='o')
plt.grid()
plt.xlabel('Tau (s)')
plt.ylabel('ADEV')
plt.title('AllanTools 2016.11 - now with Confidence Intervals!')
# just to check plot the intervals as dots also
plt.plot(taus, [ci[0] for ci in cis],'r.',label="Lower CI")
plt.plot(taus, [ci[1] for ci in cis],'g.',label="Upper CI")
plt.legend(loc='best')
plt.show()
