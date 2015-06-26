import allantools as allan
import numpy
import math
import matplotlib.pyplot as plt

import sys
sys.path.append("..")
sys.path.append("../..") # hack to import from parent directory
# remove if you have allantools installed in your python path

def read_datafile(filename):
    p=[]
    with open(filename) as f:
        for line in f:
            
            if line.startswith("#"): # skip comments
                pass
            else:
                line = line.split()
                p.append( float(line[1]) )
    return p

def to_fractional(flist,f0):
    out=[]
    for f in flist:
        out.append( f/float(f0) - 1.0 )
    return out

fname = "2015-06-26_ia_ocxo.txt"
f10MHz = read_datafile(fname)

"""
with open("ocxo.txt",'w') as f:
    for fn in f10MHz:
        f.write( '%.15f\n' % fn )
"""

f = to_fractional(f10MHz, 10e6 ) # convert to fractional frequency
my_taus = numpy.logspace(1,5,40) # log-spaced tau values from 10s and upwards
rate = 1/float(10) # data collected with 10s gate time

(oadev_taus,oadev_devs,oadev_errs,ns)  = allan.oadev(f, rate, my_taus)
(mdev_taus,mdev_devs,mdev_errs,ns)  = allan.mdev(f, rate, my_taus)
(hdev_taus,hdev_devs,hdev_errs,ns)  = allan.hdev(f, rate, my_taus)
(ohdev_taus,ohdev_devs,ohdev_errs,ns)  = allan.ohdev(f, rate, my_taus)

# now compare to results produced by Stable32
#import testutils
#(taus32,devs32,ns32) = testutils.read_stable32('rb_10s_stable32_oadev.txt', rate)


plt.subplot(111, xscale="log", yscale="log") 

plt.errorbar(oadev_taus, oadev_devs, yerr=oadev_errs, label='OADEV') 
plt.errorbar(mdev_taus, mdev_devs, yerr=mdev_errs, label='MDEV') 
plt.errorbar(hdev_taus, hdev_devs, yerr=hdev_errs, label='HDEV') 
plt.errorbar(ohdev_taus, ohdev_devs, yerr=ohdev_errs, label='OHDEV') 

#plt.errorbar(taus32, devs32, yerr=[d/math.sqrt(n) for (d,n) in zip(devs32,ns32)]) 

plt.xlabel('Taus (s)')
plt.ylabel('ADEV')


plt.legend()
plt.grid()
plt.show()
