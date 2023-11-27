"""
  NBS14 test for allantools (https://github.com/aewallin/allantools)

  nbs14 datasets are from http://www.ieee-uffc.org/frequency-control/learning-riley.asp
  
  Stable32 was used to calculate the deviations we compare against.

  The small dataset and deviations are from
  http://www.ieee-uffc.org/frequency-control/learning-riley.asp
  http://www.wriley.com/paper1ht.htm
  
  see also:
  NIST Special Publication 1065
  Handbook of Frequency Stability Analysis
  http://tf.nist.gov/general/pdf/2220.pdf
  around page 107
"""

import math
import time
import sys
print(sys.version)
import pytest
import numpy

import allantools as allan

# 1000 point deviations from:
# http://www.ieee-uffc.org/frequency-control/learning-riley.asp    Table III
# http://www.wriley.com/paper1ht.htm
# http://tf.nist.gov/general/pdf/2220.pdf  page 108
nbs14_1000_devs = [ [2.922319e-01, 9.965736e-02, 3.897804e-02],  # 0 ADEV 1, 10, 100 
                    [2.922319e-01, 9.159953e-02, 3.241343e-02],  # 1 OADEV
                    [2.922319e-01, 6.172376e-02, 2.170921e-02],  # 2 MDEV 
                    #[2.922319e-01, 9.172131e-02, 3.501795e-02],  # TOTDEV, http://www.ieee-uffc.org/frequency-control/learning-riley.asp
                    # "Calculated using bias-corrected reflected method from endpoint-matched phase data"
                    
                    [2.922319e-01, 9.134743e-02, 3.406530e-02],    # 3 TOTDEV, http://tf.nist.gov/general/pdf/2220.pdf page 108
                    # "Calculated using doubly reflected TOTVAR method"
                    
                    [2.943883e-01, 1.052754e-01, 3.910860e-02],   # 4 HDEV
                    [1.687202e-01, 3.563623e-01, 1.253382e-00],   # 5 TDEV
                    [2.943883e-01, 9.581083e-02, 3.237638e-02],   # 6 OHDEV
                    [2.884664e-01, 9.296352e-02, 3.206656e-02],   # 7 standard deviation,  sample (not population)
                    [2.943883e-01, 9.614787e-02, 3.058103e-02],   # 8 HTOTDEV
                    #[2.418528e-01, 6.499161e-02, 2.287774e-02],   # 9 MTOTDEV (from published table, WITH bias correction)
                    [2.0664e-01, 5.5529e-02, 1.9547e-02], # MTOTDEV (from Stable32 v1.60 decade run, NO bias correction)
                    #[1.396338e-01, 3.752293e-01, 1.320847e-00],    # 10 TTOTDEV (from published table, WITH bias correction)
                     [1.1930e-01, 3.2060e-01, 1.1285e+00 ],  # 10 TTOTDEV (from Stable 32 v1.60 decade run, NO bias correction)  
                    [1.0757e-01, 3.1789e-02, 5.0524e-03 ], ] # 11 THEO1 (tau= 10,100,1000, from Stable32, NO bias correction
# this generates the nbs14 1000 point frequency dataset. 
# random number generator described in 
# http://www.ieee-uffc.org/frequency-control/learning-riley.asp
# http://tf.nist.gov/general/pdf/2220.pdf   page 107
# http://www.wriley.com/tst_suit.dat
def nbs14_1000():
    """
        1000-point test dataset.
        data is fractional frequency
    """
    n = [0]*1000
    n[0] = 1234567890
    for i in range(999):
        n[i+1] = (16807*n[i]) % 2147483647
    # the first three numbers are given in the paper, so check them:
    assert( n[1] ==  395529916 and n[2] == 1209410747 and  n[3] == 633705974 )
    n = [x/float(2147483647) for x in n] # normalize so that n is in [0, 1]
    return n

nbs14_f = nbs14_1000()
nbs14_phase = allan.frequency2phase(nbs14_f, 1.0)

def check_dev(name, tau, a, b):
    print(name," tau=",tau, " ", a ," == ", b)
    assert( numpy.isclose( a, b) )

def test_oadev_rt_nbs14_1k():
    oadev_rt = allan.realtime.oadev_realtime(afs=[1,10,100],tau0=1.0)
    for x in nbs14_phase:
        oadev_rt.add_phase(x)
    for n in range(3):
        check_dev('OADEV', oadev_rt.taus()[n], oadev_rt.dev[n], nbs14_1000_devs[1][n])

def test_ohdev_rt_nbs14_1k():
    dev_rt = allan.realtime.ohdev_realtime(afs=[1,10,100],tau0=1.0)
    for x in nbs14_phase:
        dev_rt.add_phase(x)
    for n in range(3):
        check_dev('OHDEV', dev_rt.taus()[n], dev_rt.dev[n], nbs14_1000_devs[6][n])

def test_tdev_rt_nbs14_1k():
    dev_rt = allan.realtime.tdev_realtime(afs=[1,10,100],tau0=1.0)
    for x in nbs14_phase:
        dev_rt.add_phase(x)
    for n in range(3):
        check_dev('TDEV', dev_rt.taus()[n], dev_rt.dev[n], nbs14_1000_devs[5][n])
            
if __name__ == "__main__":
    test_oadev_rt_nbs14_1k()
    test_ohdev_rt_nbs14_1k()
    test_tdev_rt_nbs14_1k()
