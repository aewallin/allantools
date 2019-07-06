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
import pytest

import allantools as at
import numpy

# 10-point dataset and deviations
nbs14_phase = [ 0.00000, 103.11111, 123.22222, 157.33333, 166.44444, 48.55555,-96.33333,-2.22222, 111.88889, 0.00000 ]
nbs14_f     = [892.0,809.0,823.0,798.0,671.0,644.0,883.0,903.0,677.0]

nbs14_devs= [ (91.22945,115.8082),  # 0, ADEV(tau=1,tau=2)
              (91.22945, 85.95287), # 1, OADEV 
              (91.22945,74.78849),  # 2, MDEV
              #(91.22945,98.31100),  # TOTDEV, http://www.ieee-uffc.org/frequency-control/learning-riley.asp
              (91.22945, 93.90379), # 3, TOTDEV, http://tf.nist.gov/general/pdf/2220.pdf page 107
              (70.80608,116.7980),  # 4, HDEV
              (52.67135,86.35831),  # 5, TDEV 
              (70.80607, 85.61487), # 6, OHDEV
              #(75.50203, 75.83606),  # 7, MTOTDEV (published table, WFM bias correction)
              (6.4509e+01, 6.4794e+01), # MTOTDEV Stable32 v1.60, no bias-correction
              #(43.59112, 87.56794 ), # 8, TTOTDEV (published table, WFM bias correction)
              (3.7244e+01, 7.4818e+01), # TTOTDEV Stable32 v1.60, no bias-correction
              (70.80607, 91.16396 ), # 9 HTOTDEV (published table)
              (45.704, 81.470)] # 10 HTOTDEV Stable32 (not sure what these are!??)
              # (100.9770, 102.6039)  # Standard Deviation (sample, not population)

def test_oadev_rt_nbs14():
    oadev_rt = at.realtime.oadev_realtime(afs=[1,2],tau0=1.0)
    for x in nbs14_phase:
        oadev_rt.add_phase(x)
    print "OADEV rt ",oadev_rt.dev[0] ," == OADEV ", nbs14_devs[1][0]
    print "OADEV rt ",oadev_rt.dev[1] ," == OADEV ", nbs14_devs[1][1]
    assert( numpy.isclose( oadev_rt.dev[0], nbs14_devs[1][0] ))
    assert( numpy.isclose( oadev_rt.dev[1], nbs14_devs[1][1] ))

if __name__ == "__main__":

    #t=TestNBS14_10Point()
    #t.test_htotdev()
    pass
    # doe we need this?
    #pytest.main(["oadev_nbs10.py"])

