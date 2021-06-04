"""
  NBS14 test for allantools (https://github.com/aewallin/allantools)

  nbs14 datasets are from
  http://www.ieee-uffc.org/frequency-control/learning-riley.asp

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
# import time
# import sys
import pytest

import allantools as allan

# 10-point dataset and deviations
nbs14_phase = [0.00000, 103.11111, 123.22222, 157.33333,
               166.44444, 48.55555, -96.33333, -2.22222, 111.88889, 0.00000]
nbs14_f = [892.0, 809.0, 823.0, 798.0, 671.0, 644.0, 883.0, 903.0, 677.0]

nbs14_devs = [(91.22945, 115.8082),  # 0, ADEV(tau=1,tau=2)
              (91.22945, 85.95287),  # 1, OADEV
              (91.22945, 74.78849),  # 2, MDEV
              # (91.22945,98.31100),  # TOTDEV, http://www.ieee-uffc.org/frequency-control/learning-riley.asp
              # 3, TOTDEV, http://tf.nist.gov/general/pdf/2220.pdf page 107
              (91.22945, 93.90379),
              (70.80608, 116.7980),  # 4, HDEV
              (52.67135, 86.35831),  # 5, TDEV
              (70.80607, 85.61487),  # 6, OHDEV
              # (75.50203, 75.83606),  # 7, MTOTDEV (published table, WFM bias correction)
              # MTOTDEV Stable32 v1.60, no bias-correction
              (6.4509e+01, 6.4794e+01),
              # (43.59112, 87.56794 ), # 8, TTOTDEV (published table, WFM bias correction)
              # TTOTDEV Stable32 v1.60, no bias-correction
              (3.7244e+01, 7.4818e+01),
              (70.80607, 91.16396),  # 9 HTOTDEV (published table)
              (45.704, 81.470)]  # 10 HTOTDEV Stable32 (not sure what these are!??)
# (100.9770, 102.6039)  # Standard Deviation (sample, not population)


def check_devs(dev2, dev1, soft=False):
    rel_error = (dev2-dev1)/dev1
    tol = 1e-4
    verbose = 1

    if (abs(rel_error) < tol):
        if verbose:
            print(" OK   %0.6f \t    %0.6f \t %0.6f" % (dev1, dev2, rel_error))
        return True
    else:
        print(" ERROR   %0.6f \t %0.6f \t %0.6f" % (dev1, dev2, rel_error))
        print(" bias= %.4f" % pow(dev2/dev1, 2))
        if soft:
            return True
        return False


class TestNBS14_10Point():
    def test_adev(self):
        self.generic_test(allan.adev, nbs14_devs[0])

    def test_oadev(self):
        self.generic_test(allan.oadev, nbs14_devs[1])

    def test_mdev(self):
        self.generic_test(allan.mdev, nbs14_devs[2])

    def test_totdev(self):
        self.generic_test(allan.totdev, nbs14_devs[3])

    def test_hdev(self):
        self.generic_test(allan.hdev, nbs14_devs[4])

    def test_tdev(self):
        self.generic_test(allan.tdev, nbs14_devs[5])

    def test_ohdev(self):
        self.generic_test(allan.ohdev, nbs14_devs[6])

    def test_mtotdev(self):
        self.generic_test(allan.mtotdev, nbs14_devs[7])

    def test_ttotdev(self):
        self.generic_test(allan.ttotdev, nbs14_devs[8])

    def test_htotdev(self):
        # NOTE:
        # for tau=1, ohdev() is used instead of htotdev()
        # for tau=2 no bias correction is applied by allantools
        #           instead we correct the reference value
        #           by multiplication with sqrt(0.995)
        htotdev1, htotdev2 = nbs14_devs[9]
        self.generic_test(allan.htotdev, (htotdev1, math.sqrt(0.995)*htotdev2))

    # commented out since the meaning of the reference values is unknown
    # @pytest.mark.fails
    # def test_htotdev2(self):
    #    # compare against Stable32 values - fixme!
    #    self.generic_test( allan.htotdev, nbs14_devs[10] )

    def generic_test(self, function, ref_devs):
        taus = [1, 2]
        (taus1, adevs1, aerrs1, ns1) = function(
            nbs14_phase, rate=1.0, taus=taus)
        (taus2, adevs2, aerrs2, ns2) = function(nbs14_f, rate=1.0,
                                                data_type="freq", taus=taus)
        assert(check_devs(adevs1[0], ref_devs[0]))  # tau=1 from phase data
        assert(check_devs(adevs1[1], ref_devs[1]))  # tau=2 from phase data
        assert(check_devs(adevs2[0], ref_devs[0]))  # tau=1 from frequency data
        assert(check_devs(adevs2[1], ref_devs[1]))  # tau=2 from frequency data


if __name__ == "__main__":
    # t=TestNBS14_10Point()
    # t.test_htotdev()
    pytest.main(["test_nbs14_10point.py"])
