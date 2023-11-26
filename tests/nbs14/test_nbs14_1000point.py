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
# import math
# import time
import sys
import pytest
import numpy as np
import allantools as allan

sys.path.append("..")
sys.path.append("../..")  # hack to import from parent directory
# remove if you have allantools installed in your python path

# 1000 point deviations from:
# http://www.ieee-uffc.org/frequency-control/learning-riley.asp    Table III
# http://www.wriley.com/paper1ht.htm
# http://tf.nist.gov/general/pdf/2220.pdf  page 108
nbs14_1000_devs = [[2.922319e-01, 9.965736e-02, 3.897804e-02],  # 0 ADEV 1, 10, 100
                   [2.922319e-01, 9.159953e-02, 3.241343e-02],  # 1 OADEV
                   [2.922319e-01, 6.172376e-02, 2.170921e-02],  # 2 MDEV
                   # [2.922319e-01, 9.172131e-02, 3.501795e-02],  # TOTDEV, http://www.ieee-uffc.org/frequency-control/learning-riley.asp
                   # "Calculated using bias-corrected reflected method from endpoint-matched phase data"

                   # 3 TOTDEV, http://tf.nist.gov/general/pdf/2220.pdf page 108
                   [2.922319e-01, 9.134743e-02, 3.406530e-02],
                   # "Calculated using doubly reflected TOTVAR method"

                   [2.943883e-01, 1.052754e-01, 3.910860e-02],   # 4 HDEV
                   [1.687202e-01, 3.563623e-01, 1.253382e-00],   # 5 TDEV
                   [2.943883e-01, 9.581083e-02, 3.237638e-02],   # 6 OHDEV
                   [2.884664e-01, 9.296352e-02, 3.206656e-02],   # 7 standard deviation,  sample (not population)
                   [2.943883e-01, 9.614787e-02, 3.058103e-02],   # 8 HTOTDEV
                   # [2.418528e-01, 6.499161e-02, 2.287774e-02], # 9 MTOTDEV (from published table, WITH bias correction)
                   # MTOTDEV (from Stable32 v1.60 decade run, NO bias correction)
                   [2.0664e-01, 5.5529e-02, 1.9547e-02],
                   # [1.396338e-01, 3.752293e-01, 1.320847e-00],    # 10 TTOTDEV (from published table, WITH bias correction)
                   [1.1930e-01, 3.2060e-01, 1.1285e+00],    # 10 TTOTDEV (from Stable 32 v1.60 decade run, NO bias correction)
                   [1.0757e-01, 3.1789e-02, 5.0524e-03], ]  # 11 THEO1 (tau= 10,100,1000, from Stable32, NO bias correction

# computed with 'Pdev' from
# SigmaTheta v3.0-29-g053a Wed Oct 11 12:59:22 2023 +0000 
# available at https://gitlab.com/fm-ltfb/SigmaTheta.git
pdev_taus = [1,2,4,8,16,32,64,128,256]
pdev_devs = [ 2.9223187810675200e-01,   2.1445233564252639e-01,  1.5618112158618463e-01,
              1.1709745745448434e-01,   6.9029585189839343e-02,  4.9749707730398392e-02,
              3.8947417330713739e-02,   3.0862392741372108e-02,  1.2447414341332683e-02]
              
def nbs14_1000():
    """
    this generates the nbs14 1000 point frequency dataset.
    random number generator described in
    http://www.ieee-uffc.org/frequency-control/learning-riley.asp
    http://tf.nist.gov/general/pdf/2220.pdf   page 107
    http://www.wriley.com/tst_suit.dat
    """
    n = [0]*1000
    n[0] = 1234567890
    for i in range(999):
        n[i+1] = (16807*n[i]) % 2147483647
    # the first three numbers are given in the paper, so check them:
    assert(n[1] == 395529916 and n[2] == 1209410747 and n[3] == 633705974)
    n = [x/float(2147483647) for x in n]  # normalize so that n is in [0, 1]
    return n


fdata = nbs14_1000()
pdata = allan.frequency2phase(fdata, 1.0)


class TestNBS14_1000Point():
    def test_adev(self):
        # test with frequency data
        self.nbs14_tester(allan.adev, None, fdata, nbs14_1000_devs[0])
        self.nbs14_tester(allan.adev, pdata, None,
                          nbs14_1000_devs[0])  # test with phase data

    def test_oadev(self):
        self.nbs14_tester(allan.oadev, None, fdata, nbs14_1000_devs[1])
        self.nbs14_tester(allan.oadev, pdata, None, nbs14_1000_devs[1])

    def test_mdev(self):
        self.nbs14_tester(allan.mdev, None, fdata, nbs14_1000_devs[2])
        self.nbs14_tester(allan.mdev, pdata, None, nbs14_1000_devs[2])

    def test_totdev(self):
        self.nbs14_tester(allan.totdev, None, fdata, nbs14_1000_devs[3])
        self.nbs14_tester(allan.totdev, pdata, None, nbs14_1000_devs[3])

    def test_hdev(self):
        self.nbs14_tester(allan.hdev, None, fdata, nbs14_1000_devs[4])
        self.nbs14_tester(allan.hdev, pdata, None, nbs14_1000_devs[4])

    def test_tdev(self):
        self.nbs14_tester(allan.tdev, None, fdata, nbs14_1000_devs[5])
        self.nbs14_tester(allan.tdev, pdata, None, nbs14_1000_devs[5])

    def test_ohdev(self):
        self.nbs14_tester(allan.ohdev, None, fdata, nbs14_1000_devs[6])
        self.nbs14_tester(allan.ohdev, pdata, None, nbs14_1000_devs[6])

    def test_pdev(self):
        #d = np.loadtxt('pdev_nbs14.txt')
        #mytaus = d[:-1,0]
        #pdevs = d[:-1,1]
        #print(mytaus,pdevs)
        self.nbs14_tester(allan.pdev, None, fdata, pdev_devs, taus = pdev_taus)
        self.nbs14_tester(allan.pdev, pdata, None, pdev_devs, taus = pdev_taus)

    def notest_mtotdev(self):  # very slow, disable for now
        self.nbs14_tester(allan.mtotdev, None, fdata,
                          nbs14_1000_devs[9], soft=False)
        self.nbs14_tester(allan.mtotdev, pdata, None,
                          nbs14_1000_devs[9], soft=False)

    def notest_ttotdev(self):  # very slow, disable for now
        self.nbs14_tester(allan.ttotdev, None, fdata,
                          nbs14_1000_devs[10], soft=False)
        self.nbs14_tester(allan.ttotdev, pdata, None,
                          nbs14_1000_devs[10], soft=False)

    def test_theo1(self):
        self.nbs14_tester(allan.theo1, None, fdata,
                          nbs14_1000_devs[11], taus=[10, 100, 1000])
        self.nbs14_tester(allan.theo1, pdata, None,
                          nbs14_1000_devs[11], taus=[10, 100, 1000])

    def nbs14_tester(self, function, pdata, fdata, correct_devs, taus=[1, 10, 100], soft=False):
        rate = 1.0

        if pdata is not None and fdata is None:
            (taus2, devs, deverrs, ns) = function(pdata,
                                                  rate=rate, taus=taus)
        elif pdata is None and fdata is not None:
            (taus2, devs, deverrs, ns) = function(fdata,
                                                  data_type="freq",
                                                  rate=rate, taus=taus)
        else:
            raise Exception("Ambiguous input data")

        for i in range(3):
            if soft:
                print(' calculated ', devs[i])
                print('      table ', correct_devs[i])
            else:
                assert(self.check_devs(devs[i], correct_devs[i]))

    def check_devs(self, dev2, dev1, soft=False):
        rel_error = (dev2-dev1)/dev1
        tol = 1e-4
        verbose = 1

        if (abs(rel_error) < tol):
            if verbose:
                print(" OK   %0.6f \t    %0.6f \t %0.6f" %
                      (dev1, dev2, rel_error))
            return True
        else:
            print(" ERROR   %0.6f \t %0.6f \t %0.6f" % (dev1, dev2, rel_error))
            print("  bias corr %.4f" % pow(dev2/dev1, 2))
            if soft:
                return True
            return False


if __name__ == "__main__":
    t = TestNBS14_1000Point()
    t.test_pdev()
    pytest.main(["test_nbs14_1000point.py"])
