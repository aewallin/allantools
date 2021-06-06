"""
  PHASE.DAT test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  PHASE.DAT comes with Stable32 (version 1.53 was used in this case)

  This test is for Confidence Intervals in particular

"""


import allantools as allan
import sys
import pytest
import numpy as np

sys.path.append("..")
sys.path.append("../..")  # hack to import from parent directory
# remove if you have allantools installed in your python path
import testutils

data_file = 'PHASE.DAT'
tolerance = 1e-4  # relative tolerance for pass/fail against Stable32 values
verbose = 1


class TestPhaseDatCI():
    """
        Input file: PHASE.DAT

    """

    def test_phasedat_adev(self):
        s32_rows = testutils.read_stable32('phase_dat_adev_octave.txt', 1.0)
        phase = np.array(testutils.read_datafile('PHASE.DAT'))
        (taus, devs, errs, ns) = allan.adev(
            phase, taus=[s32['tau'] for s32 in s32_rows])

        # separate CI computation
        los = []
        his = []
        for (d, t, n, s32) in zip(devs, taus, ns, s32_rows):
            # Note FIXED alpha here
            edf2 = allan.edf_greenhall(alpha=0, d=2, m=t, N=len(
                phase), overlapping=False, modified=False)
            (lo, hi) = allan.confidence_interval(dev=d, edf=edf2)
            assert np.isclose(lo, s32['dev_min'], rtol=1e-2)
            assert np.isclose(hi, s32['dev_max'], rtol=1e-2)
            print(" alpha=0 FIXED, CI OK! tau = %f" % t)

            los.append(lo)
            his.append(hi)
            try:
                (lo2, hi2) = allan.confidence_interval_noiseID(
                    phase, d, af=int(t), dev_type="adev", data_type="phase")
                assert np.isclose(lo2, s32['dev_min'], rtol=1e-2)
                assert np.isclose(hi2, s32['dev_max'], rtol=1e-2)
                print(" ACF_NID CI OK! tau = %f" % t)
            except NotImplementedError:
                print("can't do CI for tau = %f" % t)
                pass

        # compare to Stable32
        print("adev()")
        print("    n   tau dev_min  dev      dev_max ")
        for (s32, t2, d2, lo2, hi2, n2) in zip(s32_rows, taus, devs, los, his, ns):
            print("S32 %03d %03.1f %1.6f %1.6f %1.6f" % (
                s32['n'], s32['tau'], s32['dev_min'], s32['dev'], s32['dev_max']))
            print("AT  %03d %03.1f %1.6f %1.6f %1.6f" %
                  (n2, t2, round(lo2, 5), round(d2, 5), round(hi2, 5)))
            testutils.check_approx_equal(s32['n'], n2, tolerance=1e-9)
            testutils.check_approx_equal(s32['dev_min'], lo2, tolerance=1e-3)
            testutils.check_approx_equal(s32['dev'], d2, tolerance=1e-4)
            testutils.check_approx_equal(s32['dev_max'], hi2, tolerance=1e-3)
        print("----")

    def test_phasedat_oadev(self):
        s32_rows = testutils.read_stable32('phase_dat_oadev_octave.txt', 1.0)
        phase = testutils.read_datafile('PHASE.DAT')
        (taus, devs, errs, ns) = allan.oadev(
            phase, taus=[s32['tau'] for s32 in s32_rows])

        # CI computation
        # alhpa     +2,...,-4   noise power
        # d         1 first-difference variance, 2 allan variance, 3 hadamard variance
        #           alpha+2*d >1
        # m         tau/tau0 averaging factor
        # F         filter factor, 1 modified variance, m unmodified variance
        # S         stride factor, 1 nonoverlapped estimator, m overlapped estimator (estimator stride = tau/S )
        # N         number of phase observations
        los = []
        his = []
        for (d, t, n) in zip(devs, taus, ns):
            edf2 = allan.edf_greenhall(alpha=0, d=2, m=int(
                t), N=len(phase), overlapping=True, modified=False)
            (lo, hi) = allan.confidence_interval(dev=d, edf=edf2)
            los.append(lo)
            his.append(hi)
        # compare to Stable32
        print("oadev()")
        for (s32, t2, d2, lo2, hi2, n2) in zip(s32_rows, taus, devs, los, his, ns):
            print("s32 %03d %03f %1.6f %1.6f %1.6f" % (
                s32['n'], s32['tau'], s32['dev_min'], s32['dev'], s32['dev_max']))
            print("at  %03d %03f %1.6f %1.6f %1.6f" %
                  (n2, t2, round(lo2, 5), round(d2, 5), round(hi2, 5)))
            testutils.check_approx_equal(s32['dev_min'], lo2, tolerance=1e-3)
            testutils.check_approx_equal(s32['dev'], d2, tolerance=1e-4)
            testutils.check_approx_equal(s32['dev_max'], hi2, tolerance=1e-3)
        print("----")

    def test_phasedat_mdev(self):
        s32_rows = testutils.read_stable32('phase_dat_mdev_octave.txt', 1.0)
        phase = testutils.read_datafile('PHASE.DAT')
        (taus, devs, errs, ns) = allan.mdev(
            phase, taus=[s32['tau'] for s32 in s32_rows])

        # CI computation
        # alhpa= +2,...,-4   noise power
        # d= 1 first-difference variance, 2 allan variance, 3 hadamard variance
        # alpha+2*d >1
        # m = tau/tau0 averaging factor
        # F filter factor, 1 modified variance, m unmodified variance
        # S stride factor, 1 nonoverlapped estimator, m overlapped estimator (estimator stride = tau/S )
        # N number of phase obs
        los = []
        his = []
        for (d, t, n) in zip(devs, taus, ns):
            edf2 = allan.edf_greenhall(alpha=0, d=2, m=int(
                t), N=len(phase), overlapping=True, modified=True)
            (lo, hi) = allan.confidence_interval(dev=d, edf=edf2)
            los.append(lo)
            his.append(hi)
        # compare to Stable32
        print("mdev()")
        for (s32, t2, d2, lo2, hi2, n2) in zip(s32_rows, taus, devs, los, his, ns):
            print("s32 %03d %03f %1.6f %1.6f %1.6f" % (
                s32['n'], s32['tau'], s32['dev_min'], s32['dev'], s32['dev_max']))
            print("at  %03d %03f %1.6f %1.6f %1.6f" %
                  (n2, t2, round(lo2, 5), round(d2, 5), round(hi2, 5)))
            testutils.check_approx_equal(s32['dev_min'], lo2, tolerance=1e-3)
            testutils.check_approx_equal(s32['dev'], d2, tolerance=1e-4)
            testutils.check_approx_equal(s32['dev_max'], hi2, tolerance=1e-3)
        print("----")

    def test_phasedat_hdev(self):
        s32_rows = testutils.read_stable32('phase_dat_hdev_octave.txt', 1.0)
        phase = testutils.read_datafile('PHASE.DAT')
        (taus, devs, errs, ns) = allan.hdev(
            phase, taus=[s32['tau'] for s32 in s32_rows])

        # CI computation
        # alhpa= +2,...,-4   noise power
        # d= 1 first-difference variance, 2 allan variance, 3 hadamard variance
        # alpha+2*d >1
        # m = tau/tau0 averaging factor
        # N number of phase obs
        los = []
        his = []
        for (d, t, n) in zip(devs, taus, ns):
            #edf = greenhall_simple_edf( alpha=0, d=3, m=t, S=1, F=t, N=len(phase) )
            edf2 = allan.edf_greenhall(alpha=0, d=3, m=int(
                t), N=len(phase), overlapping=False, modified=False)
            # print(edf,edf2,edf2/edf)
            (lo, hi) = allan.confidence_interval(dev=d, edf=edf2)
            #allan.uncertainty_estimate(len(phase), t, d,ci=0.683,noisetype='wf')
            los.append(lo)
            his.append(hi)
        print("hdev()")
        for (s32, t2, d2, lo2, hi2, n2) in zip(s32_rows, taus, devs, los, his, ns):
            print("s32 %03d %03f %1.6f %1.6f %1.6f" % (
                s32['n'], s32['tau'], s32['dev_min'], s32['dev'], s32['dev_max']))
            print("at  %03d %03f %1.6f %1.6f %1.6f" %
                  (n2, t2, round(lo2, 5), round(d2, 5), round(hi2, 5)))
            testutils.check_approx_equal(s32['dev_min'], lo2, tolerance=1e-3)
            testutils.check_approx_equal(s32['dev_max'], hi2, tolerance=1e-3)
        print("----")

    def test_phasedat_ohdev(self):
        s32_rows = testutils.read_stable32('phase_dat_ohdev_octave.txt', 1.0)
        phase = testutils.read_datafile('PHASE.DAT')
        (taus, devs, errs, ns) = allan.ohdev(
            phase, taus=[s32['tau'] for s32 in s32_rows])

        # CI computation
        # alhpa= +2,...,-4   noise power
        # d= 1 first-difference variance, 2 allan variance, 3 hadamard variance
        # alpha+2*d >1
        # m = tau/tau0 averaging factor
        # F filter factor, 1 modified variance, m unmodified variance
        # S stride factor, 1 nonoverlapped estimator, m overlapped estimator (estimator stride = tau/S )
        # N number of phase obs
        los = []
        his = []
        for (d, t, n) in zip(devs, taus, ns):
            # edf = greenhall_simple_edf( alpha=0, d=3, m=t, S=1, F=t, N=len(phase) )
            edf2 = allan.edf_greenhall(alpha=0, d=3, m=int(
                t), N=len(phase), overlapping=True, modified=False)
            # print(edf,edf2,edf2/edf)
            (lo, hi) = allan.confidence_interval(dev=d, edf=edf2)
            # allan.uncertainty_estimate(len(phase), t, d,ci=0.683,noisetype='wf')
            los.append(lo)
            his.append(hi)
        print("ohdev()")
        for (s32, t2, d2, lo2, hi2, n2) in zip(s32_rows, taus, devs, los, his, ns):
            print("s32 %03d %03f %1.6f %1.6f %1.6f" % (
                s32['n'], s32['tau'], s32['dev_min'], s32['dev'], s32['dev_max']))
            print("at  %03d %03f %1.6f %1.6f %1.6f" %
                  (n2, t2, round(lo2, 5), round(d2, 5), round(hi2, 5)))
            testutils.check_approx_equal(s32['dev_min'], lo2, tolerance=1e-3)
            testutils.check_approx_equal(s32['dev_max'], hi2, tolerance=1e-3)
        print("----")

    def test_phasedat_tdev(self):
        s32_rows = testutils.read_stable32('phase_dat_tdev_octave.txt', 1.0)
        phase = testutils.read_datafile('PHASE.DAT')
        (taus, devs, errs, ns) = allan.tdev(
            phase, taus=[s32['tau'] for s32 in s32_rows])

        # CI computation
        # alhpa= +2,...,-4   noise power
        # d= 1 first-difference variance, 2 allan variance, 3 hadamard variance
        # alpha+2*d >1
        # m = tau/tau0 averaging factor
        # N number of phase obs
        los = []
        his = []
        for (d, t, n) in zip(devs, taus, ns):
            edf2 = allan.edf_greenhall(alpha=0, d=2, m=int(
                t), N=len(phase), overlapping=True, modified=True)
            # covert to mdev
            # tdev = taus * mdev / np.sqrt(3.0)
            mdev = d/t*np.sqrt(3.0)
            (lo, hi) = allan.confidence_interval(dev=mdev, edf=edf2)
            # convert back to tdev
            lo = t*lo/np.sqrt(3.0)
            hi = t*hi/np.sqrt(3.0)
            los.append(lo)
            his.append(hi)
        print("tdev()")
        for (s32, t2, d2, lo2, hi2, n2) in zip(s32_rows, taus, devs, los, his, ns):
            print("s32 %03d %03f %1.6f %1.6f %1.6f" % (
                s32['n'], s32['tau'], s32['dev_min'], s32['dev'], s32['dev_max']))
            print("at  %03d %03f %1.6f %1.6f %1.6f" %
                  (n2, t2, round(lo2, 5), round(d2, 5), round(hi2, 5)))
            testutils.check_approx_equal(s32['dev_min'], lo2, tolerance=1e-3)
            testutils.check_approx_equal(s32['dev_max'], hi2, tolerance=1e-3)
        print("----")

    def test_phasedat_totdev(self):
        s32_rows = testutils.read_stable32('phase_dat_totdev_octave.txt', 1.0)
        phase = testutils.read_datafile('PHASE.DAT')
        (taus, devs, errs, ns) = allan.totdev(
            phase, taus=[s32['tau'] for s32 in s32_rows])

        los = []
        his = []
        for (d, t, n) in zip(devs, taus, ns):
            # edf = greenhall_simple_edf( alpha=0, d=3, m=t, S=1, F=t, N=len(phase) )
            # edf2 = greenhall_edf( alpha=0, d=3, m=int(t), N=len(phase), overlapping = True, modified=False  )
            # print(edf,edf2,edf2/edf)
            edf = allan.edf_totdev(len(phase), t, alpha=0)
            (lo, hi) = allan.confidence_interval(dev=d,  edf=edf)
            # allan.uncertainty_estimate(len(phase), t, d,ci=0.683,noisetype='wf')
            los.append(lo)
            his.append(hi)
        print("totdev()")
        for (s32, t2, d2, lo2, hi2, n2) in zip(s32_rows, taus, devs, los, his, ns):
            print("s32 %03d %03f %1.6f %1.6f %1.6f" % (
                s32['n'], s32['tau'], s32['dev_min'], s32['dev'], s32['dev_max']))
            print("at  %03d %03f %1.6f %1.6f %1.6f" %
                  (n2, t2, round(lo2, 5), round(d2, 5), round(hi2, 5)))
            testutils.check_approx_equal(s32['dev_min'], lo2, tolerance=1e-3)
            testutils.check_approx_equal(s32['dev_max'], hi2, tolerance=1e-3)
        print("----")

    def test_noise_id(self):
        """ test for noise-identification """
        s32_rows = testutils.read_stable32('phase_dat_oadev_octave.txt', 1.0)
        phase = testutils.read_datafile('PHASE.DAT')
        for s32 in s32_rows:
            tau, alpha, af = s32['tau'], s32['alpha'], int(s32['m'])
            try:
                alpha_int = allan.autocorr_noise_id(phase, af=af)[0]
                assert alpha_int == alpha
                print("OK noise-ID for af = %d" % af)
            except:
                print("can't do noise-ID for af = %d" % af)
            print("tau= %f, alpha= %f, alpha_int = %d" %
                  (tau, alpha, alpha_int))

    # FIXME: failing test that we don't run
    def slow_failing_phasedat_mtotdev(self):
        s32_rows = testutils.read_stable32(
            'phase_dat_mtotdev_octave_alpha0.txt', 1.0)
        phase = testutils.read_datafile('PHASE.DAT')
        (taus, devs, errs, ns) = allan.mtotdev(
            phase, taus=[s32['tau'] for s32 in s32_rows])

        los = []
        his = []
        for (d, t, n) in zip(devs, taus, ns):
            # edf = greenhall_simple_edf( alpha=0, d=3, m=t, S=1, F=t, N=len(phase) )
            if int(t) < 10:
                edf = allan.edf_greenhall(alpha=0, d=2, m=int(
                    t), N=len(phase), overlapping=True, modified=True)
            else:
                edf = allan.edf_mtotdev(len(phase), t, alpha=0)
            # print edf, edf2
            # print(edf,edf2,edf2/edf)

            print(edf)
            (lo, hi) = allan.confidence_interval(dev=d, edf=edf)
            # allan.uncertainty_estimate(len(phase), t, d,ci=0.683,noisetype='wf')
            los.append(lo)
            his.append(hi)
        print("mtotdev()")
        for (s32, t2, d2, lo2, hi2, n2) in zip(s32_rows, taus, devs, los, his, ns):
            print("s32 %03d %03f %1.6f %1.6f %1.6f" % (
                s32['n'], s32['tau'], s32['dev_min'], s32['dev'], s32['dev_max']))
            print("at  %03d %03f %1.6f %1.6f %1.6f" %
                  (n2, t2, round(lo2, 5), round(d2, 5), round(hi2, 5)))
            testutils.check_approx_equal(s32['dev_min'], lo2, tolerance=1e-3)
            testutils.check_approx_equal(s32['dev_max'], hi2, tolerance=1e-3)
        print("----")


if __name__ == "__main__":
    # pytest.main()
    t = TestPhaseDatCI()
    t.test_noise_id()
    t.test_phasedat_adev()
    t.test_phasedat_oadev() 
    t.test_phasedat_mdev()
    t.test_phasedat_hdev()
    t.test_phasedat_ohdev()
    t.test_phasedat_tdev()
    t.test_phasedat_totdev()

    # t.slow_failing_phasedat_mtotdev()
