"""
  PHASE.DAT test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  PHASE.DAT comes with Stable32 (version 1.53 was used in this case)

"""

import math
import sys
import pytest
import numpy as np
import scipy.stats # for ch2.ppf()
import scipy.special # for binom()

sys.path.append("..")
sys.path.append("../..") # hack to import from parent directory
# remove if you have allantools installed in your python path

import allantools as allan
import testutils

import os

data_file = 'PHASE.DAT'
tolerance = 1e-4
verbose = 1

def change_to_test_dir():
    # hack to run script from its own directory
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)

def read_datafile(filename):
    p = []
    if filename[-2:]=='gz':
        with gzip.open(filename,mode='rt') as f:
            for line in f:
                if not line.startswith("#"):  # skip comments
                    p.append(float(line))
    else:
        with open(filename) as f:
            for line in f:
                if not line.startswith("#"):  # skip comments
                    p.append(float(line))
    return p

# read a result-file, produced by copy/paste from Stable32
# note: header-lines need to be manually commented-out with "#"
def read_resultfile(filename):
    rows = []
    row = []
    with open(filename) as f:
        for line in f:
            if not line.startswith("#"): # skip comments
                row = []
                l2 = line.split(" ")
                l2 = [_f for _f in l2 if _f]
                for n in range(len(l2)):
                    row.append(float(l2[n]))
                rows.append(row)
    return rows


def read_stable32(resultfile, datarate):
    devresults = read_resultfile(resultfile)
    print("Read ", len(devresults), " rows from ", resultfile)
    taus = []
    devs = []
    devs_lo=[]
    devs_hi=[]
    ns = []
    
    # parse textfile produced by Stable32
    for row in devresults:
        if len(row) == 7:  # typical ADEV result file has 7 columns of data
            tau_n = row[0]  # tau in number of datapoints
            tau_s = row[1]  # tau in seconds
            taus.append(tau_n * (1 / float(datarate)))
            n = row[2]  # n averages
            lo = row[4]
            dev = row[5]  # deviation
            hi = row[6]
            devs.append(dev)
            devs_lo.append(lo)
            devs_hi.append(hi)
            ns.append(n)
    return (taus, devs, devs_lo, devs_hi, ns)

def approx_equal(a1,a2,tolerance=1e-5, assertions=True):
    rel_error = (a2 - a1) / a1
    if abs(rel_error) > tolerance:
        print("ERROR ", a1,a2,rel_error)
        if assertions:
            assert( abs(rel_error) < tolerance )
    
class TestPhaseDatCI():
    def test_phasedat_adev(self):
        (s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns) = read_stable32( 'phase_dat_adev_octave.txt' , 1.0 )
        phase = read_datafile('PHASE.DAT')
        (taus,devs,errs,ns) = allan.adev(phase, taus=s32_taus)
        
        # CI computation
        los=[]
        his=[]
        for (d,t, n) in zip(devs, taus, ns):
            #edf = greenhall_simple_edf( alpha=0, d=2, m=t, S=1, F=t, N=len(phase) )
            edf2 = allan.edf_greenhall( alpha=0, d=2, m=t, N=len(phase), overlapping = False, modified=False )
            #print(edf,edf2,edf2/edf)
            (lo,hi) = allan.confidence_intervals( dev=d, ci=0.68268949213708585, edf=edf2 )  # 0.68268949213708585
            #allan.uncertainty_estimate(len(phase), t, d,ci=0.683,noisetype='wf')
            los.append(lo)
            his.append(hi)
        
        #print greenhall_simple_edf( alpha=0, d=2, m=1, S=1, F=1, N=len(phase) )
        #print confidence_intervals( dev
        #   tau N       edf         chi2_l
        #   1   999     782.030     821.358 , 742.689
        #  2.851498,  2.922319  2.998720
        
        for (t1,d1,lo1,hi1, n1, t2, d2, lo2, hi2, n2) in zip(s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns, taus, devs, los, his, ns):
            print("s32 %03d %03f %1.6f %1.6f %1.6f" % (n1, t1, lo1, d1, hi1))
            print("at  %03d %03f %1.6f %1.6f %1.6f" % (n2, t2, round(lo2,5), round(d2,5), round(hi2,5) ))
            approx_equal(lo1, lo2,tolerance=1e-3)
            approx_equal(hi1, hi2,tolerance=1e-3)
            print("----")
    
    
    def test_phasedat_oadev(self):
        (s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns) = read_stable32( 'phase_dat_oadev_octave.txt' , 1.0 )
        phase = read_datafile('PHASE.DAT')
        (taus,devs,errs,ns) = allan.oadev(phase, taus=s32_taus)
        
        # CI computation
        # alhpa= +2,...,-4   noise power
        # d= 1 first-difference variance, 2 allan variance, 3 hadamard variance
        # alpha+2*d >1
        # m = tau/tau0 averaging factor
        # F filter factor, 1 modified variance, m unmodified variance
        # S stride factor, 1 nonoverlapped estimator, m overlapped estimator (estimator stride = tau/S )
        # N number of phase obs
        los=[]
        his=[]
        for (d,t, n) in zip(devs, taus, ns):
            #edf = greenhall_simple_edf( alpha=0, d=2, m=t, S=1, F=t, N=len(phase) )
            edf2 = allan.edf_greenhall( alpha=0, d=2, m=int(t), N=len(phase), overlapping = True, modified=False  )
            #print(edf,edf2,edf2/edf)
            (lo,hi) = allan.confidence_intervals( dev=d, ci=0.68268949213708585, edf=edf2 )  # 0.68268949213708585
            #allan.uncertainty_estimate(len(phase), t, d,ci=0.683,noisetype='wf')
            los.append(lo)
            his.append(hi)
        
        #print greenhall_simple_edf( alpha=0, d=2, m=1, S=1, F=1, N=len(phase) )
        #print confidence_intervals( dev
        #   tau N       edf         chi2    chi2    dev_lo      dev         dev_hi
        #   1   999     782.030     821.358 742.689 2.8515e-01  2.9223e-01  2.9987e-01
        #   2   997     540.681     573.374 507.975 1.9520e-01  2.0102e-01  2.0738e-01
        
        for (t1,d1,lo1,hi1, n1, t2, d2, lo2, hi2, n2) in zip(s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns, taus, devs, los, his, ns):
            print("s32 %03d %03f %1.6f %1.6f %1.6f" % (n1, t1, lo1, d1, hi1))
            print("at  %03d %03f %1.6f %1.6f %1.6f" % (n2, t2, round(lo2,5), round(d2,5), round(hi2,5) ))
            approx_equal(lo1, lo2,tolerance=1e-3)
            approx_equal(hi1, hi2,tolerance=1e-3)
            print("----")
            
    
    def test_phasedat_mdev(self):
        (s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns) = read_stable32( 'phase_dat_mdev_octave.txt' , 1.0 )
        phase = read_datafile('PHASE.DAT')
        (taus,devs,errs,ns) = allan.mdev(phase, taus=s32_taus)
        
        # CI computation
        # alhpa= +2,...,-4   noise power
        # d= 1 first-difference variance, 2 allan variance, 3 hadamard variance
        # alpha+2*d >1
        # m = tau/tau0 averaging factor
        # F filter factor, 1 modified variance, m unmodified variance
        # S stride factor, 1 nonoverlapped estimator, m overlapped estimator (estimator stride = tau/S )
        # N number of phase obs
        los=[]
        his=[]
        for (d,t, n) in zip(devs, taus, ns):
            #edf = greenhall_simple_edf( alpha=0, d=2, m=t, S=1, F=t, N=len(phase) )
            edf2 = allan.edf_greenhall( alpha=0, d=2, m=int(t), N=len(phase), overlapping = True, modified=True  )
            #print(edf,edf2,edf2/edf)
            (lo,hi) = allan.confidence_intervals( dev=d, ci=0.68268949213708585, edf=edf2 )  # 0.68268949213708585
            #allan.uncertainty_estimate(len(phase), t, d,ci=0.683,noisetype='wf')
            los.append(lo)
            his.append(hi)
        
        #print greenhall_simple_edf( alpha=0, d=2, m=1, S=1, F=1, N=len(phase) )
        #print confidence_intervals( dev
        #   tau N       edf         chi2    chi2    dev_lo      dev         dev_hi
        #   1   999     782.030     821.358 742.689 2.8515e-01  2.9223e-01  2.9987e-01
        #   2   997     540.681     573.374 507.975 1.9520e-01  2.0102e-01  2.0738e-01
        
        for (t1,d1,lo1,hi1, n1, t2, d2, lo2, hi2, n2) in zip(s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns, taus, devs, los, his, ns):
            print("s32 %03d %03f %1.6f %1.6f %1.6f" % (n1, t1, lo1, d1, hi1))
            print("at  %03d %03f %1.6f %1.6f %1.6f" % (n2, t2, round(lo2,5), round(d2,5), round(hi2,5) ))
            approx_equal(lo1, lo2,tolerance=1e-3)
            approx_equal(hi1, hi2,tolerance=1e-3)
            print("----")

    def test_phasedat_hdev(self):
        (s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns) = read_stable32( 'phase_dat_hdev_octave.txt' , 1.0 )
        phase = read_datafile('PHASE.DAT')
        (taus,devs,errs,ns) = allan.hdev(phase, taus=s32_taus)
        
        # CI computation
        # alhpa= +2,...,-4   noise power
        # d= 1 first-difference variance, 2 allan variance, 3 hadamard variance
        # alpha+2*d >1
        # m = tau/tau0 averaging factor
        # N number of phase obs
        los=[]
        his=[]
        for (d,t, n) in zip(devs, taus, ns):
            #edf = greenhall_simple_edf( alpha=0, d=3, m=t, S=1, F=t, N=len(phase) )
            edf2 = allan.edf_greenhall( alpha=0, d=3, m=int(t), N=len(phase), overlapping = False, modified=False  )
            #print(edf,edf2,edf2/edf)
            (lo,hi) = allan.confidence_intervals( dev=d, ci=0.68268949213708585, edf=edf2 )  # 0.68268949213708585
            #allan.uncertainty_estimate(len(phase), t, d,ci=0.683,noisetype='wf')
            los.append(lo)
            his.append(hi)
        
        for (t1,d1,lo1,hi1, n1, t2, d2, lo2, hi2, n2) in zip(s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns, taus, devs, los, his, ns):
            print("s32 %03d %03f %1.6f %1.6f %1.6f" % (n1, t1, lo1, d1, hi1))
            print("at  %03d %03f %1.6f %1.6f %1.6f" % (n2, t2, round(lo2,5), round(d2,5), round(hi2,5) ))
            approx_equal(lo1, lo2,tolerance=1e-3)
            approx_equal(hi1, hi2,tolerance=1e-3)
            print("----")
            
    def test_phasedat_ohdev(self):
        (s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns) = read_stable32( 'phase_dat_ohdev_octave.txt' , 1.0 )
        phase = read_datafile('PHASE.DAT')
        (taus,devs,errs,ns) = allan.ohdev(phase, taus=s32_taus)
        
        # CI computation
        # alhpa= +2,...,-4   noise power
        # d= 1 first-difference variance, 2 allan variance, 3 hadamard variance
        # alpha+2*d >1
        # m = tau/tau0 averaging factor
        # F filter factor, 1 modified variance, m unmodified variance
        # S stride factor, 1 nonoverlapped estimator, m overlapped estimator (estimator stride = tau/S )
        # N number of phase obs
        los=[]
        his=[]
        for (d,t, n) in zip(devs, taus, ns):
            #edf = greenhall_simple_edf( alpha=0, d=3, m=t, S=1, F=t, N=len(phase) )
            edf2 = allan.edf_greenhall( alpha=0, d=3, m=int(t), N=len(phase), overlapping = True, modified=False  )
            #print(edf,edf2,edf2/edf)
            (lo,hi) = allan.confidence_intervals( dev=d, ci=0.68268949213708585, edf=edf2 )  # 0.68268949213708585
            #allan.uncertainty_estimate(len(phase), t, d,ci=0.683,noisetype='wf')
            los.append(lo)
            his.append(hi)

        for (t1,d1,lo1,hi1, n1, t2, d2, lo2, hi2, n2) in zip(s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns, taus, devs, los, his, ns):
            print("s32 %03d %03f %1.6f %1.6f %1.6f" % (n1, t1, lo1, d1, hi1))
            print("at  %03d %03f %1.6f %1.6f %1.6f" % (n2, t2, round(lo2,5), round(d2,5), round(hi2,5) ))
            approx_equal(lo1, lo2,tolerance=1e-3)
            approx_equal(hi1, hi2,tolerance=1e-3)
            print("----")
            
    def test_phasedat_tdev(self):
        (s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns) = read_stable32( 'phase_dat_tdev_octave.txt' , 1.0 )
        phase = read_datafile('PHASE.DAT')
        (taus,devs,errs,ns) = allan.tdev(phase, taus=s32_taus)
        
        # CI computation
        # alhpa= +2,...,-4   noise power
        # d= 1 first-difference variance, 2 allan variance, 3 hadamard variance
        # alpha+2*d >1
        # m = tau/tau0 averaging factor
        # N number of phase obs
        los=[]
        his=[]
        for (d,t, n) in zip(devs, taus, ns):
            #edf = greenhall_simple_edf( alpha=0, d=2, m=t, S=1, F=t, N=len(phase) )
            edf2 = allan.edf_greenhall( alpha=0, d=2, m=int(t), N=len(phase), overlapping = True, modified=True  )
            
            # covert to mdev
            # taus * md / np.sqrt(3.0)
            mdev = d/t*np.sqrt(3.0)
            
            (lo,hi) = allan.confidence_intervals( dev=mdev, ci=0.68268949213708585, edf=edf2 )  # 0.68268949213708585
            
            # convert back to tdev
            lo = t*lo/np.sqrt(3.0)
            hi = t*hi/np.sqrt(3.0)
            
            los.append(lo)
            his.append(hi)
        
        for (t1,d1,lo1,hi1, n1, t2, d2, lo2, hi2, n2) in zip(s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns, taus, devs, los, his, ns):
            print("s32 %03d %03f %1.6f %1.6f %1.6f" % (n1, t1, lo1, d1, hi1))
            print("at  %03d %03f %1.6f %1.6f %1.6f" % (n2, t2, round(lo2,5), round(d2,5), round(hi2,5) ))
            approx_equal(lo1, lo2,tolerance=1e-3)
            approx_equal(hi1, hi2,tolerance=1e-3)
            print("----")

    def test_phasedat_totdev(self):
        (s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns) = read_stable32( 'phase_dat_totdev_octave.txt' , 1.0 )
        phase = read_datafile('PHASE.DAT')
        (taus,devs,errs,ns) = allan.totdev(phase, taus=s32_taus)
        

        los=[]
        his=[]
        for (d,t, n) in zip(devs, taus, ns):
            #edf = greenhall_simple_edf( alpha=0, d=3, m=t, S=1, F=t, N=len(phase) )
            #edf2 = greenhall_edf( alpha=0, d=3, m=int(t), N=len(phase), overlapping = True, modified=False  )
            #print(edf,edf2,edf2/edf)
            edf = allan.edf_totdev(len(phase), t, alpha=0)
            (lo,hi) = allan.confidence_intervals( dev=d, ci=0.68268949213708585, edf=edf )  # 0.68268949213708585
            #allan.uncertainty_estimate(len(phase), t, d,ci=0.683,noisetype='wf')
            los.append(lo)
            his.append(hi)

        for (t1,d1,lo1,hi1, n1, t2, d2, lo2, hi2, n2) in zip(s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns, taus, devs, los, his, ns):
            print("s32 %03d %03f %1.6f %1.6f %1.6f" % (n1, t1, lo1, d1, hi1))
            print("at  %03d %03f %1.6f %1.6f %1.6f" % (n2, t2, round(lo2,7), round(d2,7), round(hi2,7) ))
            approx_equal(lo1, lo2,tolerance=1e-3)
            approx_equal(hi1, hi2,tolerance=1e-3)
            print("----")
            
    def test_phasedat_mtotdev(self):
        (s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns) = read_stable32( 'phase_dat_mtotdev_octave_alpha0.txt' , 1.0 )
        phase = read_datafile('PHASE.DAT')
        (taus,devs,errs,ns) = allan.mtotdev(phase, taus=s32_taus)
        

        los=[]
        his=[]
        for (d,t, n) in zip(devs, taus, ns):
            #edf = greenhall_simple_edf( alpha=0, d=3, m=t, S=1, F=t, N=len(phase) )
            if int(t)<10:
                edf = allan.edf_greenhall( alpha=0, d=2, m=int(t), N=len(phase), overlapping = True, modified=True  )
            else:
                edf = allan.edf_mtotdev(len(phase), t, alpha=0)
            #print edf, edf2
            #print(edf,edf2,edf2/edf)
            
            print(edf)
            (lo,hi) = allan.confidence_intervals( dev=d, ci=0.68268949213708585, edf=edf )  # 0.68268949213708585
            #allan.uncertainty_estimate(len(phase), t, d,ci=0.683,noisetype='wf')
            los.append(lo)
            his.append(hi)

        for (t1,d1,lo1,hi1, n1, t2, d2, lo2, hi2, n2) in zip(s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns, taus, devs, los, his, ns):
            print("s32 %03d %03f %1.6f %1.6f %1.6f" % (n1, t1, lo1, d1, hi1))
            print("at  %03d %03f %1.6f %1.6f %1.6f" % (n2, t2, round(lo2,7), round(d2,7), round(hi2,7) ))
            approx_equal(lo1, lo2,tolerance=1e-3, assertions=False)
            approx_equal(hi1, hi2,tolerance=1e-3, assertions=False)
            print("----")

    
if __name__ == "__main__":
    #pytest.main()
    t = TestPhaseDatCI()
    t.test_phasedat_adev()
    t.test_phasedat_oadev() 
    t.test_phasedat_mdev()
    t.test_phasedat_tdev()
    t.test_phasedat_hdev()
    t.test_phasedat_ohdev()
    t.test_phasedat_totdev()
    #t.test_phasedat_mtotdev()
    
