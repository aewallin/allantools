"""
  PHASE.DAT test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  PHASE.DAT comes with Stable32 (version 1.53 was used in this case)

"""

import math
import sys
import pytest
import numpy as np
import scipy.stats

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

# greenhall ci computation

# alhpa= +2,...,-4   noise power
# d= 1 first-difference variance, 2 allaan variance, 3 hadamard variance
# alpha+2*d >1
# m = tau/tau0 averaging factor
# F filter factor, 1 modified variance, m unmodified variance
# S stride factor, 1 nonoverlapped estimator, m overlapped estimator (estimator stride = tau/S )
# N number of phase obs
J_max = 100

# this is Eqn (7) from Greenhall2004
def greenhall_sw(t, alpha):
    if alpha==2:
            return -np.abs(t)
    elif alpha==1:
        if t == 0:
            return 0
        else:
            return pow(t,2)*np.log( np.abs(t) )
    elif alpha==0:
        return np.abs( pow(t,3) )
    elif alpha==-1:
        if t == 0:
            return 0
        else:
            return pow(t,4)*np.log( np.abs(t) )
    elif alpha==-2:
        return np.abs( pow(t,5) )
    elif alpha==-3:
        if t == 0:
            return 0
        else:
            return pow(t,6)*np.log( np.abs(t) )
    elif alpha==-4:
        return np.abs( pow(t,7) )
        
    assert( 0 ) # ERROR

# this is Eqn (8) from Greenhall2004
def greenhall_sx(t, F, alpha):
    a = 2*greenhall_sw(t, alpha)
    b = greenhall_sw( t-1.0/F, alpha)
    c = greenhall_sw( t+1.0/F, alpha)
    
    # F=Inf case!?
    return pow(F,2)*(a-b-c)

# this is Eqn (9) from Greenhall2004
def greenhall_sz(t, F, alpha, d):
    if d==1:
        a = 2*greenhall_sx(t, F, alpha)
        b = greenhall_sx(t-1, F, alpha)
        c = greenhall_sx(t+1, F, alpha)
        return a-b-c
    elif d==2:
        a = 6*greenhall_sx(t, F, alpha)
        b = 4*greenhall_sx(t-1, F, alpha)
        c = 4*greenhall_sx(t+1, F, alpha)
        d = greenhall_sx(t-2, F, alpha)
        e = greenhall_sx(t+2, F, alpha)
        return a-b-c+d+e
    elif d==3:
        a = 20*greenhall_sx(t, F, alpha)
        b = 15*greenhall_sx(t-1, F, alpha)
        c = 15*greenhall_sx(t+1, F, alpha)
        d = 6*greenhall_sx(t-2, F, alpha)
        e = 6*greenhall_sx(t+2, F, alpha)
        f = greenhall_sx(t-3, F, alpha)
        g = greenhall_sx(t+3, F, alpha)
        return a-b-c+d+e-f-g
    
    assert( 0 ) # ERROR

def greenhall_BasicSum(J, M, S, F, alpha, d):
    first = pow( greenhall_sz(0, F, alpha, d), 2)
    second = (1-J/M)*pow( greenhall_sz( J/S, F, alpha, d ), 2)
    third = 0
    for j in range(1,J):
        third += 2*(1-j/M)*pow( greenhall_sz(j/S, F, alpha, d), 2)
    return first+second+third
    

def greenhall_simple_edf(alpha, d, m, S, F, N):
    L = m/F+m*d # length of filter applied to phase samples
    M = 1 + np.floor( S*(N-L) / m )
    J = min( M, (d+1)*S)
    inv_edf = (1.0/(pow(greenhall_sz(0,F,alpha,d),2)*M))*greenhall_BasicSum(J, M, S, F, alpha, d)
    return 1.0/inv_edf

def confidence_intervals(dev, ci, edf):
    ci_l = min(np.abs(ci), np.abs((ci-1))) / 2
    ci_h = 1 - ci_l
    
    chi2_l = scipy.stats.chi2.ppf(ci_l, edf)
    chi2_h = scipy.stats.chi2.ppf(ci_h, edf)
    
    # note these are variances, not deviations
    variance = dev*dev
    var_l = float(edf) * variance / chi2_h  # NIST SP1065 eqn (45) 
    var_h = float(edf) * variance / chi2_l
    return (np.sqrt(var_l), np.sqrt(var_h))

class TestPhaseDatCI():
    def test_phasedat_adev(self):
        (s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns) = read_stable32( 'phase_dat_adev_octave.txt' , 1.0 )
        phase = read_datafile('PHASE.DAT')
        (taus,devs,errs,ns) = allan.adev(phase, taus=s32_taus)
        los=[]
        his=[]
        for (d,t, n) in zip(devs, taus, ns):
            edf = greenhall_simple_edf( alpha=0, d=2, m=t, S=1, F=t, N=len(phase) )
            
            (lo,hi) = confidence_intervals( d, 0.683, edf ) 
            #allan.uncertainty_estimate(len(phase), t, d,ci=0.683,noisetype='wf')
            los.append(lo)
            his.append(hi)
        
        #print greenhall_simple_edf( alpha=0, d=2, m=1, S=1, F=1, N=len(phase) )
        #print confidence_intervals( dev
        #   tau N       edf         chi2_l
        #   1   999     782.030     821.358 , 742.689
        #  2.851498,  2.922319  2.998720
        
        for (t1,d1,lo1,hi1, n1, t2, d2, lo2, hi2, n2) in zip(s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns, taus, devs, los, his, ns):
            print(n1, t1, lo1, d1, hi1)
            print(n2, t2, lo2, d2, hi2)

            print("----")
        print(s32_devs_lo)
        print(los)
        
        #print taus, s32_taus
        
if __name__ == "__main__":
    #pytest.main()
    t = TestPhaseDatCI()
    t.test_phasedat_adev()
