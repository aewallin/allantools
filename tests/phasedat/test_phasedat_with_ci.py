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

# greenhall ci computation




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
    if F==float('inf'):
        return greenhall_sw(t,alpha+2)
    a = 2*greenhall_sw(t, alpha)
    b = greenhall_sw( t-1.0/float(F), alpha)
    c = greenhall_sw( t+1.0/float(F), alpha)
    
    return pow(F,2)*(a-b-c)

# this is Eqn (9) from Greenhall2004
def greenhall_sz(t, F, alpha, d):
    if d==1:
        a = 2*greenhall_sx(t, F, alpha)
        b = greenhall_sx(t-1.0, F, alpha)
        c = greenhall_sx(t+1.0, F, alpha)
        return a-b-c
    elif d==2:
        a = 6*greenhall_sx(t, F, alpha)
        b = 4*greenhall_sx(t-1.0, F, alpha)
        c = 4*greenhall_sx(t+1.0, F, alpha)
        dd = greenhall_sx(t-2.0, F, alpha)
        e = greenhall_sx(t+2.0, F, alpha)
        return a-b-c+dd+e
    elif d==3:
        a = 20.0*greenhall_sx(t    , F, alpha)
        b = 15.0*greenhall_sx(t-1.0, F, alpha)
        c = 15.0*greenhall_sx(t+1.0, F, alpha)
        dd = 6.0*greenhall_sx(t-2.0, F, alpha)
        e =  6.0*greenhall_sx(t+2.0, F, alpha)
        f =      greenhall_sx(t-3.0, F, alpha)
        g =      greenhall_sx(t+3.0, F, alpha)
        return a-b-c+dd+e-f-g
    
    assert( 0 ) # ERROR

# this is Eqn (10) from Greenhall2004
def greenhall_BasicSum(J, M, S, F, alpha, d):
    first = pow( greenhall_sz(0, F, alpha, d), 2)
    second = (1-float(J)/float(M))*pow( greenhall_sz( float(J)/float(S), F, alpha, d ), 2)
    third = 0
    for j in range(1,J):
        third += 2*(1.0-float(j)/float(M))*pow( greenhall_sz(float(j)/float(S), F, alpha, d), 2)
    return first+second+third
    
# this is Eqn (13) from Greenhall2004
def greenhall_simple_edf(alpha, d, m, S, F, N):
    L = m/F+m*d # length of filter applied to phase samples
    M = 1 + np.floor( S*(N-L) / m )
    J = min( M, (d+1)*S)
    inv_edf = (1.0/(pow(greenhall_sz(0,F,alpha,d),2)*M))*greenhall_BasicSum(J, M, S, F, alpha, d)
    return 1.0/inv_edf

# alhpa= +2,...,-4   noise power
# d= 1 first-difference variance, 2 Allan variance, 3 hadamard variance
# alpha+2*d >1
# m = tau/tau0 averaging factor
# F filter factor, 1 modified variance, m unmodified variance
# S stride factor, 1 nonoverlapped estimator, m overlapped estimator (estimator stride = tau/S )
# N number of phase obs

def greenhall_edf(alpha, d, m, N, overlapping = False, modified = False, verbose=True):
    if modified:
        F=1
    else:
        F=int(m)
    if overlapping:
        S=int(m)
    else:
        S=1
    assert( alpha+2*d > 1.0 )
    L = m/F+m*d # length of filter applied to phase samples
    M = 1 + np.floor( S*(N-L) / m )
    J = min( M, (d+1)*S)
    J_max = 100
    r = M/S
    if int(F)==1 and modified: # case 1, modified variances, all alpha
        if J <= J_max:
            inv_edf = (1.0/(pow( greenhall_sz(0,1,alpha,d),2)*M))*greenhall_BasicSum(J, M, S, 1, alpha, d)
            if verbose:
                print( "case 1.1 edf= %3f" % float(1.0/inv_edf) )
            return 1.0/inv_edf
        elif r>d+1:
            (a0, a1) = greenhall_table1(alpha, d)
            inv_edf = (1.0/r)*(a0-a1/r)
            if verbose:
                print( "case 1.2 edf= %3f" % float(1.0/inv_edf) )
            return 1.0/inv_edf
        else:
            m_prime = J_max/r
            inv_edf = (1.0/(pow( greenhall_sz(0,F,alpha,d),2)*J_max))*greenhall_BasicSum(J_max, J_max, m_prime, 1, alpha, d)
            if verbose:
                print( "case 1.3 edf= %3f" % float(1.0/inv_edf) )
            return 1.0/inv_edf
    elif int(F)==int(m) and int(alpha)<= 0 and not modified: # case 2, unmodified variances, alpha <= 0
        if J <= J_max:
            if m*(d+1)<= J_max:
                m_prime = m
                variant = "a"
            else:
                m_prime = float('inf')
                variant = "b"
            
            inv_edf = (1.0/(pow( greenhall_sz(0,m_prime,alpha,d),2)*M))*greenhall_BasicSum(J, M, S, m_prime, alpha, d)
            if verbose:
                print( "case 2.1%s edf= %3f" % (variant, float(1.0/inv_edf)) )
            return 1.0/inv_edf
        elif r>d+1:
            (a0, a1) = greenhall_table2(alpha, d)
            inv_edf = (1.0/r)*(a0-a1/r)
            if verbose:
                print( "case 2.2 edf= %3f" % float(1.0/inv_edf) )
            return 1.0/inv_edf
        else:
            m_prime = J_max/r
            inv_edf = (1.0/(pow( greenhall_sz(0,float('inf'),alpha,d),2)*J_max))*greenhall_BasicSum(J_max, J_max, m_prime, float('inf'), alpha, d)
            if verbose:
                print( "case 2.3 edf= %3f" % float(1.0/inv_edf) )
            return 1.0/inv_edf
    elif int(F)==int(m) and int(alpha)==1 and not modified: # case 3, unmodified variances, alpha=1
        if J <= J_max:
            inv_edf = (1.0/(pow( greenhall_sz(0,m,1,d),2)*M))*greenhall_BasicSum(J, M, S, m, 1, d) # note: m<1e6 to avoid roundoff
            return 1.0/inv_edf
        elif r>d+1:
            (a0, a1) = greenhall_table2(alpha, d)
            (b0, b1) = greenhall_table3(alpha, d)
            inv_edf = (1.0/(pow(b0+b1*np.log(m),2)*r))*(a0-a1/r)
            return 1.0/inv_edf
        else:
            m_prime=J_max/r
            (b0, b1) = greenhall_table3(alpha, d)
            inv_edf = (1.0/(pow(b0+b1*np.log(m),2)*J_max))*greenhall_BasicSum(J_max, J_max, m_prime, m_prime, 1, d)
            return 1.0/inv_edf
    elif int(F)==int(m) and int(alpha)==2 and not modified: # case 4, unmodified variances, alpha=2
        K = np.ceil(r)
        if K <= d:
            pass
        else:
            a0 = scipy.special.binom(4*d, 2*d) / pow(scipy.special.binom(2*d,d), 2)
            a1 = d/2.0
            inv_edf = (1.0/M)*(a0-a1/r)
            return 1.0/inv_edf
            
    print("greenhall_edf() no matching case!")
    assert(0) # ERROR
            
def greenhall_table3(alpha, d):
    assert(alpha==1)
    idx = d-1
    table=[ (6.0,4.0), (15.23,12.0), (47.8,40.0) ]
    return table3[idx]
                
def greenhall_table2(alpha, d):
    row_idx = -alpha+2 # map 2-> row0 and -4-> row6
    col_idx = d-1
    table2 = [ [ (3.0/2.0, 1.0/2.0) , (35.0/18.0, 1.0)  , (231.0/100.0, 3.0/2.0) ], # alpha=+2
               [ (78.6,25.2 )       , (790.0,410.0)     , (9950.0,6520.0) ],
               [ (2.0/3.0,1.0/6.0)  , (2.0/3.0,1.0/3.0) , (7.0/9.0,1.0/2.0) ], # alpha=0
               [ (-1,-1)            , (0.852,0.375)     , (0.997,0.617) ], # -1
               [ (-1,-1)            , (1.079,0.368)     , (1.033,0.607) ], #-2
               [ (-1,-1)            , (-1,-1)           , (1.053,0.553) ], #-3
               [ (-1,-1)            , (-1,-1)           , (1.302,0.535)], # alpha=-4
            ]
    #print("table2 = ", table2[row_idx][col_idx])
    return table2[row_idx][col_idx]

def greenhall_table1(alpha, d):
    row_idx = -alpha+2 # map 2-> row0 and -4-> row6
    col_idx = d-1
    table1 = [ [ (2.0/3.0, 1.0/3.0) , (7.0/9.0, 1.0/2.0)    , (22.0/25.0, 2.0/3.0) ], # alpha=+2
               [ (0.840,0.345)      , (0.997,0.616)         , (1.141,0.843) ],
               [ (1.079,0.368)      , (1.033,0.607)         , (1.184,0.848) ],
               [ (-1,-1)            , (1.048,0.534)         , (1.180,0.816) ], # -1
               [ (-1,-1)            , (1.302,0.535)         , (1.175,0.777) ], #-2
               [ (-1,-1)            , (-1,-1)               , (1.194,0.703) ], #-3
               [ (-1,-1)            , (-1,-1)               , (1.489,0.702) ], # alpha=-4
            ]
    #print("table1 = ", table1[row_idx][col_idx])
    return table1[row_idx][col_idx]
    
def confidence_intervals(dev, ci, edf):
    # for 1-sigma standard error set
    # ci = scipy.special.erf(1/math.sqrt(2))
    #    = 0.68268949213708585
    
    ci_l = min(np.abs(ci), np.abs((ci-1))) / 2
    ci_h = 1 - ci_l
    
    chi2_l = scipy.stats.chi2.ppf(ci_l, edf)
    chi2_h = scipy.stats.chi2.ppf(ci_h, edf)
    
    # note these are variances, not deviations
    variance = dev*dev
    var_l = float(edf) * variance / chi2_h  # NIST SP1065 eqn (45) 
    var_h = float(edf) * variance / chi2_l
    #print edf, chi2_h,chi2_l
    return (np.sqrt(var_l), np.sqrt(var_h))

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
            edf2 = greenhall_edf( alpha=0, d=2, m=t, N=len(phase), overlapping = False, modified=False )
            #print(edf,edf2,edf2/edf)
            (lo,hi) = confidence_intervals( dev=d, ci=0.68268949213708585, edf=edf2 )  # 0.68268949213708585
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
            edf2 = greenhall_edf( alpha=0, d=2, m=int(t), N=len(phase), overlapping = True, modified=False  )
            #print(edf,edf2,edf2/edf)
            (lo,hi) = confidence_intervals( dev=d, ci=0.68268949213708585, edf=edf2 )  # 0.68268949213708585
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
            edf2 = greenhall_edf( alpha=0, d=2, m=int(t), N=len(phase), overlapping = True, modified=True  )
            #print(edf,edf2,edf2/edf)
            (lo,hi) = confidence_intervals( dev=d, ci=0.68268949213708585, edf=edf2 )  # 0.68268949213708585
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
        # F filter factor, 1 modified variance, m unmodified variance
        # S stride factor, 1 nonoverlapped estimator, m overlapped estimator (estimator stride = tau/S )
        # N number of phase obs
        los=[]
        his=[]
        for (d,t, n) in zip(devs, taus, ns):
            #edf = greenhall_simple_edf( alpha=0, d=3, m=t, S=1, F=t, N=len(phase) )
            edf2 = greenhall_edf( alpha=0, d=3, m=int(t), N=len(phase), overlapping = False, modified=False  )
            #print(edf,edf2,edf2/edf)
            (lo,hi) = confidence_intervals( dev=d, ci=0.68268949213708585, edf=edf2 )  # 0.68268949213708585
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
            
    def test_phasedat_tdev(self):
        (s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns) = read_stable32( 'phase_dat_tdev_octave.txt' , 1.0 )
        phase = read_datafile('PHASE.DAT')
        (taus,devs,errs,ns) = allan.tdev(phase, taus=s32_taus)
        
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
            edf2 = greenhall_edf( alpha=0, d=2, m=int(t), N=len(phase), overlapping = True, modified=True  )
            
            # covert to mdev
            # taus * md / np.sqrt(3.0)
            mdev = d/t*np.sqrt(3.0)
            
            (lo,hi) = confidence_intervals( dev=mdev, ci=0.68268949213708585, edf=edf2 )  # 0.68268949213708585
            
            # convert back to tdev
            lo = t*lo/np.sqrt(3.0)
            hi = t*hi/np.sqrt(3.0)
            
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
            
if __name__ == "__main__":
    #pytest.main()
    t = TestPhaseDatCI()
    t.test_phasedat_adev()
    t.test_phasedat_oadev() 
    t.test_phasedat_mdev()
    t.test_phasedat_hdev()
    t.test_phasedat_tdev()
