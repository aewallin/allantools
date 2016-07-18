"""
  PHASE.DAT test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  PHASE.DAT comes with Stable32 (version 1.53 was used in this case)

"""

import math
import sys
import pytest

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

class TestPhaseDatCI():
    def test_phasedat_adev(self):
        (s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns) = read_stable32( 'phase_dat_adev_octave.txt' , 1.0 )
        phase = read_datafile('PHASE.DAT')
        (taus,devs,errs,ns) = allan.adev(phase, taus=s32_taus)
        los=[]
        his=[]
        for (d,t, n) in zip(devs, taus, ns):
            (lo,hi) = allan.uncertainty_estimate(len(phase), t, d,ci=0.683,noisetype='wf')
            los.append(lo)
            his.append(hi)
        
        #   tau N       edf         chi2_l
        #   1   999     782.030     821.358 , 742.689
        #  2.851498,  2.922319  2.998720
        for (t1,d1,lo1,hi1, n1, t2, d2, lo2, hi2, n2) in zip(s32_taus, s32_devs, s32_devs_lo, s32_devs_hi, s32_ns, taus, devs, los, his, ns):
            print n1, t1, lo1, d1, hi1
            print n2, t2, lo2, d2, hi2

            print "----"
        print s32_devs_lo
        print los
        #print taus, s32_taus
        
if __name__ == "__main__":
    #pytest.main()
    t = TestPhaseDatCI()
    t.test_phasedat_adev()
