import allantools as allan
import numpy
import math
import matplotlib.pyplot as plt

import sys
sys.path.append("..")
sys.path.append("../..") # hack to import from parent directory
# remove if you have allantools installed in your python path

""" this file for testing noise-ID on the ocxo dataset """

def read_datafile(filename,column=1):
    p=[]
    with open(filename) as f:
        for line in f:
            
            if line.startswith("#"): # skip comments
                pass
            else:
                line = line.split()
                p.append( float(line[column]) )
    return p
    
def read_Stable32_resultfile(filename):
    data=[]
    keys=[ "AF", "Tau", "N", "Alpha", "min", "sigma", "max" ]
    with open(filename) as f:
        for line in f:
            
            if line.startswith("#"): # skip comments
                pass
            else:
                line = line.split()
                row={}
                for column in range(7):
                    row[ keys[column] ] = float(line[column])
                print row
                data.append(row)
    return data
    
def to_fractional(flist,f0):
    out=[]
    for f in flist:
        out.append( f/float(f0) - 1.0 )
    return out

# Rn ratio= modified allan variance to allan variance
# used for MVAR, TVAR, MTOT
# but not useful at AF=1

def average(arr, n):
    end =  n * int(len(arr)/n)
    return numpy.mean(arr[:end].reshape(-1, n), 1)
    
def B1_noise_id(x,rate, m):
    """ Barnes B1 noise identification function 
        ratio of N-sample variance to 2-sample variance (AVAR)
        
        http://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote375.pdf
        
    """
    y = allan.phase2frequency(x, rate=rate)
    
    (taus, adevs, errs, ns)  =  allan.oadev(x,rate=rate, taus=[ m*1.0/rate ])
    #print adevs
    averaged_y = average( y, m )
    print len(averaged_y)
    n_var = numpy.var( averaged_y )
    
    avar = pow( adevs[0],2 )
    #print avar
    #print v
    ratio = n_var / avar
    print n_var, avar, "ratio = ",ratio, numpy.round(ratio)
    return numpy.round(ratio)

def autocorr(x, lag=1):
    #print x
    c = numpy.corrcoef( numpy.array(x[:-lag]), numpy.array(x[lag:]) )
    #print c
    return c[0,1]
    
def autocorr_noise_id(x, data_type="phase", dmin=0, dmax=2):
    # for Hadamard, set dmax=3
    done = False
    d = 0 # number of differentiations
    #dmin = 0
    #dmax = 3
    while not done:
        r1 = autocorr(x)
        rho = r1/(1.0+r1)
        if d >= dmin and ( rho < 0.25 or d >= dmax ):
            p = -2*(rho+d)
            #print r1
            #assert r1 < 0
            #assert r1 > -1.0/2.0
            alpha = p+2.0
            alpha_int = int( -1*numpy.round(2*rho) - 2*d )+2 
            #print "d=",d,"alpha=",p+2
            return alpha_int, alpha,r1, d
        else:
            x = numpy.diff(x)
            d = d + 1
            
    print "autocorr= ",c

fname = "ocxo_frequency.txt"
rate = 1/float(10) # data collected with 10s gate time

f10MHz = read_datafile(fname,column=0)

#s32results = read_Stable32_resultfile("stable32_oadev_alltau.txt")
s32results = read_Stable32_resultfile("tdev_octave.txt")

f = to_fractional(f10MHz, 10e6 ) # convert to fractional frequency
x = allan.frequency2phase(f, rate) # phase in seconds

my_taus = numpy.logspace(1,5,40) # log-spaced tau values from 10s and upwards

for row in s32results:
    m = row['AF']
    alpha = row['Alpha']
    #noise_id = B1_noise_id(x,rate,m)
    x_decimated = x[0:len(x):int(m)]
    # print "len ",len(x_decimated) # should check that length is about 32 or longer-ish.
    est_int, est,r1,d = autocorr_noise_id(x_decimated)
    #if len(x_decimated) < 32:
    #    est = autocorr_noise_id(x_decimated)
    print m,len(x_decimated),alpha,est_int,alpha==est_int, est,r1,d

(oadev_taus,oadev_devs,oadev_errs,ns)  = allan.oadev(f, data_type='freq', rate=rate, taus='octave')
#(mdev_taus,mdev_devs,mdev_errs,ns)  = allan.mdev(f, data_type='freq', rate=rate, taus=my_taus)
#(hdev_taus,hdev_devs,hdev_errs,ns)  = allan.hdev(f, data_type='freq', rate=rate, taus=my_taus)
#(ohdev_taus,ohdev_devs,ohdev_errs,ns)  = allan.ohdev(f, data_type='freq', rate=rate, taus=my_taus)

plt.subplot(111, xscale="log", yscale="log") 

plt.errorbar(oadev_taus, oadev_devs, yerr=oadev_errs, label='OADEV') 
#plt.errorbar(mdev_taus, mdev_devs, yerr=mdev_errs, label='MDEV') 
#plt.errorbar(hdev_taus, hdev_devs, yerr=hdev_errs, label='HDEV') 
#plt.errorbar(ohdev_taus, ohdev_devs, yerr=ohdev_errs, label='OHDEV') 

for row in s32results:
    plt.plot( 10*row['Tau'], row['sigma'], 'o')
    plt.text( 10*row['Tau'], 1.5*row['sigma'], '%d' % row['Alpha'])
    
plt.xlabel('Taus (s)')
plt.ylabel('ADEV')

plt.xlim((8,1e5))
plt.ylim((1e-12,1e-9))

plt.legend()
plt.grid()
plt.show()
