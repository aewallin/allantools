import allantools as allan
import numpy
import math
import matplotlib.pyplot as plt
import gzip
import sys

sys.path.append("..")
sys.path.append("../..") # hack to import from parent directory
# remove if you have allantools installed in your python path

def read_datafile(filename,column=1):
    p=[]
    with gzip.open(filename) as f:
        for line in f:
            
            if line.startswith("#"): # skip comments
                pass
            else:
                line = line.split()
                p.append( float(line[column]) )
    return p

def to_fractional(flist,f0):
    out=[]
    for f in flist:
        out.append( f/float(f0) - 1.0 )
    return out

fname = "gps_1pps_phase_data.txt.gz"
phase = read_datafile(fname,column=0)
print(  "%d values read: %.3f hours" % (len(phase), len(phase)/3600.0 ) )

#f = to_fractional(f10MHz, 10e6 ) # convert to fractional frequency
my_taus = numpy.logspace(1,6,60) # log-spaced tau values from 10s and upwards
rate = 1/float(1.0)

(adev_taus,adev_devs,adev_errs,ns)  = allan.adev(phase, rate=rate, taus=my_taus)
 
(oadev_taus,oadev_devs,oadev_errs,ns)  = allan.oadev(phase, rate=rate, taus=my_taus)

(hdev_taus,hdev_devs,hdev_errs,ns)  = allan.hdev(phase, rate=rate, taus=my_taus)
(ohdev_taus,ohdev_devs,ohdev_errs,ns)  = allan.ohdev(phase, rate=rate, taus=my_taus)

(mdev_taus,mdev_devs,mdev_errs,ns)  = allan.mdev(phase, rate=rate, taus=my_taus)

(totdev_taus,totdev_devs,totdev_errs,ns)  = allan.totdev(phase, rate=rate, taus=my_taus)

(tie_taus,tie_devs,tie_errs,ns)  = allan.tierms(phase, rate=rate, taus=my_taus)
#(mtie_taus,mtie_devs,mtie_errs,ns)  = allan.mtie(phase=phase, rate=rate, taus=my_taus)

(tdev_taus,tdev_devs,tdev_errs,ns)  = allan.tdev(phase, rate=rate, taus=my_taus)
(tdev2_taus,tdev2_devs,tdev2_errs,ns2)  = allan.tdev(allan.phase2frequency(phase,1.0), data_type='freq', rate=rate, taus=my_taus)

plt.subplot(111, xscale="log", yscale="log") 


plt.errorbar(adev_taus, adev_devs, yerr=adev_errs, label='ADEV') 
plt.errorbar(oadev_taus, oadev_devs, yerr=oadev_errs, label='OADEV') 
plt.errorbar(mdev_taus, mdev_devs, yerr=mdev_errs, label='MDEV') 
plt.errorbar(hdev_taus, hdev_devs, yerr=hdev_errs, label='HDEV') 
plt.errorbar(ohdev_taus, ohdev_devs, yerr=ohdev_errs, label='OHDEV') 
plt.errorbar(tdev_taus, tdev_devs, yerr=tdev_errs, label='TDEV(phase)') 
plt.errorbar(tdev2_taus, tdev2_devs, yerr=tdev2_errs, label='TDEV(frequency)')

plt.errorbar(totdev_taus, totdev_devs, yerr=totdev_errs, label='TOTDEV') 
plt.errorbar(tie_taus, tie_devs, yerr=tie_errs, label='TIERMS') 
#plt.errorbar(mtie_taus, mtie_devs, yerr=mtie_errs, label='MTIE') 

labels = [(60,'1 m'), (60*60,'1 h'), (12*60*60,'12 h'), (24*60*60,'1 d')]
for l in labels:
    plt.text( l[0], 1e-13, l[1] )
plt.xlabel('Taus (s)')
plt.ylabel('ADEV')


plt.legend(framealpha=0.5)
plt.grid()
plt.show()
