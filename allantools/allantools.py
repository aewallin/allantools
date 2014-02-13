#
# Allan deviation tools
# Anders Wallin (anders.e.e.wallin "at" gmail.com)
# v1.0 2014 January
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import math
import numpy


# Implemented functions:
# adev() and adev_phase()     Allan deviation
# oadev() and oadev_phase()   Overlapping Allan deviation
# mdev() and mdev_phase()     Modified Allan deviation
# tdev() and tdev_phase()     Time deviation ( modified variance scaled by (t^2/3) )
# hdev() and hdev_phase()     Hadamard deviation
# ohdev() and ohdev_phase()   Overlapping Hadamard deviation
# totdev() and totdev_phase() Total deviation
# mtie() and mtie_phase()     Maximum Time Interval Error
# tierms() and tierms_phase() Time Interval Error RMS

# to do:
# Modified Total
# The modified total variance, MTOT, is total version of the modified Allan variance.  
# It is defined for phase data as:
#                          1           N-3m+1  1  N+3m-1
# Mod s^2 total(t) = ----------------- sum     -- sum      [0zi*(m)]^2
#                     2m^2t0^2(N-3m+1) n=1     6m i=n-3m
# where the 0zi*(m) terms are the phase averages from a triply-extended 
# sequence created by uninverted even reflection at each end, 
# and the prefix 0 denotes that the linear trend has been removed.

# Time Total (modified total variance scaled by (t^2/3) )
# Hadamard Total


# References
# http://www.wriley.com/paper4ht.htm
# http://en.wikipedia.org/wiki/Allan_variance
#
# code see e.g.:
# http://www.mathworks.com/matlabcentral/fileexchange/26659-allan-v3-0
# http://www.mathworks.com/matlabcentral/fileexchange/26637-allanmodified
# http://www.leapsecond.com/tools/adev_lib.c

#
# Allan deviation of phase data
# Inputs:
# 	phase = list of phase measurements in seconds
# 	rate  = sample rate of data, i.e. interval between phase measurements is 1/rate
# 	taus  = list of tau-values for ADEV computation
# Output (tau_out, adev, adeverr, n)
# 	tau_out = list of tau-values for which ADEV was computed
# 	adev    = list of Allan deviations
# 	adeverr = list of estimated errors of allan deviations
# 	n       = list of number of pairs in allan computation. standard error is adeverr = adev/sqrt(n)



# Time deviation of phase data
def tdev_phase(phase, rate, taus):
	(taus2, md, mde, ns)  = mdev_phase(phase, rate, taus)
	td = [ t*m / math.sqrt(3.0) for (t,m) in zip(taus2,md)]
	tde = [x/math.sqrt(n) for (x,n) in zip(td,ns)]
	return (taus2, td, tde, ns)

# Time deviation of fractional frequency data
# http://en.wikipedia.org/wiki/Time_deviation
def tdev(data, rate, taus):
	phase = frequency2phase(data,rate)
	return tdev_phase(phase,rate,taus) 


# Modified Allan deviation of phase data
#
#                    1             N-3m+1     j+m-1
# Mod s2y(t) = ------------------  sum      { sum   [x(i+2m) - 2x(i+m) + x(i) ]  }**2
#              2m**2 t**2 (N-3m+1) j=1        i=j    
# 
# see http://www.leapsecond.com/tools/adev_lib.c
def mdev_phase(data,rate,taus):
	(ms,taus_used) = tau_m(data,rate,taus)
	taus = []
	md = []
	mderr = []
	ns = []
	for m in ms:
		s = 0
		v = 0
		n = 0
		tau = m/float(rate)
		i=0
		while (i+2*m) < len(data) and i < m:
			v = v + data[i+2*m] - 2*data[i+m] + data[i]
			i=i+1
		s = s + v*v
		n = n + 1

		i = 0
		while (i + 3*m) < len(data):
			v = v + data[i + 3*m] - 3*data[i + 2*m] + 3*data[i + m] - data[i]
			s = s + v * v
			n = n+1
			i = i+1
		s = s / float( 2.0 * m * m * tau * tau * n )
		#assert( n == len(data)-3*m+1 ) # n is the normalization (N-3m+1) before the sums
		s = math.sqrt( s ) 
		md.append( s ) 
		mderr.append( s / math.sqrt(n) )
		ns.append(n)
	return remove_small_ns(taus_used, md, mderr, ns)

# modified Allan deviation, fractional frequency data
def mdev(freqdata, rate, taus):
	phase = frequency2phase(freqdata,rate)
	return mdev_phase(phase,rate,taus) 

# pre-processing of the tau-list given by the user
def tau_m(data,rate,taus):
	if rate == 0:
		print "Warning! rate==0"
	rate = float(rate)
	n = len(data)
	m=[]
	for tau in taus:
		if tau > 0 and tau < (1/float(rate))*float(len(data)): # tau should be in [0, len(data)/rate]
			m.append(  int( math.floor( float(tau*rate) )) )  # m is tau in units of datapoints
	m = list(set(m)) # this removes duplicates
	m.sort() # sort from small tau to large tau
	#print "tau_m: ",m
	if len(m)==0:
		print "Warning: sanity-check on tau failed!"
		print "   len(data)=",len(data)," rate=",rate,"taus= ",taus
	taus2 = [x/float(rate) for x in m]
	return (m,taus2)

# Allan deviation
# data is a time-series of fractional frequency
# rate is the samples/s in the time-series
# taus is a list of tau-values for which we compute ADEV
def adev(data, rate, taus):
	phase = frequency2phase(data,rate)
	return adev_phase(phase,rate,taus) 

def adev_phase(data,rate,taus):
	(m,taus_used) = tau_m(data,rate,taus)
	ad    = []
	ade   = []
	adn   = []
	for mj in m:  #loop through each tau value m(j)
		(dev,deverr,n) = calc_adev_phase(data,rate,mj,mj)
		ad.append( dev ) 
		ade.append( deverr )
		adn.append( n )
	return remove_small_ns(taus_used, ad, ade, adn) # tau, adev, adeverror, naverages

def calc_adev_phase(data,rate,mj,stride):
	s=0
	n=0
	count = len(data)
	i=0
	while ((i+2*mj) < count):
		v = data[i+2*mj] - 2*data[i+mj] + data[i]
		s = s + v*v
		i = i + stride
		n=n+1
	s = s/float(2.0)
	dev = 0
	deverr = 0
	if not n==0:
		dev = math.sqrt( s /float(n)) / float(mj*(1/float(rate)))
		deverr = dev/math.sqrt(n)
	return (dev,deverr,n) 

# if n is small (==1), reject the result
def remove_small_ns(taus,devs,deverrs,ns):
	o_taus=[]
	o_dev = []
	o_err = []
	o_n = []
	for (t,d,e,n) in zip(taus,devs,deverrs,ns):
		if n > 1:
			o_taus.append(t)
			o_dev.append(d)
			o_err.append(e)
			o_n.append(n)
	return (o_taus, o_dev, o_err, o_n)

# overlapping Allan deviation of phase data
def oadev_phase(data, rate, taus):
	(m,taus_used) = tau_m(data,rate,taus)
	ad    = []
	ade   = []
	adn   = []
	for mj in m:
		(dev,deverr,n) = calc_adev_phase(data,rate,mj,1) # stride=1 for overlapping ADEV
		ad.append( dev ) 
		ade.append( deverr )
		adn.append( n )
	return remove_small_ns(taus_used, ad, ade, adn) # tau, adev, adeverror, naverages

# overlapping Allan deviation
def oadev(freqdata, rate, taus):
	phase = frequency2phase(freqdata,rate)
	return oadev_phase(phase,rate,taus) 

def phase2freqeuncy( phasedata, rate):
	freqdata = [ x*float(rate) for x in numpy.diff( phasedata ) ]
	return freqdata 

def frequency2phase( freqdata, rate):
	# integrrate freqdata to phase
	phase=[0]*( len(freqdata)+1 ) # initialize to zero
	dt = 1/float(rate)
	for i in range( len(freqdata)+1 ):
		if i == 0:
			phase[i] = 0
		else:
			phase[i] = phase[i-1] + freqdata[i-1]*dt
	return phase

# Overlapping Hadamard deviation
def ohdev(freqdata, rate, taus):
	phase= frequency2phase(freqdata, rate)
	return ohdev_phase(phase,rate,taus) 

# Overlapping Hadamard deviation of phase data
def ohdev_phase(data,rate,taus):
	rate = float(rate)
	(m,taus_used) = tau_m(data,rate,taus)
	hdevs = [] 
	hdeverrs = []
	ns = []
	for mj in m:
		(h,e,n) = hdev_phase_calc(data,rate,mj,1) # stride = 1
		hdevs.append(h)
		hdeverrs.append( e )
		ns.append(n)
	return remove_small_ns(taus_used, hdevs, hdeverrs, ns)
    
# Hadamard deviation
def hdev(freqdata, rate, taus):
	phase= frequency2phase(freqdata, rate)
	return hdev_phase(phase,rate,taus) 

# Hadamard deviation of phase data
def hdev_phase(data,rate,taus):
	rate = float(rate)
	(m,taus_used) = tau_m(data,rate,taus)
	hdevs = [] 
	hdeverrs = []
	ns = []
	for mj in m:
		(h,e, n ) = hdev_phase_calc(data,rate,mj,mj) # stride = mj
		hdevs.append(h)
		hdeverrs.append( e )
		ns.append(n)

	return remove_small_ns(taus_used, hdevs, hdeverrs, ns)

# http://www.leapsecond.com/tools/adev_lib.c
def hdev_phase_calc(data,rate,mj, stride):
	s = 0
	n = 0
	i = 0
	tau0 = 1/float(rate)
	while (i + 3*mj) < len(data):
		v = data[i + 3*mj] - 3 * data[i + 2*mj] + 3 * data[i + mj] - data[i]
		s = s + v * v
		n = n + 1 
		i = i + stride
	s = s/6.0

	if n == 0:
		n=1
	h = math.sqrt( s / float(n)) / float(tau0*mj)
	e = h/math.sqrt(n)
	return (h,e,n)

def totdev(freqdata,rate,taus):
	phasedata = frequency2phase( freqdata, rate)
	return totdev_phase(phasedata,rate,taus)

# See:
# David A. Howe, 
# The total deviation approach to long-term characterization
# of frequency stability, 
# IEEE tr. UFFC vol 47 no 5 (2000)
# 
#                  1         N-1
# totvar(t) = ------------  sum   [ x*(i-m) - 2x*(i)+x*(i+m) ]**2
#             2 t**2 (N-2)  i=2
# where x* is a new dataset with 'reflected' data at start/end
def totdev_phase(data,rate,taus):
	rate = float(rate)
	(m,taus_used) = tau_m(data,rate,taus)
	n = len(data)
	
	# totdev requires a new dataset 
	x1=[]
	for j in range(n-2):
		x1.append(float(2.0)*data[0]-data[j+1])
	x1.reverse()
	
	x2=[]
	for j in range(n-2):
		x2.append(float(2.0)*data[len(data)-1]-data[len(data)-1-(j+1)])

	x=[]
	for d in x1: # reflected data at start
		x.append(d)
	for d in data: # original data
		x.append(d)
	for d in x2: # reflected data at end
		x.append(d)

	# original dataset is now in the middle of the new dataset
	assert( data[0] == x[len(x1)] )
	devs=[]
	deverrs=[]
	ns=[]
	for mj in m:
		dev = 0
		ncount=0
		for i in range(1,len(data)-1):
			i = i + len(x1) # i=0 corresponds to start of original data
			dev = dev + pow( x[i-mj] - 2*x[i] + x[i+mj] , 2 )
			ncount=ncount+1
		dev = dev /float(2* pow(mj/rate, 2) * (n-2))
		dev = math.sqrt(dev)
		devs.append(dev)
		deverrs.append(dev/math.sqrt(ncount))
		ns.append(ncount) 

	return remove_small_ns(taus_used, devs, deverrs, ns)

def mtie(freqdata,rate,taus):
	phasedata = frequency2phase( freqdata, rate)
	return mtie_phase(phasedata,rate,taus)
    
# maximum time interval error
# this seems to correspond to Stable32 setting "Fast(u)"
# Stable32 also has "Decade" and "Octave" modes where the dataset is extended somehow?
def mtie_phase(phase, rate, taus):
	rate = float(rate)
	(m,taus_used) = tau_m(phase,rate,taus)
	n = len(phase)
	devs=[]
	deverrs=[]
	ns=[]
	for mj in m:
		dev=0
		ncount=0

		win_max=0
		win_min=0
		# we start by finding the max/min in the first window
		# when we slide the window forward by one step, one of these things happen:
		# - neither the new sample that enters the window, nor the old sample that leaves the window 
		#    changes the max/min of the window
		# - the new sample that enters the window is the new max/min of the window
		# - the old sample that leaves the window was the max/min of the window
		#    this is the most computationally expensive case, as we now have to search 
		#    the entire window for a new max/min
		for i in range(0,len(phase)-mj): # slide the start of the window over the dataset
			if i==0: # initialize window max/min
				win_max = max( phase[0:0+mj+1] ) # max and min in the first window
				win_min = min( phase[0:0+mj+1] )
			else:
				newsample = phase[i+mj] # the new sample that enters the window
				if (newsample > win_max):
					win_max=newsample
				elif (newsample < win_min):
					win_min = newsample
				oldindex = i-1
				if oldindex < 0:
					oldindex = 0
				oldsample = phase[i-1] # the old sample we throw away
				if (win_max == oldsample):
					win_max = max(phase[i:i+mj+1]) # must search for a new maximum
				if (win_min == oldsample):
					win_min = min(phase[i:i+mj+1]) # must search for new minimum sample
					
				#phases = phase[i:i+mj+1] # data window of length mj
            
			tie = win_max - win_min # largest error in this window
			if tie>dev:
				dev = tie
			ncount = ncount + 1
            
		devs.append(dev)
		deverrs.append(dev/math.sqrt(ncount))
		ns.append(ncount)

	return remove_small_ns(taus_used, devs, deverrs, ns)

def tierms(freqdata,rate,taus):
	phasedata = frequency2phase( freqdata, rate)
	return tierms_phase(phasedata,rate,taus)
    
# TIE rms
def tierms_phase(phase, rate, taus):
    rate = float(rate)
    (m,taus_used) = tau_m(phase,rate,taus)
    count = len(phase)
    devs=[]
    deverrs=[]
    ns=[]
    for mj in m:
        dev=0
        ncount=0
        tie = []
        for i in range( count-mj ):
            phases = [ phase[i], phase[i+mj] ] # pair of phases at distance mj from eachother
            tie.append( max(phases) - min(phases) ) # phase error
            ncount = ncount + 1
        # RMS of tie vector
        tie = [ pow(x,2) for x in tie ] # square
        tie =  numpy.mean( tie ) # mean
        tie = math.sqrt( tie ) # root
        devs.append(tie)
        deverrs.append(dev/math.sqrt(ncount))
        ns.append(ncount)
        
    return remove_small_ns(taus_used, devs, deverrs, ns)

# Three Cornered Hat Method
#
# Three clocks with unknown variances sa^2, sb^2, sc^3
# Three pairwise measurements give variances:
# sab^2, sbc^2, sca^2
# Assuming covariances are zero, we get:
# sa^2 = 0.5*( sab^2 + sca^2 - sbc^2 )
# (and cyclic permutations for sb and sc)
def three_cornered_hat_phase(phasedataAB, phasedataBC, phasedataCA, rate, taus, function):
	(tauAB,devAB,errAB,nsAB) = function(phasedataAB, rate, taus)
	(tauBC,devBC,errBC,nsBC) = function(phasedataBC, rate, taus)
	(tauCA,devCA,errCA,nsCA) = function(phasedataCA, rate, taus)
	varAB = [x*x for x in devAB]
	varBC = [x*x for x in devBC]
	varCA = [x*x for x in devCA]
	assert( len(varAB) == len(varBC) == len(varCA) )
	varA = [ 0.5*(ab + ca - bc) for (ab,bc,ca) in zip(varAB,varBC,varCA) ]
	devA = []
	for va in varA:
		try:
			devA.append( math.sqrt( va ) )
		except:
			devA.append( 0 )
	return ( tauAB, devA )

if __name__ == "__main__":
	print "Nothing to see here."

