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
import matplotlib.pyplot as plt # only for plotting, not required for calculations

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
def adev_phase(phase, rate, taus):
	freqdata = [ x*rate for x in numpy.diff( phase ) ] # convert from phase to fractional frequency
	return adev(freqdata,rate,taus)

# overlapping Allan deviation of phase data
def oadev_phase(phase, rate, taus):
	freqdata = [ x*rate for x in numpy.diff( phase ) ] # convert from phase to fractional frequency
	return oadev(freqdata,rate,taus)

# Time deviation of phase data
def tdev_phase(phase, rate, taus):
	freqdata = [ x*rate for x in numpy.diff( phase ) ] # convert from phase to fractional frequency
	return tdev(freqdata,rate,taus)

# Time deviation of fractional frequency data
# http://en.wikipedia.org/wiki/Time_deviation
def tdev(data, rate, taus):
	(taus2, md, mde, ns)  = mdev(data, rate, taus)
	td = [ t*m / math.sqrt(3.0) for (t,m) in zip(taus2,md)]
	tde = [x/math.sqrt(n) for (x,n) in zip(td,ns)]
	return (taus2, td, tde, ns)

# Modified Allan deviation of phase data
#
#                    1             N-3m+1     j+m-1
# Mod s2y(t) = ------------------  sum      { sum   [x(i+2m) - 2x(i+m) + x(i) ]  }**2
#              2m**2 t**2 (N-3m+1) j=1        i=j    
# 
# see http://www.leapsecond.com/tools/adev_lib.c
def mdev_phase(data,rate,taus):
	ms = tau_m(data,rate,taus)
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
		assert( n == len(data)-3*m+1 ) # n is the normalization (N-3m+1) before the sums
		s = math.sqrt( s ) 
		md.append( s ) 
		taus.append(tau)
		mderr.append( s / math.sqrt(n) )
		ns.append(n)
	return (taus, md, mderr, ns)

# modified Allan deviation, fractional frequency data
def mdev(freqdata, rate, taus):
	# integrrate freqdata to phase
	phase=[0]*( len(freqdata)+1 ) # initialize to zero
	for i in range( len(freqdata)+1 ):
		if i == 0:
			phase[i] = 0
		else:
			phase[i] = phase[i-1] + freqdata[i-1]*(1/float(rate))
	return mdev_phase(phase,rate,taus) 

# pre-processing of the tau-list given by the user
def tau_m(data,rate,taus):
	if rate == 0:
		print "Warning! rate==0"
	rate = float(rate)
	n = len(data)
	m=[]
	for tau in taus:
		if tau < (1/float(rate))*float(len(data))/3: # max tau corresponds to half-way point of dataseries
			m.append(  int( math.floor( float(tau*rate) )) )  # m is tau in units of datapoints
	m = list(set(m)) # this removes duplicates
	m.sort() # sort from small tau to large tau
	return m

# Allan deviation
# data is a time-series of fractional frequency
# rate is the samples/s in the time-series
# taus is a list of tau-values for which we compute ADEV	
def adev(data,rate, taus):
	m = tau_m(data,rate,taus)
	n = len(data)
	ad    = []
	ade   = []
	adn   = []
	for mj in m:  #loop through each tau value m(j)
		# D=zeros(1,n-m(j)+1);
		# D(1)=sum( data.freq(1:m(j)) ) / m(j);
		# for i=2:n-m(j)+1
		#     D(i)=D(i-1)+( data.freq(i+m(j)-1)-data.freq(i-1) )/m(j);
		# end
		D = [0]*( n-mj+1 ) # initialize with zeros
		D[0] = sum( data[ 0:(mj-1) ] ) / float( mj )
		for i in range(1, n-mj+1):
			D[i] = D[i-1] +( data[i-1+mj ] - data[i-1] )/ float( mj )
		allan = math.sqrt( 0.5*numpy.mean(  pow( numpy.diff( D[ 0:n-mj+1:mj ] ) , 2) ) )
		ad.append( allan ) 
		ad_n =  len( numpy.diff( D[ 0:n-mj+1:mj ] ) )
		ade.append( allan / math.sqrt( ad_n ) )
		adn.append( ad_n )
	return ([x/float(rate) for x in m], ad, ade, adn) # tau, adev, adeverror

# overlapping Allan deviation
def oadev(data,rate, taus):
	rate = float(rate)
	m = tau_m(data,rate,taus)
	n = len(data)
	odev    = []
	odeverr = []
	odevn = []
	for mj in m:
		D = [0]*( n-mj+1 ) # initialize with zeros
		D[0] = sum( data[ 0:(mj-1) ] ) / float( mj )
		for i in range(1, n-mj+1):
			D[i] = D[i-1] +( data[i-1+mj ] - data[i-1] )/ float( mj )

		z1 = D[ mj:n-mj+1  ]
		z2 = D[ 0:n-2*mj+1 ]
		assert( len(z1) == len(z2) )
		# first pair: 0 with m[j]
		# 1nd pair  : 1 with m[j]+1
		#...
		# last pair : n-2*m[j] with n-m[j]
		
		u=sum( [ (zz1-zz2)**2 for (zz1,zz2) in zip(z1,z2) ]);
		o = math.sqrt( u/ ( 2*float(n+1-2*mj) ) );
		odev.append( o )
		odeverr.append( o/ math.sqrt( len(z1) ) )
		odevn.append( len(z1) )
	return ([x/float(rate) for x in m], odev, odeverr, odevn)

# Hadamard deviation
def hdev(freqdata, rate, taus):
	# integrrate freqdata to phase
	phase=[0]*( len(freqdata)+1 ) # initialize to zero
	for i in range( len(freqdata)+1 ):
		if i == 0:
			phase[i] = 0
		else:
			phase[i] = phase[i-1] + freqdata[i-1]*(1/float(rate))
	return hdev_phase(phase,rate,taus) 

# Hadamard deviation of phase data
def hdev_phase(data,rate,taus):
	rate = float(rate)
	m = tau_m(data,rate,taus)
	n = len(data)

	hdevs = [] 
	hdeverrs = []
	ns = []
	for mj in m:
		(h,n ) = hdev_phase_calc(data,mj,mj)
		hdevs.append(h)
		hdeverrs.append( h/math.sqrt(n) )
		ns.append(n)
	taus = [x/float(rate) for x in m]
	return (taus, hdevs, hdeverrs, ns)

# http://www.leapsecond.com/tools/adev_lib.c
def hdev_phase_calc(data,mj, stride):
	s = 0
	n = 0
	i = 0
	while (i + 3*mj) < len(data):
		v = data[i + 3*mj] - 3 * data[i + 2*mj] + 3 * data[i + mj] - data[i]
		s = s + v * v
		n = n + 1 
		i = i + stride
	s = s/6.0
	h = math.sqrt( s / float(n)) / float(mj)
	return (h,n)

def plotallan(plt,y,rate,taus, style):
	(t2, ad, ade,adn) = oadev(y,rate,taus)
	plt.loglog(t2, ad, style)

def plotallan_phase(plt,y,rate,taus, style):
	(t2, ad, ade,adn) = oadev_phase(y,rate,taus)
	plt.loglog(t2, ad, style)

# plot a line with the slope alpha
def plotline(plt, alpha, taus,style):
	y = [ pow(tt,alpha) for tt in taus]
	plt.loglog(taus,y,style)


# pink noise generator
# from http://pydoc.net/Python/lmj.sound/0.1.1/lmj.sound.noise/
def iterpink(depth=20):
    '''Generate a sequence of samples of pink noise.

    Based on the Voss-McCartney algorithm, discussion and code examples at
    http://www.firstpr.com.au/dsp/pink-noise/

    depth: Use this many samples of white noise to calculate the output. A
      higher  number is slower to run, but renders low frequencies with more
      correct power spectra.

    Generates a never-ending sequence of floating-point values. Any continuous
    set of these samples will tend to have a 1/f power spectrum.
    '''
    values = numpy.random.randn(depth)
    smooth = numpy.random.randn(depth)
    source = numpy.random.randn(depth)
    sum = values.sum()
    i = 0
    while True:
        yield sum + smooth[i]

        # advance the index by 1. if the index wraps, generate noise to use in
        # the calculations, but do not update any of the pink noise values.
        i += 1
        if i == depth:
            i = 0
            smooth = numpy.random.randn(depth)
            source = numpy.random.randn(depth)
            continue

        # count trailing zeros in i
        c = 0
        while not (i >> c) & 1:
            c += 1

        # replace value c with a new source element
        sum += source[i] - values[c]
        values[c] = source[i]

def pinknoise(N):
	a=[]
	s = iterpink(80)
	for n in range(N):
		a.append( s.next() )
	
	return a

if __name__ == "__main__":
	
	# we test ADEV etc. by calculations on synthetic data
	# with known slopes of ADEV
	
	t = numpy.logspace( 0 ,4,50)
	plt.subplot(111, xscale="log", yscale="log")
	N=10000 
	
	# Colors: http://en.wikipedia.org/wiki/Colors_of_noise
	
	# pink frequency noise => constant ADEV 
	freq_pink = pinknoise(N)
	phase_p = numpy.cumsum( pinknoise(N) ) # integrate to get phase, color??
	plotallan_phase(plt, phase_p, 1, t, 'co')
	plotallan(plt, freq_pink, 1, t, 'c.')
	plotline(plt, 0, t ,'c')
	
	# white phase noise => 1/tau ADEV
	phase_white = numpy.random.randn(N)
	plotallan_phase(plt, phase_white, 1, t, 'ro')
	freq_w = numpy.diff( numpy.random.randn(N) ) # diff to get frequency, "Violet noise"
	plotallan(plt, freq_w, 1, t, 'r.')
	plotline(plt, -1, t ,'r')
	
	# white frequency modulation => 1/sqrt(tau) ADEV
	freq_white = numpy.random.randn(N)
	phase_rw = numpy.cumsum( numpy.random.randn(N) ) # integrate to get Brownian, or random walk phase
	plotallan(plt, freq_white, 1, t, 'b.')
	plotallan_phase(plt, phase_rw, 1, t, 'bo')
	plotline(plt, -0.5, t, 'b')
	
	# Brownian a.k.a random walk  frequency => sqrt(tau) ADEV
	freq_rw = numpy.cumsum( numpy.random.randn(N) )
	phase_rw_rw =numpy.cumsum( numpy.cumsum( numpy.random.randn(N) ) ) # integrate to get  phase
	plotallan(plt, freq_rw, 1, t, 'm.')
	plotallan_phase(plt, phase_rw_rw, 1, t, 'mo')
	plotline(plt, +0.5, t ,'m')
	
	plt.grid()
	plt.show()
