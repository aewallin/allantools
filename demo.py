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


import numpy
import matplotlib.pyplot as plt # only for plotting, not required for calculations

import allantools
import noise

def plotallan(plt,y,rate,taus, style):
	(t2, ad, ade,adn) = allantools.oadev(y,rate,taus)
	plt.loglog(t2, ad, style)

def plotallan_phase(plt,y,rate,taus, style):
	(t2, ad, ade,adn) = allantools.oadev_phase(y,rate,taus)
	plt.loglog(t2, ad, style)

# plot a line with the slope alpha
def plotline(plt, alpha, taus,style):
	y = [ pow(tt,alpha) for tt in taus]
	plt.loglog(taus,y,style)
    
if __name__ == "__main__":
	print "allatools tests with various noise types"
	# we test ADEV etc. by calculations on synthetic data
	# with known slopes of ADEV
	
	t = numpy.logspace( 0 ,3,50) # tau values from 1 to 1000
	plt.subplot(111, xscale="log", yscale="log")
	N=10000 
	
	# Colors: http://en.wikipedia.org/wiki/Colors_of_noise
	
	# pink frequency noise => constant ADEV 
	freq_pink = noise.pink(N)
	phase_p = numpy.cumsum( noise.pink(N) ) # integrate to get phase, color??
	plotallan_phase(plt, phase_p, 1, t, 'co')
	plotallan(plt, freq_pink, 1, t, 'c.')
	plotline(plt, 0, t ,'c')
	print "Pink frequency noise"
    
	# white phase noise => 1/tau ADEV
	phase_white = noise.white(N)
	plotallan_phase(plt, phase_white, 1, t, 'ro')
	freq_w = noise.violet(N) # diff to get frequency, "Violet noise"
	plotallan(plt, freq_w, 1, t, 'r.')
	plotline(plt, -1, t ,'r')
	print "White phase noise"
    
	# white frequency modulation => 1/sqrt(tau) ADEV
	freq_white = noise.white(N)
	phase_rw = noise.brown(N) # integrate to get Brownian, or random walk phase
	plotallan(plt, freq_white, 1, t, 'b.')
	plotallan_phase(plt, phase_rw, 1, t, 'bo')
	plotline(plt, -0.5, t, 'b')
	print "White frequency noise"
    
	# Brownian a.k.a random walk  frequency => sqrt(tau) ADEV
	freq_rw = noise.brown(N) 
	phase_rw_rw =numpy.cumsum( noise.brown(N) ) # integrate to get  phase
	plotallan(plt, freq_rw, 1, t, 'm.')
	plotallan_phase(plt, phase_rw_rw, 1, t, 'mo')
	plotline(plt, +0.5, t ,'m')
	print "Random Walk frequency noise"
    
	plt.grid()
	plt.show()
