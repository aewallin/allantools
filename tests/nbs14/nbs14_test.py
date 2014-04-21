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
import time
import sys
#sys.path.append("../..") # hack to import from parent directory
# remove if you have allantools installed in your python path

import allantools as allan

# This program tests allantools (https://github.com/aewallin/allantools)
# with datasets that have known deviations.
#
# nbs14 datasets are from http://www.ieee-uffc.org/frequency-control/learning-riley.asp
# PHASE.DAT dataset comes with Stable32, and Stable32 was used to calculate deviation
# results in files such as phade_dat_adev.txt
#
# More datasets along with known deviations (or calculated with known-good programs) are welcome.
#
# results from allantools seem correct, they agree to within 4 to 6 digits of precision with other ADEV tools.
#

# small dataset and deviations from
# http://www.ieee-uffc.org/frequency-control/learning-riley.asp
# http://www.wriley.com/paper1ht.htm
nbs14_phase = [ 0.00000, 103.11111, 123.22222, 157.33333, 166.44444, 48.55555,-96.33333,-2.22222, 111.88889, 0.00000 ]
nbs14_f     = [892,809,823,798,671,644,883,903,677]
nbs14_devs= [ (91.22945,115.8082),  # ADEV(tau=1,tau=2)
              (91.22945, 85.95287), # OADEV 
              (91.22945,74.78849),  # MDEV
              (91.22945,98.31100),  # TOTDEV
              (70.80608,116.7980),  # HDEV
              (52.67135,86.35831) ] # TDEV 
              
# 1000 point dataset and deviations from
# http://www.ieee-uffc.org/frequency-control/learning-riley.asp
# http://www.wriley.com/paper1ht.htm
nbs14_1000_devs = [ [2.922319e-01, 9.965736e-02, 3.897804e-02],  # ADEV 1, 10, 100 
                    [2.922319e-01, 9.159953e-02, 3.241343e-02],  # OADEV
                    [2.922319e-01, 6.172376e-02, 2.170921e-02],  # MDEV 
                    [2.922319e-01, 9.172131e-02, 3.501795e-02],  # TOTDEV
                    [2.943883e-01, 1.052754e-01, 3.910860e-02],  # HDEV
                    [1.687202e-01, 3.563623e-01, 1.253382e-00] ] # TDEV

# this generates the nbs14 1000 point dataset.
# random number generator described in 
# http://www.ieee-uffc.org/frequency-control/learning-riley.asp
def nbs14_1000():
	n = [0]*1000
	n[0] = 1234567890
	for i in range(999):
		n[i+1] = (16807*n[i]) % 2147483647
	# the first three numbers are given in the paper, so check them:
	assert( n[1] ==  395529916 and n[2] == 1209410747 and  n[3] == 633705974 )
	n = [x/float(2147483647) for x in n] 
	return n

def nbs14_tester( function, fdata, correct_devs ):
	rate=1.0
	taus =[1, 10, 100]
	(taus2, devs, deverrs, ns) = function( fdata, rate, taus)
	for i in range(3):
		assert( check_devs( devs[i], correct_devs[i] ) )

# TODO: this does not work!!
def nbs14_totdev_test():
	rate=1.0
	taus =[1, 10, 100]
	fdata = nbs14_1000()
	correct_devs = nbs14_1000_devs[3] 
	(taus2, devs, deverrs, ns) = allan.totdev( fdata, rate, taus)
	(taus2, adevs, adeverrs, ans) = allan.adev( fdata, rate, taus)
	for i in range(3):
		totdev_corrected = devs[i]
		if i != 0:
			# bias removal is used in the published results
			a = 0.481 # flicker frequency-noise
			#a = 0.750 # flicker frequency-noise
			ratio = pow(devs[i],2)/pow(adevs[i],2)
			print ratio-1
			print -a*taus2[i]/((len(fdata)+1)*(1/float(rate)))
			ratio_corrected = ratio*( 1-a* taus2[i]/((len(fdata)+1)*(1/float(rate))) ) # WRONG!?!
			totdev_corrected = ratio_corrected * pow(adevs[i],2)
			totdev_corrected = math.sqrt( totdev_corrected )
			print totdev_corrected, devs[i], correct_devs[i]
		assert( check_devs( totdev_corrected, correct_devs[i] ) )

def nbs14_1000_test():
	fdata = nbs14_1000()

	nbs14_tester( allan.adev, fdata, nbs14_1000_devs[0] )
	print "nbs14_1000 adev"

	nbs14_tester( allan.oadev, fdata, nbs14_1000_devs[1] )
	print "nbs14_1000 oadev"

	nbs14_tester( allan.mdev, fdata, nbs14_1000_devs[2] )
	print "nbs14_1000 mdev"

	#print "nbs13 totdev" # this test does not work, becaus we don't know how to do bias correction
	#nbs14_totdev_test()

	nbs14_tester( allan.hdev, fdata, nbs14_1000_devs[4] )
	print "nbs14_1000 hdev"

	nbs14_tester( allan.tdev, fdata, nbs14_1000_devs[5] )
	print "nbs14_1000 tdev"
	
	#########################################################
	# now we test the same data, calling the _phase functions
	pdata = allan.frequency2phase(fdata, 1.0)
	
	nbs14_tester( allan.adev_phase, pdata, nbs14_1000_devs[0] )
	print "nbs14_1000 adev_phase"

	nbs14_tester( allan.oadev_phase, pdata, nbs14_1000_devs[1] )
	print "nbs14_1000 oadev_phase"

	nbs14_tester( allan.mdev_phase, pdata, nbs14_1000_devs[2] )
	print "nbs14_1000 mdev_phase"

	#print "nbs13 totdev" # this test does not work, becaus we don't know how to do bias correction
	#nbs14_totdev_test()

	nbs14_tester( allan.hdev_phase, pdata, nbs14_1000_devs[4] )
	print "nbs14_1000 hdev_phase"

	nbs14_tester( allan.tdev_phase, pdata, nbs14_1000_devs[5] )
	print "nbs14_1000 tdev_phase"
	
	print "nbs14_1000 all tests OK"

def check_devs(dev2, dev1):
	rel_error = (dev2-dev1)/dev1
	tol = 1e-6
	verbose = 0

	if ( abs(rel_error) < tol ):
		if verbose:
			print "OK   %0.6f \t 	%0.6f \t %0.6f" % (dev1,dev2, rel_error)
		return True
	else:
		print "ERROR   %0.6f \t %0.6f \t %0.6f" % (dev1,dev2, rel_error)
		return False

def nbs14_test():
	taus = [1, 2]
	devs = []
	tol = 1e-4
	
	# first tests that call the _phase functions
	(taus2,adevs2,aerrs2,ns2) = allan.adev_phase( nbs14_phase, 1.0, taus)
	adevs = nbs14_devs[0]
	assert( check_devs( adevs2[0], adevs[0] ) )
	assert( check_devs( adevs2[1], adevs[1] ) )
	print "nbs14 adev"
	(taus2,adevs2,aerrs2,ns2) = allan.oadev_phase( nbs14_phase, 1.0, taus)
	oadevs = nbs14_devs[1]
	assert( check_devs( adevs2[0], oadevs[0] ) )
	assert( check_devs( adevs2[1], oadevs[1] ) )
	print "nbs14 oadev"
	(taus2,adevs2,aerrs2,ns2) = allan.mdev_phase( nbs14_phase, 1.0, taus)
	mdevs = nbs14_devs[2]
	assert( check_devs( adevs2[0], mdevs[0] ) )
	assert( check_devs( adevs2[1], mdevs[1] ) )
	print "nbs14 mdev"
	(taus2,adevs2,aerrs2,ns2) = allan.hdev_phase( nbs14_phase, 1.0, taus)
	hdevs = nbs14_devs[4]
	assert( check_devs( adevs2[0], hdevs[0] ) )
	assert( check_devs( adevs2[1], hdevs[1] ) )
	print "nbs14 hdev"
	(taus2,adevs2,aerrs2,ns2) = allan.tdev_phase( nbs14_phase, 1.0, taus)
	tdevs = nbs14_devs[5]
	assert( check_devs( adevs2[0], tdevs[0] ) )
	assert( check_devs( adevs2[1], tdevs[1] ) )
	print "nbs14 tdev"

	# then the same tests for frequency data
	f_fract = nbs14_f
	(taus2,adevs2,aerrs2,ns2) = allan.adev( f_fract, 1.0, taus)
	adevs = nbs14_devs[0]
	assert( check_devs( adevs2[0], adevs[0] ) )
	assert( check_devs( adevs2[1], adevs[1] ) )
	print "nbs14 freqdata adev"

	(taus2,adevs2,aerrs2,ns2) = allan.oadev( f_fract, 1.0, taus)
	oadevs = nbs14_devs[1]
	assert( check_devs( adevs2[0], oadevs[0] ) )
	assert( check_devs( adevs2[1], oadevs[1] ) )
	print "nbs14 freqdata oadev"

	(taus2,adevs2,aerrs2,ns2) = allan.mdev( f_fract, 1.0, taus)
	mdevs = nbs14_devs[2]
	assert( check_devs( adevs2[0], mdevs[0] ) )
	assert( check_devs( adevs2[1], mdevs[1] ) )
	print "nbs14 freqdata mdev"

	(taus2,adevs2,aerrs2,ns2) = allan.hdev( f_fract, 1.0, taus)
	hdevs = nbs14_devs[4]
	assert( check_devs( adevs2[0], hdevs[0] ) )
	assert( check_devs( adevs2[1], hdevs[1] ) )
	print "nbs14 freqdata hdev"

	(taus2,adevs2,aerrs2,ns2) = allan.tdev( f_fract, 1.0, taus)
	tdevs = nbs14_devs[5]
	assert( check_devs( adevs2[0], tdevs[0] ) )
	assert( check_devs( adevs2[1], tdevs[1] ) )
	print "nbs14 freqdata tdev"

	print "nbs14 all test OK"


# read a simple data-file with phase or frequency numbers on each line
def read_datfile(filename):
	p=[]
	with open(filename) as f:
		for line in f:
			if line.startswith("#"): # skip comments
				pass
			else:
				p.append( float(line) )
	return p

# read a result-file, produced by copy/paste from Stable32
# note: header-lines need to be manually commented-out with "#"
def read_resultfile(filename):
	rows=[]
	row=[]
	with open(filename) as f:
		for line in f:
			if line.startswith("#"):
				pass
			else:
				row = []
				l2 = line.split(" ")
				l2 = filter(None, l2)
				for n in range(len(l2)):
					row.append( float(l2[n]) )
				rows.append(row)
	return rows

def run():
	nbs14_test()
	nbs14_1000_test()


if __name__ == "__main__":
	run()
