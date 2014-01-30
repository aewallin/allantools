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

import allantools as allan

# This program tests allantools (https://github.com/aewallin/allantools)
# with datasets that have known deviations.
#
# nbs14 datasets are from http://www.ieee-uffc.org/frequency-control/learning-riley.asp
# PHASE.DAT dataset comes with Timelab32, and Timelab32 was used to calculate deviation
# results in files such as phade_dat_adev.txt
#
# More datasets along with known deviations (or calculated with known-good programs) are welcome.
#
# results from allantools seem correct, they agree to within 4 to 6 digits of precision with other ADEV tools.
#

# small dataset from
# http://www.ieee-uffc.org/frequency-control/learning-riley.asp
nbs14_phase = [ 0.00000, 103.11111, 123.22222, 157.33333, 166.44444, 48.55555,-96.33333,-2.22222, 111.88889, 0.00000 ]
nbs14_f     = [892,809,823,798,671,644,883,903,677]
nbs14_devs= [ (91.22945,115.8082),  # ADEV(tau=1,tau=2)
              (91.22945, 85.95287), # OADEV 
              (91.22945,74.78849),  # MDEV
              (91.22945,98.31100),  # TOTDEV
              (70.80608,116.7980),  # HDEV
              (52.67135,86.35831) ] # TDEV 

nbs14_1000_devs = [ [2.922319e-01, 9.965736e-02, 3.897804e-02],  # ADEV 1, 10, 100 
                    [2.922319e-01, 9.159953e-02, 3.241343e-02],  # OADEV
                    [2.922319e-01, 6.172376e-02, 2.170921e-02],  # MDEV 
                    [2.922319e-01, 9.172131e-02, 3.501795e-02],  # TOTDEV
                    [2.943883e-01, 1.052754e-01, 3.910860e-02],  # HDEV
                    [1.687202e-01, 3.563623e-01, 1.253382e-00] ] # TDEV

# random number generator described in http://www.ieee-uffc.org/frequency-control/learning-riley.asp
def nbs14_1000():
	n = [0]*1000
	n[0] = 1234567890
	for i in range(999):
		n[i+1] = (16807*n[i]) % 2147483647
	#print n[0:3]
	assert( n[1] ==  395529916 and n[2] == 1209410747 and  n[3] == 633705974 )
	n = [x/float(2147483647) for x in n] 
	return n

def nbs14_tester( function, fdata, correct_devs ):
	rate=1.0
	taus =[1, 10, 100]
	(taus2, devs, deverrs, ns) = function( fdata, rate, taus)
	for i in range(3):
		assert( check_devs( devs[i], correct_devs[i] ) )
	
def nbs14_1000_test():
	fdata = nbs14_1000()
	nbs14_tester( allan.adev, fdata, nbs14_1000_devs[0] )
	nbs14_tester( allan.oadev, fdata, nbs14_1000_devs[1] )
	nbs14_tester( allan.mdev, fdata, nbs14_1000_devs[2] )
	nbs14_tester( allan.hdev, fdata, nbs14_1000_devs[4] )
	nbs14_tester( allan.tdev, fdata, nbs14_1000_devs[5] )

	print "nbs14_1000 test OK"
	
def check_devs(dev2, dev1):
	rel_error = (dev2-dev1)/dev1
	tol = 1e-6
	verbose = 1
	if verbose:
		print "   %0.6f \t 	%0.6f \t %0.6f" % (dev1,dev2, rel_error)
	if ( abs(rel_error) < tol ):
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

	(taus2,adevs2,aerrs2,ns2) = allan.oadev_phase( nbs14_phase, 1.0, taus)
	oadevs = nbs14_devs[1]
	assert( check_devs( adevs2[0], oadevs[0] ) )
	assert( check_devs( adevs2[1], oadevs[1] ) )

	(taus2,adevs2,aerrs2,ns2) = allan.mdev_phase( nbs14_phase, 1.0, taus)
	mdevs = nbs14_devs[2]
	assert( check_devs( adevs2[0], mdevs[0] ) )
	assert( check_devs( adevs2[1], mdevs[1] ) )

	(taus2,adevs2,aerrs2,ns2) = allan.hdev_phase( nbs14_phase, 1.0, taus)
	hdevs = nbs14_devs[4]
	assert( check_devs( adevs2[0], hdevs[0] ) )
	assert( check_devs( adevs2[1], hdevs[1] ) )

	(taus2,adevs2,aerrs2,ns2) = allan.tdev_phase( nbs14_phase, 1.0, taus)
	tdevs = nbs14_devs[5]
	assert( check_devs( adevs2[0], tdevs[0] ) )
	assert( check_devs( adevs2[1], tdevs[1] ) )
	
	# then the same tests for frequency data
	f_fract = nbs14_f
	(taus2,adevs2,aerrs2,ns2) = allan.adev( f_fract, 1.0, taus)
	adevs = nbs14_devs[0]
	assert( check_devs( adevs2[0], adevs[0] ) )
	assert( check_devs( adevs2[1], adevs[1] ) )

	(taus2,adevs2,aerrs2,ns2) = allan.oadev( f_fract, 1.0, taus)
	oadevs = nbs14_devs[1]
	assert( check_devs( adevs2[0], oadevs[0] ) )
	assert( check_devs( adevs2[1], oadevs[1] ) )

	(taus2,adevs2,aerrs2,ns2) = allan.mdev( f_fract, 1.0, taus)
	mdevs = nbs14_devs[2]
	assert( check_devs( adevs2[0], mdevs[0] ) )
	assert( check_devs( adevs2[1], mdevs[1] ) )

	(taus2,adevs2,aerrs2,ns2) = allan.hdev( f_fract, 1.0, taus)
	hdevs = nbs14_devs[4]
	assert( check_devs( adevs2[0], hdevs[0] ) )
	assert( check_devs( adevs2[1], hdevs[1] ) )

	(taus2,adevs2,aerrs2,ns2) = allan.tdev( f_fract, 1.0, taus)
	tdevs = nbs14_devs[5]
	assert( check_devs( adevs2[0], tdevs[0] ) )
	assert( check_devs( adevs2[1], tdevs[1] ) )
	
	
	print "nbs14 test OK"


# read a simple data-file with phase or frequency numbers on each line
def read_datfile(filename):
	p=[]
	with open(filename) as f:
		for line in f:
			if line.startswith("#"):
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
				for n in range(7):
					row.append( float(l2[n]) )
				rows.append(row)
	return rows

# test a deviation function by:
# - running the function on the datafile
# - reading the correct answers from the resultfile
# - checking that tau, n, and dev are correct 
def test( function, datafile, datainterval, resultfile, verbose=0):
	phase = read_datfile(datafile)
	print "Read ",len(phase)," phase values from ", datafile
	devresults = read_resultfile(resultfile)
	print "Read ",len(devresults)," rows from ", resultfile
	taus = []
	devs = []
	ns   = []
	for row in devresults:
		tau_n = row[0] # tau in number of datapoints
		tau_s = row[1] # tau in seconds
		taus.append( tau_s )
		n = row[2] # n averages
		a = row[5]
		devs.append(a)
		ns.append(n)
		
	(taus2,devs2,errs2,ns2) = function(phase, datainterval, taus)
	for (t1,a1,n1,t2,a2,n2) in zip(taus,devs,ns, taus2,devs2,ns2):
		try:
			assert( t1 == t2 )
		except:
			print "ERROR t1=", t1, " t2=", t2
		try:
			assert( n1 == n2 )
		except:
			print "ERROR n1=", n1, " n2=",n2," at t1=", t1, " t2=", t2
		
		rel_error = (a2-a1)/a1
		tol = 1e-4
		try:
			assert( abs(rel_error) < tol )
			if verbose:
				print "OK %d %d &d \t %0.6f \t %0.6f \t %0.6f" % (t1,n1,a1,a2, rel_error)
		except:
			print "ERROR %d  %0.6f \t %0.6f \t %0.6f" % (t1,a1,a2, rel_error)
	print "test of function ",function, " Done."


if __name__ == "__main__":
	
	nbs14_1000_test()
	
	nbs14_test()
	
	data_file = 'PHASE.DAT'
	adev_result = 'phase_dat_adev.txt'
	oadev_result = 'phase_dat_oadev.txt'
	mdev_result = 'phase_dat_mdev.txt'
	tdev_result = 'phase_dat_tdev.txt'
	hdev_result = 'phase_dat_hdev.txt'
	
	verbose = 0
	test( allan.adev_phase, data_file, 1.0, adev_result , verbose)
	test( allan.oadev_phase, data_file, 1.0, oadev_result, verbose )
	test( allan.mdev_phase, data_file, 1.0, mdev_result, verbose )
	test( allan.tdev_phase, data_file, 1.0, tdev_result, verbose )
	test( allan.hdev_phase, data_file, 1.0, hdev_result, verbose )
