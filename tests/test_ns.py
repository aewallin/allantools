#!/usr/bin/python

import sys
sys.path.append("..")

import allantools as allan
from allantools import noise
import numpy

def _test( function, data, rate, taus):

	(taus2,devs2,errs2,ns2) = function(data, rate, taus)
	assert( len(taus2) == len(devs2) )
	assert( len(taus2) == len(errs2) )
	assert( len(taus2) == len(ns2) )
	for n in ns2:
		if n<= 1:
			print("test of ", function, " failed: ", n)
		assert( n > 1 ) # n should be 2 or more for each tau
	print("test_ns of function ",function, " OK.")
	
def run():
	
	# this test asks for results at unreasonable tau-values
	# either zero, not an integer multiple of the data-interval
	# or too large, given the length of the dataset
	N = 500
	rate = 1.0
	phase_white = noise.white(N)
	taus_try = numpy.logspace(0,4,4000) # try insane tau values
	_test( allan.adev_phase, phase_white, rate, taus_try)
	_test( allan.oadev_phase, phase_white, rate, taus_try)
	_test( allan.mdev_phase, phase_white, rate, taus_try)
	_test( allan.tdev_phase, phase_white, rate, taus_try)
	_test( allan.hdev_phase, phase_white, rate, taus_try)
	_test( allan.ohdev_phase, phase_white, rate, taus_try)
	_test( allan.totdev_phase, phase_white, rate, taus_try)
	_test( allan.mtie_phase, phase_white, rate, taus_try)
	_test( allan.tierms_phase, phase_white, rate, taus_try)

if __name__ == "__main__":
	run()
