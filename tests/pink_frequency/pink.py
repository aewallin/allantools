import math
import sys
sys.path.append("..")
sys.path.append("../..") # hack to import from parent directory
# remove if you have allantools installed in your python path

import allantools as allan
import testutils

import os

def run():
	# hack to run script from its own directory
	abspath = os.path.abspath(__file__)
	dname = os.path.dirname(abspath)
	os.chdir(dname)

	data_file = 'pink_frequency.txt'
	adev_result = 'adev.txt'
	oadev_result = 'oadev.txt'
	mdev_result = 'mdev.txt'
	tdev_result = 'tdev.txt'
	hdev_result = 'hdev.txt'
	ohdev_result = 'ohdev.txt'
	totdev_result = 'totdev.txt'
	mtie_result = 'mtie.txt'
	tierms_result = 'tierms.txt'
	verbose = 0
	
	tolerance = 1e-4
	rate = 1/float(42.0) # stable32 runs were done with this data-interval
	testutils.test( allan.adev, data_file, rate, adev_result , verbose, tolerance)
	#testutils.test( allan.oadev_phase, data_file, 1.0, oadev_result, verbose, tolerance )
	#testutils.test( allan.mdev_phase, data_file, 1.0, mdev_result, verbose, tolerance )
	#testutils.test( allan.tdev_phase, data_file, 1.0, tdev_result, verbose, tolerance )
	#testutils.test( allan.hdev_phase, data_file, 1.0, hdev_result, verbose, tolerance )
	#testutils.test( allan.ohdev_phase, data_file, 1.0, ohdev_result, verbose, tolerance )
	#testutils.test( allan.totdev_phase, data_file, 1.0, totdev_result, verbose, tolerance )
	#testutils.test( allan.mtie_phase, data_file, 1.0, mtie_result, verbose, tolerance )
	#testutils.test( allan.tierms_phase, data_file, 1.0, tierms_result, verbose, tolerance )
	

if __name__ == "__main__":
	run()


