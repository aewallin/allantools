"""
  Pink frequency noise test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  the pink_frequency.txt was generated with noise.py, for documentation see that file.

"""
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
	totdev_result = 'totdev_alpha0.txt'
	mtie_result = 'mtie.txt'
	tierms_result = 'tierms.txt'
	#verbose = 1
	
	tolerance = 1e-4
	rate = 1/float(42.0) # stable32 runs were done with this data-interval
	testutils.test( allan.adev, data_file, rate, adev_result , 0, tolerance)
	testutils.test( allan.oadev, data_file, rate, oadev_result, 0, tolerance )
	testutils.test( allan.mdev, data_file, rate, mdev_result, 0, tolerance )
	testutils.test( allan.tdev, data_file, rate, tdev_result, 0, tolerance )
	testutils.test( allan.hdev, data_file, rate, hdev_result, 0, tolerance )
	testutils.test( allan.ohdev, data_file, rate, ohdev_result, 0, tolerance )
	testutils.test( allan.totdev, data_file, rate, totdev_result, 0, tolerance )

if __name__ == "__main__":
	run()


