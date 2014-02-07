# phase in s for a 5071A Cs clock against a H-maser
# 1PPS time-interval counter measurements
# AW2014-02-07
#
# first datapoint          1391174210 2014-01-31 13:16:50 UTC 
# last datapoint           1391731199 2014-02-06 23:59:59 UTC
# 556990 datapoints in total

import math
import sys
sys.path.append("..")
sys.path.append("../..") # hack to import from parent directory
# remove if you have allantools installed in your python path

import allantools as allan
import testutils

import os
import time

def print_elapsed(start):
	end = time.clock()
	print " %.2f s"% ( end-start )
	return time.clock()
	
def run():
	# hack to run script from its own directory
	abspath = os.path.abspath(__file__)
	dname = os.path.dirname(abspath)
	os.chdir(dname)

	data_file = '5071A_phase.txt'
	adev_result = 'adev_decade.txt'
	oadev_result = 'oadev_decade.txt'
	mdev_result = 'mdev_decade.txt'
	tdev_result = 'tdev_decade.txt'
	hdev_result = 'hdev_decade.txt'
	ohdev_result = 'ohdev_decade.txt'
	totdev_result = 'totdev_decade.txt'
	mtie_result = 'mtie_fast.txt'
	tierms_result = 'tierms_decade.txt'
	
	verbose = 1
	
	tolerance = 1e-4
	rate = 1/float(1.0) # stable32 runs were done with this data-interval
	start0 = time.clock()
	start = print_elapsed(time.clock())
	
	testutils.test_row_by_row( allan.adev_phase, data_file, rate, adev_result , verbose, tolerance) # 0.6 s
	start = print_elapsed(start)
	testutils.test_row_by_row( allan.oadev_phase, data_file, rate, oadev_result, verbose, tolerance ) # 3.8 s
	start = print_elapsed(start)
	testutils.test_row_by_row( allan.mdev_phase, data_file, rate, mdev_result, verbose, tolerance ) # 5.5 s
	start = print_elapsed(start)
	testutils.test_row_by_row( allan.tdev_phase, data_file, rate, tdev_result, verbose, tolerance ) # 5.5 s
	start = print_elapsed(start)
	testutils.test_row_by_row( allan.hdev_phase, data_file, rate, hdev_result, verbose, tolerance ) # 0.85 s
	start = print_elapsed(start)
	testutils.test_row_by_row( allan.ohdev_phase, data_file, rate, ohdev_result, verbose, tolerance ) # 5 s
	start = print_elapsed(start)
	testutils.test_row_by_row( allan.tierms_phase, data_file, rate, tierms_result, verbose, tolerance ) # 5.9 s
	start = print_elapsed(start)
	testutils.test_row_by_row( allan.totdev_phase, data_file, rate, totdev_result, verbose, tolerance ) # 13 s
	start = print_elapsed(start)
	testutils.test_row_by_row( allan.mtie_phase, data_file, rate, mtie_result, verbose, tolerance ) # 13 s
	start = print_elapsed(start)

	print " all tests took %.2f s"% ( time.clock()-start0 ) # i7 CPU: ca 40s (without MTIE)
	
if __name__ == "__main__":
	run()


