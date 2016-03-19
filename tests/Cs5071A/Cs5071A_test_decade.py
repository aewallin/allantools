"""
  Test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  The 5071A_phase.txt is a dataset collected with a time-interval-counter
  between 1 pulse-per-second outputs from a 5071A Cs clock against a H-maser
 
  first datapoint          1391174210 2014-01-31 13:16:50 UTC 
  last datapoint           1391731199 2014-02-06 23:59:59 UTC
  556990 datapoints in total

  This test uses log-spaced tau-values (i.e. 1, 2, 4, 10, etc.)
  The CS5071A_test_all.py test much more tau-values (1,2,3,4, etc.) but is slower.

  AW2014-02-07
"""



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
	print(" %.2f s"% ( end-start ))
	return time.clock()
	
def run():
	# hack to run script from its own directory
	abspath = os.path.abspath(__file__)
	dname = os.path.dirname(abspath)
	os.chdir(dname)

	data_file = '5071A_phase.txt.gz' # input data file
    
	adev_result = 'adev_decade.txt'
	oadev_result = 'oadev_decade.txt'
	mdev_result = 'mdev_decade.txt'
	tdev_result = 'tdev_decade.txt'
	hdev_result = 'hdev_decade.txt'
	ohdev_result = 'ohdev_decade.txt'
	totdev_result = 'totdev_decade.txt'
	mtie_result = 'mtie_fast.txt'
	tierms_result = 'tierms_decade.txt'
	
	print("runtime ~15 seconds 2016-03-12 on an i7 CPU")
    
	verbose = 1
	tolerance = 1e-4
	rate = 1/float(1.0) # stable32 runs were done with this data-interval
	start0 = time.clock()
	start = print_elapsed(time.clock())
	
	testutils.test_row_by_row( allan.adev, data_file, rate, adev_result , verbose, tolerance) # 0.6 s
	start = print_elapsed(start)
	testutils.test_row_by_row( allan.oadev, data_file, rate, oadev_result, verbose, tolerance ) # 3.8 s
	start = print_elapsed(start)
	testutils.test_row_by_row( allan.mdev, data_file, rate, mdev_result, verbose, tolerance ) # 5.5 s
	start = print_elapsed(start)
	testutils.test_row_by_row( allan.tdev, data_file, rate, tdev_result, verbose, tolerance ) # 5.5 s
	start = print_elapsed(start)
	testutils.test_row_by_row( allan.hdev, data_file, rate, hdev_result, verbose, tolerance ) # 0.85 s
	start = print_elapsed(start)
	testutils.test_row_by_row( allan.ohdev, data_file, rate, ohdev_result, verbose, tolerance ) # 5 s
	start = print_elapsed(start)
	testutils.test_row_by_row( allan.tierms, data_file, rate, tierms_result, verbose, tolerance ) # 5.9 s
	start = print_elapsed(start)
	testutils.test_row_by_row( allan.totdev, data_file, rate, totdev_result, verbose, tolerance ) # 13 s
	start = print_elapsed(start)
	#testutils.test_row_by_row( allan.mtie, data_file, rate, mtie_result, verbose, tolerance ) # 13 s
	#start = print_elapsed(start)

	print(" Cs5071A_decade tests took %.2f s"% ( time.clock()-start0 )) 
    # 2014-08-31 running time without MTIE
    # Laptop: i7-3537U CPU @ 2.00GHz
	# Cs5071A_decade tests took 20.34 s
if __name__ == "__main__":
	run()


