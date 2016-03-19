"""
  Test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  AW2015-03-29
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

    data_file = 'tic_phase.txt'
    adev_result = 'tic_adev.txt'
    oadev_result = 'tic_oadev.txt'
    mdev_result = 'tic_mdev.txt'
    tdev_result = 'tic_tdev.txt'
    hdev_result = 'tic_hdev.txt'
    ohdev_result = 'tic_ohdev.txt'
    totdev_result = 'tic_totdev.txt'
    #mtie_result = 'mtie_fast.txt'
    tierms_result = 'tic_tierms.txt'
    
    verbose = 1
    
    tolerance = 1e-4
    rate = 1/float(1.0) # stable32 runs were done with this data-interval
    start0 = time.clock()
    start = print_elapsed(time.clock())
    
    testutils.test_row_by_row( allan.adev, data_file, rate, adev_result , verbose, tolerance)
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.oadev, data_file, rate, oadev_result, verbose, tolerance )
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.mdev, data_file, rate, mdev_result, verbose, tolerance )
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.tdev, data_file, rate, tdev_result, verbose, tolerance )
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.hdev, data_file, rate, hdev_result, verbose, tolerance ) # 0.85 s
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.ohdev, data_file, rate, ohdev_result, verbose, tolerance ) # 5 s
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.tierms, data_file, rate, tierms_result, verbose, tolerance ) # 5.9 s
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.totdev, data_file, rate, totdev_result, verbose, tolerance ) # 13 s
    start = print_elapsed(start)
    
    #testutils.test_row_by_row( allan.mtie_phase, data_file, rate, mtie_result, verbose, tolerance ) # 13 s
    #start = print_elapsed(start)

    print(" TIC tests took %.2f s" % ( time.clock()-start0 )) 
    # 2015-03-29 running time without MTIE
    # Laptop: i7-3537U CPU @ 2.00GHz
    # tests took 2.38 s
    
if __name__ == "__main__":
    run()


