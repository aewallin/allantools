"""
  Test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  AW2015-06-26
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
    print " %.2f s"% ( end-start )
    return time.clock()
    
def run():
    # hack to run script from its own directory
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)

    data_file = 'ocxo_frequency.txt'
    
    adev_result = 'stable32_adev_alltau.txt'
    oadev_result = 'stable32_oadev_alltau.txt'
    mdev_result = 'stable32_mdev_alltau.txt'
    tdev_result = 'stable32_tdev_alltau.txt'
    hdev_result = 'stable32_hdev_alltau.txt'
    ohdev_result = 'stable32_ohdev_alltau.txt'
    totdev_result = 'stable32_totdev_alltau.txt'
    
    verbose = 1
    tolerance = 1e-4 # relative tolerance
    rate = 1/float(1.0) # stable32 runs were done with this data-interval

    start0 = time.clock()
    start = time.clock()
    
    testutils.test_row_by_row( allan.adev, data_file, rate, adev_result , verbose, tolerance, normalize=True)
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.oadev, data_file, rate, oadev_result, verbose, tolerance, normalize=True )
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.mdev, data_file, rate, mdev_result, verbose, tolerance, normalize=True)
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.tdev, data_file, rate, tdev_result, verbose, tolerance, normalize=True )
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.hdev, data_file, rate, hdev_result, verbose, tolerance, normalize=True )
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.ohdev, data_file, rate, ohdev_result, verbose, tolerance, normalize=True  )
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.totdev, data_file, rate, totdev_result, verbose, tolerance, normalize=True  )
    start = print_elapsed(start)
    
    print " OCXO tests took %.2f s" % ( time.clock()-start0 ) 
    
if __name__ == "__main__":
    run()


