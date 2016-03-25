"""
  PHASE.DAT test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  PHASE.DAT comes with Stable32 (version 1.53 was used in this case)

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
    print(dname)

    data_file = 'PHASE.DAT'
    htotdev_result = 'phase_dat_htotdev_octave.txt'
    
    ttotdev_result1 = 'phase_dat_ttotdev_octave.txt'
    ttotdev_result2 = 'phase_dat_ttotdev_decade.txt'
    ttotdev_result3 = 'phase_dat_ttotdev_all.txt'
        
    mtotdev_result1 = 'phase_dat_mtotdev_decade.txt'
    mtotdev_result2 = 'phase_dat_mtotdev_octave.txt'
    mtotdev_result3 = 'phase_dat_mtotdev.txt'

    adev_result = 'phase_dat_adev.txt'
    oadev_result = 'phase_dat_oadev.txt'
    mdev_result = 'phase_dat_mdev.txt'
    tdev_result = 'phase_dat_tdev.txt'
    hdev_result = 'phase_dat_hdev.txt'
    ohdev_result = 'phase_dat_ohdev.txt'
    totdev_result = 'phase_dat_totdev.txt'
    mtie_result = 'phase_dat_mtie.txt'
    tierms_result = 'phase_dat_tierms.txt'
    verbose = 1
    
    tolerance = 1e-4
    testutils.test_row_by_row( allan.htotdev, data_file, 1.0, htotdev_result , verbose=verbose, tolerance=tolerance)
    
    testutils.test_row_by_row( allan.ttotdev, data_file, 1.0, ttotdev_result2 , verbose=verbose, tolerance=tolerance)
    testutils.test_row_by_row( allan.ttotdev, data_file, 1.0, ttotdev_result1 , verbose=verbose, tolerance=tolerance)

    testutils.test_row_by_row( allan.mtotdev, data_file, 1.0, mtotdev_result1 , verbose=verbose, tolerance=tolerance)
    testutils.test_row_by_row( allan.mtotdev, data_file, 1.0, mtotdev_result2 , verbose=verbose, tolerance=tolerance)
    testutils.test( allan.adev, data_file, 1.0, adev_result , verbose=verbose, tolerance=tolerance)
    testutils.test( allan.oadev, data_file, 1.0, oadev_result, verbose=verbose, tolerance=tolerance )
    testutils.test( allan.mdev, data_file, 1.0, mdev_result, verbose=verbose, tolerance=tolerance )
    testutils.test( allan.tdev, data_file, 1.0, tdev_result, verbose=verbose, tolerance=tolerance )
    testutils.test( allan.hdev, data_file, 1.0, hdev_result, verbose=verbose, tolerance=tolerance )
    testutils.test( allan.ohdev, data_file, 1.0, ohdev_result, verbose=verbose, tolerance=tolerance )
    testutils.test( allan.totdev, data_file, 1.0, totdev_result, verbose=verbose, tolerance=tolerance )
    testutils.test( allan.mtie, data_file, 1.0, mtie_result, verbose=verbose, tolerance=tolerance )
    testutils.test( allan.tierms, data_file, 1.0, tierms_result, verbose=verbose, tolerance=tolerance )
    

if __name__ == "__main__":
    run()


