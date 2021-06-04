"""
  Test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  The 5071A_phase.txt is a dataset collected with a time-interval-counter
  between 1 pulse-per-second outputs from a 5071A Cs clock against a H-maser

  first datapoint          1391174210 2014-01-31 13:16:50 UTC
  last datapoint           1391731199 2014-02-06 23:59:59 UTC
  556990 datapoints in total

  This test is quite slow, for a fater test see Cs5071A_test_decade.py

  AW2014-02-07
"""

import os
import time
import testutils
import allantools as allan
# import math
import sys
sys.path.append("..")
sys.path.append("../..")  # hack to import from parent directory
# remove if you have allantools installed in your python path


def print_elapsed(label, start, start0):
    end = time.clock()
    print("%s test done in %.2f s, elapsed= %.2f min" %
          (label, end-start, (end-start0)/60))
    return time.clock()


def run():
    # hack to run script from its own directory
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)

    data_file = '5071A_phase.txt.gz'  # input data

    adev_result = 'adev_all.txt'
    oadev_result = 'oadev_all.txt'
    mdev_result = 'mdev_all.txt'
    tdev_result = 'tdev_all.txt'
    hdev_result = 'hdev_all.txt'
    ohdev_result = 'ohdev_all.txt'
    totdev_result = 'totdev_all.txt'
    mtie_result = 'mtie_all.txt'
    tierms_result = 'tierms_all.txt'
    verbose = 1

    tolerance = 1e-4
    rate = 1/float(1.0)  # stable32 runs were done with this data-interval

    print("WARNING: this test takes a long time to run!!")
    print(" ~12 minutes as of 2014-02-07 on an i7 CPU")
    print("")
    verbose = 0
    start0 = time.clock()
    start = print_elapsed('Starting..', time.clock(), start0)

    testutils.test_row_by_row(allan.adev, data_file,
                              rate, adev_result, verbose, tolerance)   # 1.34 s
    start = print_elapsed('ADEV', start, start0)

    testutils.test_row_by_row(allan.hdev, data_file,
                              rate, hdev_result, verbose, tolerance)  # 1.9 s
    start = print_elapsed('HDEV', start, start0)

    testutils.test_row_by_row(allan.mtie, data_file,
                              rate, mtie_result, verbose, tolerance)  # 13 s
    start = print_elapsed('MTIE', start, start0)

    testutils.test_row_by_row(allan.oadev, data_file,
                              rate, oadev_result, verbose, tolerance)  # 63 s
    start = print_elapsed('OADEV', start, start0)

    testutils.test_row_by_row(allan.ohdev, data_file,
                              rate, ohdev_result, verbose, tolerance)  # 88 s
    start = print_elapsed('OHDEV', start, start0)

    testutils.test_row_by_row(allan.mdev, data_file,
                              rate, mdev_result, verbose, tolerance)  # 98 s
    start = print_elapsed('MDEV', start, start0)

    testutils.test_row_by_row(allan.tdev, data_file,
                              rate, tdev_result, verbose, tolerance)  # 99 s
    start = print_elapsed('TDEV', start, start0)

    testutils.test_row_by_row(allan.tierms, data_file,
                              rate, tierms_result, verbose, tolerance)  # 117 s
    start = print_elapsed('TIERMS', start, start0)

    testutils.test_row_by_row(allan.totdev, data_file,
                              rate, totdev_result, verbose, tolerance)  # 245 s
    start = print_elapsed('TOTDEV', start, start0)

    end = time.clock()
    # total time 470s on i7 CPU (without mtie!)
    print("Tests done in %2.3f s" % (end-start0))


if __name__ == "__main__":
    run()
