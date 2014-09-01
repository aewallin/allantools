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
	print dname

	data_file = 'PHASE.DAT'
	adev_result = 'phase_dat_adev.txt'
	oadev_result = 'phase_dat_oadev.txt'
	mdev_result = 'phase_dat_mdev.txt'
	tdev_result = 'phase_dat_tdev.txt'
	hdev_result = 'phase_dat_hdev.txt'
	ohdev_result = 'phase_dat_ohdev.txt'
	totdev_result = 'phase_dat_totdev.txt'
	mtie_result = 'phase_dat_mtie.txt'
	tierms_result = 'phase_dat_tierms.txt'
	verbose = 0
	
	tolerance = 1e-4
	testutils.test( allan.adev_phase, data_file, 1.0, adev_result , verbose, tolerance)
	testutils.test( allan.oadev_phase, data_file, 1.0, oadev_result, verbose, tolerance )
	testutils.test( allan.mdev_phase, data_file, 1.0, mdev_result, verbose, tolerance )
	testutils.test( allan.tdev_phase, data_file, 1.0, tdev_result, verbose, tolerance )
	testutils.test( allan.hdev_phase, data_file, 1.0, hdev_result, verbose, tolerance )
	testutils.test( allan.ohdev_phase, data_file, 1.0, ohdev_result, verbose, tolerance )
	testutils.test( allan.totdev_phase, data_file, 1.0, totdev_result, verbose, tolerance )
	testutils.test( allan.mtie_phase, data_file, 1.0, mtie_result, verbose, tolerance )
	testutils.test( allan.tierms_phase, data_file, 1.0, tierms_result, verbose, tolerance )
	

if __name__ == "__main__":
	run()


