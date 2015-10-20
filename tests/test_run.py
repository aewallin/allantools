"""
  top-level test-script for allantools - for use with py.test framework
  https://github.com/aewallin/allantools

  this files runs all the tests defined in the /tests/ subdirectories
  
  More datasets along with known deviations (or calculated with 
  known-good programs) are welcome.
  
  results from allantools seem correct, they agree to within 4 to 6 digits 
  of precision with other ADEV tools.
"""

from nbs14 import nbs14_test
from phasedat import phase_dat_test
from pink_frequency import pink
from Cs5071A import Cs5071A_test_decade
from Keysight53230A_ti_noise_floor import TIC_test
from ocxo import ocxo_test
from test_ns import run as test_ns_run

def test_ocxo():
    # high-stability OCXO timebase on HP instrument
    ocxo_test.run()

def test_TIC(): 
    # 53230A counter noise floor dataset
    TIC_test.run()

def test_nbs14():
    # NBS14 test data with published deviations
    nbs14_test.run()

def test_phase_dat():
    # phase.dat from Stable32
    phase_dat_test.run() 

def test_pink():
    # synthetic pink frequency noise
    pink.run() 

def test_Cs5971A():
    # HP 5071A Cs-clock measured against H-maser
    Cs5071A_test_decade.run() 

def test_ns():
    # sanity-checks for tau values
    test_ns_run() 
