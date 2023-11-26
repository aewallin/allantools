#!/usr/bin/python
"""
  top-level test-script for allantools
  https://github.com/aewallin/allantools

  this files runs all the tests defined in the /tests/ subdirectories
  
  More datasets along with known deviations (or calculated with 
  known-good programs) are welcome.
  
  results from allantools seem correct, they agree to within 4 to 6 digits 
  of precision with other ADEV tools.
"""

import time
import sys
import pytest

if __name__ == "__main__":
    print("running all AllanTools tests. May take up to 200 seconds.")
    start = time.time()
    pytest.main()
    end = time.time()
    print("-------------------------")
    print("All tests done in %2.3f s" % (end-start)) 
    # 2014-08-31 running time without MTIE on laptop with i7-3537U CPU @ 2.00GHz
    # All tests done in 19.841 s
    # 2015-05-05 17.373 s
    # 2016-04-08: 80 tests in 174 s
    print("if we came this far without assertions or errors all is OK!")
