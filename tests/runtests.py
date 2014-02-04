

# this files runs all the tests defined in the subdirectories

import time
import sys
sys.path.append("..") # allows importing allantools from parent dir


from nbs14 import nbs14_test
from phasedat import phase_dat_test
from pink_frequency import pink

if __name__ == "__main__":
	
	start = time.clock()
	
	nbs14_test.run()
	phase_dat_test.run()
	pink.run()
	
	end = time.clock()
	print "Tests done in %2.3f s" % (end-start)
