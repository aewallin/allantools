#from allantools.dataset import Dataset
#from allantools import noise
import allantools as at
import pytest
import numpy as np

"""
    unit-tests for allantools.Noise() Kasdin & Walter
    noise-generator
"""

@pytest.fixture
def noisegen():
    return at.Noise()


def test_timeseries_length(noisegen):
    """
        check that the time-series is of correct length
    """
    for n in range(2,20):
        nr = pow(2,n)
        noisegen.set_input(nr=nr)
        noisegen.generateNoise()
        print nr
        print len( noisegen.time_series )
        assert( len( noisegen.time_series ) == nr )

def test_adev(noisegen):
    qd_list = [2e-20, 3e-15, 5e-10, 6e-9, 7e-6]
    b_list  = [0, -1 , -2, -2, -2]
    for (qd,b) in zip(qd_list, b_list):
        #qd=2e-20
        #b=0
        for tau in [1,2,3,4,5]:
            noisegen.set_input(nr=pow(2,14), qd=qd , b=b)
            noisegen.generateNoise()
            (taus,devs,errs,ns)=at.adev( noisegen.time_series, taus=[tau], rate=1.0 )
            
            adev_calculated = devs[0]
            adev_predicted = noisegen.adev(tau0=1.0, tau=tau)
            #print taus,devs
            print adev_calculated, adev_predicted, adev_calculated/adev_predicted 
            assert np.isclose( adev_calculated, adev_predicted, rtol=1e-1, atol=0)

if __name__ == "__main__":
    #test_timeseries_length( at.Noise() )
    #test_adev( at.Noise() )

    pytest.main()
