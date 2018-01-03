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
    """
        check that time-series has the ADEV that we expect
    """
    qd_list = [2e-20, 3e-15, 5e-10, 6e-9, 7e-6]
    b_list  = [0, -1, -2 , -3, -4]
    tau_list = [1,2,3,4,5, 20, 30]
    for b in b_list:
        for qd in qd_list:
            #qd=2e-20
            #b=0
            for tau in tau_list:
                noisegen.set_input(nr=pow(2,14), qd=qd , b=b)
                noisegen.generateNoise()
                (taus,devs,errs,ns)=at.adev( noisegen.time_series, taus=[tau], rate=1.0 )
                
                adev_calculated = devs[0]
                adev_predicted = noisegen.adev(tau0=1.0, tau=tau)
                #print taus,devs
                print b, tau, qd, adev_calculated, adev_predicted, adev_calculated/adev_predicted 
                assert np.isclose( adev_calculated, adev_predicted, rtol=3e-1, atol=0)
                # NOTE high relative tolarence here !!

def test_mdev(noisegen):
    """
        check that time-series has the MDEV that we expect
    """
    qd_list = [2e-20, 3e-15, 5e-10, 6e-9, 7e-6]
    b_list  = [0, -1, -3, -4, -2 ]
    tau_list = [1,2,3,4,5, 20, 30]
    for b in b_list:
        for qd in qd_list:
            #qd=2e-20
            #b=0
            for tau in tau_list:
                noisegen.set_input(nr=pow(2,14), qd=qd , b=b)
                noisegen.generateNoise()
                (taus,devs,errs,ns)=at.mdev( noisegen.time_series, taus=[tau], rate=1.0 )
                
                mdev_calculated = devs[0]
                mdev_predicted = noisegen.mdev(tau0=1.0, tau=tau)
                #print taus,devs
                print b, tau, qd, mdev_calculated, mdev_predicted, mdev_calculated/mdev_predicted 
                assert np.isclose( mdev_calculated, mdev_predicted, rtol=3e-1, atol=0)
                # NOTE high relative tolarence here !!
                
if __name__ == "__main__":
    #test_timeseries_length( at.Noise() )
    #test_adev( at.Noise() )
    #test_mdev( at.Noise() )
    pytest.main()
