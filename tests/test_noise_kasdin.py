#from allantools.dataset import Dataset
#from allantools import noise
import allantools as at
import pytest


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


if __name__ == "__main__":
    #test_timeseries_length( at.Noise() )
    pytest.main()
