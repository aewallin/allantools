.. fits2hdf documentation master file, created by
   sphinx-quickstart on Fri May 22 16:29:56 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Getting started
===============

Installation
------------

To install, you first need to clone the directory from github::

    git clone https://github.com/awallin/allantools

and then run::

    python setup.py install
    
from the command line. 

Basic usage
------------

::

    import allantools
    import pylab as plt
    t = numpy.logspace(0, 3, 50)                            # tau values from 1 to 1000
    white_noise_data = allantools.noise.white(10000)        # Generate some data
    (t2, ad, ade, adn) = allantools.oadev(freqyency=y, rate=rate, taus=taus)    # Compute the overlapping ADEV
    plt.loglog(t2, ad)                                      # Plot the results
    plt.show()


Further examples are given in the ``examples`` directory of the package. For more exciting API musings you
can `read the API <api.html>`_.
