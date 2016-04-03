.. fits2hdf documentation master file, created by
   sphinx-quickstart on Fri May 22 16:29:56 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Getting started
===============

Installation
------------

from Github
~~~~~~~~~~~

To install, you first need to clone the directory from github::

    $ git clone https://github.com/awallin/allantools

and then run::

    $ python setup.py install
    
from the command line. 

from PyPi
~~~~~~~~~

To install the latest released versions allantools from PyPi::

    $ pip install allantools
    
Basic usage
------------

For a general description see the `API <api.html>`_.


Minimal example, phase data
~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can call allantools with only one parameter.
This is suitable for time-interval measurements at 1 Hz, for example
from a time-interval-counter measuring the 1PPS output of two clocks.

::

    import allantools
    import pylab as plt
    x = allantools.noise.white(10000)        # Generate some phase data, in seconds.
    (taus, adevs, errors, ns) = allantools.oadev(x)
    # when only one input parameter is given, phase data in seconds is assummed
    # when no rate parameter is given, rate=1.0 is the default
    # when no taus parameter is given, taus='octave' is the default

Frequency data example
~~~~~~~~~~~~~~~~~~~~~~

Note that allantools assumes nondimensional frequency data input.
Normalization, by e.g. dividing all datapoints with the average frequency, is left to the user.

::

    import allantools
    import pylab as plt
    t = numpy.logspace(0, 3, 50)             # tau values from 1 to 1000
    y = allantools.noise.white(10000)        # Generate some frequency data
    r = 12.3  # sample rate in Hz of the input data
    (t2, ad, ade, adn) = allantools.oadev(freqyency=y, rate=r, taus=t)    # Compute the overlapping ADEV
    plt.loglog(t2, ad)                                      # Plot the results
    plt.show()


Further examples are given in the ``examples`` directory of the package. For more exciting API musings you
can `read the API <api.html>`_.
