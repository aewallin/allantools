Getting started
===============

Installation
------------

from Github
~~~~~~~~~~~

To clone and install from github::

    $ git clone https://github.com/aewallin/allantools
    $ python setup.py install
    
This installs the latest development version of allantools. Cloning the
repository gives you examples, tests, test-data, and iPython notebooks also.

from PyPi
~~~~~~~~~

To install the latest released version of allantools from PyPi::

    $ pip install allantools
    
No examples, tests, or test-data is installed from the PyPi package.
    
Basic usage
------------

For a general description see the `API <api.html>`_.


Minimal example, phase data
~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can call allantools with only one parameter - an array of phase data.
This is suitable for time-interval measurements at 1 Hz, for example
from a time-interval-counter measuring the 1PPS output of two clocks.

::

    >>> import allantools
    >>> import pylab as plt
    >>> x = allantools.noise.white(10000)        # Generate some phase data, in seconds.
    >>> (taus, adevs, errors, ns) = allantools.oadev(x)

when only one input parameter is given, phase data in seconds is assummed
when no rate parameter is given, rate=1.0 is the default
when no taus parameter is given, taus='octave' is the default

Frequency data example
~~~~~~~~~~~~~~~~~~~~~~

Note that allantools assumes nondimensional frequency data input.
Normalization, by e.g. dividing all datapoints with the average 
frequency, is left to the user.

::

    >>> import allantools
    >>> import pylab as plt
    >>> import numpy as np
    >>> t = np.logspace(0, 3, 50)  # tau values from 1 to 1000
    >>> y = allantools.noise.white(10000)  # Generate some frequency data
    >>> r = 12.3  # sample rate in Hz of the input data
    >>> (t2, ad, ade, adn) = allantools.oadev(y, rate=r, data_type="freq", taus=t)  # Compute the overlapping ADEV
    >>> fig = plt.loglog(t2, ad) # Plot the results
    >>> # plt.show()

variations on taus parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The taus parameter can be given in a number of ways:

::

    >>> (t, d, e, n) = allantools.adev(x) # omitted, defaults to taus='octave'
    >>> (t, d, e, n) = allantools.adev(x, taus=[]) # empty list, defaults to taus='octave'
    >>> (t, d, e, n) = allantools.adev(x, taus='all') # 1, 2, 3, 4, 5, 6, 7, ...
    >>> (t, d, e, n) = allantools.adev(x, taus='octave') # 1, 2, 4, 8, 16, ...
    >>> (t, d, e, n) = allantools.adev(x, taus='decade') # 1, 2, 4, 10, 20, 40, 100, ...
    >>> my_taus=[1,5,15,28]
    >>> (t, d, e, n) = allantools.adev(x, taus=my_taus) # python list
    >>> (t, d, e, n) = allantools.adev(x, taus=np.array(my_taus)) # numpy array
    
    
Further examples are given in the ``examples`` directory of the package. For more exciting API musings you
can `read the API <api.html>`_.
