

AllanTools
==========

.. image:: https://badge.fury.io/py/AllanTools.svg
    :target: https://badge.fury.io/py/AllanTools 
.. image:: https://img.shields.io/conda/vn/conda-forge/allantools.svg
    :target: https://anaconda.org/conda-forge/allantools
.. image:: https://img.shields.io/conda/dn/conda-forge/allantools.svg
    :target: https://anaconda.org/conda-forge/allantools
.. image:: https://github.com/aewallin/allantools/actions/workflows/python-pytest.yml/badge.svg
    :target: https://github.com/aewallin/allantools/actions/workflows/python-pytest.yml
.. image:: https://github.com/aewallin/allantools/actions/workflows/python-flake8.yml/badge.svg
    :target: https://github.com/aewallin/allantools/actions/workflows/python-flake8.yml
    :alt: flake8 Status
.. image:: https://readthedocs.org/projects/allantools/badge/?version=latest
    :target: https://allantools.readthedocs.io/en/latest/
.. image:: https://coveralls.io/repos/github/aewallin/allantools/badge.svg
    :target: https://coveralls.io/github/aewallin/allantools

A python library for calculating Allan deviation and related 
time & frequency statistics. `LGPL v3+ license <https://www.gnu.org/licenses/lgpl.html>`_.

* Development at https://github.com/aewallin/allantools
* Installation package at https://pypi.python.org/pypi/AllanTools
* Discussion group at https://groups.google.com/d/forum/allantools
* Documentation available at https://allantools.readthedocs.org
 

Input data should be evenly spaced observations of either fractional frequency,
or phase in seconds. Deviations are calculated for given tau values in seconds.

=====================================   ====================================================   ====================================================
Function                                Description                                            Comment
=====================================   ====================================================   ====================================================
``adev()``                              Allan deviation                                        Classic - use only if required - relatively poor confidence.
``oadev()``                             Overlapping Allan deviation                            General purpose - most widely used - first choice
``mdev()``                              Modified Allan deviation                               Used to distinguish between White and Flicker Phase Modulation.
``tdev()``                              Time deviation                                         Based on modified Allan variance.
``hdev()``                              Hadamard deviation                                     Rejects frequency drift, and handles divergent noise.
``ohdev()``                             Overlapping Hadamard deviation                         Better confidence than normal Hadamard.
``pdev()``                              Parabolic deviation                                    Estimate uncertainty of Omega-counter data
``totdev()``                            Total deviation                                        Better confidence at long averages for Allan deviation.
``mtotdev()``                           Modified total deviation                               Modified Total deviation. Better confidence at long averages for modified Allan
``ttotdev()``                           Time total deviation
``htotdev()``                           Hadamard total deviation
``theo1()``                             Theo1 deviation                                        Theo1 is a two-sample variance with improved confidence and extended averaging factor range.
``mtie()``                              Maximum Time Interval Error
``tierms()``                            Time Interval Error RMS
``gradev()``                            Gap resistant overlapping Allan deviation
``gcodev()``                            Groslambert Covariance                                 Improved three-corner-hat analysis
=====================================   ====================================================   ====================================================

Noise generators for creating synthetic datasets are also included:

* violet noise with f^2 PSD
* white noise with f^0 PSD
* pink noise with f^-1 PSD
* Brownian or random walk noise with f^-2 PSD 

More details on available statistics and noise generators : `full list of available functions <https://allantools.readthedocs.io/en/latest/functions.html>`_  

see /tests for tests that compare allantools output to other 
(e.g. Stable32) programs. More test data, benchmarks, ipython notebooks, 
and comparisons to known-good algorithms are welcome!

Installation 
------------


Install from pypi::
    
    pip install allantools

Latest version + examples, tests, test data, iPython notebooks : clone from github, then install ::  

    python setup.py install

(see `python setup.py --help install` for install options)

These commands should be run as root for system-wide installation, or 
you can use the `--user` option to install for your account only. 
Exact command names may vary depending on your OS / package manager / target python version.

Basic usage 
-----------

Minimal example, phase data
~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can call allantools with only one parameter - an array of phase data.
This is suitable for time-interval measurements at 1 Hz, for example
from a time-interval-counter measuring the 1PPS output of two clocks.

::

    >>> import allantools
    >>> x = allantools.noise.white(10000)        # Generate some phase data, in seconds.
    >>> (taus, adevs, errors, ns) = allantools.oadev(x)

when only one input parameter is given, phase data in seconds is assumed
when no rate parameter is given, rate=1.0 is the default
when no taus parameter is given, taus='octave' is the default

Frequency data example
~~~~~~~~~~~~~~~~~~~~~~

Note that allantools assumes non-dimensional frequency data input.
Normalization, by e.g. dividing all data points with the average 
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


*New in 2016.11* : simple top-level `API <api.html>`_, using dedicated classes for data handling and plotting.

::

    import allantools # https://github.com/aewallin/allantools/
    import numpy as np

    # Compute a deviation using the Dataset class
    a = allantools.Dataset(data=np.random.rand(1000))
    a.compute("mdev")

    # New in 2019.7 : write results to file
    a.write_results("output.dat")

    # Plot it using the Plot class
    b = allantools.Plot()
    # New in 2019.7 : additional keyword arguments are passed to 
    # matplotlib.pyplot.plot()
    b.plot(a, errorbars=True, grid=True)
    # You can override defaults before "show" if needed
    b.ax.set_xlabel("Tau (s)")
    b.show()


Jupyter notebooks with examples 
-------------------------------

Jupyter notebooks are interactive python scripts, embedded in a browser, 
allowing you to manipulate data and display plots like easily. For guidance 
on installing jupyter, please refer to https://jupyter.org/install.

See /examples for some examples in notebook format.

github formats the notebooks into nice web-pages, for example 

* https://github.com/aewallin/allantools/blob/master/examples/noise-color-demo.ipynb
* https://github.com/aewallin/allantools/blob/master/examples/three-cornered-hat-demo.ipynb


Authors 
-------
* Anders E.E. Wallin, anders.e.e.wallin "at" gmail.com , https://github.com/aewallin
* Danny Price, https://github.com/telegraphic 
* Cantwell G. Carson, carsonc "at" gmail.com 
* Frédéric Meynadier, https://github.com/fmeynadier
* Yan Xie, https://github.com/yxie-git
* Erik Benkler, https://github.com/EBenkler
