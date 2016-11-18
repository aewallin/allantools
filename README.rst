AllanTools
==========

.. image:: https://badge.fury.io/py/allantools.svg
    :target: https://badge.fury.io/py/allantools
.. image:: https://travis-ci.org/aewallin/allantools.svg?branch=master
    :target: https://travis-ci.org/aewallin/allantools
.. image:: http://readthedocs.org/projects/allantools/badge/?version=latest
    :target: http://allantools.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
.. image:: https://coveralls.io/repos/github/aewallin/allantools/badge.svg?branch=master 
    :target: https://coveralls.io/github/aewallin/allantools?branch=master 

A python library for calculating Allan deviation and related 
time & frequency statistics. GPL v3+ license.

Developed at https://github.com/aewallin/allantools and also available 
on PyPi at https://pypi.python.org/pypi/AllanTools

Input data should be evenly spaced observations of either fractional frequency,
or phase in seconds. Deviations are calculated for given tau values in seconds.

These statistics are currently included:

* adev()    Allan deviation
* oadev()   overlapping Allan deviation,
* mdev()    modified Allan deviation,
* tdev()    Time deviation
* hdev()    Hadamard deviation
* ohdev()   overlapping Hadamard deviation
* totdev()  total Allan deviation
* mtie()    Maximum time interval error
* tierms()  Time interval error RMS
* mtotdev() Modified total deviation
* ttotdev() Time total deviation
* htotdev() Hadamard total deviation
* theo1()   Thêo1 deviation

Noise generators for creating synthetic datasets are also included:

* violet noise with f^2 PSD
* white noise with f^0 PSD
* pink noise with f^-1 PSD
* Brownian or random walk noise with f^-2 PSD 


see /tests for tests that compare allantools output to other 
(e.g. Stable32) programs. More test data, benchmarks, ipython notebooks, 
and comparisons to known-good algorithms are welcome!

Documentation
=============
See /docs for documentation in sphinx format. On Ubuntu this requires 
the **python-sphinx** and **python-numpydoc** packages.
html/pdf documentation using sphinx can be built locally with::

    /docs$ make html
    /docs$ make latexpdf

this generates html documentation in docs/_build/html and pdf 
documentation in docs/_build/latex.

The sphinx documentation is also auto-generated online

* http://allantools.readthedocs.org

IPython notebooks with examples 
=============================== 
See /examples for some examples in IPython notebook format.


github formats the notebooks into nice web-pages, for example 

* https://github.com/aewallin/allantools/blob/master/examples/noise-color-demo.ipynb
* https://github.com/aewallin/allantools/blob/master/examples/three-cornered-hat-demo.ipynb

todo: add here a very short guide on how to get started with ipython

Authors 
======= 
* Anders E.E. Wallin, anders.e.e.wallin "at" gmail.com 
* Danny Price, https://github.com/telegraphic 
* Cantwell G. Carson, carsonc "at" gmail.com 
* Frédéric Meynadier, https://github.com/fmeynadier

Installation 
============


clone from github, then install with::  

    sudo python setup.py install    

(see `python setup.py --help install` for install options)

or download from pypi::
    
    sudo pip install allantools


Usage 
=====

New in 2016.11 : simple top-level API, using dedicated classes for data handling and plotting.

::

    import allantools # https://github.com/aewallin/allantools/
    import numpy as np

    # Compute a deviation using the Dataset class
    a = allantools.Dataset(data=np.random.rand(1000))
    a.compute("mdev")

    # Plot it using the Plot class
    b = allantools.Plot()
    b.plot(a, errorbars=True, grid=True)
    # You can override defaults before "show" if needed
    b.ax.set_xlabel("Tau (s)")
    b.show()

Lower-level access to the algorithms is still possible :

::

    import allantools # https://github.com/aewallin/allantools/ 
    rate = 1/float(data_interval_in_s) # data rate in Hz 
    taus = [1,2,4,8,16] #  tau-values in seconds
    # fractional frequency data
    (taus_used, adev, adeverror, adev_n) = allantools.adev(fract_freqdata, data_type='freq', rate=rate, taus=taus)
    # phase data
    (taus_used, adev, adeverror, adev_n) = allantools.adev(phasedata, data_type='phase', rate=rate, taus=taus)

    # notes:
    #  - taus_used may differ from taus, if taus has a non-integer multiples 
    #  of data_interval - adeverror assumes 1/sqrt(adev_n) errors

Tests
=====

The tests compare the output of allantools to other programs such
as Stable32. Tests may be run using py.test (http://pytest.org).
Slow tests are marked 'slow' and tests failing because of a known
reason are marked 'fails'. To run all tests::
    
    $ py.test

To exclude known failing tests::

    $ py.test -m "not fails" --durations=10

To exclude tests that run slowly::

    $ py.test -m "not slow" --durations=10

To exclude both (note option change)::

    $ py.test -k "not (slow or fails)" --durations=10

To run the above command without installing the package::

    $ python setup.py test --addopts "-k 'not (fails or slow)'"

Test coverage may be obtained with the 
(https://pypi.python.org/pypi/coverage) module::

    coverage run --source allantools setup.py test --addopts "-k 'not (fails or slow)'"
    coverage report # Reports on standard output 
    coverage html # Writes annotated source code as html in ./htmlcov/

On Ubuntu this requires packages **python-pytest** and 
**python-coverage**.

Testing on multiple python versions can be done with tox (https://testrun.org/tox)

    $ tox

Tests run continuously on Travis-CI at https://travis-ci.org/aewallin/allantools



