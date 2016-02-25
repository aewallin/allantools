allantools
==========

A python library for calculating Allan deviation and related time & frequency statistics. GPL v3+ license.

Developed here https://github.com/aewallin/allantools, but also available on PyPi at https://pypi.python.org/pypi/AllanTools

Input data should be evenly spaced observations of either fractional frequency,
or phase in seconds. Deviations are calculated for given tau values in seconds.

These statistics are currently included:
* ADEV    Allan deviation
* OADEV   overlapping Allan deviation,
* MDEV    modified Allan deviation,
* TDEV    Time deviation
* HDEV    Hadamard deviation
* OHDEV   overlapping Hadamard deviation
* TOTDEV  Total deviation
* MTIE    Maximum time interval error
* TIERMS  Time interval error RMS

Noise generators for creating synthetic datasets are also included:
* violet noise with f^2 PSD
* white noise with f^0 PSD
* pink noise with f^-1 PSD
* Brownian or random walk noise with f^-2 PSD 


see /tests for tests that compare allantools output to other (e.g. Stable32) programs.
More test data, benchmarks, ipython notebooks, and comparisons to known-good algorithms are welcome.

Documentation
=============
See /docs for documentation in sphinx format. On Ubuntu this requires the 'python-sphinx' and 'python-numpydoc' packages.
html documentation using sphinx can be built locally with

> make html

this generates html documentation in /docs/_build/html

The sphinx documentation is also auto-generated online

* http://allantools.readthedocs.org
* http://allantools-aewallin.readthedocs.org

IPython notebooks with examples 
=============================== 
See /examples for some examples in IPython notebook format.


github formats the notebooks into nice web-pages, for example * 
https://github.com/aewallin/allantools/blob/master/examples/noise-color-demo.ipynb


todo: add here a very short guide on how to get started with ipython


Authors 
======= 
* Anders E.E. Wallin, anders.e.e.wallin "at" gmail.com 
* Danny Price, https://github.com/telegraphic 
* Cantwell G. Carson, carsonc "at" gmail.com 
* Frédéric Meynadier, https://github.com/fmeynadier

Installation 
============


clone from github, or download from pypi.

> sudo python setup.py install

Usage 
=====

```python 
import allantools # https://github.com/aewallin/allantools/ 
rate = 1/float(data_interval) # data rate in Hz 
taus = [1,2,4,8,16] #  tau-values in seconds
# fractional frequency data
(taus_used, adev, adeverror, adev_n) = allantools.adev(fract_freqdata, rate, taus)
# phase data
(taus_used, adev, adeverror, adev_n) = allantools.adev_phase(phasedata, rate, taus)

# notes:
#  - taus_used may differ from taus, if taus has a non-integer multiples 
#  of data_interval - adeverror assumes 1/sqrt(adev_n) errors
```

Development 
===========

To do: 
* Stable32-style plots using matplotlib 
* Modified Total variance 
* Time Total (modified total variance scaled by (t^2/3) ) 
* Hadamard Total 
* Tests for different noise types according to IEEE 1139, include power-spectral-density calculations 
* Conversion between phase noise and Allan variance 
* The mtie_phase_fast approach to MTIE, using a binary tree (see BREGNI reference)

Make sure your patch does not break any of the tests, and does not 
significantly reduce the readability of the code.

Tests may be run using [py.test](http://pytest.org) (automatically finds 
tests/test_run.py) Test coverage may be obtained with the 
[coverage](https://pypi.python.org/pypi/coverage) module :

``` coverage run --source allantools setup.py test coverage report # 
Reports on standard output coverage html # Writes annotated source code 
as html in ./htmlcov/ ```

On Ubuntu this requires packages **python-pytest** and 
**python-coverage**.

References 
========== 

http://en.wikipedia.org/wiki/Allan_variance

1139-2008 - IEEE Standard Definitions of Physical Quantities for 
Fundamental Frequency and Time Metrology - Random Instabilities 
http://dx.doi.org/10.1109/IEEESTD.2008.4797525

F. Vernotte, "Variance Measurements", 2011 IFCS & EFTF
http://www.ieee-uffc.org/frequency-control/learning/pdf/Vernotte-Varience_Measurements.pdf

S. Stein, Frequency and Time - Their Measurement and Characterization. 
Precision Frequency Control Vol 2, 1985, pp 191-416. 
http://tf.boulder.nist.gov/general/pdf/666.pdf

W.J.Riley, "THE CALCULATION OF TIME DOMAIN FREQUENCY STABILITY" 
http://www.wriley.com/paper1ht.htm

Tom Van Baak http://www.leapsecond.com/tools/adev_lib.c

Fabian Czerwinski, Matlab code
http://www.mathworks.com/matlabcentral/fileexchange/26659-allan-v3-0

M. A. Hopcroft, Matlab code
http://www.mathworks.com/matlabcentral/fileexchange/26637-allanmodified

SESIA I., GALLEANI L., TAVELLA P., Application of the Dynamic Allan Variance 
for the Characterization of Space Clock Behavior, 
http://dx.doi.org/10.1109/TAES.2011.5751232
       
S. BREGNI, Fast Algorithms for TVAR and MTIE Computation in Characterization of
Network Synchronization Performance. 
http://home.deib.polimi.it/bregni/papers/cscc2001_fastalgo.pdf

David A. Howe, The total deviation approach to long-term characterization
of frequency stability, IEEE tr. UFFC vol 47 no 5 (2000)
http://dx.doi.org/10.1109/58.869040

Ilaria Sesia and Patrizia Tavella, Estimating the Allan variance in the 
presence of long periods of missing data and outliers.
2008 Metrologia 45 S134 http://dx.doi.org/10.1088/0026-1394/45/6/S19
