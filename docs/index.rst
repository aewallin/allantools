.. fits2hdf documentation master file, created by
   sphinx-quickstart on Fri May 22 16:29:56 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

allantools documentation
========================

.. toctree::
   :maxdepth: 2
   
   getting_started
   api
..
   
About allantools
----------------

A python library for calculating Allan deviation and related time & frequency statistics. GPL v3+ license.

Developed here: 

    https://github.com/aewallin/allantools, 

but also available on PyPi at: 

    https://pypi.python.org/pypi/AllanTools

Input data should be evenly spaced observations of either fractional frequency, or phase in seconds. Deviations are calculated for given tau values in seconds.

These statistics are currently included:

* ADEV Allan deviation
* OADEV overlapping Allan deviation
* GRADEV Gap-resistant overlapping allan deviation
* MDEV modified Allan deviation
* TDEV Time deviation
* HDEV Hadamard deviation
* OHDEV overlapping Hadamard deviation
* TOTDEV Total deviation
* MTIE Maximum time interval error
* TIERMS Time interval error RMS

see ``/tests`` for tests that compare allantools output to other (e.g. Stable32) programs.
More test data, benchmarks, and comparisons to known-good algorithms are welcome.

Authors
-------

* Anders E.E. Wallin, anders.e.e.wallin "at" gmail.com
* Danny Price, https://github.com/telegraphic
* Cantwell G. Carson, carsonc "at" gmail.com

References
----------

http://en.wikipedia.org/wiki/Allan_variance

1139-2008 - IEEE Standard Definitions of Physical Quantities for Fundamental Frequency and Time Metrology - Random Instabilities http://dx.doi.org/10.1109/IEEESTD.2008.4797525

S. Stein, Frequency and Time - Their Measurement and Characterization. Precision Frequency Control Vol 2, 1985, pp 191-416. http://tf.boulder.nist.gov/general/pdf/666.pdf

W.J.Riley, "THE CALCULATION OF TIME DOMAIN FREQUENCY STABILITY" http://www.wriley.com/paper1ht.htm

Tom Van Baak http://www.leapsecond.com/tools/adev_lib.c

Fabian Czerwinski, Matlab code http://www.mathworks.com/matlabcentral/fileexchange/26659-allan-v3-0

M. A. Hopcroft, Matlab code http://www.mathworks.com/matlabcentral/fileexchange/26637-allanmodified

SESIA I., GALLEANI L., TAVELLA P., Application of the Dynamic Allan Variance for the Characterization of Space Clock Behavior, http://dx.doi.org/10.1109/TAES.2011.5751232

S. BREGNI, Fast Algorithms for TVAR and MTIE Computation in Characterization of Network Synchronization Performance. http://home.deib.polimi.it/bregni/papers/cscc2001_fastalgo.pdf

David A. Howe, The total deviation approach to long-term characterization of frequency stability, IEEE tr. UFFC vol 47 no 5 (2000) http://dx.doi.org/10.1109/58.869040

Ilaria Sesia and Patrizia Tavella, Estimating the Allan variance in the presence of long periods of missing data and outliers. 2008 Metrologia 45 S134 http://dx.doi.org/10.1088/0026-1394/45/6/S19
