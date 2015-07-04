.. fits2hdf documentation master file, created by
   sphinx-quickstart on Fri May 22 16:29:56 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

allantools documentation
========================

.. toctree::
   :maxdepth: 2
   
   about
   getting_started
   api
   ...
   
About allantools
----------------

A small GPL v3+ Python library for calculating Allan deviation and related statistics.

Developed here: 

    https://github.com/aewallin/allantools, 

but also available on PyPi at: 

    https://pypi.python.org/pypi/AllanTools

Input data should be evenly spaced observations of either fractional frequency,
or phase in seconds. Deviations are calculated for given tau values in seconds.