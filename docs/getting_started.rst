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

To use ``allantools`` ... to ... ::

    from fits2hdf.io.fitsio import read_fits
    from fits2hdf.io.hdfio import export_hdf
    a = read_fits('my_file.fits')
    export_hdf(a, 'my_file.hdf')

and to convert the other way::

    from fits2hdf.io.fitsio import export_fits
    from fits2hdf.io.hdfio import read_hdf
    a = read_hdf('my_file.hdf')
    export_hdf(a, 'my_file.fits')

For more exciting API musings you can `read the API <api.html>`_.
