Development 
===========

To-do list
----------

Here follows an un-ordered to do list:

* Statistics and core algorithms

    * The mtie_phase_fast approach to MTIE, using a binary tree (see BREGNI reference)
    * TheoH
    * Confidence intervals based on identified noise-type and equivalent degrees of freedom.
    * Bias corrections for biased statistics (totdev, mtotdev, htotdev, theo1)
    
* Improve documentation
* Improve packaging for PyPi and/or other packaging systems (PPA for Ubuntu/Debian?)
* Stable32-style plots using matplotlib 
* Tests for different noise types according to IEEE 1139, include power-spectral-density calculations 
* Conversion between phase noise and Allan variance 
* Phase noise calculations and plots
* Comparison to other libraries such as GPSTk

Make sure your patch does not break any of the tests, and does not 
significantly reduce the readability of the code.

Notes for Pypi
--------------

Creating a source distribution

    python setup.py sdist

Testing the source distribution. The install takes a long time while 
compiling numpy and scipy.

::

    $ virtualenv tmp
    $ tmp/bin/pip install dist/AllanTools-2016.2.tar.gz 
    $ tmp/bin/python
    >>> import allantools

Registering, uploading and testing  source distribution to PyPi test server
(requries a ~/.pypirc with username and password)

::

    $ python setup.py register -r test
    $ python setup.py sdist upload -r test
    $ pip install -i https://testpypi.python.org/pypi AllanTools

Registering and uploading to PyPi

::

    $ python setup.py register
    $ python setup.py sdist upload

