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
    * Multi-threading for faster processing of (very) large datasets
    * Faster algorithms for mtotdev() and htotdev() which are currently very slow
    
* Improve documentation
* Improve packaging for PyPi and/or other packaging systems: Conda?, PPA for Ubuntu/Debian?
* Tests for different noise types according to IEEE 1139, include power-spectral-density calculations 
* Conversion between phase noise and Allan variance 
* Phase noise calculations and plots
* Comparison to other libraries such as GPSTk

Make sure your patch does not break any of the tests, and does not 
significantly reduce the readability of the code.

Documentation generation
------------------------
See /docs for documentation in sphinx format. On Ubuntu this requires 
the **python-sphinx** and **python-numpydoc** packages.
html/pdf documentation using sphinx can be built locally with::

    /docs$ make html
    /docs$ make latexpdf

this generates html documentation in docs/_build/html and pdf 
documentation in docs/_build/latex.

The sphinx documentation is also auto-generated online

* https://allantools.readthedocs.org

Tests
-----

The tests compare the output of allantools to other programs such
as Stable32. Tests may be run using py.test (http://pytest.org). 
Package managers may install it with different binary names (e.g. pytest-3 
for the python3 version on debian).


Slow tests are marked 'slow' and tests failing because of a known
reason are marked 'fails'. To run all tests::
    
    $ py.test

To exclude known failing tests::

    $ py.test -m "not fails" --durations=10

To exclude tests that run slowly::

    $ py.test -m "not slow" --durations=10

To exclude both (note option change) and also check docstrings is ReST files ::

    $ py.test -k "not (slow or fails)" --durations=10 --doctest-glob='*.rst'

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

Notes for Pypi
--------------

Creating a source distribution

::

    python setup.py sdist

This creates a package in dist/
Testing the source distribution. The install takes a long time while 
compiling numpy and scipy.

::

    $ virtualenv tmp
    $ tmp/bin/pip install dist/AllanTools-2016.2.tar.gz 
    $ tmp/bin/python
    >>> import allantools

Uploading to PyPi.
(requries a ~/.pypirc with username and password)

::

    $ twine upload dist/*

Check result at https://pypi.python.org/pypi/AllanTools
