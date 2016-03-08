#!/usr/bin/env python

from setuptools import setup
# import numpy

setup(name='AllanTools',
      version='2016.02',
      description='Allan deviation and related time/frequency statistics',
      author='Anders Wallin',
      author_email='anders.e.e.wallin@gmail.com',
      url='https://github.com/aewallin/allantools',
      license='GPLv3+',
      packages=['allantools',],
      requires=['numpy', 'scipy'],
      #include_dirs=[numpy.get_include()],
      setup_requires=['pytest-runner'],
      tests_require=['pytest'],
      long_description="""Given phase or fractional frequency data this package calculates:
                        Allan deviation, overlapping Allan deviation, modified Allan deviation,
                        Hadamard deviation, overlapping Hadamard deviation, time deviation,
                        total deviation, MTIE, TIE-RMS. Synthetic noise data generators are also included."""
     )

