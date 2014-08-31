#!/usr/bin/env python

from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(name='AllanTools',
      version='0.23',
      description='Allan deviation and related time/frequency statistics',
      author='Anders Wallin',
      author_email='anders.e.e.wallin@gmail.com',
      url='https://github.com/aewallin/allantools',
      license='GPLv3+',
      packages=['allantools',],
      requires=['numpy'],
      include_dirs=[numpy.get_include()],
      long_description="""Given phase or fractional frequency data this package calculates:
                        Allan deviation, overlapping Allan deviation, modified Allan deviation
                        Hadamard deviation, overlapping Hadamard deviation, time deviation,
                        total deviation, MTIE, TIE-rms"""
     )

