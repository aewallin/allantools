#!/usr/bin/env python

from distutils.core import setup

setup(name='AllanTools',
      version='0.21',
      description='Allan deviation and relates time/frequency statistics',
      author='Anders Wallin',
      author_email='anders.e.e.wallin@gmail.com',
      url='https://github.com/aewallin/allantools',
      license='GPLv3+',
      packages=['allantools',],
      long_description ='Given phase or fractional frequency data this package calculates: ' \
                        'Allan deviation, overlapping Allan deviation, modified Allan deviation, ' \
                        'Hadamard deviation, overlapping Hadamard deviation, time deviation, ' \
                        'total deviation, MTIE, TIE-rms'
     )
