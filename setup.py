#!/usr/bin/env python

from setuptools import setup
from setuptools.command.test import test as TestCommand
# import numpy
import sys

class PyTest(TestCommand):
    """ Setup tests with py.test framework """
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)


setup(name='AllanTools',
      version='2016.02',
      description='Allan deviation and related time/frequency statistics',
      author='Anders Wallin',
      author_email='anders.e.e.wallin@gmail.com',
      url='https://github.com/aewallin/allantools',
      license='GPLv3+',
      packages=['allantools',],
      requires=['numpy'],
      #include_dirs=[numpy.get_include()],
      tests_require=['pytest'],
      cmdclass={'test': PyTest},
      long_description="""Given phase or fractional frequency data this package calculates:
                        Allan deviation, overlapping Allan deviation, modified Allan deviation
                        Hadamard deviation, overlapping Hadamard deviation, time deviation,
                        total deviation, MTIE, TIE-rms"""
     )

