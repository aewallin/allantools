#!/usr/bin/env python

from setuptools import setup
import json
import os

pkginfo_path = os.path.join(os.path.dirname(__file__),
                            'allantools',
                            'allantools_info.json')
pkginfo = json.load(open(pkginfo_path))

setup(name=pkginfo['name'],
      version=pkginfo['version'],
      description=pkginfo['description'],
      author=pkginfo['main_author'],
      author_email=pkginfo['main_author_email'],
      url=pkginfo['url'],
      license=pkginfo['license'],
      packages=['allantools', ],
      package_data={'allantools': ['allantools_info.json']},
      install_requires=['numpy', 'scipy'],
      requires=['numpy', 'scipy'],
      # include_dirs=[numpy.get_include()],
      setup_requires=['pytest-runner'],
      tests_require=['pytest', 'numpy'],
      long_description=open('README.rst', 'r').read(),
      )
