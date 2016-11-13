#!/bin/sh

# run tests and generate html-report in /htmlcov
# works at least on Ubuntu 16.04LTS
# AW2016-11-13

python-coverage run --source allantools setup.py test
python-coverage html
python-coverage report
