#!/usr/bin/python

import sys
sys.path.append("..")
import allantools as at
import pytest


def test_plot():
    ds = at.Dataset(data=at.noise.white(1000), rate=1.234)
    ds.compute("adev")
    p = at.Plot(no_display=True)
    p.plot(ds, errorbars=True)
    # p.show()  # can't show() in test

if __name__ == "__main__":
    test_plot()
