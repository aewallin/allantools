#!/usr/bin/python
import allantools as at
import numpy as np
import os.path

# Taken from examples
f = 1
x = np.loadtxt(os.path.join(os.path.dirname(__file__), "x.dat"))
y = np.loadtxt(os.path.join(os.path.dirname(__file__), "y.dat"))
y_gap = np.loadtxt(os.path.join(os.path.dirname(__file__), "y_gap.dat"))


def test_gradev():
    x_ax, y_ax, (el, eh), ns = at.gradev(
        y, data_type='phase', rate=f, taus=x)
    x_ax_gap, y_ax_gap, (el_gap, eh_gap), ns_gap = at.gradev(
        y_gap, data_type='phase', rate=f, taus=x)
    for i in range(len(x_ax)):
        i_gap = np.where(x_ax_gap == x_ax[i])
        if len(i_gap[0]) > 0:
            # Seems like a loose condition
            assert(np.log(y_ax[i]/y_ax_gap[i_gap[0][0]]) < 1)

if __name__ == "__main__":
    test_gradev()
