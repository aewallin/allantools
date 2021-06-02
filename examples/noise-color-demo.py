"""
 Allan deviation tools
 Anders Wallin (anders.e.e.wallin "at" gmail.com)
 v1.0 2014 January

 This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


import numpy
# only for plotting, not required for calculations
import matplotlib.pyplot as plt
import allantools
from allantools import noise


def plotallan(plt, y, rate, taus, style, label=""):
    (t2, ad, ade, adn) = allantools.mdev(
        y, data_type='freq', rate=rate, taus=taus)
    plt.loglog(t2, ad, style, label=label)


def plotallan_phase(plt, y, rate, taus, style, label="", alpha=1.0):
    (t2, ad, ade, adn) = allantools.mdev(
        y, data_type='phase', rate=rate, taus=taus)
    plt.loglog(t2, ad, style, label=label, alpha=alpha)


def plotline(plt, alpha, taus, style, label=""):
    """ plot a line with the slope alpha """
    y = [pow(tt, alpha) for tt in taus]
    plt.loglog(taus, y, style, label=label)


if __name__ == "__main__":
    print("allatools tests with various noise types")
    # we test ADEV etc. by calculations on synthetic data
    # with known slopes of ADEV

    t = [tt for tt in numpy.logspace(0, 3, 50)]  # tau values from 1 to 1000
    plt.subplot(111, xscale="log", yscale="log")
    N = 100000

    # Colors: http://en.wikipedia.org/wiki/Colors_of_noise

    # Brownian a.k.a random walk  frequency => sqrt(tau) ADEV
    print("Random Walk frequency noise - should have sqrt(tau) ADEV")
    freq_rw = noise.brown(N)
    phase_rw_rw = numpy.cumsum(noise.brown(N))  # integrate to get  phase
    plotallan(plt, freq_rw, 1, t, 'm.')
    plotallan_phase(plt, phase_rw_rw, 1, t, 'mo',
                    label='random walk frequency')
    plotline(plt, +0.5, t, 'm', label="f^(+1/2)")

    # pink frequency noise => constant ADEV
    print("Pink frequency noise - should have constant ADEV")
    freq_pink = noise.pink(N)
    phase_p = numpy.cumsum(noise.pink(N))  # integrate to get phase, color??
    plotallan_phase(plt, phase_p, 1, t, 'co',
                    label="pink/flicker frequency noise")
    plotallan(plt, freq_pink, 1, t, 'c.')
    plotline(plt, 0, t, 'c', label="f^0")

    # white frequency modulation => 1/sqrt(tau) ADEV
    print("White frequency noise - should have 1/sqrt(tau) ADEV")
    freq_white = noise.white(N)
    # integrate to get Brownian, or random walk phase
    phase_rw = noise.brown(N)
    plotallan(plt, freq_white, 1, t, 'b.')
    plotallan_phase(plt, phase_rw, 1, t, 'bo',
                    label='random walk phase a.k.a. white frequency noise')
    plotline(plt, -0.5, t, 'b', label="f^(-1/2)")

    # pink phase noise => 1/tau ADEV and MDEV
    print("Pink phase noise - should tau^(-3/2) MDEV")
    phase_pink = noise.pink(N)
    plotallan_phase(plt, phase_pink, 1, t, 'ko',
                    label="pink/flicker phase noise")
    plotline(plt, -1, t, 'k', label="f^(-1)")

    # white phase noise => 1/tau ADEV  and tau^(-3/2) MDEV
    print("White phase noise - should have 1/tau ADEV")
    phase_white = noise.white(N)
    plotallan_phase(plt, phase_white, 1, t, 'ro', label="white phase noise")
    freq_w = noise.violet(N)  # diff to get frequency, "Violet noise"
    plotallan(plt, freq_w, 1, t, 'r.')
    plotline(plt, -1.5, t, 'r', label="f^(-3/2)")

    plt.title('allantools noise type demo')
    plt.xlabel('Tau')
    plt.ylabel('Modified Allan deviation')
    print("Done.")
    plt.legend(loc="lower left", framealpha=0.7)
    plt.grid()
    plt.show()
