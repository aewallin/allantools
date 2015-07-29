"""
Allan deviation tools
Anders Wallin (anders.e.e.wallin "at" gmail.com)
v1.0 2014 January

Functions for generating noise
AW 2014

See: http://en.wikipedia.org/wiki/Colors_of_noise

This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy
# import math


def white(N):
    """ white noise has constant PSD (up to the nyquist frequency) """
    return numpy.random.randn(N)


def violet(N):
    """ violet noise with f^2 PSD """
    return numpy.diff(numpy.random.randn(N))


def brown(N):
    """ Brownian or random walk (diffusion) noise with 1/f^2 PSD
        not really a color... rather Brownian or random-walk
    """
    return numpy.cumsum(numpy.random.randn(N))


def pink(N, depth=80):
    """ N-length vector with (approximate) pink noise
    pink noise has 1/f PSD """
    a = []
    s = iterpink(depth)
    for n in range(N):
        a.append(next(s))
    return a


def iterpink(depth=20):
    """Generate a sequence of samples of pink noise.

    pink noise generator
    from http://pydoc.net/Python/lmj.sound/0.1.1/lmj.sound.noise/

    Based on the Voss-McCartney algorithm, discussion and code examples at
    http://www.firstpr.com.au/dsp/pink-noise/

    depth: Use this many samples of white noise to calculate the output. A
      higher  number is slower to run, but renders low frequencies with more
      correct power spectra.

    Generates a never-ending sequence of floating-point values. Any continuous
    set of these samples will tend to have a 1/f power spectrum.
    """
    values = numpy.random.randn(depth)
    smooth = numpy.random.randn(depth)
    source = numpy.random.randn(depth)
    sumvals = values.sum()
    i = 0
    while True:
        yield sumvals + smooth[i]

        # advance the index by 1. if the index wraps, generate noise to use in
        # the calculations, but do not update any of the pink noise values.
        i += 1
        if i == depth:
            i = 0
            smooth = numpy.random.randn(depth)
            source = numpy.random.randn(depth)
            continue

        # count trailing zeros in i
        c = 0
        while not (i >> c) & 1:
            c += 1

        # replace value c with a new source element
        sumvals += source[i] - values[c]
        values[c] = source[i]


def pinknoise_to_file(N=10000, filename="pinknoise_frequency.txt"):
    n = pink(N, depth=80)
    f0 = 10e6  # imagined carrier, 10 MHz
    n = [x + f0 for x in n]  # add noise and carrier
    fr = [(x - f0) / float(f0) for x in n]  # fractional frequency
    with open(filename, 'w') as f:
        f.write("# simulated oscillator with pink frequency noise\n")
        f.write("# fractional frequency data \n")
        f.write("# number of datapoints %d\n" % len(fr))
        for ff in fr:
            f.write("%0.6g\n" % ff)
        #print sum(n)/float(len(n))


if __name__ == "__main__":
    #pinknoise_to_file()
    print("Nothing to see here.")
