"""
Allantools plotting utilities

**Authors:** Frederic Meynadier (frederic.meynadier "at" gmail.com),
    Mike DePalatis (http://mike.depalatis.net)

Version history
---------------

**2019.07**
- Initial release

License
-------

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


class Plot(object):
    """ A class for plotting data once computed by Allantools

    :Example:
        ::

            import allantools
            import numpy as np
            a = allantools.Dataset(data=np.random.rand(1000))
            a.compute("mdev")
            b = allantools.Plot()
            b.plot(a)
            b.show()

    Uses matplotlib. self.fig and self.ax stores the return values of
    matplotlib.pyplot.subplots(). plot() sets various defaults, but you
    can change them by using standard matplotlib method on self.fig and self.ax
    """
    def __init__(self, no_display=False):
        """ set ``no_display`` to ``True`` when we don't have an X-window
        (e.g. for tests)
        """
        try:
            import matplotlib
            if no_display:
                matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            self.plt = plt
        except ImportError:
            raise RuntimeError("Matplotlib is required for plotting")
        self.fig, self.ax = plt.subplots()
        self.ax.set_xscale("log")
        self.ax.set_yscale("log")

    def plot(self, atDataset,
             errorbars=False,
             grid=False,
             **kwargs):
        """ Use matplotlib methods for plotting

        Additional keywords arguments are passed to
        :py:func:`matplotlib.pyplot.plot`.

        Parameters
        ----------
        atDataset : allantools.Dataset()
            a dataset with computed data
        errorbars : boolean
            Plot errorbars. Defaults to False
        grid : boolean
            Plot grid. Defaults to False

        """
        if errorbars:
            self.ax.errorbar(atDataset.out["taus"],
                             atDataset.out["stat"],
                             yerr=atDataset.out["stat_err"],
                             **kwargs)
        else:
            self.ax.plot(atDataset.out["taus"],
                         atDataset.out["stat"],
                         **kwargs)
        self.ax.set_xlabel("Tau")
        self.ax.set_ylabel(atDataset.out["stat_id"])

        if grid:
            self.ax.grid(True, which="minor", ls="-", color='0.65')
            self.ax.grid(True, which="major", ls="-", color='0.25')
        else:
            self.ax.grid(False)

    def show(self):
        """Calls matplotlib.pyplot.show()

        Keeping this separated from ``plot()`` allows to tweak display before
        rendering
        """
        self.plt.show()

    def save(self, f):
        """Save figure to file
        """
        self.plt.savefig(f)
