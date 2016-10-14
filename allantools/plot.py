"""
Allantools plotting utilities

**Authors:** Frédéric Meynadier (frederic.meynadier "at" gmail.com),
    Mike DePalatis (http://mike.depalatis.net)

Version history
---------------

**unreleased**
- Initial commit

License
-------

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
import matplotlib.pyplot as plt


class atPlot():
    def __init__(self):
        self.fig, self.ax = plt.subplots()
        self.ax.set_xscale("log")
        self.ax.set_yscale("log")

    def plot(self, atDataset):
        self.ax.plot(atDataset.out["taus"],
                     atDataset.out["stat"])

    def show(self):
        plt.show()
