"""
Allantools frontend object

**Author:** Frédéric Meynadier (frederic.meynadier "at" gmail.com)

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


class AllantoolsFrontend():
    """ Frontend class for Allantools

    """
    def __init__(self):
        self.inp = {"data": None,
                    "rate": None,
                    "data_type": None,
                    "taus": None}
        self.out = {"taus": None,
                    "stat": None,
                    "stat_err": None,
                    "stat_n": None,
                    "stat_unc": None}

    def set_input(self, data,
                  rate=1.0, data_type="phase", taus=None):
        self.inp["data"] = data
        self.inp["rate"] = rate
        self.inp["data_type"] = data_type
        self.inp["taus"] = taus

    def compute(self, function):
        (self.out["taus"],
         self.out["stat"],
         self.out["stat_err"],
         self.out["stat_n"]) = function(self.inp["data"],
                                        rate=self.inp["rate"],
                                        data_type=self.inp["data_type"],
                                        taus=self.inp["taus"])
        return self.out
