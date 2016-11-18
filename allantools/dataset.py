"""
Allantools dataset object

**Authors:** Frederic Meynadier (frederic.meynadier "at" gmail.com),
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

from . import allantools


class Dataset():
    """ Dataset class for Allantools

    :Example:
        ::

            import numpy as np
            # Load random data
            a = allantools.Dataset(data=np.random.rand(1000))
            # compute mdev
            a.compute("mdev")
            print(a.out["stat"])

    compute() returns the result of the computation and also stores it in the
    object's ``out`` member.

    """
    #: input data Dict, will be initialized by __init__()
    inp = {"data": None,
           "rate": None,
           "data_type": None,
           "taus": None}

    #: output data Dict, to be populated by compute()
    out = {"taus": None,
           "stat": None,
           "stat_err": None,
           "stat_n": None,
           "stat_unc": None,
           "stat_id": None}

    def __init__(self, data=None, rate=1.0, data_type="phase", taus=None):
        """ Initialize object with input data

        Parameters
        ----------
        data: np.array
            Input data. Provide either phase or frequency (fractional,
            adimensional)
        rate: float
            The sampling rate for data, in Hz. Defaults to 1.0
        data_typ: {'phase', 'freq'}
            Data type, i.e. phase or frequency. Defaults to "phase".
        taus: np.array
            Array of tau values, in seconds, for which to compute statistic.
            Optionally set taus=["all"|"octave"|"decade"] for automatic
            calculation of taus list

        Returns
        -------
        Dataset()
            A Dataset() instance

        """
        self.inp["data"] = data
        self.inp["rate"] = rate
        self.inp["data_type"] = data_type
        self.inp["taus"] = taus

    def set_input(self, data,
                  rate=1.0, data_type="phase", taus=None):
        """ Optionnal method if you chose not to set inputs on init

        Parameters
        ----------
        data: np.array
            Input data. Provide either phase or frequency (fractional,
            adimensional)
        rate: float
            The sampling rate for data, in Hz. Defaults to 1.0
        data_typ: {'phase', 'freq'}
            Data type, i.e. phase or frequency. Defaults to "phase".
        taus: np.array
            Array of tau values, in seconds, for which to compute statistic.
            Optionally set taus=["all"|"octave"|"decade"] for automatic
        """
        self.inp["data"] = data
        self.inp["rate"] = rate
        self.inp["data_type"] = data_type
        self.inp["taus"] = taus

    def compute(self, function):
        """Evaluate the passed function with the supplied data.

        Stores result in self.out.

        Parameters
        ----------
        function: str
            Name of the :mod:`allantools` function to evaluate

        Returns
        -------
        result: dict
            The results of the calculation.

        """
        try:
            func = getattr(allantools, function)
        except AttributeError:
            raise AttributeError("function must be defined in allantools")

        whitelisted = ["theo1", "mtie", "tierms"]

        if function[-3:] != "dev" and function not in whitelisted:
            # this should probably raise a custom exception type so
            # it's easier to distinguish from other bad things
            raise RuntimeError("function must be one of the 'dev' functions")
        result = func(self.inp["data"], rate=self.inp["rate"],
                      data_type=self.inp["data_type"], taus=self.inp["taus"])
        keys = ["taus", "stat", "stat_err", "stat_n"]
        result = {key: result[i] for i, key in enumerate(keys)}
        self.out = result.copy()
        self.out["stat_id"] = function
        return result
