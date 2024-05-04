"""
Allantools Noise object

**Authors:** Julia Leute (julia.leute "at" gmail.com)
    Anders Wallin (anders.e.e.wallin "at" gmail.com)

Version history
---------------

**2019.07**
- Initial release, see also https://github.com/jleute/colorednoise

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

import numpy as np


class Noise(object):
    """ Generate discrete colored noise

    Python / Numpy implementation of [Kasdin1992]_
    Kasdin, N.J., Walter, T., "Discrete simulation of power law noise [for
    oscillator stability evaluation]," Frequency Control Symposium, 1992.
    46th., Proceedings of the 1992 IEEE, pp.274,283, 27-29 May 1992
    http://dx.doi.org/10.1109/FREQ.1992.270003

    Parameters
    ----------
    nr: integer
        length of generated time-series
        must be power of two
    qd: float
        discrete variance
    b: float
        noise type

    +----+--------------------------------------------+
    | b  |  noise type                                |
    +====+============================================+
    | 0  | White Phase Modulation (WPM)               |
    +----+--------------------------------------------+
    | -1 | Flicker Phase Modulation (FPM)             |
    +----+--------------------------------------------+
    | -2 | White Frequency Modulation (WFM)           |
    +----+--------------------------------------------+
    | -3 | Flicker Frequency Modulation (FFM)         |
    +----+--------------------------------------------+
    | -4 | Random Walk Frequency Modulation (RWFM)    |
    +----+--------------------------------------------+

    Returns
    -------
    Noise()
        A Noise() instance

    :Example:
        ::

            import numpy as np
            noise = allantools.Noise(nr=2*8, qd=1.0e-20, b=-1)
            noise.generateNoise()
            print noise.time_series

    """

    def __init__(self, nr=2, qd=1, b=0):
        """ Initialize object with input data
        """
        self.nr = nr
        self.qd = qd
        self.b = b
        self.time_series = np.array([])

    def set_input(self, nr=2, qd=1, b=0):
        """ Set inputs after initialization

        Parameters
        ----------
        nr: integer
            length of generated time-series
            number must be power of two
        qd: float
            discrete variance
        b: float
            noise type

            +----+--------------------------------------------+
            | b  |  noise type                                |
            +====+============================================+
            | 0  | White Phase Modulation (WPM)               |
            +----+--------------------------------------------+
            | -1 | Flicker Phase Modulation (FPM)             |
            +----+--------------------------------------------+
            | -2 | White Frequency Modulation (WFM)           |
            +----+--------------------------------------------+
            | -3 | Flicker Frequency Modulation (FFM)         |
            +----+--------------------------------------------+
            | -4 | Random Walk Frequency Modulation (RWFM)    |
            +----+--------------------------------------------+

        """
        self.nr = nr
        self.qd = qd
        self.b = b

    def generateNoise(self):
        """ Generate noise time series based on input parameters

        Returns
        -------
        time_series: np.array
            Time series with colored noise.
            len(time_series) == nr

        """
        # Fill wfb array with white noise based on given discrete variance
        wfb = np.zeros(self.nr*2)
        wfb[:self.nr] = np.random.normal(0, np.sqrt(self.qd), self.nr)
        # Generate the hfb coefficients based on the noise type
        mhb = -self.b/2.0
        hfb = np.zeros(self.nr*2)
        hfb = np.zeros(self.nr*2)
        hfb[0] = 1.0
        indices = np.arange(self.nr-1)
        hfb[1:self.nr] = (mhb+indices)/(indices+1.0)
        hfb[:self.nr] = np.multiply.accumulate(hfb[:self.nr])
        # Perform discrete Fourier transform of wfb and hfb time series
        wfb_fft = np.fft.rfft(wfb)
        hfb_fft = np.fft.rfft(hfb)
        # Perform inverse Fourier transform of the product of wfb and hfb FFTs
        time_series = np.fft.irfft(wfb_fft*hfb_fft)[:self.nr]
        self.time_series = time_series

    def phase_psd_from_qd(self, tau0=1.0):
        """ return phase power spectral density coefficient :math:`g_b`
            for noise-type defined by (qd, b, tau0)
            where tau0 is the interval between data points

            Colored noise generated with (qd, b, tau0) parameters will
            show a phase power spectral density of
            :math:`S_x(f) = Phase_{PSD}(f) = g_b  f^b`

            [Kasdin1992]_ eqn (39)
        """
        return self.qd*2.0*pow(2.0*np.pi, self.b)*pow(tau0, self.b+1.0)

    def frequency_psd_from_qd(self, tau0=1.0):
        """ return frequency power spectral density coefficient :math:`h_a`
            for the noise type defined by (qd, b, tau0)

            Colored noise generated with (qd, b, tau0) parameters will
            show a frequency power spectral density of

            .. math::

                S_y(f) = Frequency_{PSD}(f) = h_a  f^a

            where the slope :math:`a` comes from the phase PSD slope :math:`b`:
            :math:`a = b + 2`

            [Kasdin1992]_ eqn (39)
        """
        a = self.b + 2.0
        return self.qd*2.0*pow(2.0*np.pi, a)*pow(tau0, a-1.0)

    def adev(self, tau0, tau):
        """ return predicted ADEV of noise-type at given tau

        """
        prefactor = self.adev_from_qd(tau0=tau0, tau=tau)
        c = self.c_avar()
        avar = pow(prefactor, 2)*pow(tau, c)
        return np.sqrt(avar)

    def mdev(self, tau0, tau):
        """ return predicted MDEV of noise-type at given tau

        """
        prefactor = self.mdev_from_qd(tau0=tau0, tau=tau)
        c = self.c_mvar()
        mvar = pow(prefactor, 2)*pow(tau, c)
        return np.sqrt(mvar)

    def c_avar(self):
        """ return tau exponent "c" for noise type.
            AVAR = prefactor * h_a * tau^c
        """
        if self.b == -4:
            return 1.0
        elif self.b == -3:
            return 0.0
        elif self.b == -2:
            return -1.0
        elif self.b == -1:
            return -2.0
        elif self.b == 0:
            return -2.0

    def c_mvar(self):
        """ return tau exponent "c" for noise type.
            MVAR = prefactor * h_a * tau^c
        """
        if self.b == -4:
            return 1.0
        elif self.b == -3:
            return 0.0
        elif self.b == -2:
            return -1.0
        elif self.b == -1:
            return -2.0
        elif self.b == 0:
            return -3.0

    def adev_from_qd(self, tau0=1.0, tau=1.0):
        """ prefactor for Allan deviation for noise type defined by (qd, b, tau0)

        Colored noise generated with (qd, b, tau0) parameters will
        show an Allan variance of:

        .. math::

            AVAR = prefactor \\cdot h_a \\cdot \\tau^c

        where :math:`a = b + 2` is the slope of the frequency PSD.
        and :math:`h_a` is the frequency PSD prefactor :math:`S_y(f) = h_a  f^a`

        The relation between a, b, c is:

        +---+---+---------+----------+
        | a | b | c(AVAR) | c(MVAR)  |
        +===+===+=========+==========+
        |-2 |-4 |  1      | 1        |
        +---+---+---------+----------+
        |-1 |-3 |  0      | 0        |
        +---+---+---------+----------+
        | 0 |-2 | -1      | -1       |
        +---+---+---------+----------+
        |+1 |-1 | -2      | -2       |
        +---+---+---------+----------+
        |+2 | 0 | -2      | -3       |
        +---+---+---------+----------+

        Coefficients from [Dawkins2007]_.

        [Vernotte2015]_ Table I

        """
        # g_b = self.phase_psd_from_qd(tau0)
        h_b = self.frequency_psd_from_qd(tau0)
        f_h = 0.5/tau0
        if self.b == 0:  # WPM
            coeff = 3.0*f_h / (4.0*pow(np.pi, 2))  # E, White PM, tau^-1
        elif self.b == -1:
            # D, Flicker PM, tau^-1
            gamma = 0.57721566490153286060651209  # https://en.wikipedia.org/wiki/Euler%27s_constant
            coeff = (3*gamma-np.log(2)+3*np.log(2.0*np.pi*f_h*tau))/(4.0*pow(np.pi, 2))
        elif self.b == -2:  # C, white FM,  1/sqrt(tau)
            coeff = 0.5
        elif self.b == -3:
            coeff = 2*np.log(2)  # B, flicker FM,  constant ADEV
        elif self.b == -4:
            coeff = 2.0*pow(np.pi, 2)/3.0  # A, RW FM, sqrt(tau)

        # return np.sqrt(coeff*g_b*pow(2.0*np.pi, 2))
        return np.sqrt(coeff*h_b)

    def mdev_from_qd(self, tau0=1.0, tau=1.0):
        # FIXME: tau is unused here - can we remove it?
        """ prefactor for Modified Allan deviation for noise
            type defined by (qd, b, tau0)

            Colored noise generated with (qd, b, tau0) parameters will
            show an Modified Allan variance of:

            .. math::

                MVAR = prefactor \\cdot h_a \\cdot \\tau^c

            where :math:`a = b + 2` is the slope of the frequency PSD.
            and :math:`h_a` is the frequency PSD prefactor :math:`S_y(f) = h_a  f^a`

        """
        # g_b = self.phase_psd_from_qd(tau0)   # not used
        h_b = self.frequency_psd_from_qd(tau0)
        # f_h = 0.5/tau0 #unused!?
        if self.b == 0:
            coeff = 3.0/(8.0*pow(np.pi, 2))  # E, White PM, tau^-{3/2}
        elif self.b == -1:
            # D, Flicker PM, tau^-1
            coeff = (24.0*np.log(2)-9.0*np.log(3))/8.0/pow(np.pi, 2)
        elif self.b == -2:
            # C, white FM,  1/sqrt(tau)
            coeff = 0.25
        elif self.b == -3:
            # B, flicker FM,  constant MDEV
            # coeff = 2.0*np.log(3.0*pow(3.0, 11.0/16.0)/4.0)
            coeff = (27.0*np.log(3)-32.0*np.log(2))/8.0/pow(np.pi, 2)  # Vernotte Table I
            coeff = (27.0/20.0)*np.log(2)  # Benkler2015 Table 1
        elif self.b == -4:
            # A, RW FM, sqrt(tau)
            coeff = 11.0/20.0*pow(np.pi, 2)

        return np.sqrt(coeff*h_b)

    def pdev_from_qd(self, tau0=1.0, tau=1.0):
        # FIXME: tau is unused here - can we remove it?
        """ prefactor for Parabolic Allan deviation for noise
            type defined by (qd, b, tau0)

            Colored noise generated with (qd, b, tau0) parameters will
            show an Parabolic Allan variance of:

            .. math::

                PVAR = prefactor \\cdot h_a \\cdot \\tau^c

            where :math:`a = b + 2` is the slope of the frequency PSD.
            and :math:`h_a` is the frequency PSD prefactor :math:`S_y(f) = h_a  f^a`

        """
        g_b = self.phase_psd_from_qd(tau0)
        # f_h = 0.5/tau0 #unused!?
        if self.b == 0:  # WPM, tau^(-3/2)
            coeff = 3.0/(2.0*pow(np.pi, 2))  # Vernotte Table I
        elif self.b == -1:
            # D, Flicker PM, tau^-1
            coeff = (3.0*np.log(16)-1.0)/2.0/pow(np.pi, 2)
        elif self.b == -2:  # C, white FM,  1/sqrt(tau)
            coeff = 3.0/5.0
        elif self.b == -3:
            # B, flicker FM,  constant xDEV
            coeff = 2.0*(7.0-np.log(16.0))/5.0
        elif self.b == -4:
            # A, RW FM, sqrt(tau)
            coeff = 26.0/35.0*pow(np.pi, 2)

        return np.sqrt(coeff*g_b*pow(2.0*np.pi, 2))


# end of file noise_kasdin.py
