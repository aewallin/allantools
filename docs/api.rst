.. py:currentmodule:: allantools
.. py:module: allantools

API
===

Implemented functions
---------------------

=====================================   ====================================================
Function                                Description
=====================================   ====================================================
``adev()``                              Allan deviation
``oadev()``                             Overlapping Allan deviation
``mdev()``                              Modified Allan deviation
``tdev()``                              Time deviation
``hdev()``                              Hadamard deviation
``ohdev()``                             Overlapping Hadamard deviation
``totdev()``                            Total deviation
``mtotdev()``                           Modified total deviation
``ttotdev()``                           Time total deviation
``htotdev()``                           Hadamard total deviation
``mtie()``                              Maximum Time Interval Error
``tierms()``                            Time Interval Error RMS
``gradev()``                            Gap resistant overlapping Allan deviation
``theo1()``                             Theo1 deviation
``uncertainty_estimate()``              Determine the uncertainty of a given two-sample variance estimate
``three_cornered_hat_phase()``          Three Cornered Hat Method
=====================================   ====================================================


References
-----------

* http://www.wriley.com/paper4ht.htm
* http://en.wikipedia.org/wiki/Allan_variance
* http://tf.nist.gov/general/pdf/2220.pdf
* http://www.stable32.com/Handbook.pdf

for code see e.g.:
    
* http://www.mathworks.com/matlabcentral/fileexchange/26659-allan-v3-0
* http://www.mathworks.com/matlabcentral/fileexchange/26637-allanmodified
* http://www.leapsecond.com/tools/adev_lib.c


General usage
--------------

The deviation functions are generally of the form::

    (tau_out, adev, adeverr, n) = allantools.adev(frequency=y, phase=x, rate=r, taus=t)

*Inputs:*

    * **phase** = list of phase measurements in seconds, e.g. from a time-interval-counter
    * **frequency** = list of fractional frequency measurements (nondimensional), e.g. from a frequency-counter
    * **rate**  = sample rate of data in Hz , i.e. interval between phase measurements is 1/rate seconds.
    * **taus**  = list of tau-values for ADEV computation. The keywords "all", "octave", or "decade" can also be used.

*Outputs*

    * **tau_out** = list of tau-values for which deviations were computed
    * **adev**    = list of ADEV (or another statistic) deviations
    * **adeverr** = list of estimated errors of allan deviations. some functions instead return  a confidence interval (**err_l**, **err_h**)
    * **n**       = list of number of pairs in allan computation. standard error is adeverr = adev/sqrt(n)
    
Function listing
================

Statistics
----------

.. autofunction:: adev
.. autofunction:: oadev
.. autofunction:: mdev
.. autofunction:: hdev
.. autofunction:: ohdev
.. autofunction:: tdev
.. autofunction:: totdev
.. autofunction:: mtotdev
.. autofunction:: ttotdev
.. autofunction:: htotdev
.. autofunction:: theo1
.. autofunction:: mtie
.. autofunction:: tierms

Utilities
---------

.. autofunction:: frequency2phase
.. autofunction:: phase2frequency
.. autofunction:: phase2radians
.. autofunction:: uncertainty_estimate
.. autofunction:: three_cornered_hat_phase
