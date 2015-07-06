.. py:currentmodule:: allantools
.. py:module: allantools

API
===

Implemented functions
---------------------

=====================================   ====================================================
Function                                Notes
=====================================   ====================================================
``adev()`` and ``adev_phase()``         Allan deviation
``oadev()`` and ``oadev_phase()``       Overlapping Allan deviation
``mdev()`` and ``mdev_phase()``         Modified Allan deviation
``tdev()`` and ``tdev_phase()``         Time deviation ( modified variance scaled by (t^2/3))
``hdev()`` and ``hdev_phase()``         Hadamard deviation
``ohdev()`` and ``ohdev_phase()``       Overlapping Hadamard deviation
``totdev()`` and ``totdev_phase()``     Total deviation
``mtie()`` and ``mtie_phase()``         Maximum Time Interval Error
``tierms()`` and ``tierms_phase()``     Time Interval Error RMS
``gradev()`` and ``gradev_phase()``     Gap resistant overlapping Allan deviation
``uncertainty_estimate()``              Determine the uncertainty of a given two-sample variance estimate
``three_cornered_hat_phase()``          Apply Three Cornered Hat Method (see examples)
=====================================   ====================================================

To implement
------------

* Modified Total. The modified total variance, MTOT, is total version of the modified Allan variance.

It is defined for phase data as::
            
                                1           N-3m+1  1  N+3m-1
       Mod s^2 total(t) = ----------------- sum     -- sum      [0zi*(m)]^2
                           2m^2t0^2(N-3m+1) n=1     6m i=n-3m
                       
where the ``0zi*(m)`` terms are the phase averages from a triply-extended
sequence created by uninverted even reflection at each end,
and the prefix 0 denotes that the linear trend has been removed.

* Time Total (modified total variance scaled by (t^2/3) )
* Hadamard Total


References
-----------

* http://www.wriley.com/paper4ht.htm
* http://en.wikipedia.org/wiki/Allan_variance

for code see e.g.:
    
* http://www.mathworks.com/matlabcentral/fileexchange/26659-allan-v3-0
* http://www.mathworks.com/matlabcentral/fileexchange/26637-allanmodified
* http://www.leapsecond.com/tools/adev_lib.c


General usage:
--------------

*Inputs (phase data):*

    * **phase** = list of phase measurements in seconds, e.g. from a time-interval-counter
    * **rate**  = sample rate of data, i.e. interval between phase measurements is 1/rate
    * **taus**  = list of tau-values for ADEV computation
    
*Inputs (frequency data):*

    * **data** = list of fractional frequency measurements (nondimensional), e.g. from a frequency-counter
    * **rate**  = sample rate of data, i.e. gate time of a zero-dead-time counter is 1/rate
    * **taus**  = list of tau-values for ADEV computation
    
*Output (tau_out, adev, adeverr, n)*

    * **tau_out** = list of tau-values for which deviations were computed
    * **adev**    = list of ADEV (or another statistic) deviations
    * **adeverr** = list of estimated errors of allan deviations
    * **n**       = list of number of pairs in allan computation. standard error is adeverr = adev/sqrt(n)
    
Function listing
================

.. autofunction:: adev
.. autofunction:: adev_phase
::
    
                  1      
     s2y(t) = --------- sum [x(i+2) - 2x(i+1) + x(i) ]^2
               2*tau^2  
     
.. autofunction:: oadev
.. autofunction:: oadev_phase
.. autofunction:: mdev
.. autofunction:: mdev_phase
               
::
    
                   N-3m+1     j+m-1
    Mod s2y(t) = ------------------  sum      { sum   [x(i+2m) - 2x(i+m) + x(i) ]  }**2
                 2m**2 t**2 (N-3m+1) j=1        i=j
     
.. autofunction:: hdev
.. autofunction:: hdev_phase
.. autofunction:: ohdev
.. autofunction:: ohdev_phase
.. autofunction:: tdev
.. autofunction:: tdev_phase                
.. autofunction:: totdev
.. autofunction:: totdev_phase

::
    
                     1         N-1
    totvar(t) = ------------  sum   [ x*(i-m) - 2x*(i)+x*(i+m) ]**2
                2 t**2 (N-2)  i=2
                
.. autofunction:: mtie
.. autofunction:: mtie_phase
.. autofunction:: tierms
.. autofunction:: tierms_phase
.. autofunction:: frequency2phase
.. autofunction:: uncertainty_estimate
.. autofunction:: three_cornered_hat_phase