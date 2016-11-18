.. py:currentmodule:: allantools
.. py:module: allantools
.. py:module: allantools.noise

API
===

The Dataset() class
-------------------

**New in version 2016.11**

This class allows simple data handling.

.. autoclass:: Dataset
    :members:

    .. automethod:: __init__


The Plot() class
-------------------

**New in version 2016.11**

This class allows simple data plotting.


.. autoclass:: Plot
    :members:
    
    .. automethod:: __init__

Implemented statistics functions
--------------------------------

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
``theo1()``                             Theo1 deviation
``mtie()``                              Maximum Time Interval Error
``tierms()``                            Time Interval Error RMS
``gradev()``                            Gap resistant overlapping Allan deviation
=====================================   ====================================================

Low-level access to the statistics functions
--------------------------------------------

The deviation functions are generally of the form::

    (tau_out, adev, adeverr, n) = allantools.adev(data, rate=1.0, data_type="phase", taus=None)

*Inputs:*

    * **data** = list of phase measurements in seconds, or list of fractional frequency measurements (nondimensional)
    * **rate**  = sample rate of data in Hz , i.e. interval between phase measurements is 1/rate seconds.
    * **data_type=** = either "phase" or "freq"
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

Noise Generation
----------------

.. autofunction:: allantools.noise.white
.. autofunction:: allantools.noise.brown
.. autofunction:: allantools.noise.violet
.. autofunction:: allantools.noise.pink

Utilities
---------

.. autofunction:: frequency2phase
.. autofunction:: phase2frequency
.. autofunction:: phase2radians
.. autofunction:: edf_simple
.. autofunction:: edf_greenhall
.. autofunction:: edf_totdev
.. autofunction:: edf_mtotdev
.. autofunction:: three_cornered_hat_phase
