allantools
==========

A small Python library for calculating Allan deviation and related statistics.

Input data should be evenly spaced observations of either fractional frequency,
or phase in seconds. Deviations are calculated for given tau values in seconds.

Usage:
```python
import allantools
rate = 1/float(data_interval) # in Hz
taus = [1,2,4,8,16] # in seconds
# fractional frequency data
(taus_used, adev, adeverror, adev_n) = allantools.adev(fract_freqdata, rate, taus)
# phase data
(taus_used, adev, adeverror, adev_n) = allantools.adev_phase(phasedata, rate, taus)
```

These statistics are currently included:
* ADEV, Allan deviation
* OADEV, overlapping Allan deviation,
* MDEV, modified Allan deviation,
* TDEV, time deviation
* HDEV, Hadamard deviation

see /tests for tests that compare allantools output to other programs.
