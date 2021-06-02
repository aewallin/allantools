"""
    This file is part of allantools https://github.com/aewallin/allantools

    Benchmark for real-time statistics computations

"""
import allantools as at
import numpy
import time


def bench(decades, n_pts, function, n_add=1e4):
    # decades = 4
    afs = numpy.floor(numpy.logspace(0, decades, n_pts*decades))
    afs = [int(xx) for xx in afs]
    # print len(afs)
    # print afs

    # realtime object
    dev_rt = function(afs=afs, tau0=1.0)

    # push initial phase-samples
    x = numpy.random.randn(pow(10, decades+2))
    for x in range(pow(10, decades)):
        dev_rt.add_phase(x)
    # print "init"

    # now time how long it takes to push 1000 points
    t = time.time()
    for x in range(pow(10, decades), int(pow(10, decades)+n_add)):
        dev_rt.add_phase(x)
    t = time.time()-t
    # print len(astr.x)
    return (len(afs), t, n_add)


print('DEV\tn_AFs\tms/point\tus/(AF*point)')
for n_pts in [2, 4, 8, 16, 32, 64, 128]:
    (n_afs, t, n_add) = bench(4, n_pts, function=at.realtime.tdev_realtime)
    # print "done"
    print('%s\t%d\t%.3f\t%.3f' %
          ('TDEV', n_afs, 1e3*t/n_add, 1e6*t/(n_afs*n_add)))

print('DEV\tn_AFs\tms/point\tus/(AF*point)')
for n_pts in [2, 4, 8, 16, 32, 64, 128]:
    (n_afs, t, n_add) = bench(4, n_pts, function=at.realtime.ohdev_realtime)
    # print "done"
    print('%s\t%d\t%.3f\t%.3f' %
          ('OHDEV', n_afs, 1e3*t/n_add, 1e6*t/(n_afs*n_add)))

print('DEV\tn_AFs\tms/point\tus/(AF*point)')
for n_pts in [2, 4, 8, 16, 32, 64, 128]:
    (n_afs, t, n_add) = bench(4, n_pts, function=at.realtime.oadev_realtime)
    # print "done"
    print('%s\t%d\t%.3f\t%.3f' %
          ('OADEV', n_afs, 1e3*t/n_add, 1e6*t/(n_afs*n_add)))
