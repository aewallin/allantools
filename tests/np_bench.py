#!/usr/bin/env python

import numpy as np
import pylab as plt
import allantools.allantools_pure_python as alt
import allantools.allantools as alp
import time


def bench_function(p_func, np_func, N = 10000, noise_func = np.random.random ):
    """ benchmark p_func against np_func
        we do this by timing the functions on an N-sized random dataset
        we then check that both funtions return the same results
        and we return the cpu-seconds for each function.
    """
    taus = [1, 3, 5, 16, 128]    
    rate = 2.1
    data = noise_func(N)
    t1 = time.time()
    o_taus, o_dev, o_err, o_n = p_func(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    o_taus_, o_dev_, o_err_, o_n_ = np_func(data, rate, taus)
    t4 = time.time()

    assert np.allclose(o_taus, o_taus_)
    assert np.allclose(o_dev, o_dev_)
    assert np.allclose(o_err, o_err_)
    t_python = (t2 - t1)
    t_numpy = (t4 - t3)
    return (N, t_python, t_numpy)

def benchmark_plot(func1, func2, name, max_data_log = 7, noise_func = np.random.random):
    """ benchmark func1 against func2 with synthetic data of increasing
        size N up to log(N)=max_data_log
        plot a log-log graph
    """
    print "\ntesting %s" % name
    print "N \t pure_python \t numpy   \t speedup "
    t_p_times = []
    t_np_times = []
    Ns = []
    mean_speedup=0
    speedup_n=0
    N_big = 1e5
    n=0
    n_max=30
    for N in np.logspace(2,max_data_log,n_max):
        (N, t_p, t_np) = bench_function(func1, func2, N, noise_func)
        print "%02d/%02d: %d \t %2.3fs \t %2.3fs \t %2.2fx " % (n+1,n_max, N, t_p, t_np, (t_p / t_np) )
        n=n+1
        t_p_times.append(t_p)
        t_np_times.append(t_np)
        Ns.append(N)
        if N>N_big:
            speedup_n=speedup_n+1
            mean_speedup = mean_speedup + t_p/t_np
    if speedup_n==0:
        speedup_n = 1
    mean_speedup = mean_speedup / speedup_n
    plt.figure()
    plt.loglog(Ns, t_p_times,'ro',label='pure Python %s' % name )
    plt.loglog(Ns, t_np_times,'bs',label='numpy %s' % name )
    plt.xlabel('Input data size')
    plt.ylabel('CPU seconds')
    speedup_txt = "mean speedup for N>%.2g : %3.1fx" % (N_big, mean_speedup)
    plt.text( 1000, 1, speedup_txt )
    plt.legend(loc='upper left')
    plt.title('allantools numpy benchmark, AW 2014-08-30')
    plt.show()
    
def brownian_noise(N):
    """ Brownian or random walk (diffusion) noise with 1/f^2 PSD
        not really a color... rather Brownian or random-walk
    """
    return np.cumsum(np.random.randn(N))

if __name__ == "__main__":
    #benchmark_plot( alt.adev, alp.adev, "ADEV",6)
    #benchmark_plot( alt.oadev, alp.oadev, "OADEV",6)
    #benchmark_plot( alt.mdev, alp.mdev, "MDEV",6)
    #benchmark_plot( alt.tdev, alp.tdev, "TDEV",6)
    #benchmark_plot( alt.hdev, alp.hdev, "HDEV",5)
    #benchmark_plot( alt.ohdev, alp.ohdev, "OHDEV",5)
    #benchmark_plot( alt.totdev, alp.totdev, "TOTDEV",5)
    benchmark_plot( alt.mtie, alp.mtie, "MTIE",7, brownian_noise)
    #benchmark_plot( alt.tierms, alp.tierms, "TIERMS",7, brownian_noise)
    pass
    


