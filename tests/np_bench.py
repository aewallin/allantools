#!/usr/bin/env python

""" 
    Benchmarking for allantools (https://github.com/aewallin/allantools)
    
    First version: AW 2014-08-31
"""

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

def benchmark_run(func1, func2, name, max_data_log = 7, noise_func = np.random.random):
    """ benchmark func1 against func2 with synthetic data of increasing
        size N up to log(N)=max_data_log
    """
    print "\nBenchmarking: %s" % name
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
        if N>N_big: # numpy times are meaningless for short datasets,
                    # since time.time() seems to have a 'noise floor' of about 2 ms (?)
            speedup_n=speedup_n+1
            mean_speedup = mean_speedup + t_p/t_np
    if speedup_n==0:
        speedup_n = 1
    mean_speedup = mean_speedup / speedup_n
    return (Ns, t_p_times, t_np_times, mean_speedup, name)

def benchmark_plot(data):
    """ log-log graph of the benchmark results """
    plt.figure()
    idx = 0
    symbols = ['o', 's', 'h', 'D','1','8','*','+','x']
    for (N, tp, tnp, speedup, name) in data:
        plt.loglog(N, tp,'r%s' % symbols[idx],label='pure Python %s' % name )
        plt.loglog(N, tnp,'b%s' % symbols[idx],label='numpy %s' % name )
        idx += 1
        speedup_txt = "%s speedup: %3.1fx" % (name, speedup)
        plt.text( 700, (0.9)/np.exp(0.3*idx), speedup_txt , fontsize=14)
    plt.xlabel('Input data size')
    plt.ylabel('CPU seconds')

    plt.legend(loc='upper left')
    plt.title('allantools numpy benchmark, AW 2014-08-31')
    plt.show()
    
def brownian_noise(N):
    """ Brownian or random walk (diffusion) noise with 1/f^2 PSD
        not really a color... rather Brownian or random-walk
    """
    return np.cumsum(np.random.randn(N))

if __name__ == "__main__":
    N_log_max = 5
    # runs on desktop computer with i7-2600K @ 3.4 GHz CPU:
    # N_log_max   seconds
    # 3           1.1
    # 4           6.8
    # 5           45.7
    # 6           359.9
    data=[]
    # run the benchmarks and store results into one list of tuples
    t0 = time.time()
    #data.append( benchmark_run( alt.adev  , alp.adev  , "ADEV"  ,N_log_max) )
    #data.append( benchmark_run( alt.oadev , alp.oadev , "OADEV" ,N_log_max) )
    #data.append( benchmark_run( alt.mdev  , alp.mdev  , "MDEV"  ,N_log_max) )
    #data.append( benchmark_run( alt.tdev  , alp.tdev  , "TDEV"  ,N_log_max) )
    #data.append( benchmark_run( alt.hdev  , alp.hdev  , "HDEV"  ,N_log_max) )
    #data.append( benchmark_run( alt.ohdev , alp.ohdev , "OHDEV" ,N_log_max) )
    #data.append( benchmark_run( alt.totdev, alp.totdev, "TOTDEV",N_log_max) )
    data.append( benchmark_run( alt.mtie  , alp.mtie  , "MTIE"  ,N_log_max) )
    #data.append( benchmark_run( alt.tierms, alp.tierms, "TIERMS",N_log_max) )
    t1 = time.time()
    print "Benchmarks done in %.1f seconds" % (t1-t0)
    # log-log plot of all benchmark data
    benchmark_plot( data )
    pass
    


