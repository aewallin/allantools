"""
 Useful collection of functions for the allantools test-suite
 
"""

import sys
import gzip
import numpy


# read a simple data-file with phase or frequency numbers on each line
def read_datafile(filename):
    p = []
    if filename[-2:]=='gz':
        with gzip.open(filename) as f:
            for line in f:
                if not line.startswith("#"):  # skip comments
                    p.append(float(line))
    else:
        with open(filename) as f:
            for line in f:
                if not line.startswith("#"):  # skip comments
                    p.append(float(line))
    return p


# read a result-file, produced by copy/paste from Stable32
# note: header-lines need to be manually commented-out with "#"
def read_resultfile(filename):
    rows = []
    row = []
    with open(filename) as f:
        for line in f:
            if not line.startswith("#"): # skip comments
                row = []
                l2 = line.split(" ")
                l2 = [_f for _f in l2 if _f]
                for n in range(len(l2)):
                    row.append(float(l2[n]))
                rows.append(row)
    return rows


def read_stable32(resultfile, datarate):
    devresults = read_resultfile(resultfile)
    print("Read ", len(devresults), " rows from ", resultfile)
    taus = []
    devs = []
    ns = []
    # parse textfile produced by Stable32
    for row in devresults:
        if len(row) == 7:  # typical ADEV result file has 7 columns of data
            tau_n = row[0]  # tau in number of datapoints
            tau_s = row[1]  # tau in seconds
            taus.append(tau_n * (1 / float(datarate)))
            n = row[2]  # n averages
            a = row[5]  # deviation
            # Note we don't read the error-columns, since they are mostly zero for 'all tau' runs of Stable32
            devs.append(a)
            ns.append(n)
        elif len(row) == 4:  # the MTIE/TIErms results are formatted slightly differently
            tau_n = row[0]  # tau in number of datapoints
            tau_s = row[1]  # tau in seconds
            taus.append(tau_n * (1 / float(datarate)))
            n = row[2]  # n averages
            a = row[3]  # MTIE or TIErms
            devs.append(a)
            ns.append(n)
    return (taus, devs, ns)


# test a deviation function by:
# - running the function on the datafile
# - reading the correct answers from the resultfile
# - checking that tau, n, and dev are correct 
def test(function, datafile, datarate, resultfile, verbose=0, tolerance=1e-4):
    # if Stable32 results were given with more digits we could decrease tolerance

    phase = read_datafile(datafile)
    print("Read ", len(phase), " phase values from ", datafile)
    (taus, devs, ns) = read_stable32(resultfile, datarate)

    # run allantools algorithm
    (taus2, devs2, errs2, ns2) = function(phase, datarate, taus)
    # check that allantools and Stable32 agree on length of DEV, Tau, and N results

    assert ( len(taus) == len(taus2) )
    assert ( len(devs) == len(devs2) )
    assert ( len(ns) == len(ns2) )

    n_errors = 0  # number of errors/problems detected
    for (t1, a1, n1, t2, a2, n2) in zip(taus, devs, ns, taus2, devs2, ns2):
        # Check that allantools and Stable32 give exactly the same Tau and N results
        errs = check_deviations((t1, a1, n1, t2, a2, n2), tolerance, verbose)
        n_errors = n_errors + errs


    assert ( n_errors == 0)  # no errors allowed!
    print("test of function ", function, " Done. Errors=", n_errors)


def to_fractional(data):
    mu = numpy.mean(data)
    return [(x-mu)/mu for x in data]

# test one tau-value (i.e. one row in the result file) at a time
def test_row_by_row(function, datafile, datarate, resultfile, verbose=0, tolerance=1e-4, normalize=False):
    # if Stable32 results were given with more digits we could decrease tolerance

    phase = read_datafile(datafile)
    if normalize: # convert frequencies in Hz to fractional frequencies
        phase = to_fractional(phase)
        
    print("Read ", len(phase), " values from ", datafile)

    (taus, devs, ns) = read_stable32(resultfile, datarate)

    if verbose:
        print("Tau N  \t DEV(Stable32) \t DEV(allantools) \t relative error")
    # run allantools algorithm, row by row
    for (tau, dev, n) in zip(taus, devs, ns):
        (taus2, devs2, errs2, ns2) = function(phase, datarate, [tau])
        check_deviations((tau, dev, n, taus2[0], devs2[0], ns2[0]), tolerance, verbose)


def check_deviations(xxx_todo_changeme, tolerance, verbose):
    # Check that allantools and Stable32 give exactly the same Tau and N results
    (t1, a1, n1, t2, a2, n2) = xxx_todo_changeme
    n_errors = 0
    try:
        assert ( t1 == t2 )
    except:
        print("ERROR tau1=", t1, " tau2=", t2)
        n_errors = n_errors + 1
    try:
        assert ( n1 == n2 )
    except:
        n_errors = n_errors + 1
        print("ERROR n1=", n1, " n2=", n2, " at t1=", t1, " t2=", t2)

    # check the DEV result, with a given relative error tolerance
    rel_error = (a2 - a1) / a1
    # tolerance = 1e-4 # if Stable32 results were given with more digits we could decrease tol
    try:
        assert ( abs(rel_error) < tolerance )
        if verbose:
            print("%d %d  \t %0.6g \t %0.6g \t %0.6f OK!" % (t1, n1, a1, a2, rel_error))
    except:
        n_errors = n_errors + 1
        print("ERROR %d  %d %0.6g \t %0.6g \t %0.6f" % (t1, n1, a1, a2, rel_error))
        n_errors = n_errors + 1

    return n_errors


if __name__ == "__main__":
    print("nothing to see here.")
    pass

