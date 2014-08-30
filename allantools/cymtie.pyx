from math import sqrt

def calc_mtie_phase(phase, mj):
    """ maximum time interval error
    this seems to correspond to Stable32 setting "Fast(u)"
    Stable32 also has "Decade" and "Octave" modes where the dataset is extended somehow?
    rate = float(rate) """

    cdef int ncount, i, oldindex
    cdef double dev, win_max, win_min
    cdef double newsample, oldsample, tie

    dev = 0
    ncount = 0
    win_max = 0
    win_min = 0

    # we start by finding the max/min in the first window
    # when we slide the window forward by one step, one of these things happen:
    # - neither the new sample that enters the window, nor the old sample that leaves the window
    #    changes the max/min of the window
    # - the new sample that enters the window is the new max/min of the window
    # - the old sample that leaves the window was the max/min of the window
    #    this is the most computationally expensive case, as we now have to search
    #    the entire window for a new max/min
    for i in range(0, len(phase) - mj):  # slide the start of the window over the dataset
        if i == 0:  # initialize window max/min
            win_max = max(phase[0:0 + mj + 1])  # max and min in the first window
            win_min = min(phase[0:0 + mj + 1])
        else:
            newsample = phase[i + mj]  # the new sample that enters the window
            if newsample > win_max:
                win_max = newsample
            elif newsample < win_min:
                win_min = newsample
            oldindex = i - 1
            if oldindex < 0:
                oldindex = 0
            oldsample = phase[i - 1]  # the old sample we throw away
            if win_max == oldsample:
                win_max = max(phase[i:i + mj + 1])  # must search for a new maximum
            if win_min == oldsample:
                win_min = min(phase[i:i + mj + 1])  # must search for new minimum sample

                #phases = phase[i:i+mj+1] # data window of length mj

        tie = win_max - win_min  # largest error in this window
        if tie > dev:
            dev = tie
        ncount += 1

    deverr = dev / sqrt(ncount)

    return dev, deverr, ncount