"""
    This file is part of AllanTools
    https://github.com/aewallin/allantools

    TDEV and MTIE masks for ITU Standards
    AW2017-2018

    ITU Primary Reference Clock
    https://www.itu.int/rec/T-REC-G.811/en

    ITU Primary Reference Time Clock
    https://www.itu.int/rec/T-REC-G.8272/en

    ITU Enhanced Primary Reference Time Clock
    https://www.itu.int/rec/T-REC-G.8272.1/en
"""


def prc_tdev(tau):
    """ Time Deviation, ITU Primary Reference Clock
        https://www.itu.int/rec/T-REC-G.811/en
        Recommendation G.811 (1997) Amendment 1 (04/16)
        Approved in 2016-04-13
    """
    if tau < 100:
        return 3e-9
    if tau < 1000:
        return 0.03*1.0e-9*tau
    else:
        return 30e-9


def prc_mtie(tau):
    """ MTIE, ITU Primary Reference Clock
        https://www.itu.int/rec/T-REC-G.811/en
        Recommendation G.811 (1997) Amendment 1 (04/16)
        Approved in 2016-04-13
    """
    if tau < 1000:
        return 1e-6 * (tau*0.275e-3 + 0.025)
    else:
        return 1e-6 * (tau*1.0e-5 + 0.29)


def eprtc_tdev(tau):
    """ Time Deviation, ITU Enhanced Primary Reference Time Clock
        https://www.itu.int/rec/T-REC-G.8272.1/en
        Recommendation G.8272.1/Y.1367.1 (2016) Amendment 1 (08/17)
        Approved in 2017-08-13

        Table 2 TDEV, page 4
    """
    if tau < 30e3:
        return 1e-9
    if tau < 300e3:
        return 3.33333e-5*1.0e-9*tau
    else:
        return 10e-9


def eprtc_mtie(tau):
    """ MTIE
        ITU Enhanced Primary Reference Time Clock
        https://www.itu.int/rec/T-REC-G.8272.1/en
        Recommendation G.8272.1/Y.1367.1 (2016) Amendment 1 (08/17)
        Approved in 2017-08-13

        Table 1, page 3
    """
    if tau <= 1:
        return 4e-9
    if tau <= 100:
        return 1e-9*(0.11114*tau+3.89)
    if tau <= 400000:
        return 1e-9*(0.0375e-3*tau+15)
    else:
        return 1e-9*30  # 30 ns


def prtcA_tdev(tau):
    """ Time Deviation
        ITU Primary Reference Time Clock
        https://www.itu.int/rec/T-REC-G.8272/en
        Recommendation G.8272/Y.1367 (11/18)
        Approved in 2018-11-29

        Table 3
    """
    if tau < 100:
        return 3e-9
    if tau < 1000:
        return 0.03e-9*tau
    else:
        return 30e-9


def prtcB_tdev(tau):
    """ Time Deviation
        ITU Primary Reference Time Clock
        https://www.itu.int/rec/T-REC-G.8272/en
        Recommendation G.8272/Y.1367 (11/18)
        Approved in 2018-11-29

        Table 4
    """
    if tau < 100:
        return 1e-9
    if tau < 500:
        return 0.01e-9*tau
    else:
        return 5e-9


def prtcA_mtie(tau):
    """ MTIE
        ITU Primary Reference Time Clock
        https://www.itu.int/rec/T-REC-G.8272/en

        Table 1
    """
    if tau < 273:
        return 1e-6 * (0.275e-3*tau + 0.025)
    else:
        return 1e-6 * 0.1  # 100ns


def prtcB_mtie(tau):
    """ MTIE
        ITU Primary Reference Time Clock
        https://www.itu.int/rec/T-REC-G.8272/en

        Table 2
    """
    if tau < 54.5:
        return 1e-6 * (0.275e-3*tau + 0.025)
    else:
        return 1e-6 * 0.04  # 40ns
