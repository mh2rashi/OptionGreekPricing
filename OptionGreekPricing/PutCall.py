import numpy as np
from math import e


def NonDividend_Stock(S, K, T, i=None, C=None, P=None, r=None):
    if r is None:
        r = np.log(1 + i)

    if C is None and P is None:
        return "Please input the Price of a Call (C) or Put (P) option"

    if C is None:
        C = S - K * e ** (-r * T) + P
        return C

    if P is None:
        P = C - S + K * e ** (-r * T)
        return P


def DiscreteDividend_Stock(S, K, T, delta, i=None, C=None, P=None, r=None):
    if r is None:
        r = np.log(1 + i)

    if C is None and P is None:
        return "Please input the Price of a Call (C) or Put (P) option"

    if C is None:
        C = S * e ** (-delta * T) - K * e ** (-r * T) + P
        return C

    if P is None:
        P = C - S * e ** (-delta * T) + K * e ** (-r * T)
        return P


def ContinuousDividend_Stock(S, K, T, delta, i=None, C=None, P=None, r=None):
    if r is None:
        r = np.log(1 + i)

    if C is None and P is None:
        return "Please input the Price of a Call (C) or Put (P) option"

    if C is None:
        C = S * e ** (-delta * T) - K * e ** (-r * T) + P
        return C

    if P is None:
        P = C - S * e ** (-delta * T) + K * e ** (-r * T)
        return P


def Currency_Option(x0, K, T, rf, rd, C=None, P=None):
    if C is None and P is None:
        return "Please input the Price of a Call (C) or Put (P) option"

    if C is None:
        C = x0 * e ** (-rf * T) - K * e ** (-rd * T) + P
        return C

    if P is None:
        P = C - x0 * e ** (-rf * T) + K * e ** (-rd * T)
        return P
