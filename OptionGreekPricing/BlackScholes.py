import numpy as np
from math import e
from scipy.stats import norm


def BlackScholes_Stock(Option_type, S, K, sigma, r, T, delta):
    d_1 = (np.log(S / K) + (r - delta + 0.5 * (sigma ** 2)) * T) / (sigma * np.sqrt(T))
    d_2 = d_1 - sigma * np.sqrt(T)

    N_d1 = norm.cdf(d_1)
    N_d2 = norm.cdf(d_2)

    if Option_type == 'C':
        premium = (S * e ** (-T * delta)) * (N_d1) - (K * e ** (-r * T)) * (N_d2)

    else:
        premium = (K * e ** (-r * T)) * (1 - N_d2) - (S * e ** (-T * delta)) * (1 - N_d1)

    return premium


def BlackScholes_Discrete_Stock(Option_type, S, K, sigma, D, t, D_rate, T, r):
    def PV_Divs(D, t, D_rate, T, r):
        import math
        e = math.e
        remaining_divs = int((T - t) // (D_rate))

        present_value = 0
        for i in range(remaining_divs + 1):
            present_value += D * e ** (-r * (t + i * D_rate))

        return present_value

    Prepaid_Forward_Price = S - PV_Divs(D, t, D_rate, T, r)

    d_1 = (np.log(Prepaid_Forward_Price / K) + (r + 0.5 * (sigma ** 2)) * T) / (sigma * np.sqrt(T))
    d_2 = d_1 - sigma * np.sqrt(T)

    from scipy.stats import norm

    N_d1 = norm.cdf(d_1)
    N_d2 = norm.cdf(d_2)

    if Option_type == 'C':
        premium = Prepaid_Forward_Price * (N_d1) - (K * e ** (-r * T)) * (N_d2)

    else:
        premium = (K * e ** (-r * T)) * (1 - N_d2) - Prepaid_Forward_Price * (1 - N_d1)

    return premium


def BlackScholes_Currency(Option_type, x, K, sigma, r, T, r_f):
    d_1 = (np.log(x / K) + (r - r_f + 0.5 * (sigma ** 2)) * T) / (sigma * np.sqrt(T))
    d_2 = d_1 - sigma * np.sqrt(T)

    from scipy.stats import norm

    N_d1 = norm.cdf(d_1)
    N_d2 = norm.cdf(d_2)

    if Option_type == 'C':
        premium = x * (e ** -r_f * T) * (N_d1) - (K * e ** (-r * T)) * (N_d2)

    else:
        premium = (K * e ** (-r * T)) * (1 - N_d2) - x * (e ** -r_f * T) * (1 - N_d1)

    return premium


def BlackScholes_Futures(Option_type, F, K, sigma, r, T):
    d_1 = (np.log(F / K) + (0.5 * (sigma ** 2)) * T) / (sigma * np.sqrt(T))
    d_2 = d_1 - sigma * np.sqrt(T)

    from scipy.stats import norm

    N_d1 = norm.cdf(d_1)
    N_d2 = norm.cdf(d_2)

    if Option_type == 'C':
        premium = F * (e ** -r * T) * (N_d1) - (K * e ** (-r * T)) * (N_d2)

    else:
        premium = (K * e ** (-r * T)) * (1 - N_d2) - F * (e ** -r * T) * (1 - N_d1)

    return premium
