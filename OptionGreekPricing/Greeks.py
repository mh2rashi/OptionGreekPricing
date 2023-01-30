import numpy as np
from math import e
from scipy.stats import norm


def Greek_Delta(Option_type, S, K, sigma, r, T, delta):
    d_1 = (np.log(S / K) + (r - delta + 0.5 * (sigma ** 2)) * T) / (sigma * np.sqrt(T))
    d_2 = d_1 - sigma * np.sqrt(T)

    from scipy.stats import norm

    N_d1 = norm.cdf(d_1)
    N_d2 = norm.cdf(d_2)

    if Option_type == 'C':
        Delta = e ** (-T * delta) * (N_d1)

    else:
        Delta = e ** (-T * delta) * (1 - N_d1)

    return Delta


def Greek_Gamma(S, K, sigma, r, T, delta):
    d_1 = (np.log(S / K) + (r - delta + 0.5 * (sigma ** 2)) * T) / (sigma * np.sqrt(T))
    d_2 = d_1 - sigma * np.sqrt(T)

    from scipy.stats import norm

    Gamma = e ** (-T * delta) * norm.pdf(d_1, 0, 1) / (S * sigma * np.sqrt(T))

    return Gamma


def Greek_Vega(S, K, sigma, r, T, delta):
    d_1 = (np.log(S / K) + (r - delta + 0.5 * (sigma ** 2)) * T) / (sigma * np.sqrt(T))
    d_2 = d_1 - sigma * np.sqrt(T)

    Vega = S * e ** (-T * delta) * norm.pdf(d_1, 0, 1) * np.sqrt(T)

    return Vega


def Greek_Theta(Option_type, S, K, sigma, r, T, delta):
    d_1 = (np.log(S / K) + (r - delta + 0.5 * (sigma ** 2)) * T) / (sigma * np.sqrt(T))
    d_2 = d_1 - sigma * np.sqrt(T)

    if Option_type == 'C':
        theta_call = -S * norm.pdf(d_1, 0, 1) * sigma / (2 * np.sqrt(T)) - r * (K * e ** (-r * T)) * norm.cdf(d_2)
        return theta_call

    else:
        theta_put = -S * norm.pdf(d_1, 0, 1) * sigma / (2 * np.sqrt(T)) + r * (K * e ** (-r * T)) * norm.cdf(-d_2)
        return theta_put


def Greek_Rho(Option_type, S, K, sigma, r, T, delta):
    d_1 = (np.log(S / K) + (r - delta + 0.5 * (sigma ** 2)) * T) / (sigma * np.sqrt(T))
    d_2 = d_1 - sigma * np.sqrt(T)

    if Option_type == 'C':
        rho_call = T * (K * e ** (-r * T)) * norm.cdf(d_2)
        return rho_call

    else:
        rho_put = -T * (K * e ** (-r * T)) * norm.cdf(-d_2)
        return rho_put


def Greek_Psi(Option_type, S, K, sigma, r, T, delta):
    d_1 = (np.log(S / K) + (r - delta + 0.5 * (sigma ** 2)) * T) / (sigma * np.sqrt(T))
    d_2 = d_1 - sigma * np.sqrt(T)

    if Option_type == 'C':
        psi_call = -T * (S * e ** (-delta * T)) * norm.cdf(d_1)
        return psi_call

    else:
        psi_put = T * (S * e ** (-delta * T)) * norm.cdf(-d_1)
        return psi_put
