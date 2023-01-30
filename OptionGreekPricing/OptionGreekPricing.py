import numpy as np
from math import e, factorial, sqrt
from scipy.stats import norm
from random import gauss


class PutCall:
    def __init__(self):
        pass

    @staticmethod
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

    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def Currency_Option(x0, K, T, rf, rd, C=None, P=None):

        if C is None and P is None:
            return "Please input the Price of a Call (C) or Put (P) option"

        if C is None:
            C = x0 * e ** (-rf * T) - K * e ** (-rd * T) + P
            return C

        if P is None:
            P = C - x0 * e ** (-rf * T) + K * e ** (-rd * T)
            return P


class BinomialTrees:
    pass

    @staticmethod
    def NonDivStock_OnePeriod(Option_type, K, u=None, d=None, p=None, S=None, S_u=None, S_d=None, i=None, r=None):
        if r is None:
            r = np.log(1 + i)

        if u is None and d is None:
            p = (S / (e ** -r) - S_d) / (S_u - S_d)

        if S_u is None and S_d is None:
            p = (1 / (e ** -r) - d) / (u - d)
            S_u = S * u
            S_d = S * d

        if Option_type == 'C':
            Option_premium = (e ** -r) * ((max(S_u - K, 0)) * p + (1 - p) * (max(S_d - K, 0)))
            return Option_premium
        else:
            Option_premium = (e ** -r) * ((max(K - S_u, 0)) * p + (1 - p) * (max(K - S_d, 0)))
            return Option_premium

    @staticmethod
    def NonDivStock_MultiPeriod_European(Option_type, K, S, u, d, n, t, r=None, i=None):

        def fact(y):
            return factorial(y)

        if r is None:
            r = np.log(1 + i)

        p = (e ** (r * t / n) - d) / (u - d)

        option_premium = 0

        if Option_type == 'P':
            for x in range(n + 1):
                S_t = (u ** x) * (d ** (n - x)) * S
                payoff = max(K - S_t, 0)
                expected_value = fact(n) / (fact(n - x) * fact(x)) * (p ** x) * ((1 - p) ** (n - x)) * payoff
                option_premium += expected_value * e ** (-r * t)

        if Option_type == 'C':
            for x in range(n + 1):
                S_t = (u ** x) * (d ** (n - x)) * S
                payoff = max(S_t - K, 0)
                expected_value = fact(n) / (fact(n - x) * fact(x)) * (p ** x) * ((1 - p) ** (n - x)) * payoff
                option_premium += expected_value * e ** (-r * t)

        return option_premium

    @staticmethod
    def NonDivStock_MultiPeriod_American(Option_type, K, S, u, d, n, t, r=None, i=None):

        def fact(y):
            return factorial(y)

        if r is None:
            r = np.log(1 + i)

        p = ((e ** (r * t / n)) - d) / (u - d)

        option_premium = [0 for i in range(n + 1)]

        if Option_type == 'P':
            for x in range(n + 1):
                S_t = (u ** x) * (d ** (n - x)) * S
                payoff = max(K - S_t, 0)
                option_premium[x] = payoff

            for i in range(n - 1, -1, -1):
                for j in range(0, i + 1):
                    S_t = (u ** j) * (d ** (i - j)) * S
                    option_premium[j] = (option_premium[j + 1] * (p) + (1 - p) * option_premium[j]) * e ** (-r * t / n)
                    option_premium[j] = max(option_premium[j], K - S_t)

        if Option_type == 'C':
            for x in range(n + 1):
                S_t = (u ** x) * (d ** (n - x)) * S
                payoff = max(S_t - K, 0)
                option_premium[x] = payoff

            for i in range(n - 1, -1, -1):
                for j in range(0, i + 1):
                    S_t = (u ** j) * (d ** (i - j)) * S
                    option_premium[j] = (option_premium[j + 1] * (p) + (1 - p) * option_premium[j]) * e ** (-r * t / n)
                    option_premium[j] = max(option_premium[j], S_t - K)

        return option_premium[0]

    @staticmethod
    def Currency_MultiPeriod_European(Option_type, K, E, u, d, n, t, r_f, r_d):

        def fact(y):
            return factorial(y)

        p = (e ** ((r_d - r_f) * t / n) - d) / (u - d)

        option_premium = 0

        if Option_type == 'P':
            for x in range(n + 1):
                E_t = (u ** x) * (d ** (n - x)) * E
                payoff = max(K - E_t, 0)
                expected_value = fact(n) / (fact(n - x) * fact(x)) * (p ** x) * ((1 - p) ** (n - x)) * payoff
                option_premium += expected_value * e ** (-r_d * t)

        if Option_type == 'C':
            for x in range(n + 1):
                E_t = (u ** x) * (d ** (n - x)) * E
                payoff = max(E_t - K, 0)
                expected_value = fact(n) / (fact(n - x) * fact(x)) * (p ** x) * ((1 - p) ** (n - x)) * payoff
                option_premium += expected_value * e ** (-r_d * t)

        return option_premium

    @staticmethod
    def Currency_MultiPeriod_American(Option_type, K, E, u, d, n, t, r_f, r_d):

        def fact(y):
            return factorial(y)

        p = (e ** ((r_d - r_f) * t / n) - d) / (u - d)

        option_premium = [0 for i in range(n + 1)]

        if Option_type == 'P':
            for x in range(n + 1):
                E_t = (u ** x) * (d ** (n - x)) * E
                payoff = max(K - E_t, 0)
                option_premium[x] = payoff

            for i in range(n - 1, -1, -1):
                for j in range(0, i + 1):
                    E_t = (u ** j) * (d ** (i - j)) * E
                    option_premium[j] = (option_premium[j + 1] * (p) + (1 - p) * option_premium[j]) * e ** (
                            -r_d * t / n)
                    option_premium[j] = max(option_premium[j], K - E_t)

        if Option_type == 'C':
            for x in range(n + 1):
                E_t = (u ** x) * (d ** (n - x)) * E
                payoff = max(E_t - K, 0)
                option_premium[x] = payoff

            for i in range(n - 1, -1, -1):
                for j in range(0, i + 1):
                    E_t = (u ** j) * (d ** (i - j)) * E
                    option_premium[j] = (option_premium[j + 1] * (p) + (1 - p) * option_premium[j]) * e ** (
                            -r_d * t / n)
                    option_premium[j] = max(option_premium[j], E_t - K)

        return option_premium[0]

    @staticmethod
    def Futures_MultiPeriod_American(Option_type, K, F, u, d, n, t, i=None, r=None):

        def fact(y):
            return factorial(y)

        if r is None:
            r = np.log(1 + i)

        p = (1 - d) / (u - d)

        option_premium = [0 for i in range(n + 1)]

        if Option_type == 'P':
            for x in range(n + 1):
                F_t = (u ** x) * (d ** (n - x)) * F
                payoff = max(K - F_t, 0)
                option_premium[x] = payoff

            for i in range(n - 1, -1, -1):
                for j in range(0, i + 1):
                    F_t = (u ** j) * (d ** (i - j)) * F
                    option_premium[j] = (option_premium[j + 1] * (p) + (1 - p) * option_premium[j]) * e ** (-r * t / n)
                    option_premium[j] = max(option_premium[j], K - F_t)

        if Option_type == 'C':
            for x in range(n + 1):
                F_t = (u ** x) * (d ** (n - x)) * F
                payoff = max(F_t - K, 0)
                option_premium[x] = payoff

            for i in range(n - 1, -1, -1):
                for j in range(0, i + 1):
                    F_t = (u ** j) * (d ** (i - j)) * F
                    option_premium[j] = (option_premium[j + 1] * (p) + (1 - p) * option_premium[j]) * e ** (-r * t / n)
                    option_premium[j] = max(option_premium[j], F_t - K)

        return option_premium[0]


class BlackScholes:
    pass

    @staticmethod
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

    @staticmethod
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

    @staticmethod
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

    @staticmethod
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


class Greeks:
    pass

    @staticmethod
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

    @staticmethod
    def Greek_Gamma(S, K, sigma, r, T, delta):

        d_1 = (np.log(S / K) + (r - delta + 0.5 * (sigma ** 2)) * T) / (sigma * np.sqrt(T))
        d_2 = d_1 - sigma * np.sqrt(T)

        from scipy.stats import norm

        Gamma = e ** (-T * delta) * norm.pdf(d_1, 0, 1) / (S * sigma * np.sqrt(T))

        return Gamma

    @staticmethod
    def Greek_Vega(S, K, sigma, r, T, delta):

        d_1 = (np.log(S / K) + (r - delta + 0.5 * (sigma ** 2)) * T) / (sigma * np.sqrt(T))
        d_2 = d_1 - sigma * np.sqrt(T)

        Vega = S * e ** (-T * delta) * norm.pdf(d_1, 0, 1) * np.sqrt(T)

        return Vega

    @staticmethod
    def Greek_Theta(Option_type, S, K, sigma, r, T, delta):

        d_1 = (np.log(S / K) + (r - delta + 0.5 * (sigma ** 2)) * T) / (sigma * np.sqrt(T))
        d_2 = d_1 - sigma * np.sqrt(T)

        if Option_type == 'C':
            theta_call = -S * norm.pdf(d_1, 0, 1) * sigma / (2 * np.sqrt(T)) - r * (K * e ** (-r * T)) * norm.cdf(d_2)
            return theta_call

        else:
            theta_put = -S * norm.pdf(d_1, 0, 1) * sigma / (2 * np.sqrt(T)) + r * (K * e ** (-r * T)) * norm.cdf(-d_2)
            return theta_put

    @staticmethod
    def Greek_Rho(Option_type, S, K, sigma, r, T, delta):

        d_1 = (np.log(S / K) + (r - delta + 0.5 * (sigma ** 2)) * T) / (sigma * np.sqrt(T))
        d_2 = d_1 - sigma * np.sqrt(T)

        if Option_type == 'C':
            rho_call = T * (K * e ** (-r * T)) * norm.cdf(d_2)
            return rho_call

        else:
            rho_put = -T * (K * e ** (-r * T)) * norm.cdf(-d_2)
            return rho_put

    @staticmethod
    def Greek_Psi(Option_type, S, K, sigma, r, T, delta):

        d_1 = (np.log(S / K) + (r - delta + 0.5 * (sigma ** 2)) * T) / (sigma * np.sqrt(T))
        d_2 = d_1 - sigma * np.sqrt(T)

        if Option_type == 'C':
            psi_call = -T * (S * e ** (-delta * T)) * norm.cdf(d_1)
            return psi_call

        else:
            psi_put = T * (S * e ** (-delta * T)) * norm.cdf(-d_1)
            return psi_put


class MonteCarlo:
    pass

    @staticmethod
    def MC_riskneutral_European(Option_type, S, K, sigma, r, T, delta, simulations=90000):

        # Generating stock prices
        def generate_asset_price(S, r, sigma, T, delta):
            return S * e ** ((r - delta - 0.5 * (sigma ** 2)) * T + sigma * sqrt(T) * gauss(0, 1))

        # Option payoff
        def option_payoff(Option_type, S_T, K):

            if Option_type == 'C':
                payoff = max(S_T - K, 0.0)

            else:
                payoff = max(K - S_T, 0.0)

            return payoff

        # Simulation
        sum_payoffs = 0

        for i in range(simulations):
            S_T = generate_asset_price(S, r, sigma, T, delta)
            sum_payoffs += option_payoff(Option_type, S_T, K)

        price = ((sum_payoffs) / float(simulations)) * e ** (-r * T)

        return price
