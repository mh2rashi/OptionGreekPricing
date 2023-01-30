import numpy as np
from math import e, factorial


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
