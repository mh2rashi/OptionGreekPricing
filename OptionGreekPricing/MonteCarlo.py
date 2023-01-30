from math import e, sqrt
from random import gauss


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
