from OptionGreekPricing import PutCall
from OptionGreekPricing import BinomialTrees
from OptionGreekPricing import BlackScholes
from OptionGreekPricing import Greeks
from OptionGreekPricing import MonteCarlo

# Put Call
# print(PutCall.NonDividend_Stock(S=40, K=45, T=0.75, C=2.84, r=0.05))
# print(PutCall.DiscreteDividend_Stock(40, 50, 1, 0.02, C=2.34, r=0.08))
# print(PutCall.ContinuousDividend_Stock(S=40, K=50, T=1, delta=0.02, C=2.34, r=0.08))
# print(PutCall.Currency_Option(x0=1.4, K=1.5, T=0.75, rf=0.08, rd=0.05, C=.0223))

# Binomial Trees
# print(BinomialTrees.NonDivStock_OnePeriod(Option_type='C', K=55, S=50, S_u=60, S_d=40, r=0.05))
# print(BinomialTrees.NonDivStock_MultiPeriod_European(Option_type='P', K=160, S=150, u=1.3, d=0.7, n=2, t=0.5, r=0.06))
# print(BinomialTrees.NonDivStock_MultiPeriod_American(Option_type='P', K=160, S=150, u=1.3, d=0.7, n=2, t=0.5, r=0.06))
# print(BinomialTrees.NonDivStock_MultiPeriod_American(Option_type='P', K=100, S=100, u=1.1, d=1/1.1, n=3, t=1, r=0.06))
# print(BinomialTrees.Currency_MultiPeriod_European(Option_type='C', K=1.25, E=1.15, u=1.053903, d=0.95361, n=2, t=0.5, r_f=0.04, r_d=0.05))
# print(BinomialTrees.Currency_MultiPeriod_American(Option_type='C', K=1.25, E=1.15, u=1.053903, d=0.95361, n=2, t=0.5, r_f=0.04, r_d=0.05))
# print(BinomialTrees.Futures_MultiPeriod_American(Option_type='C', K=60, F=60, u=1.090463, d=0.917042, n=3, t=3/12, r=0.05))
# print(BinomialTrees.Futures_MultiPeriod_American(Option_type='P', K=640, F=650, u=1.193365, d=0.837967, n=2, t=1, r=0.04))

# Black Scholes
# print(BlackScholes.BlackScholes_Stock("C", 50, 52, 0.4, 0.08, 3/12, 0.04))
# print(BlackScholes.BlackScholes_Stock("P", 50, 52, 0.4, 0.08, 3/12, 0.04))
# print(BlackScholes.BlackScholes_Discrete_Stock('P',42, 40, 0.3, 0.75, 3/12, 3/12, 6/12, 0.04))
# print(BlackScholes.BlackScholes_Currency("C",0.009, 0.010, 0.05, 0.04, 1, 0.02))
# print(BlackScholes.BlackScholes_Futures('C', 10, 10, 0.25, 0.04, 1))

# Greeks
# print(Greeks.Greek_Delta('C', 45, 43, 0.2, 0.05, 3/12, 0.02))
# print(Greeks.Greek_Delta('P', 45, 43, 0.2, 0.05, 3/12, 0.02))
# print(Greeks.Greek_Gamma(45, 43, 0.2, 0.05, 3/12, 0.02))
# print(Greeks.Greek_Vega(45, 43, 0.2, 0.05, 3/12, 0.02))
# print(Greeks.Greek_Theta('C',45, 43, 0.2, 0.05, 3/12, 0.02))
# print(Greeks.Greek_Rho('C',45, 43, 0.2, 0.05, 3/12, 0.02))
# print(Greeks.Greek_Psi('C',45, 43, 0.2, 0.05, 3/12, 0.02))

# Monte Carlo
# print(MonteCarlo.MC_riskneutral_European('C', 857.29 , 860, 0.2076 , 0.0014 , 9/12, 0))







