# OptionGreekPricing

## Table of Contents

- [Overview](#overview)
- [Setup](#setup)
- [Usage](#usage)


## Overview

Current Version: **1.0.3**

A package written in Python with equations to calculate option Prices and Greeks  using the following methods: 

1. Put-Call Parity. The `PutCall` object contains the following equations:
   - Option price for a non-dividend stock.
   - Option price for a discrete-dividend stock.
   - Option price for a continuous-dividend stock.
   - Option price for a currency option

2. Binomial Trees. The `BinomialTress` object contains the following equations:
   - Option price for a one period non-dividend stock.
   - Option price for a multi-period non-dividend stock european option.
   - Option price for a multi-period non-dividend stock american option.
   - Option price for a multi-period currency exchange rate european option.
   - Option price for a multi-period currency exchange rate american option.
   - Option price for a multi-period futures contract american option.
   
3. Black-Scholes. The `BlackScholes` object contains the following equations:
   - Option price for a common stock.
   - Option price for a common stock with discrete dividends.
   - Option price for a currency exchange rate.
   - Option price for a futures contract.
   
4. Greeks. The `Greeks` object contains the following equations for both call and put options:
   - Delta.
   - Gamma.
   - Vega.
   - Theta.
   - Rho.
   - Psi.
   
5. Monte Carlo. The `MonteCarlo` object contains the following equation:
   - Option price for a risk-neutral european option.
   
## Setup

**Setup - PyPi Install:**

The project can be found at PyPI, if you'd like to view the project please use this
[link](https://pypi.org/project/python-trading-robot/). To **install** the library,
run the following command from the terminal.

```bash
pip install option-greek-pricing
```


## Usage

For more detailed examples and explanations, please look at the following file: [Option Greeks & Pricing.ipynb](https://github.com/mh2rashi/OptionGreekPricing/blob/main/Option%20Greeks%20%26%20Pricing.ipynb)
