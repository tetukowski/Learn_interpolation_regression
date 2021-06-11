#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 15:58:22 2021

@author: tetuko
"""

import matplotlib.pyplot as plt
import numpy as np


def LinFit(x, y, sigmy, n, iopt):
    """
    Determines the coefficients a and b of a linear model P(x;a,b) = a * x + b
    and the associated uncertainties, sigma and sigmb, for a set of observed
    data points by "Chi-square" regression

    x[]- x-coordinates of observed data, i = 1,..,n
    y[]- y-coordinates of observed data
    sigmy[] - standard deviations of observed data
    n - number of observed data points
    iopt - iopt == 0 - initializes sigmy[i] = 1 (equal weights)
         - iopt != 0 - uses the received values of sigmy[i]
    a, b
    - parameters of the linear model (output)
    sigma - uncertainties associated with the model parameters (output)
    sigmb
    chi2 - value of Chi-square merit function for the output parameters
"""

    if (iopt == 0):
        for i in range(0, n):
            sigmy[i] = 1e0  # iopt = 0: equall weights

    s = sx = sy = sxx = sxy = 0e0                               # prepare sums
    for i in range(0, n):
        f = 1e0/(sigmy[i]*sigmy[i])
        s += f
        sx += x[i] * f
        sxx += x[i] * x[i] * f
        sy += y[i] * f
        sxy += x[i] * y[i] * f

    f = 1e0/(s*sxx - sx*sx)
    a = (s*sxy - sx*sy) * f
    sigma = (s*f)**0.5           # model parameters
    b = (sy*sxx - sx*sxy) * f
    sigmb = (sxx*f)**0.5         # and uncertainties

    chi2 = 0e0                   # value of Chi-square function

    for i in range(0, n):
        f = (y[i] - a*x[i] - b)/sigmy[i]
        chi2 += f*f

    return (a, b, sigma, sigmb, chi2)


x = [1, 2, 3, 4, 5]
y = [0.8, 2.1, 2.8, 4.0, 4.4]
sigmy = []
for yi in y:
    sigmy.append(yi*0.15)

a, b, sigma, sigmb, chi2 = LinFit(x, y, sigmy, len(x), 0)
print(f'a={a}, b={b}, chi2={chi2}')

fig, ax = plt.subplots(dpi=200)
ax.plot(x, y, ls='none', marker='o', markersize=6, mfc='none',
        label='data')
ax.errorbar(x, y, yerr=sigmy, color = 'k', ls='none')

line_x = np.linspace(x[0], x[-1])
ax.plot(line_x, line_x*a+b, ls='solid',
        label='linear regression')

