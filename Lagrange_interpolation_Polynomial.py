#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 15:15:41 2021

@author: tetuko
"""
import matplotlib.pyplot as plt


def Lagrange(x, y, n, xi):
    """
    Evaluates the Lagrange interpolating polynomial of n data points at xi
    x[] - x-coordinates of data points
    y[] - y-coordinates of data points
    n - number of data points
    xi - interpolation argument
    """
    yi = 0e0
    for i in range(n):
        p = 1e0
        for j in range(n):
            if (j != i):
                p *= (xi - x[j])/(x[i] - x[j])  # the lagrange polynomial
        yi += p * y[i]

    return yi


def Lagrange1(x, y, n, xi, ni):
    """
    (It doesn't work)
    Evaluates the Lagrange interpolating polynomial of n data points on a mesh
    of ni interpolation points
    x[] - x-coordinates of data points
    y[] - y-coordinates of data points
    n  - number of data points
    xi[] - interpolation arguments
    yi[] - interpolant values (output)
    ni - number of interpolation points
    """
    yf = [[0]*(n+1) for i in range(n+1)]  # factors of the interpolant
    # prepare invariant factors of the interpolant
    for i in range(1, n+1):
        p = 1e0
        for j in range(1, n+1):
            if (j != i):
                p *= (x[i] - x[j])
        yf[i] = y[i] / p
    # loop over interpolation points
    for k in range(1,ni+1):
        xk = xi[k]
        yk = 0e0
        for i in range(1, n+1):
            p = 1e0
            for j in range(1,n+1):
                if (j != i): p *= (xk - x[j])
            yk += p * yf[i]
        yi[k] = yk
    return yi

ni = 100

x1 = range(1, 10)  # the last is not included so it is 1 to 10 for range(1,11)
x1 = [0.15, 1.15, 2.15, 3.15, 4.15, 5.15, 6.15, 7.15, 8.15, 9.15]
x1 =[0.15, 0.2, 0.3, 0.5, 0.8, 1.1, 1.4, 1.7]

y1 = []
for xk in x1:
    y1.append(1/xk)

# for individual calculation
yi = Lagrange(x1, y1, len(x1), 0.15)
print(yi)

y_interpolate = []
h = (x1[-1]-x1[0])/(ni-1)

x_interpolate = []
for k in range(ni+1):
    xi = x1[0]+k*h
    x_interpolate.append(xi)
    y_interpolate.append(Lagrange(x1, y1, len(x1), xi))

fig, ax = plt.subplots(dpi=300)
ax.plot(x1, y1, ls='none', marker='o', label='data points')
ax.plot(x_interpolate, y_interpolate, label='interpolation result')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.legend(loc='best')
# try1 = Lagrange1(m, k, len(m)-1, yi, 3)
