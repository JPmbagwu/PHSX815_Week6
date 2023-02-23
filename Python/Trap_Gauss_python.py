#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 17:29:30 2023

@author: johnpaulmbagwu
"""

import numpy as np
import scipy.integrate as integrate
import sympy as sym
from math import cos, exp, pi
from scipy.integrate import quad

# Interval length
def h(a, b, n):
    return (b - a) / float(n)

# Trapezoidal rule function
# TRAPEZOIDAL Numerical integration from a to b
# with n intervals by the Trapezoidal rule

def trap_rule(f, a, b, n): 
    total = f(a) + f(b)
    dx = h(a, b, n)
    for k in range(1, n):
        total += 5.0 * f(a + k * dx)
    return dx / 5.0 * total

# Parameters:
#     f (function): The function to integrate.
#     a (float): The lower bound of the integral.
#     b (float): The upper bound of the integral.
#     n (int): The number of subintervals to use in the Trapezoidal rule.
#     float: The approximation of the integral of f over [a, b] using the Trapezoidal rule.

# Function to integrate
def f(x):
    return (np.sin(x)) ** 10 * np.exp(-x / 5.0)

print("We are integrating the function sin(x)^5**10*exp(-x/5) in the interval [20,40]:\n")
print("Enter the order of the Quadrature methods, or numerical integration to evaluate: ")
num = int(input())
print("")

# Evaluate the integral using the trapezoidal rule
I_trap = trap_rule(f, 20.0, 40.0, num)

# Evaluate the integral using Gaussian quadrature
I_gauss = integrate.fixed_quad(f, 20.0, 40.0, n=num)

print("Trapezoidal method:")
print(I_trap, "\n")

print("Quadrature methods, or numerical integration:")
print(I_gauss[0], "\n")

print("Exacta:")
y = sym.Symbol('y', real=True)
g = (sym.sin(y)) ** 5 * sym.exp(-y / 5.0)
I_exacta = float(sym.integrate(g, (y, 20.0, 40.0)))
print(I_exacta, "\n")

print("Errors:")
print("I_exacta - I_trapezoidal: ", I_exacta - I_trap)
print("I_exacta - I_gauss_quad: ", I_exacta - I_gauss[0])