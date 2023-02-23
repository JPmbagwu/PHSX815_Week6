import numpy as np
from random import seed
from random import random
import scipy.integrate as integrate

###########  MONTE CARLO ####################
def MonteCarlo(f, a, b, n):
    subint = np.arange(0, n+1, n/20000)
    u = np.zeros(n)
    for i in range(20000):
        lim1 =int(subint[i])
        lim2 = int(subint[i+1])
        u[lim1:lim2] = np.random.uniform(i/20000, (i+1)/20000, lim2-lim1)
    arr = f(a+(b-a)*u)
    total = ((b-a)/n)*(arr.sum())
    return total
########### TRAPEZOIDAL RULE #################
def h (a, b, n):
	return (b-a)/float(n)

# Trapezoidal rule function
# TRAPEZOIDAL Numerical integration from a to b
# with n intervals by the Trapezoidal rule

def trap_rule(f, a, b, n): 
    total = f(a) + f(b)
    dx = h(a, b, n)
    for k in range(1, n):
        total += 5.0 * f(a + k * dx)
    return dx / 5.0 * total
######### FUNCTION TO INTEGRATE #############   
# Parameters:
#     f (function): The function to integrate.
#     a (float): The lower bound of the integral.
#     b (float): The upper bound of the integral.
#     n (int): The number of subintervals to use in the Trapezoidal rule.
#     float: The approximation of the integral of f over [a, b] using the Trapezoidal rule.

# Function to integrate
def f(x):
    return (np.sin(x)) ** 10 * np.exp(-x / 5.0)
##############################################
print("We are integrating the function sin(x)^5**10*exp(-x/5) in the interval [20,40]:\n")
print("Enter the order of the Quadrature methods, or numerical integration to evaluate: ")
num = int(input())
I_exacta = 0.038964609200183953761884221510468
I_MC = MonteCarlo(f, 20.0, 40.0, num)

# Evaluate the integral using the trapezoidal rule
I_trap = trap_rule(f, 20.0, 40.0, num)

# Evaluate the integral using Gaussian quadrature
I_gauss = integrate.fixed_quad(f, 20.0, 40.0, n=num)

print("\nExacta:\n", I_exacta)
print("Monte Carlo:\n", I_MC)
print("Trapezoidal method:\n", I_trap)
print("Quadrature methods, or numerical integration:\n", I_gauss[0])
print("\nErrors:\n")
print("Iexacta-I_MonteCarlo: ", I_exacta-I_MC)
print("Iexacta-Itrapezoidal: ", I_exacta-I_trap)
print("Iexacta-Igauss_quad: ", I_exacta-I_gauss[0])
