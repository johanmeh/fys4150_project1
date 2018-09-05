#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 16:18:42 2018

@author: maritkollstuen
"""


"""
Solution to task 1 b) of assignment 1 in FYS4150

Function that uses gaussian elimination (forward and backward substitution) to 
create a numerical solution to the second derivative of the function u(x) and 
compare it to the analytical solution. 

Input: 
    n = size of the quadratic matrix
    a = value of the elements of the below-leading diagonal
    b = value of the elements of the leading diagonal
    c = value of the elements of the above-leading diagonal
"""


import numpy as np
import matplotlib.pyplot as plt

def f(x): # Source term of the differential equation 
    return 100*np.exp(-10*x)


def u(x): # Closed - form solution of the differential equation 
    return 1-(1-np.exp(-10))*x - np.exp(-10*x)


def gauss_elim (a, b, c, n):
    x = np.linspace(0, 1, (n+2))
    h = 1/(n+1)
    f_vec = h**2*f(x)
    
    f_vec[0] = 0
    f_vec[-1]= 0
     
    
    b_vec = np.ones(n)*b
    
    # Forward substitution 
    for i in range(0,n-1):
        coef = a/b_vec[i]
        b_vec[i+1] = b_vec[i+1] - coef*c
        f_vec[i+2] = f_vec[i+2] - coef*f_vec[i+1]
    
    # Backward substitution 
    for i in range(n,0,-1):
        f_vec[i-1] = f_vec[i-1] - (c/b_vec[i-1]) * f_vec[i]
    
    v_vec = np.zeros(n+2)   
    
    for i in range(1, n+1):
        v_vec[i] = f_vec[i]/b_vec[i-1]
        
    return (v_vec)

def plot_stuff(a, b, c, n):
    x = np.linspace(0, 1, (n+2))
    u_exact = u(x)
    
    v  = gauss_elim(a, b, c, n)
    plt.figure()
    plt.plot(x, u_exact, color='blue', label = 'u - analytic')
    plt.plot(x, v, color='green', label = 'v - numeric, n = {}'.format(n))
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.savefig('size_{}.png'.format(n))
    

if __name__ == "__main__":
    
    for n in [10, 100, 1000]:
        plot_stuff(-1, 2, -1, n)

