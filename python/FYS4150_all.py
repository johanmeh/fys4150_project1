#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 14:05:59 2018

@author: maritkollstuen
"""


"""
Solution to task 1 c) of assignment 1 in FYS4150 (missing CPU time!)

Contains also:
    - solution to 1 d) (missing a log of the relative error + plots are strange!)
    - Solution to 1 4) (Plots are strange!)


Function that uses gaussian elimination (forward and backward substitution) to 
create a numerical solution to the second derivative of the function u(x) and 
compare it to the analytical solution. 

Input: 
    n = size of the quadratic matrix
    a = value of the elements of the above- and below-leading diagonal
    b = value of the elements of the leading diagonal
"""


import numpy as np
import matplotlib.pyplot as plt
import time

def f(x): # Source term of the differential equation 
    return 100*np.exp(-10*x)


def u(x): # Closed - form solution of the differential equation 
    return 1-(1-np.exp(-10))*x - np.exp(-10*x)


def gauss_elim_general (a, b, c, n):
    """
    Yields v_vec for the general matrix (Task 1b)
    
    """
    x = np.linspace(0, 1, (n+2))
    h = 1/(n+1)
    f_vec = h**2*f(x)
    
    f_vec[0] = 0
    f_vec[-1]= 0
     
    
    b_vec = np.ones(n)*b
    
    t1 = time.time()
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
    t2 = time.time()    
    print('Gaussian elimination (general algorithm), n = {}: {}'.format(n, (t2 -t1))) 
    return (v_vec)

def gauss_elim_special (a, b, n): 
    """
    Yields v_vec for our special matrix (Task 1c)
    
    """
    x = np.linspace(0, 1, (n+2))
    h = 1/(n+1)
    f_vec = h**2*f(x)
    
    f_vec[0] = 0
    f_vec[-1]= 0
     
    
    b_vec = np.ones(n)*b
    
    t1 = time.time()
    # Forward substitution 
    for i in range(0,n-1):
        coef = a/b_vec[i]
        b_vec[i+1] = b_vec[i+1] - coef*a
        f_vec[i+2] = f_vec[i+2] - coef*f_vec[i+1]
    
    # Backward substitution 
    for i in range(n,0,-1):
        f_vec[i-1] = f_vec[i-1] - (a/b_vec[i-1]) * f_vec[i]
    
    v_vec = np.zeros(n+2)   
    
    for i in range(1, n+1):
        v_vec[i] = f_vec[i]/b_vec[i-1]
    t2 = time.time()
    print('Gaussian elimination (our special matrix), n = {}: {}'.format(n, (t2 -t1)))    
    return (v_vec)




def relative_error(v ,n):
    x = np.linspace(0, 1, (n+2))
    u_ex = u(x)
    
    epsilon = np.zeros(len(u_ex))
    for i in range(len(x)):
        epsilon[i] = np.log10(np.abs((v[i] - u_ex[i])/u_ex[i]))
    return epsilon, np.nanmax(np.abs(epsilon))



def make_matrix(n, a, b, c): # Make tridiagonal matrix for LU-decomposition
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                A[i, j] = b
        for j in range(n-1):
            if i == (j+1):
                A[j, i] = c
                A[i, j] = a
    return A


def lu_dec(A, n): 
    """
    LU - deconposistion
    """
    x = np.linspace(0, 1, (n+2))
    h = 1/(n+1)
    f_vec = h**2*f(x)
    t1 = time.time()
    v = np.linalg.solve(A, f_vec[1:-1])
    t2 = time.time()
    print('LU-decomposition, n= {}: {}'.format(n, (t2 -t1)))
    return v, f_vec


def plot_stuff(v, v_lu, n):
    x = np.linspace(0, 1, (n+2))
    u_exact = u(x)
    
    
    # Plot exact versus numeric result and result from LU-decomposition

    plt.figure()
    plt.plot(x, u_exact, color='blue', label = 'u - analytic', Linestyle = ':')
    plt.plot(x, v, color='green', label = 'v - numeric, n = {}'.format(n), Linestyle ='--')
    plt.plot(x[1:-1], v_lu, color='red', label = 'v - LU-decomposition, n = {}'.format(n), Linestyle ='-.')
    plt.legend(fontsize = 17)
    plt.xlabel('x', fontsize = 17)
    plt.ylabel('u(x)', fontsize = 17)
    


if __name__ == "__main__":
    
    
    
    for n in [10, 100, 1000, 10**4]: # Don't try LU-decomposition with higher n than this!  
        v_gauss_ge = gauss_elim_general(-1, 2, -1, n)
        v_gauss_sp = gauss_elim_special(-1, 2, n) 
        A = make_matrix(n, -1, 2, -1)
        v_lu, f_vec = lu_dec(A, n)
        plot_stuff(v_gauss_sp, v_lu, n)
    

    sizes = [10, 100, 1000, 10**4, 10**5]#, 10**6, 10**7] # sizes of matrices   
    h = np.zeros(len(sizes))
    max_e = np.zeros(len(sizes))
    for n in sizes:
        i = sizes.index(n)
        h[i] = 1/(n+1)
        v_gauss_sp = gauss_elim_special(-1, 2, n) 
        epsilon, max_e[i] = relative_error(v_gauss_sp, n)
    


    plt.figure()
    plt.semilogx(h, max_e)
    plt.xlabel('Stepsize, h', fontsize = 17)
    plt.ylabel('Max relative error', fontsize = 17)
    
 
    





