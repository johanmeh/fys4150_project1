#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 16:23:26 2018

@author: maritkollstuen
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 08:30:41 2018

@author: maritkollstuen

"""

"""
Solution to task 1 c) of assignment 1 in FYS4150

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

def f(x): # Source term of the differential equation 
    return 100*np.exp(-10*x)


def u(x): # Closed - form solution of the differential equation 
    return 1-(1-np.exp(-10))*x - np.exp(-10*x)


def gauss_elim (a, b, n):
    x = np.linspace(0, 1, (n+2))
    h = 1/(n+1)
    f_vec = h**2*f(x)
    
    f_vec[0] = 0
    f_vec[-1]= 0
     
    
    b_vec = np.ones(n)*b
    
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
        
    return (v_vec)

def plot_stuff(a, b, n):
    x = np.linspace(0, 1, (n+2))
    u_exact = u(x)
    
    v  = gauss_elim(a, b, n)
    plt.figure()
    plt.plot(x, u_exact, color='blue', label = 'u - analytic')
    plt.plot(x, v, color='green', label = 'v - numeric, n = {}'.format(n))
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.savefig('size_{}.png'.format(n))
    

if __name__ == "__main__":
    
    for n in [10, 100, 1000]:
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 16:23:26 2018

@author: maritkollstuen
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 08:30:41 2018

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
import timeit

def f(x): # Source term of the differential equation 
    return 100*np.exp(-10*x)


def u(x): # Closed - form solution of the differential equation 
    return 1-(1-np.exp(-10))*x - np.exp(-10*x)


def gauss_elim (a, b, n):
    x = np.linspace(0, 1, (n+2))
    h = 1/(n+1)
    f_vec = h**2*f(x)
    
    f_vec[0] = 0
    f_vec[-1]= 0
     
    
    b_vec = np.ones(n)*b
    
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
        
    return (v_vec)

def relative_error(v ,n):
    x = np.linspace(0, 1, (n+2))
    u_ex = u(x)
    
    epsilon = np.zeros(len(u_ex))
    for i in range(len(x)):
        epsilon[i] = np.log10(np.abs((v[i] - u_ex[i])/v[i]))
    return epsilon



def make_matrix(n, a, b, c): # Make tridiagonal matrix for LU-decomposition
    A = np.zeros((n+2, n+2))
    for i in range(n+2):
        for j in range(n+2):
            if i == j:
                A[i, j] = a
        for j in range(n+1):
            if i == (j+1):
                A[j, i] = c
                A[i, j] = b
    return A


def lu_dec(A, n): # LU- decomposition
    x = np.linspace(0, 1, (n+4))
    f_vec = f(x)
    
    v = np.linalg.solve(A, f_vec[1:-1])
    return v


def plot_stuff(v, n, epsilon, v_lu):
    x = np.linspace(0, 1, (n+2))
    u_exact = u(x)
    
    
    # Plot exact versus numeric result and result from LU-decomposition
 #   v  = gauss_elim(a, b, n)
    plt.figure()
    plt.plot(x, u_exact, color='blue', label = 'u - analytic')
    plt.plot(x, v, color='green', label = 'v - numeric, n = {}'.format(n))
    plt.plot(x, v_lu, color='red', label = 'v - LU-decomposition, n = {}'.format(n))
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('u(x)')
    
    plt.figure()
    plt.plot(x, epsilon, color= 'red', label = 'Relative error')
    plt.xlabel('x')
    plt.ylabel('Epsilon')
 #   plt.savefig('size_{}.png'.format(n))
    


if __name__ == "__main__":
    
    
    
    for n in [10, 100, 1000]:
        x = gauss_elim(-1, 2, n)
        epsilon = relative_error(x, n)
        A = make_matrix(n, -1, 2, -1)
        v_lu = lu_dec(A, n)
        plot_stuff(x, n, epsilon, v_lu)
        
 #        timeit.timeit(x, number = 100)        plot_stuff(-1, 2, n)

