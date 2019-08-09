from numpy import r_
from math import sin, cos, pi
from random import uniform
import numpy as np

def annulus_example(R= 1.5, d= .5, n= 100):
    '''
    Creates a point cloud in R^2 of randomly sampled points from an 
    annulus with inner radius R and outer radius R+d.
    
    Parameters
    ----------
    R: float() - inner radius of the annulus
    d: float() - thickness of the annulus
    n: int() - number of points
    
    
    Output
    ------
    np.array of dimensions (n,2)
    '''
    def rand_annulus_pt():
        r = uniform(0,1)
        th = uniform(0,2*pi)
        return (R+r*d) * cos(th), (R+r*d) * sin(th)
    return np.array([(rand_annulus_pt()) for x in range(n)])

def annulus_variable_d_example(R= 1.5, d= .5, n= 100):
    '''
    Creates a point cloud in R^2 of randomly sampled points from an 
    annulus with inner radius R and outer radius between (R+d, R).
    
    Parameters
    ----------
    R: float() - inner radius of the annulus
    d: float() - thickness of the annulus
    n: int() - number of points
    
    
    Output
    ------
    np.array of dimensions (n,2)
    '''
    annulus = np.zeros((n, 2))
    idx = 0
    while (idx < n):
        x = uniform(-R-d, R+d)
        y = uniform(-R-d, R+d)
        if ((x*x + y*y) / ((R+d)*(R+d)) < 1) & (x*x/((R-d)*(R-d)) + y*y/(R*R) > 1):
            annulus[idx, :] = [x, y]
            idx = idx + 1
    return annulus

def annulus2_example(R= 1.5, d= .5, n= 100):
    '''
    Creates a point cloud in R^2 of randomly sampled points from two overlapping 
    annuli with inner radius R and outer radius R+d.
    
    Parameters
    ----------
    R: float() - inner radius of the annulus
    d: float() - thickness of the annulus
    n: int() - number of points
    
    
    Output
    ------
    np.array of dimensions (n,2)
    '''
    def rand_annulus_pt():
        r = uniform(0,1)
        th = uniform(0,2*pi)
        return (R+r*d) * cos(th), (R+r*d) * sin(th)
    return r_[np.array([(rand_annulus_pt()) for x in range(n//2)]) + [[-R,0]],
             np.array([(rand_annulus_pt()) for x in range(n//2)]) + [[R,0]]]

def annulus_bar_example(R= 1.5, d= .5, n= 100):
    '''
    Creates a point cloud in R^2 of randomly sampled points from an 
    annulus with inner radius R and outer radius R+d with an bar of thickness d cutting it in half.
    
    Parameters
    ----------
    R: float() - inner radius of the annulus
    d: float() - thickness of the annulus
    n: int() - number of points
    
    
    Output
    ------
    np.array of dimensions (n,2)
    '''

    def rand_annulus_pt():
        r = uniform(0,1)
        th = uniform(0,2*pi)
        return (R+r*d) * cos(th), (R+r*d) * sin(th)

    def vert_bar_pt():
        xx = uniform(-d, d)
        yy = uniform(-R, R)
        return xx, yy

    return r_[np.array([(rand_annulus_pt()) for x in range(2*n//3)]), 
                 np.array([(vert_bar_pt()) for x in range(n//3)])]
    