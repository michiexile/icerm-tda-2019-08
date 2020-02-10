from numpy import r_, sqrt
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

def annulus_example_AP(R= 1.5, d= .5, n= 100):
    '''
    Creates a point cloud in R^2 of randomly sampled points from an 
    annulus with inner radius R and outer radius R+d. Area-preserveing sampling
    
    Parameters
    ----------
    R: float() - inner radius of the annulus
    d: float() - thickness of the annulus
    n: int() - number of points
    
    
    Output
    ------
    np.array of dimensions (n,2)
    '''
    def psi(z1,z2):
        den = d*(2*R+d)
        num = -R+sqrt(den*z1+R**2)
        return num/den,z2
    
    def rand_annulus_pt():
        r = uniform(0,1)
        th = uniform(0,1)
        r, th = psi(r,th)
        return (R+r*d) * cos(th*2*pi), (R+r*d) * sin(th*2*pi)
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
    
def lorenz_example():
    '''
    TO DO: need to add number of points or something of the sort
    '''
    def f(state, t):
        (x,y,z) = state
    return (sigma*(y-x), x*(rho-z)-y, x*y-beta*z)

    ts = arange(0,40,0.025)
    states = scipy.integrate.odeint(f, (1.,1.,1.), ts)
    return states[-800:,:]

def torus_example(R1= 1, R2= .3, n= 100):
    '''
    R1 #donut hole
    R2 #thickness
    n # number of points
    '''
    def rand_torus_pt():
        th = uniform(0,2*pi)
        ph = uniform(0,2*pi)
        return (R1+R2*cos(th))*cos(ph), (R1+R2*cos(th))*sin(ph),R2*sin(th)
    return np.array([(rand_torus_pt()) for x in range(n)])

def pinched_torus_example(R= 1.5, n= 100):
    '''
    Creates a point cloud in R^3 of randomly sampled points from an
    piched torus with radius R.
    Parameters
    ----------
    R: float() - radius of the annulus
    n: int() - number of points
    Output
    ------
    np.array of dimensions (n,3)
    '''
    def rand_S1_pt():
        #https://en.wikipedia.org/wiki/Pinched_torus
        x = uniform(0,2*pi)
        y = uniform(0,2*pi)
        g_xy = 2 + sin(x/2)*cos(y)
        return R*g_xy*cos(x), R*g_xy*sin(x), R*sin(x/2)*sin(y)
    return np.array([(rand_S1_pt()) for x in range(n)])


def double_pinched_torus_example(R= 1.5, n= 100, a=.5, b=.5, c=0, d=0):
    '''
    Creates a point cloud in R^3 of randomly sampled points from an
    piched torus with radius R.
    Parameters
    ----------
    R: float() - radius of the annulus
    n: int() - number of points
    Output
    ------
    np.array of dimensions (n,3)
    '''
    def rand_S1_pt(a=a,b=b,c=c,d=d):
        #https://demonstrations.wolfram.com/DupinCyclides/
        #https://en.wikipedia.org/wiki/Dupin_cyclide
        u = uniform(-pi,pi)#phi
        v = uniform(-pi,pi)#psi
        #if c*c!=(a*a-b*b):
        #   print('ERROR in parameters')
     
        X_coord = d*(c-a*cos(u)*cos(v))+b*b*cos(u)
        X_coord = X_coord/(c-a*cos(u)*cos(v))

        Y_coord = b*sin(u)*(a-d*cos(v))
        Y_coord = Y_coord/(c-a*cos(u)*cos(v))

        Z_coord = b*sin(v)*(c*cos(u)-d)
        Z_coord = Z_coord/(c-a*cos(u)*cos(v))
    
        return X_coord, Y_coord, Z_coord
    return np.array([(rand_S1_pt()) for x in range(n)])


def klein_bottle_example_4D(R=0.9, P=1.2, epsilon=0.01, n= 100):
    '''
    Creates a point cloud in R^4 of randomly sampled points from an
    Kelin bottle with radius R/P and epsilon bumped in R^4.
    Parameters
    ----------
    R: float() - radius (aspect ratio 1) of the immersed torus
    P: float() - radius (aspect ratio 2) of the immersed torus
    epsilon: float() - small bumping parameter
    n: int() - number of points
    Output
    ------
    np.array of dimensions (n,4)
    '''
    def rand_TR_pt():
        #https://en.wikipedia.org/wiki/Klein_bottle
        theta = uniform(0,2*pi)
        v = uniform(0,2*pi)
        return R*( cos(theta/2)*cos(v)-sin(theta/2)*sin(2*v) ) , R*( sin(theta/2)*cos(v)+cos(theta/2)*sin(2*v) ), P*cos(theta)*( 1+epsilon*sin(v) ), P*sin(theta)*( 1+epsilon*sin(v) )
    return np.array([(rand_TR_pt()) for x in range(n)])

def klein_bottle_example_Fig8(r=2.1, n= 100):
    if r <= 2:
        raise Exception('r should >2.')
    '''
    Creates a point cloud in R^3 of randomly sampled points from an
    Klein bottle with radius r and immersed like a figure 8 in R^3.
    Parameters
    ----------
    r: float() - radius of the immersed 8.
    n: int() - number of points
    Output
    ------
    np.array of dimensions (n,3)
    '''
    def rand_F8_pt():
        #https://en.wikipedia.org/wiki/Klein_bottle
        theta = uniform(0,2*pi)
        v = uniform(0,2*pi)
        return ( r + cos(theta/2)*sin(v)-sin(theta/2)*sin(2*v) )*cos(theta), ( r + cos(theta/2)*sin(v)-sin(theta/2)*sin(2*v) )*sin(theta), sin(theta/2)*sin(v)+cos(theta/2)*sin(2*v)
    return np.array([(rand_F8_pt()) for x in range(n)])
