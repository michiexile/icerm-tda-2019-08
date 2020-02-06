import numpy as np
from numpy import *

def cost_interpolate(Z, L, F, B):
    '''
    Compute the cost function $c = (1 - \lambda ) \|f-\delta z\|_2 + \lambda \|f-\delta z\|_1$
    '''
    return (1-L) * np.sum(np.square(F - B*Z.reshape(B.shape[1],1))) + L* np.linalg.norm(B*Z.reshape(Z.shape[0],1),ord=1)

def grad_interpolate(Z, L, F, B):
    '''
    Compute the gradient of the cost function $c = (1 - \lambda ) \|f-\delta z\|_2 + \lambda \|f-\delta z\|_1$
    '''
    dz = B*Z.reshape(B.shape[1],1)
    return np.array((1-L) *2 * B.T * (dz - F) + L* B.T * np.sign(dz - F)).flatten()

#----------
def cost_1norm(Z, F, B):
    '''
    Compute the cost function $c = \|f-\delta z\|_1$
    '''
    return np.linalg.norm(B*Z.reshape(Z.shape[0],1),ord=1)

def grad_1norm(Z, F, B):
    '''
    Compute the gradient of the cost function $c = \|f-\delta z\|_1$
    '''
    dz = B*Z.reshape(B.shape[1],1)
    return np.array(B.T * np.sign(dz - F)).flatten()

#----------
#L^2 norm penalty
def cost_2norm(Z, F, B):
    '''
    Compute the cost function $c = \|f-\delta z\|_2$
    B is the boundary delta parameter;
    F is the evaluated values on a finite grid;
    Z is the coordinate mapping.
    '''
    return np.sum(np.square(F - B*Z.reshape(B.shape[1],1)))

def grad_2norm(Z, F, B):
    '''
    Compute the gradient of the cost function $c = \|f-\delta z\|_2$
    with analytic expression.
    '''
    dz = B*Z.reshape(B.shape[1],1)
    return np.array(2 * B.T * (dz - F)).flatten()

def lpnorm(Z, ord=2):
    return float_power(float_power(Z, ord).sum(), 1/ord)

def cost_Lpnorm_mvj(Z, F, B, p=2, alpha=0, q=2):
    '''
    Compute the cost function $c = \|f-\delta z\|_p + alpha*\|f-\delta z\|_q$
    '''
    return lpnorm(array(F + B * Z.reshape(-1,1)), ord=p) + \
        alpha*lpnorm(array(F + B * Z.reshape(-1,1)), ord=q)

def grad_Lpnorm_mvj(Z, F, B, p=2, alpha=0, q=2):
    def fbzlpgrad(fbz, z, p):
        #Calculate the lp-norm gradient for the quantity F-bz
        fbz=fbz.reshape(-1,1)
        tmp_scalar=float_power(abs(fbz), p-1).ravel()
        tmp_scalar=tmp_scalar.reshape(-1,1)
        #print tmp_scalar.shape
        #print B.shape
        #print diag(sign(z)).shape
        return float_power(lpnorm(fbz, p), 1-p) * ( ( tmp_scalar * (B*diag(sign(z)))))
    bz = B * Z.reshape(-1,1)
    fbz = F - bz
    fbzlpgrad_tmp1 = fbzlpgrad(fbz, Z, p)
    fbzlpgrad_tmp2 = fbzlpgrad(fbz, Z, q)
    return array(fbzlpgrad_tmp1 + alpha * fbzlpgrad_tmp2).ravel()
