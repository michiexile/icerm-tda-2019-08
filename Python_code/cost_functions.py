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
def cost_2norm(Z, F, B):
    '''
    Compute the cost function $c = \|f-\delta z\|_2$
    '''
    return np.sum(np.square(F - B*Z.reshape(B.shape[1],1)))

def grad_2norm(Z, F, B):
    '''
    Compute the gradient of the cost function $c = \|f-\delta z\|_2$
    '''
    dz = B*Z.reshape(B.shape[1],1)
    return np.array(2 * B.T * (dz - F)).flatten()

def lpnorm(Z, ord=2):
    return float_power(float_power(Z, ord).sum(), 1/ord)

def cost_Lpnorm_mvj(Z, F, B, p=2, alpha=0, q=2):
    return lpnorm(array(F + B @ Z.reshape(-1,1)), ord=p) + \
        alpha*lpnorm(array(F + B @ Z.reshape(-1,1)), ord=q)

def grad_Lpnorm_mvj(Z, F, B, p=2, alpha=0, q=2):
    def fbzlpgrad(fbz, z, p):
        return float_power(lpnorm(array(fbz), p), 1-p) * ( ( float_power(abs(fbz), p-1).ravel() @ (B@diag(sign(z)))))
    bz = B @ Z.reshape(-1,1)
    fbz = F - bz
    return array(fbzlpgrad(fbz, Z, p) + alpha * fbzlpgrad(fbz, Z, q)).ravel()