import numpy as np

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