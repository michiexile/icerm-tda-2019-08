import dionysus
import scipy as sp
import numpy as np

def coboundary_1(vr, thr):
    D = [[],[]]
    data = []
    indexing = {}
    ix = [0]*2
    for s in vr:
        if s.dimension() != 1:
            continue
        elif s.data > thr:
            break
        indexing.setdefault(s.dimension(),{})
        indexing.setdefault(s.dimension()-1,{})
        if not s in indexing[s.dimension()]:
            indexing[s.dimension()][s] = ix[s.dimension()]
            ix[s.dimension()] += 1
        for dat, k in enumerate(s.boundary()): 
            if not k in indexing[s.dimension()-1]:
                indexing[k.dimension()][k] = k[0]
                ix[k.dimension()] += 1
            D[0].append(indexing[s.dimension()][s]) #rows
            D[1].append(indexing[k.dimension()][k]) #cols
            data.append(1. if dat % 2 == 0 else -1.)
    return sp.sparse.csr_matrix((data, (D[0], D[1]))), indexing

def optimizer_inputs(vr, bars, cocycle, init_z, prime):
    bdry,indexing = coboundary_1(vr,max(bar.death for bar in bars))
    n, m = bdry.shape # edges X nodes
    #-----------------
    l2_cocycle = [0]*len(init_z) #reorganize the coordinates so they fit with the coboundary indices
    for i, coeff in enumerate(init_z):
        l2_cocycle[i] = coeff
    l2_cocycle = np.array(l2_cocycle)
    #-----------------
    f = np.zeros((n,1)) # cocycle we need to smooth out, reorganize to fit coboundary
    for c2 in cocycle:
        if c2.element<(prime//2):
            f[indexing[1][vr[c2.index]]] += c2.element
        else:
            f[indexing[1][vr[c2.index]]] += c2.element-prime  
    return l2_cocycle,f,bdry