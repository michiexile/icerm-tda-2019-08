import numpy as np
import dionysus

def coboundary(vr, thr, n_nodes=100):
    '''
    Compute coboundary matrices of the Vietoris-Rips complex with ball radius thr for dionysus 2.0
    
    Parameters
    ----------
    vietoris_rips: dionysus._dionysus.Filtration

    thr: float()

    Returns
    -------
    list of coboundary matrices [d1,d2,...]. 
    
    dict(dict()) for every dimension k, a dictionary with keys the simplices of dimension k and the relative column id in the coboundary matrix
    '''

    D = {}
    data = {}
    indexing = {}
    ix = [0]*n_nodes
    for s in vr:
        if s.dimension() == 0:
            continue
        elif s.data > thr:
            break
        D.setdefault(s.dimension(),[[],[]])
        data.setdefault(s.dimension(),[])
        indexing.setdefault(s.dimension(),{})
        indexing.setdefault(s.dimension()-1,{})
        if not s in indexing[s.dimension()]:
            indexing[s.dimension()][s] = ix[s.dimension()]
            ix[s.dimension()] += 1
        for dat, k in enumerate(s.boundary()): 
            if not k in indexing[s.dimension()-1]:
                indexing[s.dimension()-1][k] = ix[s.dimension()-1]
                ix[s.dimension()-1] += 1
            D[s.dimension()][0].append(indexing[s.dimension()][s]) #rows
            D[s.dimension()][1].append(indexing[s.dimension()-1][k]) #cols
            data[s.dimension()].append(1. if dat % 2 == 0 else -1.)
    import scipy as sp
    return [sp.sparse.csr_matrix((data[d], (D[d][0], D[d][1]))).todense() for d in range(1,max(D.keys())+1)], indexing
