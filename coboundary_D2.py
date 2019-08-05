def coboundary(vietoris_rips, thr, n_nodes=100):
    '''
    Parameters
    ----------
    vietoris_rips: dionysus._dionysus.Filtration

    thr: float()

    Returns
    -------

    list of coboundary matrices
    '''

    D = {}
    data = {}
    indexing = {}
    ix = [0]*n_nodes
    thr = 1.5
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
            if dat%2:
                data[s.dimension()].append(-1)
            else:
                data[s.dimension()].append(1)
    import scipy as sp
    return [sp.sparse.csr_matrix((data[d], (D[d][0], D[d][1]))).todense() for d in range(1,max(D.keys())+1)]
