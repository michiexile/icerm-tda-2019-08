#!/usr/bin/env python

from    cvxopt          import spmatrix, matrix
from    cvxopt.blas     import copy
from    lsqr            import lsqr
from    sys             import argv, exit
from    scipy           import polyfit, savetxt, array, c_
import  os.path
import argparse, math

def smooth(boundary_list, cocycle_list):
    dimension = max((max(d[1], d[2]) for d in boundary_list))
    dimension += 1

    # NB: D is a coboundary matrix; 1 and 2 below are transposed
    D = spmatrix([d[0] for d in boundary_list],
                 [d[2] for d in boundary_list],
                 [d[1] for d in boundary_list], (dimension, dimension))

           
    z = spmatrix([zz[0] for zz in cocycle_list],
                 [zz[1] for zz in cocycle_list],
                 [0     for zz in cocycle_list], (dimension, 1))

    v1 = D * z
    # print "D^2 is zero:", not bool(D*D)
    # print "D*z is zero:", not bool(v1)
    z = matrix(z)

    def Dfun(x,y,trans = 'N'):
        if trans == 'N':
            copy(D * x, y)
        elif trans == 'T':
            copy(D.T * x, y)
        else:
            assert False, "Unexpected trans parameter"

    tol = 1e-10
    show = False
    maxit = None
    solution = lsqr(Dfun, matrix(z), show = show, atol = tol, btol = tol, itnlim = maxit)
    
    v = z - D*solution[0]

    # print sum(v**2)
    # assert sum((D*v)**2) < tol and sum((D.T*v)**2) < tol, "Expected a harmonic cocycle"
    if not (sum((D*v)**2) < tol and sum((D.T*v)**2) < tol):
        print "Expected a harmonic cocycle:", sum((D*v)**2), sum((D.T*v)**2) 

    return solution[0], v


def vertex_values(solution, vertices):
    values = [None]*len(vertices)
    for i,v in vertices:
        values[v] = solution[i]
    return values
    

def read_list_file(filename):
    list = []
    with open(filename) as fp:
        for line in fp.xreadlines():
            if line.startswith('#'): continue
            list.append(map(int, line.split()))
    return list

def unwrap_simple(vals):
        lastP = vals[0]
        np = [lastP]
        offset = 0.0
	for p in vals[1:]:
                offset = math.floor(lastP-p+0.5)
                lastP = p+offset
                np.append(lastP)        
        
        return np


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Compute a coordinate from a cocycle")
    parser.add_argument('boundary', metavar='boundary', help='boundary file')
    parser.add_argument('cocycle', metavar='cocycle', help='cocycle file')
    parser.add_argument('vertexmap', metavar='vertexmap', help='vertexmap file')
    args = parser.parse_args()

    boundary_filename = args.boundary
    cocycle_filename = args.cocycle
    vertexmap_filename = args.vertexmap

    boundary_list = read_list_file(boundary_filename)
    cocycle_list = read_list_file(cocycle_filename)
    vertexmap_list = read_list_file(vertexmap_filename)

    solution, v = smooth(boundary_list, cocycle_list)
    values = vertex_values(solution, vertexmap_list)

    ccluw = unwrap_simple(values)
    (ar,vr) = polyfit(range(len(ccluw)),ccluw,1)
    if ar < 0: # decreasing coordinate with time
	valarr = c_[range(len(values)),[-v for v in values]]
    else:
	valarr = c_[range(len(values)),values]

    outfn = os.path.splitext(cocycle_filename)[0] + '.val'
    savetxt(outfn, valarr)
