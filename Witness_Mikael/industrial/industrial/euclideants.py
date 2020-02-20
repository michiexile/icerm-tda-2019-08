#!/usr/bin/env python

import argparse, sys
from scipy import savetxt, loadtxt
from scipy.spatial.distance import pdist, squareform

parser = argparse.ArgumentParser(description="Euclidean metric preprocessor")
parser.add_argument('inputfile', help='points')
parser.add_argument('outputfile', help='metric matrix')
parser.add_argument('-m', dest='m', default='euclidean', help='metric; from scipy.spatial.distances')
args = parser.parse_args()

pts = loadtxt(args.inputfile)
mat = squareform(pdist(pts, args.m))

# Modify metric for timeseries representation
for i in range(mat.shape[0]-1):
    mat[i,i+1] = 0
    mat[i+1,i] = 0

savetxt(args.outputfile, mat)
