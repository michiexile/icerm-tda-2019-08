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

savetxt(args.outputfile, mat)
