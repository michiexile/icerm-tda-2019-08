#!/usr/bin/env python

import argparse, sys, re, scipy

parser = argparse.ArgumentParser(description="Compute a delay embedding")
parser.add_argument('input', metavar='input', help='data file')
parser.add_argument('dimension', type=int, help='dimension increase multiplier')
parser.add_argument('delay', type=int, help='delay size, given in sample count')
parser.add_argument('output', default=None, help='output file')
args = parser.parse_args()


dim = args.dimension
delay = args.delay

mat = scipy.loadtxt(args.input)
if mat.ndim < 2:
    mat.resize((len(mat),1))
    
ret = scipy.concatenate([mat[i*delay:-(dim-i)*delay,:] for i in range(dim)],1)

if args.output:
    scipy.savetxt(args.output, ret)
else:
    print ret
