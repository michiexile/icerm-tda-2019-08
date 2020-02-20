#!/usr/bin/env python

import matplotlib
matplotlib.use('AGG')

import sys, re, scipy, glob, argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Produce a 2d PCA plot of the data")
parser.add_argument('valf', metavar='valf', help='valfile')
args = parser.parse_args()

val = scipy.loadtxt(args.valf)
ccl = scipy.mod(val[:,1],1.0)
plt.figure(figsize=(10,1))
plt.scatter(val[:,0],ccl,edgecolor='none')
plt.savefig("%s-ts.pdf" % args.valf[:-4])
plt.close()
