#!/usr/bin/env python

import matplotlib
matplotlib.use('AGG')

import sys, re, scipy, glob, argparse, os.path
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Produce a 2d PCA plot of the data")
parser.add_argument('val1', metavar='val1', help='First valfile')
parser.add_argument('val2', metavar='val2', help='First valfile')
args = parser.parse_args()

val1 = scipy.loadtxt(args.val1)
ccl1 = scipy.mod(val1[:,1],1.0)
val2 = scipy.loadtxt(args.val2)
ccl2 = scipy.mod(val2[:,1],1.0)
plt.figure(figsize=(4,4))
plt.scatter(ccl1,ccl2,edgecolor='none')
plt.savefig("%s-%s.pdf" % (args.val1[:-4], os.path.basename(args.val2)))
plt.close()
