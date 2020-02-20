#!/usr/bin/env python

import matplotlib
matplotlib.use('AGG')

import sys, re, scipy, mdp, glob, argparse
import matplotlib.pyplot as plt
from scipy import r_

parser = argparse.ArgumentParser(description="Plot a persistence diagram",
				 epilog='If no dimension is specified, all dimensions will be plotted, distinguishing by color.')
parser.add_argument('slug', metavar='slug', help='prefix for the persistence computation')
parser.add_argument('-d', metavar='d', type=int, help='dimension to plot')
args = parser.parse_args()

dgm = scipy.loadtxt("%s.dgm" % args.slug)

infval = 2*dgm[dgm!=float('inf')].max()
dgm[dgm==float('inf')] = infval
plt.figure(figsize=[4,4])
for (c,i) in enumerate(scipy.unique(dgm[:,0])):
    plt.plot(dgm[dgm[:,0]==i,1],dgm[dgm[:,0]==i,2],marker='.',markersize=10,linestyle='None')

plt.axis([0,infval,0,infval])
#plt.axis('equal')

plt.savefig('%s-dgm.pdf' % args.slug)

