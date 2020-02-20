#!/usr/bin/env python

import matplotlib
matplotlib.use('AGG')

import sys, re, scipy, mdp, glob, argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Produce a 2d PCA plot of the data")
parser.add_argument('xyf', metavar='xyf', help='xyfile')
parser.add_argument('valf', metavar='valf', help='valfile')
args = parser.parse_args()

xy = scipy.loadtxt(args.xyf)

pca = mdp.nodes.PCANode(svd=True)
pca.train(xy)
pca.stop_training()
xyp = pca.execute(xy)

plt.hsv()

val = scipy.loadtxt(args.valf)
ccl = scipy.mod(val[:,1],1.0)
plt.figure(figsize=(1.5,1.5))
plt.scatter(xyp[:,0],xyp[:,1],c=ccl,edgecolor='none',s=2)
plt.axis('off')
plt.savefig("%s.pdf" % args.valf[:-4])
plt.close()

