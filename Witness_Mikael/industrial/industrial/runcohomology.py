#!/usr/bin/env python

import argparse, sys, subprocess, glob, re, os

cocycle="./cocycle.py"
delay="./delay.py"
rips="./rips-pairwise-cohomology"
ripsdist="./rips-explicit-cohomology"
plotpersistence='./plotpersistence.py'
plotcoord='./plotcoord.py'
plotts='./plotts.py'
plotcorr='./plotcorr.py'

parser = argparse.ArgumentParser(description="Run the cohomology gait normalization pipeline")
parser.add_argument('slug', metavar='slug', help='slug to generate all filenames')
parser.add_argument('-m', dest='m', type=float, default=50, help='Vietoris-Rips parameter value for cohomology computation')
parser.add_argument('-delay', dest='delay', nargs=2, default=None, help='use delay embedding. Arguments: embedding dimension then window size.')
parser.add_argument('-dist', dest='dist', default=None, help='use this executable to generate a distance matrix for explicit cohomology computation')
parser.add_argument('-distm', dest='distm', default='euclidean', help='When using -dist, use this metric from scipy.spatial.distance')
args = parser.parse_args()

print args

xyf = "%s.xy" % args.slug
ptsf = "%s.xy" % args.slug

if args.delay:
    print "\n\tDelay embedding into %s dimensions with %s delay" % (args.delay[0],args.delay[1])
    demf = "%s.del" % args.slug
    (ut,st,cut,cst,rt) = os.times()
    subprocess.call([delay, ptsf, args.delay[0], args.delay[1], demf])
    (aut,ast,acut,acst,art) = os.times()
    print "\t\t\tDelay: user %0.5f system %0.5f wall %0.5f" % (acut-cut,acst-cst,art-rt)
    ptsf = demf

if args.dist:
    print "Computing external distance"
    xydf = '%s.xyd' % args.slug
    subprocess.call([args.dist,
		     '-m', args.distm,
		     ptsf,
		     xydf])
    p = subprocess.Popen(['wc', '-l', xydf], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
	raise IOError(err)
    xydn = int(result.strip().split()[0])
    subprocess.call([ripsdist, 
		 "-m", "%f" % args.m,
		 "-b", "%s.bdry" % args.slug,
		 "-v", "%s.vtx" % args.slug,
		 "-c", args.slug,
		 "-d", "%s.dgm" % args.slug,
		 xydf, '%d' % xydn])
else:
    subprocess.call([rips,
		 "-m", "%f" % args.m,
		 "-b", "%s.bdry" % args.slug,
		 "-v", "%s.vtx" % args.slug,
		 "-c", args.slug,
		 "-d", "%s.dgm" % args.slug,
		 ptsf])

f = open("%s.dgm" % args.slug)
for line in f:
    if "inf" in line:
	print line,
f.close()

print "\n\tPlotting persistence diagram"
subprocess.call([plotpersistence,args.slug])

print "\n\tComputing cocycles"
cclfiles = glob.glob("%s*ccl" % args.slug)
for cclf in cclfiles:
    print "%s" % cclf
    (ut,st,cut,cst,rt) = os.times()
    subprocess.call([cocycle,
		     "%s.bdry" % args.slug,
		     cclf,
		     "%s.vtx" % args.slug],stderr=sys.stderr)
    (aut,ast,acut,acst,art) = os.times()
    print "\t\t\tCocycle: user %0.5f system %0.5f wall %0.5f" % (acut-cut,acst-cst,art-rt)


print "\n\tPlotting resulting coordinates"
valfiles = glob.glob("%s*.val" % args.slug)
for valf in valfiles:
    print valf
    subprocess.call([plotcoord, ptsf, valf])
    subprocess.call([plotts, valf])

for i in range(len(valfiles)):
    for j in range(i+1,len(valfiles)):
	val1 = valfiles[i]
	val2 = valfiles[j]
	subprocess.call([plotcorr,val1,val2])

