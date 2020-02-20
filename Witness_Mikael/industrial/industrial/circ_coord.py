#!/usr/bin/env python
import os
from sys import argv, exit
import scipy
import numpy
import cocycle

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def comp(x,y):
    a = x[1]-x[0]
    b = y[1]-y[0]
    if a<b:
        return 1
    if a==b:
        return 0
    if a>b:
        return -1
    

def parse_dgm(dgm_file,max_rips):
    persistence=[]

    l1h = open(dgm_file)
    
    for l in l1h:
        (d,bi,de) = l.split()
        di = float(de)
       
        if float(d)==1:
            if numpy.isinf(di):
                persistence.append((float(bi),max_rips))
            else:
                persistence.append((float(bi),di))

        
            
    return persistence

#==============================
# a wrapper for the rips
# cohomology Dionysus code
#==============================


    

def circ_coord(max_rips,data,dgm_file,coord_file):


    exe = './rips-pairwise-cohomology'	
    cmd = "%s -m %s -d tmp_dgm '%s'" % (exe, max_rips, data)
    print cmd
    os.system(cmd)

    
    intervals  = parse_dgm('tmp_dgm',float(max_rips))
    
    
    intervals.sort(comp)
    coordinates = []
    num = 3
    if len(intervals)<num:
        num=len(intervals)

    print "Saving barcode: %s" % dgm_file
    scipy.savetxt(dgm_file,intervals)
    print "Barcode saved."
    
    i = 0
    print "i"
    coordinates = []
    for x in range(num):
        print "loop: %d" % x
        print intervals
        half_way  = (intervals[x][1]+intervals[x][0])/2
        print "halfway: %f" % half_way
        cmd = "%s -m %f -b tmp_bdy.bdry -c cocycles -v vertex.vtx '%s'" % (exe, half_way, data)
        print cmd
        os.system(cmd)

        boundary_filename = 'tmp_bdy.bdry'
        cocycle_filename = 'cocycles'+'-'+str(i)+'.ccl'  #This should be the correct cocycle
        i=i+1
        vertexmap_filename = 'vertex.vtx'
       
        boundary_list = cocycle.read_list_file(boundary_filename)
        cocycle_list =  cocycle.read_list_file(cocycle_filename)
        vertexmap_list =  cocycle.read_list_file(vertexmap_filename)

        solution, v =  cocycle.smooth(boundary_list, cocycle_list)
        values =  cocycle.vertex_values(solution, vertexmap_list)
        if coordinates == []:
            coordinates = values
        else:
            coordinates = scipy.vstack((coordinates,values))

    print "Writing coordinates"
    scipy.savetxt(coord_file,scipy.matrix(coordinates))



