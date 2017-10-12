#!/usr/bin/python
import sys,os
import math
import numpy as np
import random
from copy import deepcopy

def dot( a, m ):
    matrix = np.array(m)
    (I,J) = matrix.shape
    prod = np.zeros( (len(a) ) ) 
    try:
        for i in range(3):
            for j in range(3):
                prod[i] = prod[i] + a[j]*matrix[j][i]
    except:
        print "Dot error"
        sys.exit(1)
    return prod

def make_matrix( x,y,z,theta):
    costt = math.cos(theta)
    sintt = math.sin(theta)
    matrix = [ [ costt+(1-costt)*x*x  , (1-costt)*x*y-sintt*z, (1-costt)*x*z+sintt*y ],
               [ (1-costt)*y*x+sintt*z, costt+(1-costt)*y*y  , (1-costt)*y*z-sintt*x ],
               [ (1-costt)*z*x-sintt*y, (1-costt)*z*y+sintt*x,  costt+(1-costt)*z*z  ] ]
    matrix = np.array(matrix)
    return matrix

def build_matrix( c1, c2 , theta ):
    x1,y1,z1 = c1[:3]
    x2,y2,z2 = c2[:3]
    dx = x1-x2
    dy = y1-y2
    dz = z1-z2
    mod = math.sqrt( dx*dx + dy*dy + dz*dz )
    matrix = make_matrix( dx/mod,dy/mod,dz/mod,theta )
    return matrix

def rotate_atom( coord, c1 , matrix ):
    x,y,z = coord[:3]
    x0,y0,z0 = c1[:3]
    dx,dy,dz = x-x0,y-y0,z-z0
    nx,ny,nz = dot( (dx,dy,dz), matrix ) 
    return nx+x0,ny+y0,nz+z0
