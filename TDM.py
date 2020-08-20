#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 10:50:15 2019
Create Diagonal Matrix
Uses numba and jit conda enviroment for parallelization and speed

@author: kirahart
"""
import numpy as np
from numba import jit, f8 
import matplotlib

def CreateDiagonals(nx,deltax,deltaz,k0,sgn):
    delta    = deltaz/(4*k0*deltax**2)
    idelta   = 1j*delta
    #diagonal matrix
    a = np.ones(nx-1,dtype = 'complex')
    b = np.ones(nx,dtype = 'complex')
    c = np.ones(nx-1,dtype = 'complex')
    b[0:nx] = 1 + sgn * 2 *idelta
    #offdiagonal
    a = np.multiply(a,(idelta*sgn))
    c = a
    return(a,b,c)
    
def tridiag(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)


@jit(f8[:] (f8[:],f8[:],f8[:],f8[:]))
def TDMsolve(a, b, c, d):
    nf = len(b) # number of equations
    aa, bb, cc, dd = map(np.array, (a, b, c, d)) # copy arrays
    for i in range(1, nf):
        if bb[i-1]==0:   
            print("Danger, Will Robinson")
        m = aa[i-1]/bb[i-1]
        bb[i] = bb[i] - m*cc[i-1] 
        dd[i] = dd[i] - m*dd[i-1]
    x = bb
    x[-1] = dd[-1]/bb[-1]

    for i in range(nf-2, -1, -1):
        x[i] = (dd[i]-cc[i]*x[i+1])/bb[i]
    return(x)
    
def TDMmultiply(a,b,c,r):
    f = np.zeros(len(r),dtype = 'complex')
    f[0] = b[0]*r[0] + c[0]*r[1] 
    i = 1
    while i < (len(b)-1):
        f[i] =  a[i-1]*r[i-1]+ b[i]*r[i]+c[i]*r[i+1]
        i = i +1
    f[i] = a[i-1]*r[i-1]+ b[i]*r[i]
    return(f)
    
def CrankNicolson(dx,k0,Eold,zstep):
    nx = len(Eold)
    #make L+ and L- matricites
    ap,bp,cp = CreateDiagonals(nx,dx,zstep,k0,1)
    am,bm,cm = CreateDiagonals(nx,dx,zstep,k0,-1)
    #calculate L+.Eold
    d = TDMmultiply(ap,bp,cp,Eold)
    #use TDM solve
    x = TDMsolve(am,bm,cm,d)
    return(x)
    