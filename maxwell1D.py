#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 09:36:41 2019

1D Maxwell Solver
Yee Algorithm

includes functions to be imported to solve other problems
no variables hard coded
push changes to git repository

courant condition = CFL = c dt/dx

@author: kirahart
"""

import numpy as np
import cmath

def take1Dstep(E0,H0,CFL,n):
    i=0
    E=E0
    H=H0
    while i < n-1:
        E[i] = E[i] + CFL * (H[i+1] - H[i]) 
        H[i] = H[i] + CFL * (E[i] - E[i-1]) 
        i=i+1
        
    E[n-1] = E[0]
    H[0] = H[n-1]
    return(E,H)
        
def gaussian(la,k0,w0,x,x0):
    return( cmath.sin(k0*x)*cmath.exp(-((x - x0)/w0)**2.0))

#sets up initial gaussian fields, error in inital condition
def gaussianField_err(la, k0, w0, x0, N, CFL,RL):
    Ei = np.zeros(N,dtype=np.complex)
    Hi = np.zeros(N,dtype=np.complex)
    x0 = 0.5 * N
    i = 0
    while i < N:
        x = i 
        Ei[i] = np.multiply(1,gaussian(la,k0,w0,x,x0))
        Hi[i] = np.multiply(RL,gaussian(la,k0,w0,x,x0))
        i=i+1
    return(Ei,Hi)
    
#sets up initial gaussian fields
def gaussianField(la, k0, w0, x0, N, CFL,RL):
    Ei = np.zeros(N,dtype=np.complex)
    Hi = np.zeros(N,dtype=np.complex)
    x0 = 0.5 * N
    i = 0
    while i < N:
        xe = i 
        xh = xe
        xh = xh - 0.5 * CFL  # shift E field by 1/2c dt
        Ei[i] = np.multiply(1,gaussian(la,k0,w0,xe,x0))
        Hi[i] = np.multiply(RL,gaussian(la,k0,w0,xh,x0))
        i=i+1
    return(Ei,Hi)

def propagate1D(numsteps , Ei, Hi ,CFL):
    E = Ei
    H = Hi
    N=len(Ei)
    i =0
    while i < numsteps:
        E,H = take1Dstep(E,H,CFL,N)
        i=i+1
    return(E, H)
    
def save_steps1D(numsteps , Ei, Hi ,CFL):
    N=int(len(Ei))
    E = np.zeros((numsteps,N))
    H = np.zeros((numsteps,N))
    E[0] = Ei
    H[0] = Hi
    i =1
    while i < numsteps:
        E[i],H[i] = take1Dstep(E[i-1],H[i-1],CFL,N)
        i=i+1
    return(E, H)
    
def simple_wave(x,k0):
    return(cmath.cos(k0*x))
    
def simple_field(k0, N):
    Ei = np.zeros(N,dtype=np.complex)
    Hi = np.zeros(N,dtype=np.complex)
    for i in range(N):
        Ei[i] =  simple_wave(i,k0)
        Hi[i] = 0
    return(Ei,Hi)
    
def initial_whitenoise(N):
    Ei=np.random.uniform(-1,1,N)
    Hi=np.zeros(N)
    return(Ei,Hi)
    
    
def zero_counting(Ei,Hi,CFL,N,k0,countnum):
    E=Ei;H=Hi
    istart=0;
    start=0;count=0; i = 0
    Efold = E[0]
    while count < countnum:
        E,H = take1Dstep(E,H,CFL,N)
        if E[0]* Efold <=0:
            if start == 0:
                start = 1 ; istart= i ; count = 0
            else:
                count=count+1
        i = i+1
        Efold= E[0]
    i = i -1
    T = 2* CFL *(i-istart)/count
    omega= 2*cmath.pi/T
    YeeValue = 2.0/CFL*cmath.asin( CFL* cmath.sin(k0/2.0))
    return(omega,YeeValue)
    
def find_envelope_analytic():
    print("hi")
    
    
    
