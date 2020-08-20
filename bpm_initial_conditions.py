#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 10:14:55 2019
Inital Confitions for FFT Methods
@author: kirahart
"""

import numpy as np
from cmath import sin,cos, sqrt, exp, pi

def GaussianBeam1DRotated( x, z, w, k, angle):
    beam = np.zeros(len(x),dtype = 'complex')
    for i in range(len(x)):
        zt = +cos(angle)*z + sin(angle)*x[i];
        xt = -sin(angle)*z + cos(angle)*x[i];
        aux = 1.0 + 2.0* 1j*zt/(k*w**2);
        beam[i] = exp( -xt * xt/(w**2*aux) + 1j*k*zt)/sqrt(aux)
    return(beam)

def GaussianBeam_BPM_FFT(lmb,beamwaist,angle,LX,NX):
    k0 = 2*pi/lmb;
    dk = 2*pi/LX;
    dx = LX/NX;
    zR = pi* beamwaist**2/lmb;
    cx = dx*(np.linspace(0,NX-1,NX)-NX/2) #real space coordinates
    #define transverse wavenumbers
    kx = np.zeros(NX)
    k = 0
    while k < NX/2 + 1:
        kx[k] = dk*k
        k = k+1
    while k < NX :
        kx[k] = dk*(k-NX)
        k = k+1
    return(k0,kx,zR,cx)
    

def GaussianBeam_BPM_FFT2D(lamb,w,LX,LZ,NX,dz,z):
    k0 = 2*pi/lamb;
    dk = 2*pi/LX;
    dx = LX/NX;
    zR = pi* w**2/lamb;
    stps = LZ/dz;
    cx = dx*(np.linspace(0,NX-1,NX)-NX/2) #real space coordinates
    #define transverse wavenumbers
    kx = np.zeros(NX)
    k = 0
    while k < NX/2 + 1:
        kx[k] = dk*k
        k = k+1
    while k < NX :
        kx[k] = dk*(k-NX)
        k = k+1
    #make initial condition
    aux = 1.0 + 2.0* 1j*z/(k*w**2)
    am0 = np.zeros((NX,NX),dtype = "complex")
    for x in range(NX):
        for y in range(NX):
            val =  exp(-((cx[x]**2 + cx[y]**2)/w**2*aux)+ 1j*k*z)/sqrt(aux)
            am0[x,y]=val
    return(k0,kx,zR,stps,cx,am0)
    
def spot_of_arago(lamb,a,LX,LZ,NX,dz):
    k0 = 2*pi/lamb;
    dk = 2*pi/LX;
    dx = LX/NX;
    zR = pi* a**2/lamb;
    stps = LZ/dz;
    cx = dx*(np.linspace(0,NX-1,NX)-NX/2) #real space coordinates
    #define transverse wavenumbers
    kx = np.zeros(NX)
    #max prop angle
    k = 0
    while k < NX/2 + 1:
        kx[k] = dk*k
        k = k+1
    while k < NX :
        kx[k] = dk*(k-NX)
        k = k+1
    #make initial condition
    am0 = np.ones((NX,NX),dtype = "complex")
    for x in range(NX):
        for y in range(NX):
            val =  (cx[x]**2 + cx[y]**2)
            if val < a**2:
                am0[x,y]=val
    return(k0,kx,zR,stps,cx,am0)
    
        
def spot_of_arago_err(lamb,a,LX,LZ,NX,dz):
    k0 = 2*pi/lamb;
    dk = 2*pi/LX;
    dx = LX/NX;
    zR = pi* a**2/lamb;
    stps = LZ/dz;

    cx = dx*(np.linspace(0,NX-1,NX)-NX/2) #real space coordinates
    #definee transverse wavenumbers
    kx = np.zeros(NX)
    #max prop angle
    k = 0
    while k < NX/2 + 1:
        kx[k] = dk*k
        k = k+1
    while k < NX :
        kx[k] = dk*(k-NX)
        k = k+1
    #make initial condition
    am0 = np.ones((NX,NX),dtype = "complex")
    for x in range(0,int(NX*.55)):
        for y in range(NX):
            val =  (cx[x]**2 + cx[y]**2)
            if val < a**2:
                am0[x,y]=val
    return(k0,kx,zR,stps,cx,am0)
    
    


    