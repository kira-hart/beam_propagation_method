#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 16:56:04 2019

@author: kirahart
"""

import numpy as np
import matplotlib.pyplot as plt
from cmath import exp,sqrt
import math

def fft2D(arr):
    a_fft=np.fft.fft2(arr)
    return(np.fft.fftshift(a_fft))
    
def turkey_window(N,a):
    window = np.ones(N)
    n = int(0)
    while n < a*N/2-1:
        window[n]=0.5*(1+ math.cos(math.pi*(n * 2 /(a*N)-1)))
        n = n+1
    n =int( N*(1-a/2))
    while n < N-1:
        window[n]=0.5*(1+math.cos(math.pi*(n * 2 /(a*N)-2/a+1)))
        n = n+1
    return(window)

def fft_window(arr,l,n):
    if l >1 :
        l = 1
    if l < 0 :
        l = 0
    N = len(arr)
    window = turkey_window(N,l)
    new = np.zeros((N,n))
    for i in range(N):
        new[i]=np.multiply(arr[i],window[i])
    return(new)
    
def fft2d_bpm(propagator,am0):
    am0_f = np.fft.fft2(am0)
    new = np.dot(am0_f,propagator)
    return(np.fft.ifft2(new))
    

def boundary_gaurd(NX,LX,dx):
    bxy = np.zeros((NX,NX),dtype = "complex")
    cx = dx*(np.linspace(0,NX-1,NX)-NX/2) #coordinates
    for x in range(NX):
        for y in range(NX):
            bxy[x,y] = exp(-( (cx[x]**2 + cx[y]**2)/(LX*LX/5))**8 )
    return(bxy)
    
def linear_step(a0,propagator):
    a0fft=np.fft.fft2(a0)
    arg = np.multiply(propagator,a0fft)
    return(np.fft.ifft2(arg))
    
def propagator_paraxial(NX,kx,k0,dz):
    pxy = np.zeros((NX,NX),dtype = "complex")
    for x in range(NX):
        for y in range(NX):
            pxy[x,y] = exp(-1j*(kx[x]**2 + kx[y]**2)/(2*k0)*dz )
    return(pxy)
    
def propagator(NX,kx,k0,dz):
    pxy = np.zeros((NX,NX),dtype = "complex")
    for x in range(NX):
        for y in range(NX):
            pxy[x,y] = exp(+1j*sqrt(k0**2 - (kx[x]**2 + kx[y]**2))*dz )
    return(pxy)
    
def propagate_gaurd(am0,steps,propagator,boundary):
   am1 = am0
   for i in range(steps):
        print("on step "+str(i)+"/"+str(steps))
        am1 =  linear_step(am1,propagator)
        am1 =  np.multiply(am1,boundary)
   return(am1)
   
def longitudinal_field(arr,NX,k0,kx):
    spatial_spectrum = np.fft.fft2(arr)
    operator = np.zeros((NX,NX),dtype = "complex")
    for x in range(NX):
        for y in range(NX):
            if k0**2 > kx[x]**2 + kx[y]**2:
                operator[x,y] = -kx[x]/sqrt(k0**2 - kx[x]**2 - kx[y]**2 )
    return(np.fft.ifft2(np.multiply(operator,spatial_spectrum)))
    


    