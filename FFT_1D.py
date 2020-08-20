#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 10:19:47 2019
% 1-D beam propagator, FFT based method
@author: kirahart
"""
from cmath import exp , sqrt
import numpy as np


def paraxial_propagator(kx,k0,dz):
    px  = np.zeros(len(kx),dtype = 'complex')
    for i in range(len(kx)):
        px[i] = exp(-1j*(kx[i]*kx[i])/(2*k0)*dz + 1j*dz*k0)
    return(px)
    
def non_paraxial_propagator(kx,k0,dz):
    px  = np.zeros(len(kx),dtype = 'complex')
    for i in range(len(kx)):
        px[i] = 1j*dz*((sqrt(k0*k0 - kx[i]*kx[i]))) 
    return(px)
    
def linear_step(a0,propagator):
    a0fft=np.fft.fft(a0)
    arg = np.multiply(propagator,a0fft)
    return(np.fft.ifft(arg))