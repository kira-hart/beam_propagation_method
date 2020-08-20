#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 12:26:20 2019

@author: kirahart
"""

import numpy as np
from FFT_2D import fft2d_bpm, propagate_gaurd, linear_step,propagator,boundary_gaurd
from bpm_initial_conditions import  GaussianBeam_BPM_FFT2D ,spot_of_arago
import matplotlib.pyplot as plt

lamb = 633.0e-09;          
LX     = 5.0e-3;
NX     = 1024#4096;
dz     = 0.05;
dx     = LX/NX;
LZ     = 0.25;
beamwaist = 5e-4

k0,kx,zR,stps,cx,am0 = spot_of_arago(lamb,beamwaist,LX,LZ,NX,dz)
#GaussianBeam_BPM_FFT2D(lamb,beamwaist,LX,LZ,NX,dz)

b = boundary_gaurd(NX,LX,dx)
p = propagator(NX,kx,k0,dz)
#am1 = fft2d_bpm(p,am0)
am1 = propagate_gaurd(am0,5,p,b)


fig, axs = plt.subplots(1,2, sharex=False, sharey=False)
fig.suptitle('')
axs[0].matshow(np.real(am0))
axs[1].matshow(np.real(am1))
plt.show()



row = int(NX/2)
fig, axs = plt.subplots(4, sharex=False, sharey=False)
fig.suptitle('')
axs[0].plot(kx,np.real(p[row]))
axs[1].plot(kx,np.abs(np.fft.fft(am0[0:NX,int(NX/4)])))
axs[2].plot(cx,np.abs(am0)[row])
axs[3].plot(cx,np.abs(am1)[row])
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=10, hspace=.05)

plt.show()