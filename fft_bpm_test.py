#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 09:49:00 2019

@author: kirahart
"""

import numpy as np
from FFT_1D import paraxial_propagator, non_paraxial_propagator,linear_step
from bpm_initial_conditions import GaussianBeam1DRotated, GaussianBeam_BPM_FFT
import matplotlib.pyplot as plt

# beam parameters
lmb           = 800.0e-09;
beamwaist     = 1.0e-06;
angle         = 0.0;

# computational domain parameters
LX     = 1.0e-03;
NX     = 1024*8;

k0,kx,zR,cx = GaussianBeam_BPM_FFT(lmb,beamwaist,angle,LX,NX)

#integrating steps
nsteps   = 10;
distance = 100*zR;
dz       = 2*distance/nsteps;

am0 = GaussianBeam1DRotated( cx, -distance, beamwaist, k0, angle)
am1 = am0
propagator = non_paraxial_propagator(kx,k0,dz)

for s in range(1,nsteps):
    am1 = linear_step(am1,propagator)
    
#expected solution
amt  = GaussianBeam1DRotated(cx,distance,beamwaist,k0,angle);    
    
#plt.plot(cx, np.abs(am0))
#plt.plot(cx, np.abs(am1))
plt.plot(cx, np.abs(amt))
plt.show()

