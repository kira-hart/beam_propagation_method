#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 12:04:48 2019

@author: kirahart
"""

import numpy as np
import matplotlib.pyplot as plt
import cmath
from numba import jit, f8 
import h5py
import argparse
from collections import namedtuple
import seaborn as sns
sns.set()
plt.rcParams['axes.grid'] = False

#load other written functions
from ADI_BPM_2d import onetoteoD,get_coords, adi_propagate
from bpm_initial_conditions import GaussianBeam_BPM_FFT2D,GaussianBeam_BPM_FFT


lamb = 800e-9
beamwaist = 0.5e-3
LX = 1e-2
LY = 1e-2
LZ = 1.5
NX = 200
NY = NX
NZ = 10
dz = LZ/NZ
dx = LX/NX
dy = LY/NY

k0,kx,zR,stps,cx,A0 = GaussianBeam_BPM_FFT2D(lamb,beamwaist,LX,LZ,NX,dz ,0)

cx,cy =  get_coords(NX,NY,LX,LY)
cx = np.real(cx)
cy = np.real(cy)

fig, ax = plt.subplots(1,1)
img= plt.imshow(np.real(A0),Cmap="RdBu_r")
plt.title("Initial Beam")
fig.colorbar(img)
plt.grid(False)
plt.xlabel("NX")
plt.ylabel("NY")
plt.show()






As= adi_propagate(LX, NX, LY,NY,LZ, NZ,k0,A0)

fig, ax = plt.subplots(1,1)
img= plt.imshow(np.real(As),Cmap="RdBu_r")
fig.colorbar(img)
plt.grid(False)
plt.title("ADI Propagated Beam")
plt.xlabel("NX")
plt.ylabel("NY")
plt.show()

plt.plot(cx,np.real(A0[int(NX/2),:]),c='black',label = "$E_0$")
plt.grid(True)
plt.plot(cx,np.real(As[int(NX/2),:]),c='red',label = "$E_{y=0}$")
plt.plot(cx,np.real(As[:,int(NX/2)]),c='blue',label = "$E_{x=0}$")
plt.plot(cx,np.imag(As[int(NX/2),:]),ls = '--',c='red')
plt.plot(cx,np.imag(As[:,int(NX/2)]),ls ='--',c='blue')
plt.legend(frameon=True,fontsize = 14)
plt.show()