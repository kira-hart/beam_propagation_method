#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 15:35:55 2019
Implementation of 2D ADI Solver
@author: kirahart
"""

import numpy as np
import matplotlib.pyplot as plt
import cmath
import h5py
import argparse
from collections import namedtuple

#load other written functions
from TDM import TDMsolve, TDMmultiply
from bpm_initial_conditions import GaussianBeam_BPM_FFT2D

def get_coords(NX,NY,LX,LY):
   dimX = NX;
   dimY = NY;
   cx  = np.zeros(dimX, dtype='complex')
   cy  = np.zeros(dimY, dtype='complex')
   i=0
   while i < dimX:
       cx[i] = (LX/dimX)*(i - dimX/2.0);
       cy[i] = (LY/dimY)*(i - dimY/2.0);
       i = i+1
   return(cx,cy)
  


def setstep(dx, dy, dz,k0):
  ideltaX = (0+1j)*dz/(4.0*k0*dx*dx);
  ideltaY = (0+1j)*dz/(4.0*k0*dy*dy);
  return(ideltaX, ideltaY)

def ADI_BPM_solver(LX, NX, LY,NY,ideltaX,ideltaY):
   '''allocation'''
    # set up coordinate grid
   dimX = NX;
   dimY = NY;
   cx  = np.zeros(dimX, dtype='complex')
   cy  = np.zeros(dimY, dtype='complex')
   i=0
   while i < dimX:
       cx[i] = (LX/dimX)*(i - dimX/2.0);
       cy[i] = (LY/dimY)*(i - dimY/2.0);
       i = i+1
  
   #grid spacing
   dx = cx[1]-cx[0]
   dy = cy[1]-cy[0]
    
   #2D array holders
   num2D=6;
   Holder2D = np.zeros(dimX*dimY*num2D,dtype='complex');
   E00 = Holder2D[0*dimX*dimY:1*dimX*dimY]#current amplitude
   E12 = Holder2D[1*dimX*dimY:2*dimX*dimY]#half-step amplitude
   E11 = Holder2D[2*dimX*dimY:3*dimX*dimY] #full-step, new amplitde
   R00 = Holder2D[3*dimX*dimY:4*dimX*dimY] #first half-step right-hand-side
   R12 = Holder2D[4*dimX*dimY:5*dimX*dimY]  #second half-step right-hand-side
   
   # 1D array holders fir C- N matricies
   num1D=7
   Holder1DX = np.zeros(dimX*num1D, dtype='complex')
   Holder1DY = np.zeros(dimY*num1D, dtype='complex')

   #matrix L+ along X
   LPAX = Holder1DX[0*dimX:1*dimX]
   LPBX = Holder1DX[1*dimX:2*dimX]
   LPCX = Holder1DX[2*dimX:3*dimX]
   #matrix L- along X
   LMAX = Holder1DX[3*dimX:4*dimX]
   LMBX = Holder1DX[4*dimX:5*dimX]
   LMCX = Holder1DX[5*dimX:6*dimX]
  
   #matrix L+ along Y
   LPAY = Holder1DY[0*dimY:1*dimY]
   LPBY = Holder1DY[1*dimY:2*dimY]
   LPCY = Holder1DY[2*dimY:3*dimY]
   #matrix L- along Y
   LMAY = Holder1DY[3*dimY:4*dimY]
   LMBY = Holder1DY[4*dimY:5*dimY]
   LMCY = Holder1DY[5*dimY:6*dimY]
   #source vectors and targets
   auxX = Holder1DX[6*dimX:7*dimX]
   auxY = Holder1DY[6*dimY:7*dimY]
 
   i = 0
   while i < dimX:
       auxX[i] = 1.0 + 2.0*ideltaX
       LPAX[i]=+1.0*ideltaX
       LPBX[i]=1.0 -2.0*ideltaX
       LPCX[i] = 1.0*ideltaX
       LMAX[i]=-1.0*ideltaX
       LMBX[i]=1.0 + 2.0*ideltaX
       LMCX[i] = -1.0*ideltaX
       i = i + 1
   i = 0
   while i < dimY:
       auxY[i] =  1.0+2.0*ideltaY
       LPAY[i] = +1.0*ideltaY
       LPBY[i] = 1.0 -2.0*ideltaY
       LPCY[i] = 1.0*ideltaY
       LMAY[i] =-1.0*ideltaY
       LMBY[i] = 1.0+2.0*ideltaY
       LMCY[i] = -1.0*ideltaY
       i = i + 1
   return( Holder2D ,Holder1DX, Holder1DY )

def SetAmplitude(A0,Holder2D,dimX,dimY):
    y = 0
    while y < dimY:
        x = 0
        while x < dimX:
            Holder2D[x+y*dimX] = A0[x][y]
            x = x+1
        y = y+1
  

         
def adi_step(Holder2D ,Holder1DX, Holder1DY,dimX,dimY):
   E00 = onetoteoD(Holder2D[0*dimX*dimY:1*dimX*dimY],dimX,dimY)#current amplitude
   E12 = onetoteoD(Holder2D[1*dimX*dimY:2*dimX*dimY],dimX,dimY)#half-step amplitude
   E11 =onetoteoD( Holder2D[2*dimX*dimY:3*dimX*dimY],dimX,dimY) #full-step, new amplitde
   R00 = onetoteoD(Holder2D[3*dimX*dimY:4*dimX*dimY],dimX,dimY) #first half-step right-hand-side
   R12 = onetoteoD(Holder2D[4*dimX*dimY:5*dimX*dimY],dimX,dimY)  #second half-step right-hand-side
   
   
   LPAX = Holder1DX[0*dimX:1*dimX]
   LPBX = Holder1DX[1*dimX:2*dimX]
   LPCX = Holder1DX[2*dimX:3*dimX]
   LMAX = Holder1DX[3*dimX:4*dimX]
   LMBX = Holder1DX[4*dimX:5*dimX]
   LMCX = Holder1DX[5*dimX:6*dimX]
   LPAY = Holder1DY[0*dimY:1*dimY]
   LPBY = Holder1DY[1*dimY:2*dimY]
   LPCY = Holder1DY[2*dimY:3*dimY]
   LMAY = Holder1DY[3*dimY:4*dimY]
   LMBY = Holder1DY[4*dimY:5*dimY]
   LMCY = Holder1DY[5*dimY:6*dimY]
   auxX = Holder1DX[6*dimX:7*dimX]
   auxY = Holder1DY[6*dimY:7*dimY]
    
   # Multiply R00 = LPX * E00
   y = 0
   while y < dimY:    
       R00[y]=TDMmultiply(LPAX,LPBX,LPCX, E00[y]);
       y = y +1 
   #Solve LMY*E12 = R00
   x= 0
   while x < dimX:   
       E12[: ,x]=TDMsolve(LMAY,LMBY,LMCY, R00[: ,x])
       x = x +1 
   #replace destroyed diagonal
   y = 0
   while y < dimY:    
       LMBY[y*dimX:y*dimX+dimY]= auxY[y*dimX:y*dimX+dimY]
       y = y +1 
   x= 0
   while x < dimX:    
       R12[:, x]=TDMmultiply(LPAY,LPBY,LPCY, E12[:, x]) 
       x = x +1 
   y = 0
   while y < dimY:    
       E11[y]=TDMsolve(LMAX,LMBX,LMCX, R12[y])
       y = y +1 
   x = 0
   while x < dimX:    
       LMBX[x*dimY: x*dimY+dimX]= auxX[x*dimY: x*dimY+dimX]
       x = x +1 
   
   Enew = np.ravel(E11)
   return( Enew )
 

def onetoteoD(Holder2D,NX,NY):
    AE=np.zeros((NX,NY),dtype='complex')
    y = 0
    while y < NY:
        x = 0
        while x < NX:
            AE[x][y]= Holder2D[x+NX*y] 
            x = x+1
        y = y+1
    return(AE)
  

def adi_propagate(LX, NX, LY,NY,LZ, NZ,k0,A0):
    dx = LX/NX
    dy = LY/NY
    dz = LZ/NZ
    ideltaX, ideltaY = setstep(dx, dy, dz, k0)
    Holder2D ,Holder1DX, Holder1DY = ADI_BPM_solver(LX, NX, LY,NY,ideltaX,ideltaY)
    SetAmplitude(A0,Holder2D,NX,NY)
    for i in range(NZ):
        Enew = adi_step(Holder2D ,Holder1DX, Holder1DY,NX,NY)
        Holder2D ,Holder1DX, Holder1DY = ADI_BPM_solver(LX, NX, LY,NY,ideltaX,ideltaY)
        Holder2D[0:1*NX*NY] = Enew
    return(onetoteoD(Holder2D,NX,NY))

    