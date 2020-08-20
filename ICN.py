#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 11:40:47 2019
Creates Field Class for Iterated CN
@author: kirahart
"""
import numpy as np
from numba import jit
from TDM import TDMmultiply, TDMsolve
import matplotlib.pyplot as plt


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

class Field():
    def __init__(self,NX,NY,dx,dy,dz,kr,k0):
        self.NX    = NX
        self.NY    = NY
        self.dx    = dx
        self.dy    = dy
        self.dz    = np.abs(dz)
        self.kr    = kr
        self.k0    = k0
        self.nref  = 1.0 
        
        self.ideltaX = (0+1j)*dz/(4.0*k0*dx*dx);
        self.ideltaY = (0+1j)*dz/(4.0*k0*dy*dy);
        
        self.data2d   = np.zeros((9,NX,NY),dtype = "complex")
        '''X and Y components are stored adjacent'''
        
        self.coeff_xx = 1j*self.dz/(2.0* self.k0 *self.nref *self.dx *self.dx); 
        self.coeff_yy = 1j*self.dz/(2.0*self.k0*self.nref*self.dy*self.dy); 
        self.coeff_di = 1j*self.dz*self.k0/(2.0*self.nref);
        self.coeff_xy = 1j*self.dz/(8.0*self.k0*self.nref*self.dx*self.dy); 

        
       
    def initial_field(self, Ei,pol):
        #use True for x-pol, False for y-pol
        if pol:
            for x in range(self.NX):
                for y in range(self.NY):
                    self.data2d[0][x][y] = Ei[x][y]
                    self.data2d[7][x][y] = Ei[x][y]
        else:
            for x in range(self.NX):
                for y in range(self.NY):
                    self.data2d[1][x][y] = Ei[x][y]
                    self.data2d[8][x][y] = Ei[x][y]
   
    
    """---------Define A Matrix Application-------"""
    def ApplyAXX(self,n,src,tar):
        epsilon = np.square(n)
        T = self.data2d[tar]
        S = self.data2d[src]
       #zeroing recieving arrays 
        x = 0
        while x < self.NX:
            T[0,x] = 0.0
            T[self.NY-1,x] = 0.0
            x = x+1
        y = 0
        while y < self.NY:
            T[y,0] = 0.0
            T[y, self.NX-1] = 0.0
            y= y+1
        #set epsilons 
        y = 1;
        while y < self.NY -1:
            x = 1; 
            while x < self.NX-1:
                epsL = epsilon[x-1,y]
                epsC = epsilon[x,  y]
                epsR = epsilon[x+1,y]
                
                auxL = (epsL + epsC)/(2.0*epsC)
                auxC = (epsL + epsC)/(2.0*epsL) + (epsC + epsR)/(2.0*epsR)
                auxR = (epsR + epsC)/(2.0*epsC)
                
                T[x,y]  =  self.coeff_xx*(auxR*S[x+1,y] - auxC*S[x,y] + auxL*S[x-1,y])
                T[x,y]= T[x,y]+  self.coeff_yy*(S[x,y+1] - 2.0*S[x,y] + S[x,y-1]);
                T[x,y]= T[x,y]+  self.coeff_xx*(S[x+1,y] - S[x-1,y] );
                T[x,y] +=  self.coeff_yy*(S[x,y+1] - S[x,y-1] );
                T[x,y] +=  self.coeff_di*(epsC - self.nref*self.nref)*S[x,y]

                x=x+1
            y=y+1
        
        
    def ApplyAYY(self,n,src,tar):
        epsilon = np.square(n)
        T = self.data2d[tar+1]
        S = self.data2d[src+1]
       #zeroing recieving arrays 
        x = 0
        while x < self.NX:
            T[x,0] = 0.0
            T[x,self.NY-1] = 0.0
            x = x+1
        y = 0
        while y < self.NY:
            T[0,y] = 0.0
            T[self.NX-1,y] = 0.0
            y= y+1
        #set epsilons 
        y = 1;
        while y < self.NY -1:
            x = 1; 
            while x < self.NX-1:
                epsL = epsilon[x,y-1]
                epsC = epsilon[x,  y]
                epsR = epsilon[x,y+1]
                
                auxL = (epsL + epsC)/(2.0*epsC)
                auxC = (epsL + epsC)/(2.0*epsL) + (epsC + epsR)/(2.0*epsR)
                auxR = (epsR + epsC)/(2.0*epsC)
                
                T[x,y]=   self.coeff_xx*(S[x+1,y]- 2.0*S[x,y]+S[x-1,y])
                T[x,y]= T[x,y]+ self.coeff_yy*(auxR*S[x,y+1] - auxC*S[x,y] + auxL*S[x,y-1]);
                T[x,y]= T[x,y]+  self.coeff_xx*(S[x+1,y]- S[x-1,y]);
                T[x,y]= T[x,y]+  self.coeff_yy*( S[x,y+1] - S[x,y-1]);
                T[x,y]= T[x,y]+ self.coeff_di*(epsC - self.nref*self.nref)*S[x,y];
                
                x=x+1
            y=y+1

    def ApplyAXY(self,n,src,tar):
        epsilon = np.square(n)
        T = self.data2d[tar]
        S = self.data2d[src+1]
       #zeroing recieving arrays 
        x = 0
        while x < self.NX:
            T[x,0] = 0.0
            T[x,self.NY-1] = 0.0
            x = x+1
        y = 0
        while y < self.NY:
            T[0,y] = 0.0
            T[self.NX-1,y] = 0.0 
            y = y+1
        y = 1;
        while y < self.NY -1:
            x = 1; 
            while x < self.NX-1:
                epsPP = epsilon[x+1,y+1];
                epsPM = epsilon[x+1,y-1];
                epsMP = epsilon[x-1,y+1];
                epsMM = epsilon[x-1,y-1];
                epsP0 = epsilon[x+1,y  ];
                epsM0 = epsilon[x-1,y  ];
                
                T[x,y]= T[x,y]+((epsPP/epsP0 - 1.0)*S[x+1,y+1]-(epsMP/epsM0 - 1.0)*S[x-1,y+1] -(epsPM/epsP0 - 1.0)*S[x+1,y-1]+(epsMM/epsM0 - 1.0)*S[x-1,y-1])*self.coeff_xy;
                
                x=x+1
            y=y+1
    
    def ApplyAYX(self,n,src,tar):
        epsilon = np.square(n)
        T = self.data2d[tar+1]
        S = self.data2d[src]
       #zeroing recieving arrays 
        x = 0
        while x < self.NX:
            T[x,0] = 0.0
            T[x,self.NY-1] = 0.0
            x = x+1
        y = 0
        while y < self.NY:
            T[0,y] = 0.0
            T[self.NX-1,y] = 0.0 
            y = y+1
        y = 1;
        while y < self.NY -1:
            x = 1;
            while x < self.NX-1:
                epsPP = epsilon[x+1,y+1];
                epsPM = epsilon[x+1,y-1];
                epsMP = epsilon[x-1,y+1];
                epsMM = epsilon[x-1,y-1];
                eps0P = epsilon[x , y+1];
                eps0M = epsilon[x  ,y-1];
                
                T[x,y]= T[x,y]+((epsPP/eps0P - 1.0)*S[x+1,y+1] -(epsMP/eps0P - 1.0)*S[x-1,y+1] -(epsPM/eps0M - 1.0)*S[x+1,y-1]+(epsMM/eps0M - 1.0)*S[x-1,y-1])*self.coeff_xy;
                
                x=x+1
            y=y+1
    

    def ApplyRHS(self,n,src,tar):
        self.ApplyAXX(n,src,tar)
        self.ApplyAYY(n,src,tar)
        self.ApplyAXY(n,src,tar)
        self.ApplyAYX(n,src,tar)
        
    def Add(self, src,tar,coeff):
        TX = self.data2d[tar]
        SX = self.data2d[src]
        TY = self.data2d[tar+1]
        SY = self.data2d[src+1]
        for x in range(self.NX):
                for y in range(self.NY):
                    TX[x,y] += SX[x,y]*coeff
                    TY[x,y] += SY[x,y]*coeff
        
        
        
    '''------------Propagation Functions-------'''
    def step(self,n):
        self.ApplyRHS(n,0,2)
        self.Add(2,0,1.0)
        
        self.ApplyRHS(n,2,4)
        self.Add(4,0,0.5)
        
        self.ApplyRHS(n,4,2)
        self.Add(2,0,0.25)
        
        
    '''------------ADI Solver-------'''
    def ADI_BPM_solver(self,field,n):
       '''allocation'''
        # set up coordinate grid
       dimX = self.NX;
       dimY = self.NY;
     
       #grid spacing
       dx = self.dx
       dy = self.dy
       ideltaX= self.ideltaX
       ideltaY= self.ideltaY
       k0 = self.k0
       dz = self.dz
        
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
       
        
       E00 = onetoteoD(Holder2D[0*dimX*dimY:1*dimX*dimY],dimX,dimY)#current amplitude
       E00 = self.data2d[field]
       E12 = onetoteoD(Holder2D[1*dimX*dimY:2*dimX*dimY],dimX,dimY)#half-step amplitude
       E11 =onetoteoD( Holder2D[2*dimX*dimY:3*dimX*dimY],dimX,dimY) #full-step, new amplitde
       R00 = onetoteoD(Holder2D[3*dimX*dimY:4*dimX*dimY],dimX,dimY) #first half-step right-hand-side
       R12 = self.data2d[field]
       # Multiply R00 = LPX * E00
       '''y = 0
       while y < dimY:    
           R00[y]=TDMmultiply(LPAX,LPBX,LPCX, E00[y]);
           R00[y]=  R00[y]*np.exp(0.5j*k0*(1-n[y])*dz)
           y = y +1 '''
       #Solve LMY*E12 = R00
       x= 0
       while x < dimX:   
           E12[: ,x]=TDMsolve(LMAY,LMBY,LMCY, R00[: ,x])
           x = x +1 
       #replace destroyed diagonal
       '''y = 0
       while y < dimY:    
           LMBY[y*dimX:y*dimX+dimY]= auxY[y*dimX:y*dimX+dimY]
           y = y +1 
       x= 0
       while x < dimX:    
           R12[:, x]=TDMmultiply(LPAY,LPBY,LPCY, E12[:, x]) 
           R12[:,x]= R12[:,x]*np.exp(0.5j*k0*(1-n[:,x])*dz)
           x = x +1 '''
       y = 0
       while y < dimY:    
           E11[y]=TDMsolve(LMAX,LMBX,LMCX, R12[y])
           y = y +1 
       x = 0
       while x < dimX:    
           LMBX[x*dimY: x*dimY+dimX]= auxX[x*dimY: x*dimY+dimX]
           x = x +1 
       
       Enew = E11
       self.data2d[field] = Enew
       self.data2d[field+7] = Enew
       
    def ADI_BPM_solver_nonvec(self,field,n):
       '''allocation'''
        # set up coordinate grid
       dimX = self.NX;
       dimY = self.NY;
     
       #grid spacing
       dx = self.dx
       dy = self.dy
       ideltaX= self.ideltaX
       ideltaY= self.ideltaY
       k0 = self.k0
       dz = self.dz
        
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
       
        
       E00 = onetoteoD(Holder2D[0*dimX*dimY:1*dimX*dimY],dimX,dimY)#current amplitude
       E00 = self.data2d[field]
       E12 = onetoteoD(Holder2D[1*dimX*dimY:2*dimX*dimY],dimX,dimY)#half-step amplitude
       E11 =onetoteoD( Holder2D[2*dimX*dimY:3*dimX*dimY],dimX,dimY) #full-step, new amplitde
       R00 = onetoteoD(Holder2D[3*dimX*dimY:4*dimX*dimY],dimX,dimY) #first half-step right-hand-side
       R12 = onetoteoD(Holder2D[4*dimX*dimY:5*dimX*dimY],dimX,dimY)  #second half-step right-hand-side
       # Multiply R00 = LPX * E00
       y = 0
       while y < dimY:    
           R00[y]=TDMmultiply(LPAX,LPBX,LPCX, E00[y]);
           R00[y]=  R00[y]*np.exp(0.5j*k0*(1-n[y])*dz)
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
           R12[:,x]= R12[:,x]*np.exp(0.5j*k0*(1-n[:,x])*dz)
           x = x +1 
       y = 0
       while y < dimY:    
           E11[y]=TDMsolve(LMAX,LMBX,LMCX, R12[y])
           y = y +1 
       x = 0
       while x < dimX:    
           LMBX[x*dimY: x*dimY+dimX]= auxX[x*dimY: x*dimY+dimX]
           x = x +1 
       
       Enew = E11
       self.data2d[field] = Enew
       
     
       
       
       
    
    
    
    
    
    
    