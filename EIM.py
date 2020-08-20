#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 13:55:52 2019
Effective Index Method 
@author: kirahart
"""
import numpy as np
from numba import jit
import matplotlib.pyplot as plt

class geom:
    def __init__(self,NX,NY):
        self.NX    = NX
        self.NY    = NY
        self.n     = np.ones((NX,NY),dtype = "complex")
       
    def circle(self, R, nn):
        for x in range(self.NX):
            for y in range(self.NX):
                if np.sqrt((x-2*self.NX)**2 +(y-2*self.NY)**2) < R:
                    self.n[x][y] = nn
   
    def waveguide(self, LX, LY, Lclad,nc,nclad):
        #set up waveguide core
        xi = int(self.NX/2 - LX/2);
        xf = int(self.NX/2 + LX/2);
        yi=  int(self.NY/2 - LY/2);
        yf = int(self.NY/2 + LY/2);
        self.n[yi:yf, xi:xf] = nc
        
        #set up substrate
        yci = int(self.NY/2 + LY/2)-1;
        ycf = yci + Lclad
        self.n[yci:ycf, 0: self.NX] = nclad
        
    def block(self, xc, yc, Lx,Ly,n):
        #set up waveguide core
        xi = int(xc-Lx/2);
        xf = int(xc+Lx/2);
        yi=  int(yc - Ly/2);
        yf = int(yc + Ly/2);
        self.n[yi:yf, xi:xf] = n
        
    
    def plot(self,LX,LY,cma,cQ,showQ):
        NX = np.linspace(-LX/2,LX/2,self.NX)
        NY = np.linspace(-LY/2,LY/2,self.NY)
        if cQ:
            plt.contour(NX,NY,np.abs(self.n),cmap = cma)
        else:
            plt.contourf(NX,NY,np.abs(self.n),cmap = cma)
            plt.colorbar()
            plt.title("Refractive Index Map")
        plt.xlabel("X [$\mu m$]")
        plt.ylabel("Y [$\mu m$]")
        if showQ:
            plt.show()
        
    def applyPEC(self):
        print("PEC applied")
        