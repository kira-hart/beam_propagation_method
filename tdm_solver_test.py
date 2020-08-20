#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 14:42:33 2019

@author: kirahart
"""

import numpy as np
import matplotlib.pyplot as plt
from TDM import CreateDiagonals,tridiag, TDMsolve, TDMmultiply, CrankNicolson
from bpm_initial_conditions import GaussianBeam1DRotated
import time
from math import pi
from mpl_toolkits.axes_grid1 import make_axes_locatable


s = 100
a = np.random.rand(s-1)
b = np.random.rand(s)
c = np.random.rand(s-1)

diag = tridiag(a,b,c)
def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)

fig, (ax1) = plt.subplots(ncols=1)
img1 = ax1.imshow(diag, cmap = 'bone')
ax1.title.set_text("Random Diagonal Matrix")
colorbar(img1)
plt.tight_layout(h_pad=1)
plt.show()

din =  np.random.rand(s)

dout = TDMmultiply(a,b,c,din)

din_solved = TDMsolve(a,b,c,dout)

plt.plot(abs((din-din_solved)))
plt.title("Error in solved random $E$ using TDMsolver",fontsize=12)
plt.xlabel("index",fontsize=12)
plt.ylabel("$|E_{i}-E_{i, solved}|$",fontsize=12)
plt.show()

plt.plot(din,'.',c='black',label= "$E_{i}$")
plt.plot(dout,'-',c='black',label= "$E _{out}$")
plt.plot(din_solved,'--',c='red',label= "$E_{i,solved}$")
plt.title("Random vector retrieved using TDMSolver")
plt.xlabel("index")
plt.ylabel("E")
plt.legend()
plt.show()


#TIME TEST
TIMES=[]
grid_sizes = [5,10,100,1000,10000,100000,1000000]
for s in grid_sizes:
    start = time.time()
    for n in range(10):
        a = np.random.rand(s-1)
        b = 2* np.random.rand(s)
        c = np.random.rand(s-1)
        d = np.random.rand(s)
        mat = TDMsolve(a,b,c,d)
    end = time.time()
    print("for nx = " +str(s))
    print("time elapsed (sec) = " +str((end-start)/10))
    TIMES.append(abs(end-start)/10)
   
    
plt.plot(np.log(grid_sizes),TIMES)
plt.title('Time for TDMsolve' )
plt.xlabel('Grid Points $\log(n_x)$')
plt.ylabel('Time (sec)')
plt.show()



