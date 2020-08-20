#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 11:17:59 2019

@author: kirahart
"""

import numpy as np
import matplotlib.pyplot as plt
from TDM import CreateDiagonals,tridiag, TDMsolve, TDMmultiply, CrankNicolson
from bpm_initial_conditions import GaussianBeam1DRotated
import timeit
from math import pi
from mpl_toolkits.axes_grid1 import make_axes_locatable

#tunable parameters
lamb   = 800e-9;          
w0     = 2.0e-06;         
xrange = 100*lamb
nx     = 200

# calculated parameters
deltax  = xrange/nx
deltaz = deltax
k0     = 2*pi/lamb;  


###-------------Task 1
'''a,b,c = CreateDiagonals(nx,deltax,deltaz,k0,1)
diag = tridiag(a,b,c)
a,b,c = CreateDiagonals(nx,deltax,deltaz,k0,-1)
diag2 = tridiag(a,b,c)

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)

fig, (ax1, ax2) = plt.subplots(ncols=2)
img1 = ax1.imshow(np.imag(diag), cmap = 'bone')
ax1.title.set_text("Im{$L^{+}$}")
colorbar(img1)
img2 = ax2.imshow(np.imag(diag2), cmap = 'bone')
ax2.title.set_text("Im{$L^{-}$}")
colorbar(img2)
plt.tight_layout(h_pad=1)
plt.show()'''
  

#############------TASK 2 
'''see tdm_solve_test'''

######-----------TASK 3
'''see tdm_solve_test'''

#####-----------Task 4



####-----------Task 5
#tunable parameters
lamb   = 800e-9;          
w0     = 5.0e-06;         
xrange = 500e-06
nx     = 1000

# calculated parameters
deltax  = xrange/nx
deltaz = deltax
k0     = 2*pi/lamb;  

x = np.linspace(-nx*deltax/2,nx*deltax/2,nx)
x_um = np.multiply(x,1e6)
E0 = GaussianBeam1DRotated( x, 0, w0, k0, 0)

z = 500e-6
Ef = GaussianBeam1DRotated( x, z, w0, k0, 0)


E1 = E0
for i in range(int(z/deltaz)):
    E1 = CrankNicolson(deltax,k0,E1,deltaz)

#plt.plot(x_um,np.real(E0),'black',label = "Analytic $E_0$")
plt.plot(x_um,np.real(E1),'red',label ="CN $E_f$")
plt.plot(x_um,np.real(Ef),'black',ls = '--',label ="Analytic $E_f$")
plt.legend(frameon=True)
plt.xlabel('x $\mu m$')
plt.ylabel('Re[E]')
plt.xlim(-20,150)
plt.show()

#plt.plot(x_um,np.imag(E0),'black',label = "Analytic $E_0$")
plt.plot(x_um,np.imag(E1),'red',label ="CN $E_f$")
plt.plot(x_um,np.imag(Ef),'black',ls = '--',label ="Analytic $E_f$")
plt.legend(frameon=True)
plt.xlabel('x $\mu m$')
plt.ylabel('Im[E]')
#plt.xlim(-.01,0.01)
plt.show()


plt.plot(x_um,np.abs(E0),'black',label = "Analytic $E_0$")
plt.plot(x_um,np.abs(E1),'red',label ="CN $E_f$")
plt.plot(x_um,np.abs(Ef),'black',ls = '--',label ="Analytic $E_f$")
plt.legend(frameon=True)
plt.xlabel('x $\mu m$')
plt.ylabel('Abs[E]')
#plt.xlim(-1,20)
plt.show()


