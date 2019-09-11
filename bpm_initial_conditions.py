#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 10:14:55 2019
Inital Confitions for FFT Methods
@author: kirahart
"""

import numpy as np
import cmath
from cmath import sin,cos, sqrt, exp

def GaussianBeam1DRotated( x, z, w, k, angle):
    zt = +cos(angle)*z + sin(angle)*x;
    xt = -sin(angle)*z + cos(angle)*x;
    aux = 1.0 + 2.0* cmath.i*zt/(k*w**2);
    return(cmath.exp( -xt * xt/(w**2*aux) + cmath.i *k*zt)/cmath.sqrt(aux))
