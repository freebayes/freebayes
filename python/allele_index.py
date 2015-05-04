# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 16:36:12 2015

@author: dima
"""
import numpy as np
import math

def stirling1(x, m):
    out = 1
    for y in range(x, m+x):
        out = out * y
        # print(y, out)
    return out

def ai(*args):
    np.sort(args)
    out = 0
    for m,k in enumerate(args):
        # print(m, k)
        out += stirling1(k, m+1) / math.factorial(m+1)
    return int(out)
        

ai(0,0)
ai(0,1)
ai(1,1)
ai(0,2)
ai(1,2)
ai(2,2)

ai(0,0,1)

def inverse_ai(ind, ploidy):
    x_m = 1/2*(-1 + np.sqrt(1+ 8 *ind) )
    k_last = np.floor(x_m)
    k_last = np.floor(x_m)
    
    return x_m, ind - np.floor(x_m)
    