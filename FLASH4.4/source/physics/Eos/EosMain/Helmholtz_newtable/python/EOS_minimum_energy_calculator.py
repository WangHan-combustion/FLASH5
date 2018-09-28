#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 16:59:04 2018

@author: alexgrannan
"""

from math import *
import numpy as np
from numpy import ma
import h5py
import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm
from datetime import datetime

## Constants
avo     = 6.0221417930e23
kerg    = 1.3806504240000000e-16
me      = 9.1093821545e-28
clight  = 2.99792458e10
mecc    = me * clight * clight
h       = 6.6260689633e-27
hbar    = h/2/pi
hion    = 13.605698140
ev2erg  = 1.60217648740e-12
rt2     = 1.4142135623730951
rt3     = 1.7320508075688772
rtpi    = 1.7724538509055159
cpf0    = h/(me*clight)
cpf1    = 3.0/(8.0*pi) * cpf0**3
cpf2    = 4.0/cpf1
cpf3    = 2.0*rt3*rtpi/(rt2*cpf1)
twoth   = 2.0/3.0
fa0     = 64.0/(9.0*pi)
forpi   = 4.0 * pi
xconst  = 8*pi*sqrt(2)*(me**3)*(clight**3)/(h**3)
econst  = xconst * me * clight**2
pconst  = xconst * 2.0 * mecc / 3.0
qe      = 4.8032042712e-10
esqu    = qe*qe
ssol    = 5.6704e-5
asol    = 4.0 * ssol / clight



def EOS_minimum_energy_calculator(den,Ye,sumy):
    
    if type(den) is not np.ndarray:
        den = den*np.ones(1)
    
    if type(sumy) is not np.ndarray:
        sumy = sumy*np.ones(1)
        
    if type(Ye) is not np.ndarray:
        Ye = Ye*np.ones(1)    
    
    
    e_min = np.ones((len(den),len(Ye),len(sumy)))
     
    for j in range(0,len(den)):
        for k in range(0,len(Ye)):
            for l in range(0,len(sumy)):
                
                num_density_0 = 8*pi*(me*clight/h)**3.0/3
            
                num_density_ele = den[j]*avo*Ye[k]
            
                x0 = (num_density_ele/num_density_0)**(1.0/3.0)
            
                fermi_momentum = x0*me*clight
            
                den_trans = 8*pi*(2.5*me*clight)**3.0/3/h**3/Ye[k]/avo
                c1 = num_density_0*mecc
                
##                if (x0 >= 1.0):
#                e_min_1[j,k,l] = c1*(-8*x0**3 \
#                        + 3*sqrt(x0**2 + 1)*(2*x0**3 \
#                        + x0) - 3*log(x0 \
#                        + sqrt(1 + x0**2)))/8/den[j]
#                else:
#                    e_min_1[j,k,l] = num_density_0*mecc*(3*x0**5/10)/den[j] 
                
                c2 = pi*me**4.0*clight**5.0/3.0/h**3.0
                e_min[j,k,l] = c2*(8.0*x0**3.0*(sqrt(x0**2.0 + 1.0) - 1.0) -\
                (x0*(2.0*x0**2.0 - 3.0)*sqrt(x0**2.0 + 1.0) + 3.0*asinh(x0)))/den[j]
    
# Analytical Solution   
#    e_min = num_density_0*mecc*(-8*fermi_momentum_norm**3 \
#            + 3*sqrt(fermi_momentum_norm**2 + 1)*(2*fermi_momentum_norm**3 \
#            + fermi_momentum_norm) - 3*log(fermi_momentum_norm \
#            + sqrt(1 + fermi_momentum_norm**2)))/8/den
                            
    return (e_min)


def EOS_minimum_energy_calculator_test(den,Ye,sumy):
    
    e_min_1 = EOS_minimum_energy_calculator(den,Ye,sumy)
    
    return(e_min)
abar_lo = 1.0
abar_hi = 32.0
sumy_range = np.linspace(1.0/abar_hi, 1.0/abar_lo)

Ye_lo = 0.1
Ye_hi = 1.0
Ye_range = np.linspace(Ye_lo,Ye_hi)

dlo = log10(1.0e-15)
dhi = log10(1.0e15)
d_range=np.logspace(dlo, dhi)

e_min = EOS_minimum_energy_calculator(d_range,Ye_range,sumy_range)

#num_density_0 = 8*pi*(me*clight/h)**3.0/3
#    
#num_density_ele = d_range*Ye*avo
#    
#fermi_momentum_norm = (num_density_ele/num_density_0)**(1.0/3.0)
#
#fermi_momentum = fermi_momentum_norm*me*clight
#        
##e_min = num_density_0*mecc*(-8*fermi_momentum_norm**3 \
##        + 3*sqrt(fermi_momentum_norm**2 + 1)*(2*fermi_momentum_norm**3 \
##        + fermi_momentum_norm) - 3*log(fermi_momentum_norm \
##        + sqrt(1 + fermi_momentum_norm**2)))/8/den
#
#e_min_1 = num_density_0*mecc*(3*np.sqrt(fermi_momentum_norm**2 + 1)*(2*fermi_momentum_norm**3 \
#        + fermi_momentum_norm) - 3*np.arcsinh(fermi_momentum_norm))/8/d_range
#
#e_min_2 = num_density_0*mecc*(-8*fermi_momentum_norm**3 \
#           + 3*np.sqrt(fermi_momentum_norm**2 + 1)*(2*fermi_momentum_norm**3 \
#           + fermi_momentum_norm) - 3*np.arcsinh(fermi_momentum_norm))/8/d_range
#           
#e_min_3 = num_density_0*mecc*0.3*(fermi_momentum_norm**5)/d_range
#e_min_4 = num_density_0*mecc*0.75*(fermi_momentum_norm**4.0)/d_range
#         
#plt.figure(1)
##plt.plot(d_range_log, e_min_1,label='w/ Rest Mass')
#plt.plot(d_range, e_min_1[:,-1,0],label='emin1')
#plt.plot(d_range, e_min_2[:,-1,0],label='emin2')
##plt.plot(d_range, e_min_3,label='NR Limit')
##plt.plot(d_range, e_min_4,label='UR Limit')
#plt.yscale('log')
#plt.xscale('log')
#plt.legend()
#plt.grid()
   