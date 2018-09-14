#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 15:14:49 2018

@author: alexgrannan
"""

from math import *
import numpy as np
from numpy import ma
import h5py
import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm
from datetime import datetime

## Write Table Data
#hf = h5py.File('Entropic_EOS_Table.h5', 'w')
#
#g1 = hf.create_group('Variable_Table')
#
#g1.create_dataset('rho',data=d_mesh)
#g1.create_dataset('norm_energy',data=diff_e_mesh)
#
#g2 = hf.create_group('Entropy_Table')
#
#g2.create_dataset('S3',data=s3)
#g2.create_dataset('dS3drho',data=ds3dd)
#g2.create_dataset('dS3_2drho_2',data=ds3ddd)
#g2.create_dataset('dS3de3',data=ds3de)
#g2.create_dataset('dS3_2de3_2',data=ds3dee)
#g2.create_dataset('dS3_2de3drho',data=ds3ded)
#
#g3 = hf.create_group('Eta_Table')
#
#g3.create_dataset('eta',data=eta_totable)
#g3.create_dataset('detadrho',data=detadd_totable)
#g3.create_dataset('deta_2drho_2',data=detaddd_totable)
#g3.create_dataset('detade3',data=detade_totable)
#g3.create_dataset('deta_2de3_2',data=detadee_totable)
#g3.create_dataset('deta_2de3drho',data=detaded_totable)
#
#hf.close() 


## Read in File

filename = 'Entropic_EOS_Table.h5'
f = h5py.File(filename, 'r')

norm_energy_mesh = np.array(f.get('Variable_Table/norm_energy'))
norm_energy_lo = norm_energy_mesh[0,0,0,0]
norm_energy_hi = norm_energy_mesh[norm_energy_mesh.shape[0]-1,0,0,0]

den_mesh = np.array(f.get('Variable_Table/den'))
den_lo = den_mesh[0,0,0,0]
den_hi = den_mesh[0,den_mesh.shape[1]-1,0,0]

abar_mesh = np.array(f.get('Variable_Table/abar'))
abar_lo = abar_mesh[0,0,0,0]
abar_hi = abar_mesh[0,0,den_mesh.shape[2]-1,0]

zbar_mesh = np.array(f.get('Variable_Table/zbar'))
zbar_lo = zbar_mesh[0,0,0,0]
zbar_hi = zbar_mesh[0,0,0,den_mesh.shape[3]-1]

## Entropy and its derivatives
s = np.array(f.get('Entropy_Table/S'))

dsdd = np.array(f.get('Entropy_Table/dSdd'))
dsde = np.array(f.get('Entropy_Table/dSde'))
dsda = np.array(f.get('Entropy_Table/dSda'))
dsdz = np.array(f.get('Entropy_Table/dSdz'))

dsddd = np.array(f.get('Entropy_Table/dSddd'))
dsdde = np.array(f.get('Entropy_Table/dSdde'))
dsdda = np.array(f.get('Entropy_Table/dSdda'))
dsddz = np.array(f.get('Entropy_Table/dSddz'))
dsdee = np.array(f.get('Entropy_Table/dSdee'))
dsdea = np.array(f.get('Entropy_Table/dSdea'))
dsdez = np.array(f.get('Entropy_Table/dSdez'))
dsdaa = np.array(f.get('Entropy_Table/dSdaa'))
dsdaz = np.array(f.get('Entropy_Table/dSdaz'))
dsdzz = np.array(f.get('Entropy_Table/dSdzz'))

maxw1 = np.array(f.get('Thermo_Consist/maxw1'))
maxw2 = np.array(f.get('Thermo_Consist/maxw2'))

plt.figure(1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.contourf(den_mesh[:,:,1,4], norm_energy_mesh[:,:,1,4], 1/(dsde[:,:,1,4]),locator=ticker.LogLocator(), cmap='rainbow')
plt.colorbar(label='Entropy');
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'Normalized Energy')
plt.xlabel(r'Density (g/cc)')
plt.grid()
#plt.xlim((1e6, 1.4e6))

f.close() 


den_hydro = 1.0
E3_hydro = 2.0453130951119775E26
abar_hydro = 2.0
zbar_hydro = 0.5
Ye_hydro = zbar_hydro/abar_hydro
Emin_hydro = EOS_minimum_energy_calculator(den_hydro,abar_hydro,zbar_hydro)
norm_e_hydro = (E3_hydro - Emin_hydro)/Emin_hydro


## Find data

den_index = floor((log10(den_hydro) - log10(den_lo))*den_mesh.shape[1]/(log10(den_hi) - log10(den_lo)))
norm_e_index = floor((log10(norm_e_hydro) - log10(norm_energy_lo))*den_mesh.shape[0]/(log10(norm_energy_hi) - log10(norm_energy_lo)))
abar_index = floor((abar_hydro - abar_lo)*den_mesh.shape[2]/(abar_hi - abar_lo))
zbar_index = floor((zbar_hydro - zbar_lo)*den_mesh.shape[3]/(zbar_hi - zbar_lo))


ind_pt = [norm_e_index,den_index,abar_index,zbar_index]
print(ind_pt)

#######################  Construct Training Data  #############################

den_tr    = den_mesh[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()
norm_e_tr = norm_energy_mesh[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()
abar_tr   = abar_mesh[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()
zbar_tr   = zbar_mesh[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()
s_tr      = s[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()

dsde_tr   = dsde[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()
dsdd_tr   = dsdd[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()
dsda_tr   = dsda[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()
dsdz_tr   = dsda[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()

dsdee_tr  = dsdee[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()
dsdde_tr  = dsdde[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()
dsdea_tr  = dsdea[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()
dsdez_tr  = dsdez[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()

dsddd_tr  = dsddd[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()
dsdda_tr  = dsdda[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()
dsddz_tr  = dsddz[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()

dsdaa_tr  = dsdaa[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()
dsdaz_tr  = dsdaz[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()

dsdzz_tr  = dsdzz[norm_e_index:norm_e_index+2,den_index:den_index+2,abar_index:abar_index+2,zbar_index-2:zbar_index].ravel()

gph_sigma_f = 1.0
sigma_n = 1.0e-10
hyperparams_in=[1.0,1.0,1.0,1.0]

f_prior = np.concatenate((s_tr,dsde_tr,dsdd_tr,dsda_tr,dsdz_tr,dsdee_tr,dsdde_tr,\
                          dsdea_tr,dsdez_tr,dsddd_tr,dsdda_tr,dsddz_tr,dsdaa_tr,\
                          dsdaz_tr,dsdzz_tr))

[hyperparams_out_0,K0,fnew0,hyperparams_out_1,K1,fnew1,hyperparams_out_2,K2,fnew2]\
    =EOS_MLL_Solver(norm_e_tr,den_tr,abar_tr,zbar_tr,\
                             norm_e_hydro,den_hydro,abar_hydro,zbar_hydro,\
                             f_prior,gph_sigma_f,sigma_n,hyperparams_in)

################## Calculate Thermodynamic Quantitites #######################

s_post     = fnew2[0]

dsde_out   = fnew2[1]
dsdd_out   = fnew2[2]
dsda_out   = fnew2[3]
dsdz_out   = fnew2[4]

dsdee_out   = fnew2[5]
dsdde_out   = fnew2[6]
dsdea_out   = fnew2[7]
dsdea_out   = fnew2[8]

dsddd_out   = fnew2[9]
dsdda_out   = fnew2[10]
dsddz_out   = fnew2[11]

dsdaa_out   = fnew2[12]
dsdaz_out   = fnew2[13]

dsdzz_out   = fnew2[4]

#temp=1.0/dsde_out
#
#print('interpolated temp',temp)

#plt.figure(1)
##plt.rc('text', usetex=True)
##plt.rc('font', family='serif')
#plt.contourf(K0, cmap='rainbow')
#plt.axis('equal')
#plt.colorbar()
#plt.grid()

## for the gas
## the temperature and density exponents (c&g 9.81 9.82)
## the specific heat at constant volume (c&g 9.92)
## the third adiabatic exponent (c&g 9.93)
## the first adiabatic exponent (c&g 9.97)
## the second adiabatic exponent (c&g 9.105)
## the specific heat at constant pressure (c&g 9.98)
## and relativistic formula for the sound speed (c&g 14.29)
#
#zz        = pgas/den
#chit_gas  = temp/pgas * dpgasdt
#chid_gas  = dpgasdd/zz
#cv_gas    = degasdt
#x         = zz * chit_gas/(temp * cv_gas)
#gam3_gas  = x + 1.0
#gam1_gas  = chit_gas*x + chid_gas
#nabad_gas = x/gam1_gas
#gam2_gas  = 1.0/(1.0 - nabad_gas)
#cp_gas    = cv_gas * gam1_gas/chid_gas
#z         = 1.0 + (egas + clight*clight)/zz
#sound_gas = clight * sqrt(gam1_gas/z)
# 
# 
## for the totals
#zz    = pres/den
#chit  = temp/pres * dpresdt
#chid  = dpresdd/zz
#cv    = denerdt
#x     = zz * chit/(temp * cv)
#gam3  = x + 1.0d0
#gam1  = chit*x + chid
#nabad = x/gam1
#gam2  = 1.0d0/(1.0d0 - nabad)
#cp    = cv * gam1/chid
#z     = 1.0d0 + (ener + clight*clight)/zz
#sound = clight * sqrt(gam1/z)
# 
#
### Calculate Thermodynamic Consistency
# 
## maxwell relations; each is zero if the consistency is perfect
## delicate subtraction in very degenerate regions causes roundoff error
# 
#        dse = temp*dentrdt/denerdt - 1.0d0
# 
#        dpe = (denerdd*den**2 + temp*dpresdt)/pres - 1.0d0
# 
#        dsp = -dentrdd*(den**2/dpresdt) - 1.0d0

