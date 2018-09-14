#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 12:24:56 2018

@author: alexgrannan
"""

from math import *
import numpy as np
from numpy import ma
import h5py
import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm
from datetime import datetime

start_time = datetime.now()

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


###################################### Inputs ################################
##############################################################################
temp_0    = 1.0000E3
beta      = kerg * temp_0 / mecc
etaele_0  = -21.653832766028515  
abar      = 1.0
zbar      = 0.1
den_0     = 1.0000E-12
einput_0  = 144754322608.70731
Ye        = zbar/abar

# set the on/off switches
radmult  = 1
ionmult  = 1
ionized  = 1
elemult  = 1
coulmult = 1
potmult  = 0
mode     = 0
print('Inputs set up')
############################### Initialize Meshes ############################
##############################################################################

abar_lo = 1.0
abar_hi = 4.0
abar_range = np.linspace(abar_lo, abar_hi,5)
print('abar range set')
Ye_lo = 0.1
Ye_hi = 1.0
Ye_range = np.linspace(Ye_hi,Ye_lo,5)
print('Ye range set')
den_lo = log10(den_0)
den_hi = log10(1e10)
den_range=np.logspace(den_lo, den_hi, num=round(5*(den_hi-den_lo)))
print('density range set')


e_min = EOS_minimum_energy_calculator(den_range,abar_range,Ye_range*abar_range)
print('minimum energy calculated')
elo = log10((einput_0 - e_min[0,0,0])/e_min[0,0,0])
ehi = log10(1.0e26)
e_range_log = np.logspace(elo, ehi, num=round(5*(ehi-elo)))
print('energy range set')


eta_range  = np.zeros_like(e_range_log)
temp_range = np.zeros_like(den_range)

d_mesh, diff_e_mesh, abar_mesh, zbar_mesh = np.meshgrid(den_range, e_range_log, abar_range, Ye_range*abar_range)

## Build Arrays
t_mesh   = temp_0*np.ones_like(d_mesh)
eta_mesh = etaele_0*np.ones_like(d_mesh)

s       = np.zeros_like(d_mesh)
dsdd    = np.zeros_like(d_mesh)
dsde    = np.zeros_like(d_mesh)
dsda    = np.zeros_like(d_mesh)
dsdz    = np.zeros_like(d_mesh)
dsddd   = np.zeros_like(d_mesh)
dsdde   = np.zeros_like(d_mesh)
dsdda   = np.zeros_like(d_mesh)
dsddz   = np.zeros_like(d_mesh)
dsdee   = np.zeros_like(d_mesh)
dsdea   = np.zeros_like(d_mesh)
dsdez   = np.zeros_like(d_mesh)
dsdaa   = np.zeros_like(d_mesh)
dsdaz   = np.zeros_like(d_mesh)
dsdzz   = np.zeros_like(d_mesh)
maxw1   = np.zeros_like(d_mesh)
maxw2   = np.zeros_like(d_mesh)

eta_totable     = np.zeros_like(d_mesh)
#detadd_totable  = np.zeros_like(d_mesh)
#detaddd_totable = np.zeros_like(d_mesh)
#detade_totable  = np.zeros_like(d_mesh)
#detadee_totable = np.zeros_like(d_mesh)
#detaded_totable = np.zeros_like(d_mesh)
l1_mesh = np.zeros_like(d_mesh)
l2_mesh = np.zeros_like(d_mesh)
l3_mesh = np.zeros_like(d_mesh)
l4_mesh = np.zeros_like(d_mesh)

## Initialize Array
#t_mesh[0,0]      = temp_0
#temp_range[0]    = temp_0
#eta_totable[0,0] = etaele_0
#eta_range[0]     = etaele_0
eta_mesh[0,0,0,0] = etaele_0
print('arrays and meshes initialized')
print('looping begins')

for l in range(len(Ye_range)):

    for k in range(len(abar_range)):
        
        for i in range(len(e_range_log)):
            
            for j in range(len(den_range)):
                
    
    #print ('j=',j)

        #print ('$$')
                print ('i=',i,'of ',len(e_range_log)-1,' and j=',j,' of ',len(den_range)-1,'k=',k,'of ',len(abar_range)-1,'l=',l,'of ',len(Ye_range)-1)
                
                if ((i == 0) and (j == 0) and (k == 0) and (l == 0)):
                    #print('A')
                    continue
                
                elif ((i != 0) and (j == 0) and (k == 0) and (l == 0)):
                    
                    beta_in = kerg*t_mesh[i-1,j,k,l]/mecc
                    eta_in  = eta_mesh[i-1,j,k,l]
                elif ((i == 0) and (j != 0) and (k == 0) and (l == 0)):
                
                    beta_in = kerg*t_mesh[i,j-1,k,l]/mecc
                    eta_in  = eta_mesh[i,j-1,k,l]
                    
                elif ((i == 0) and (j == 0) and (k != 0) and (l == 0)):
                
                    beta_in = kerg*t_mesh[i,j,k-1,l]/mecc
                    eta_in  = eta_mesh[i,j,k-1,l]
                    
                elif ((i != 0) and (j != 0) and (k == 0) and (l == 0)):
                
                    beta_in = kerg*t_mesh[i-1,j-1,k,l]/mecc
                    eta_in  = eta_mesh[i-1,j-1,k,l]
                    
                elif ((i == 0) and (j != 0) and (k != 0) and (l == 0)):
                
                    beta_in = kerg*t_mesh[i,j-1,k-1,l]/mecc
                    eta_in  = eta_mesh[i,j-1,k-1,l]
                
                elif ((i != 0) and (j == 0) and (k != 0) and (l == 0)):
                
                    beta_in = kerg*t_mesh[i-1,j,k-1,l]/mecc
                    eta_in  = eta_mesh[i-1,j,k-1,l]
                
                elif ((i != 0) and (j != 0) and (k != 0) and (l == 0)):
                    beta_in = kerg*t_mesh[i-1,j-1,k-1,l]/mecc
                    eta_in  = eta_mesh[i-1,j-1,k-1,l]
                    #print('beta_in',beta_in,'eta_in',eta_in)
                    #print('B')
                elif ((i == 0) and (j == 0) and (k == 0) and (l != 0)):
                    
                    beta_in = kerg*t_mesh[i,j,k,l-1]/mecc
                    eta_in  = eta_mesh[i,j,k,l-1]
                elif ((i != 0) and (j == 0) and (k == 0) and (l != 0)):
                
                    beta_in = kerg*t_mesh[i-1,j,k,l-1]/mecc
                    eta_in  = eta_mesh[i-1,j,k,l-1]
                    
                elif ((i == 0) and (j != 0) and (k == 0) and (l != 0)):
                
                    beta_in = kerg*t_mesh[i,j-1,k,l-1]/mecc
                    eta_in  = eta_mesh[i,j-1,k,l-1]
                    
                elif ((i == 0) and (j == 0) and (k != 0) and (l != 0)):
                
                    beta_in = kerg*t_mesh[i,j,k-1,l-1]/mecc
                    eta_in  = eta_mesh[i,j,k-1,l-1]
                    
                elif ((i != 0) and (j != 0) and (k == 0) and (l != 0)):
                
                    beta_in = kerg*t_mesh[i-1,j-1,k,l-1]/mecc
                    eta_in  = eta_mesh[i-1,j-1,k,l-1]
                
                elif ((i == 0) and (j != 0) and (k != 0) and (l != 0)):
                
                    beta_in = kerg*t_mesh[i,j-1,k-1,l-1]/mecc
                    eta_in  = eta_mesh[i,j-1,k-1,l-1]
                
                elif ((i != 0) and (j == 0) and (k != 0) and (l != 0)):
                    beta_in = kerg*t_mesh[i-1,j,k-1,l-1]/mecc
                    eta_in  = eta_mesh[i-1,j,k-1,l-1]
                
                    #print('C')
                else:
                    beta_in = kerg*t_mesh[i-1,j-1,k-1,l-1]/mecc
                    eta_in  = eta_mesh[i-1,j-1,k-1,l-1]
                    #print('D')
        
                abar_in = abar_mesh[i,j,k,l]
                zbar_in = zbar_mesh[i,j,k,l]
                den_in  = d_mesh[i,j,k,l]
                e_in    = diff_e_mesh[i,j,k,l]*e_min[j,k,l] + e_min[j,k,l]
                    
                #print('den_in',den_in,'e_in',e_in,'temp_in',beta_in*mecc/kerg,'eta_in',eta_in)
                print("%2s %6.2e %2s %6.2e %2s %6.2e %2s %6.2e %2s %6.2e %2s %6.2e"% ('den_in =',den_in,'e_in = ',e_in,'temp_in=',beta_in*mecc/kerg,'eta_in',eta_in,'zbar=',zbar_in,'abar=',abar_in))
                
                
                
                [s_tmp,dsde_tmp,dsdd_tmp,dsda_tmp,dsdz_tmp,dsdee_tmp,dsded_tmp,dsdea_tmp,dsdez_tmp,dsddd_tmp,dsdda_tmp,dsddz_tmp,\
                    dsdaa_tmp,dsdaz_tmp,dsdzz_tmp,eta_out,beta_out,maxw1_tmp,maxw2_tmp]=\
                                Entropic_EOS(beta_in,den_in,eta_in,abar_in,zbar_in,e_in,\
                                radmult,ionmult,ionized,elemult,coulmult,potmult,mode)
                              
                 # Total Entropy
                temp_out = beta_out*mecc/kerg
                s[i,j,k,l]     = s_tmp
                # First deriviative of total entropy with density
                dsdd[i,j,k,l]  = dsdd_tmp
                dsde[i,j,k,l]  = dsde_tmp
                dsda[i,j,k,l]  = dsda_tmp
                dsdz[i,j,k,l]  = dsdz_tmp
                # Second deriviative of total entropy with density
                dsddd[i,j,k,l] = dsddd_tmp
                dsdde[i,j,k,l] = dsded_tmp
                dsdda[i,j,k,l] = dsdda_tmp
                dsddz[i,j,k,l] = dsddz_tmp
                dsdee[i,j,k,l] = dsdee_tmp
                dsdea[i,j,k,l] = dsdea_tmp
                dsdez[i,j,k,l] = dsdez_tmp
                dsdaa[i,j,k,l] = dsdaa_tmp
                dsdaz[i,j,k,l] = dsdaz_tmp
                dsdzz[i,j,k,l] = dsdzz_tmp
#                # First derivative of entropy with energy
#                ds3de[i,j]  = temp_out
#                # Second drivative of entropy with energy
#                ds3dee[i,j] = -1/((deiondt + deeledt + deposdt)*temp_out**2)
#                # Mixed derivative of entropy with temperature and energy
#                ds3ded[i,j] = dsded
#                print('dsded',dsded)
        
                # Chemical Potential
        #        eta_totable[i,j]  = etaele_out
        #        # First Derivative of Chemical Potential with density
        #        detadd_totable[i,j] = detadd
        #        # Second Derivative of Chemical Potential with density
        #        detaddd_totable[i,j] = detaddd
        #        # First Derivative of Chemical Potential with energy
        #        detade_totable[i,j] = 1/(deion_deta + deele_deta + depos_deta)
        #        # Second Derivative of Chemical Potential with energy
        #        detadee_totable[i,j] = 1/(deion_deta2 + deele_deta2 + depos_deta2)
        #        # Mixed Derivative of Chemical Potential
        #        detaded_totable[i,j] = 0.0
                
                t_mesh[i,j,k,l]= temp_out
#                temp_range[j] = temp_out
                eta_mesh[i,j,k,l]  = eta_out
#                eta_range[i]   = etaele_out
                maxw1[i,j,k,l] = maxw1_tmp
                maxw2[i,j,k,l] = maxw2_tmp
                print("%2s %6.2e %2s %6.2e %2s %6.2e %2s %6.2e %2s %6.2e %2s %6.2e"% ('den_in =',den_in,'e_in = ',e_in,'temp_out=',temp_out,'eta_out',eta_out,'zbar=',zbar_in,'abar=',abar_in))

                #print('den_in',den_in,'e_in',e_in,'temp_out',temp_out,'eta_out',eta_out)

## 
                
############################## Find hyperparameters ##########################
##############################################################################
                
                



#print ('i=',i,'j=',j)
#t_plot = t_mesh[0,:]
#plt.figure(1)   
#plt.plot(e_range_log,t_plot)
##plt.yscale('log')
##plt.xscale('log')
#plt.show
#
#plt.figure(2)   
#plt.plot(e_range_log,eta_mesh[0,j])
##plt.yscale('log')
#plt.xscale('log')
#plt.show
        
# The following is not strictly essential, but it will eliminate
# a warning.  Comment it out to see the warning.
# eta_mesh = ma.masked_where(eta_mesh <= 0, eta_mesh)

plt.figure(1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.contourf(d_mesh[:,:,0,0], diff_e_mesh[:,:,0,0], t_mesh[:,:,0,0], 200,locator=ticker.LogLocator(), cmap='rainbow')
plt.colorbar(label='Temp (K)');
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'Normalized Energy')
plt.xlabel(r'Density (g/cc)')
plt.grid()
#plt.xlim((1e6, 1.4e6))
#
#plt.figure(2)
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.contourf(d_mesh, diff_e_mesh,(eta_mesh), 200,cmap='rainbow')
#plt.colorbar(label=r'\eta');
#plt.xscale('log')
#plt.yscale('log')
#plt.ylabel('Normalized Energy')
#plt.xlabel('Density (g/cc)')
#plt.grid()
##plt.xlim((1e6, 1.4e6))
#
#T_degen = 1e5*(den_range**(2.0/3.0))
#
#plt.figure(3)
#plt.plot(den_range, temp_range)
#plt.plot(den_range, T_degen)
#plt.yscale('log')
#plt.xscale('log')
##
#plt.figure(4)
#plt.plot(e_range_log,eta_range)
##plt.yscale('log')
#plt.xscale('log')
##plt.xlim((1e6, 1.4e6))
##plt.xlim(1e6, 1.4e6)
#
#plt.figure(5)
#plt.plot(den_range,e_min)
#plt.yscale('log')
#plt.xscale('log')
##plt.xlim((1e6, 1.4e6))
##plt.xlim(1e6, 1.4e6)
#
#plt.figure(6)
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.contourf(d_mesh, diff_e_mesh, s3, 200,locator=ticker.LogLocator(), cmap='rainbow')
#plt.colorbar(label=r'S_{3} (erg/g/K)');
#plt.xscale('log')
#plt.yscale('log')
#plt.ylabel('Normalized Energy')
#plt.xlabel('Density (g/cc)')
#plt.grid()
##plt.savefig('destination_path.eps', format='eps', dpi=1000)
##plt.xlim((1e6, 1.4e6))
#
#
#
###################### Create Data File and save data to file ##################
################################################################################
hf = h5py.File('Entropic_EOS_Table.h5', 'w')

g1 = hf.create_group('Variable_Table')

g1.create_dataset('den',data=d_mesh)
g1.create_dataset('norm_energy',data=diff_e_mesh)
g1.create_dataset('abar',data=abar_mesh)
g1.create_dataset('zbar',data=zbar_mesh)


g2 = hf.create_group('Entropy_Table')

g2.create_dataset('S',data=s)
g2.create_dataset('dSdd',data=dsdd)
g2.create_dataset('dSde',data=dsde)
g2.create_dataset('dSda',data=dsda)
g2.create_dataset('dSdz',data=dsdz)
g2.create_dataset('dSddd',data=dsddd)
g2.create_dataset('dSdde',data=dsdde)
g2.create_dataset('dSdda',data=dsdda)
g2.create_dataset('dSddz',data=dsddz)
g2.create_dataset('dSdee',data=dsdee)
g2.create_dataset('dSdea',data=dsdea)
g2.create_dataset('dSdez',data=dsdez)
g2.create_dataset('dSdaa',data=dsdaa)
g2.create_dataset('dSdaz',data=dsdaz)
g2.create_dataset('dSdzz',data=dsdzz)

g3 = hf.create_group('Thermo_Consist')
g3.create_dataset('maxw1',data=maxw1)
g3.create_dataset('maxw2',data=maxw2)

#g3 = hf.create_group('hyperparams_0')
#
#g3.create_dataset('l1_0',data=l1_0)
#g3.create_dataset('l2_0',data=l2_0)
#g3.create_dataset('l3_0',data=l3_0)
#g3.create_dataset('l4_0',data=l4_0)
#
#g3 = hf.create_group('hyperparams_1')
#
#g3.create_dataset('l1_1',data=l1_1)
#g3.create_dataset('l2_1',data=l2_1)
#g3.create_dataset('l3_1',data=l3_1)
#g3.create_dataset('l4_1',data=l4_1)
#
#g3 = hf.create_group('hyperparams_2')
#
#g3.create_dataset('l1_2',data=l1_2)
#g3.create_dataset('l2_2',data=l2_2)
#g3.create_dataset('l3_2',data=l3_2)
#g3.create_dataset('l4_2',data=l4_2)

hf.close() 

time_elapsed = datetime.now() - start_time 

print('Time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))