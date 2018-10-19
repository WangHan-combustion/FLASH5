#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 16:39:34 2018

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

temp_var = 1.0e3
beta_var = kerg*temp_var/mecc
den_var = 1.0e-12*1.000
eta_var = -19.351247671779717
abar_var = 1.0
zbar_var = 1.0
einput_var = 249434204973.66833*1.00
# set the on/off switches
radmult  = 1
ionmult  = 1
ionized  = 1
elemult  = 1
coulmult = 1
potmult  = 0


def EntropicEOS_singlefile(beta_var,den_var,eta_var,abar_var,zbar_var,einput_var,\
                           radmult,ionmult,ionized,elemult,coulmult,potmult):
    
    
    ## Set all values to zero
    prad     = 0.0
    dpraddd  = 0.0
    dpraddt  = 0.0
    dpradda  = 0.0
    dpraddz  = 0.0
    dpradddd = 0.0
    dpradddt = 0.0
    dpraddda = 0.0
    dpradddz = 0.0
    dpraddtt = 0.0
    dpraddta = 0.0
    dpraddtz = 0.0
    dpraddaa = 0.0
    dpraddaz = 0.0
    dpraddzz = 0.0
    
    erad     = 0.0
    deraddd  = 0.0
    deraddt  = 0.0
    deradda  = 0.0
    deraddz  = 0.0
    deradddd = 0.0
    deradddt = 0.0
    deraddda = 0.0
    deradddz = 0.0
    deraddtt = 0.0
    deraddta = 0.0
    deraddtz = 0.0
    deraddaa = 0.0
    deraddaz = 0.0
    deraddzz = 0.0
     
    srad     = 0.0
    dsraddd  = 0.0
    dsraddt  = 0.0
    dsradda  = 0.0
    dsraddz  = 0.0
    dsradddd = 0.0
    dsradddt = 0.0
    dsraddda = 0.0
    dsradddz = 0.0
    dsraddtt = 0.0
    dsraddta = 0.0
    dsraddtz = 0.0
    dsraddaa = 0.0
    dsraddaz = 0.0
    dsraddzz = 0.0
     
    xni     = 0.0
    dxnidd  = 0.0
    dxnidt  = 0.0
    dxnida  = 0.0
    dxnidz  = 0.0
    dxniddd = 0.0
    dxniddt = 0.0
    dxnidda = 0.0
    dxniddz = 0.0
    dxnidtt = 0.0
    dxnidta = 0.0
    dxnidtz = 0.0
    dxnidaa = 0.0
    dxnidaz = 0.0
    dxnidzz = 0.0
    
    pion     = 0.0
    dpiondd  = 0.0
    dpiondt  = 0.0
    dpionda  = 0.0
    dpiondz  = 0.0
    dpionddd = 0.0
    dpionddt = 0.0
    dpiondda = 0.0
    dpionddz = 0.0
    dpiondtt = 0.0
    dpiondta = 0.0
    dpiondtz = 0.0
    dpiondaa = 0.0
    dpiondaz = 0.0
    dpiondzz = 0.0
    
    eion     = 0.0
    deiondd  = 0.0
    deiondt  = 0.0
    deionda  = 0.0
    deiondz  = 0.0
    deionddd = 0.0
    deionddt = 0.0
    deiondda = 0.0
    deionddz = 0.0
    deiondtt = 0.0
    deiondta = 0.0
    deiondtz = 0.0
    deiondaa = 0.0
    deiondaz = 0.0
    deiondzz = 0.0
    
    sion     = 0.0
    dsiondd  = 0.0
    dsiondt  = 0.0
    dsionda  = 0.0
    dsiondz  = 0.0
    dsionddd = 0.0
    dsionddt = 0.0
    dsiondda = 0.0
    dsionddz = 0.0
    dsiondtt = 0.0
    dsiondta = 0.0
    dsiondtz = 0.0
    dsiondaa = 0.0
    dsiondaz = 0.0
    dsiondzz = 0.0
    
    etaion   = 0.0
    detaidd  = 0.0
    detaidt  = 0.0
    detaida  = 0.0
    detaidz  = 0.0
    detaiddd = 0.0
    detaiddt = 0.0
    detaidda = 0.0
    detaiddz = 0.0
    detaidtt = 0.0
    detaidta = 0.0
    detaidtz = 0.0
    detaidaa = 0.0
    detaidaz = 0.0
    detaidzz = 0.0
    
    xne     = 0.0
    dxnedd  = 0.0
    dxnedt  = 0.0
    dxneda  = 0.0
    dxnedz  = 0.0
    dxneddd = 0.0
    dxneddt = 0.0
    dxnedda = 0.0
    dxneddz = 0.0
    dxnedtt = 0.0
    dxnedta = 0.0
    dxnedtz = 0.0
    dxnedaa = 0.0
    dxnedaz = 0.0
    dxnedzz = 0.0
     
    etaele  = 0.0
    detadd  = 0.0
    detadt  = 0.0
    detada  = 0.0
    detadz  = 0.0
    detaddd = 0.0
    detaddt = 0.0
    detadda = 0.0
    detaddz = 0.0
    detadtt = 0.0
    detadta = 0.0
    detadtz = 0.0
    detadaa = 0.0
    detadaz = 0.0
    detadzz = 0.0
     
    etapos   = 0.0
     
    xnefer     = 0.0
    dxneferdd  = 0.0
    dxneferdt  = 0.0
    dxneferda  = 0.0
    dxneferdz  = 0.0
    dxneferddd = 0.0
    dxneferddt = 0.0
    dxneferdda = 0.0
    dxneferddz = 0.0
    dxneferdtt = 0.0
    dxneferdta = 0.0
    dxneferdtz = 0.0
    dxneferdaa = 0.0
    dxneferdaz = 0.0
    dxneferdzz = 0.0
    
    xnpfer     = 0.0
    dxnpferdd  = 0.0
    dxnpferdt  = 0.0
    dxnpferda  = 0.0
    dxnpferdz  = 0.0
    dxnpferddd = 0.0
    dxnpferddt = 0.0
    dxnpferdda = 0.0
    dxnpferddz = 0.0
    dxnpferdtt = 0.0
    dxnpferdta = 0.0
    dxnpferdtz = 0.0
    dxnpferdaa = 0.0
    dxnpferdaz = 0.0
    dxnpferdzz = 0.0
    
    pele     = 0.0
    dpeledd  = 0.0
    dpeledt  = 0.0
    dpeleda  = 0.0
    dpeledz  = 0.0
    dpeleddd = 0.0
    dpeleddt = 0.0
    dpeledda = 0.0
    dpeleddz = 0.0
    dpeledtt = 0.0
    dpeledta = 0.0
    dpeledtz = 0.0
    dpeledaa = 0.0
    dpeledaz = 0.0
    dpeledzz = 0.0
     
    eele     = 0.0
    deeledd  = 0.0
    deeledt  = 0.0
    deeleda  = 0.0
    deeledz  = 0.0
    deeleddd = 0.0
    deeleddt = 0.0
    deeledda = 0.0
    deeleddz = 0.0
    deeledtt = 0.0
    deeledta = 0.0
    deeledtz = 0.0
    deeledaa = 0.0
    deeledaz = 0.0
    deeledzz = 0.0
     
    sele     = 0.0
    dseledd  = 0.0
    dseledt  = 0.0
    dseleda  = 0.0
    dseledz  = 0.0
    dseleddd = 0.0
    dseleddt = 0.0
    dseledda = 0.0
    dseleddz = 0.0
    dseledtt = 0.0
    dseledta = 0.0
    dseledtz = 0.0
    dseledaa = 0.0
    dseledaz = 0.0
    dseledzz = 0.0
    
    ppos     = 0.0
    dpposdd  = 0.0
    dpposdt  = 0.0
    dpposda  = 0.0
    dpeledz  = 0.0
    dpposddd = 0.0
    dpposddt = 0.0
    dpposdda = 0.0
    dpposddz = 0.0
    dpposdtt = 0.0
    dpposdta = 0.0
    dpposdtz = 0.0
    dpposdaa = 0.0
    dpposdaz = 0.0
    dpposdzz = 0.0
    
    epos     = 0.0
    deposdd  = 0.0
    deposdt  = 0.0
    deposda  = 0.0
    deeledz  = 0.0
    deposddd = 0.0
    deposddt = 0.0
    deposdda = 0.0
    deposddz = 0.0
    deposdtt = 0.0
    deposdta = 0.0
    deposdtz = 0.0
    deposdaa = 0.0
    deposdaz = 0.0
    deposdzz = 0.0
     
    spos     = 0.0
    dsposdd  = 0.0
    dsposdt  = 0.0
    dsposda  = 0.0
    dseledz  = 0.0
    dseleddd = 0.0
    dseleddt = 0.0
    dseledda = 0.0
    dseleddz = 0.0
    dseledtt = 0.0
    dseledta = 0.0
    dseledtz = 0.0
    dseledaa = 0.0
    dseledaz = 0.0
    dseledzz = 0.0
     
    pep     = 0.0
    dpepdd  = 0.0
    dpepdt  = 0.0
    dpepda  = 0.0
    dpepdz  = 0.0
    dpepddd = 0.0
    dpepddt = 0.0
    dpepdda = 0.0
    dpepddz = 0.0
    dpepdtt = 0.0
    dpepdta = 0.0
    dpepdtz = 0.0
    dpepdaa = 0.0
    dpepdaz = 0.0
    dpepdzz = 0.0
     
    eep     = 0.0
    deepdd  = 0.0
    deepdt  = 0.0
    deepda  = 0.0
    deepdz  = 0.0
    deepddd = 0.0
    deepddt = 0.0
    deepdda = 0.0
    deepddz = 0.0
    deepdtt = 0.0
    deepdta = 0.0
    deepdtz = 0.0
    deepdaa = 0.0
    deepdaz = 0.0
    deepdzz = 0.0
    
    sep     = 0.0
    dsepdd  = 0.0
    dsepdt  = 0.0
    dsepda  = 0.0
    dsepdz  = 0.0
    dsepddd = 0.0
    dsepddt = 0.0
    dsepdda = 0.0
    dsepddz = 0.0
    dsepdtt = 0.0
    dsepdta = 0.0
    dsepdtz = 0.0
    dsepdaa = 0.0
    dsepdaz = 0.0
    dsepdzz = 0.0
    
    eip      = 0.0
    deipdd   = 0.0
    deipdt   = 0.0
    deipda   = 0.0
    deipdz   = 0.0
    
    sip      = 0.0
    dsipdd   = 0.0
    dsipdt   = 0.0
    dsipda   = 0.0
    dsipdz   = 0.0
    
    pcoul     = 0.0
    dpcouldd  = 0.0
    dpcouldt  = 0.0
    dpcoulda  = 0.0
    dpcouldz  = 0.0
    dpcoulddd = 0.0
    dpcoulddt = 0.0
    dpcouldda = 0.0
    dpcoulddz = 0.0
    dpcouldtt = 0.0
    dpcouldta = 0.0
    dpcouldtz = 0.0
    dpcouldaa = 0.0
    dpcouldaz = 0.0
    dpcouldzz = 0.0
    
    ecoul     = 0.0
    decouldd  = 0.0
    decouldt  = 0.0
    decoulda  = 0.0
    decouldz  = 0.0
    decoulddd = 0.0
    decoulddt = 0.0
    decouldda = 0.0
    decoulddz = 0.0
    decouldtt = 0.0
    decouldta = 0.0
    decouldtz = 0.0
    decouldaa = 0.0
    decouldaz = 0.0
    decouldzz = 0.0
    
    decoul_deta=0.0
    decoul_dbeta=0.0
     
    scoul     = 0.0
    dscouldd  = 0.0
    dscouldt  = 0.0
    dscoulda  = 0.0
    dscouldz  = 0.0
    dscoulddd = 0.0
    dscoulddt = 0.0
    dscouldda = 0.0
    dscoulddz = 0.0
    dscouldtt = 0.0
    dscouldta = 0.0
    dscouldtz = 0.0
    dscouldaa = 0.0
    dscouldaz = 0.0
    dscouldzz = 0.0
     
    chit     = 0.0
    dchitdd  = 0.0
    dchitdt  = 0.0
    dchitda  = 0.0
    dchitdz  = 0.0
    
    chid     = 0.0
    dchiddd  = 0.0
    dchiddt  = 0.0
    dchidda  = 0.0
    dchiddz  = 0.0
     
    cv     = 0.0
    dcvdd  = 0.0
    dcvdt  = 0.0
    dcvda  = 0.0
    dcvdz  = 0.0
     
    cp     = 0.0
    dcpdd  = 0.0
    dcpdt  = 0.0
    dcpda  = 0.0
    dcpdz  = 0.0
     
    gam1     = 0.0
    dgam1dd  = 0.0
    dgam1dt  = 0.0
    dgam1da  = 0.0
    dgam1dz  = 0.0
     
    gam2     = 0.0
    dgam2dd  = 0.0
    dgam2dt  = 0.0
    dgam2da  = 0.0
    dgam2dz  = 0.0
    
    gam3     = 0.0
    dgam3dd  = 0.0
    dgam3dt  = 0.0
    dgam3da  = 0.0
    dgam3dz  = 0.0
     
    nabad   = 0.0
    dnabdd  = 0.0
    dnabdt  = 0.0
    dnabda  = 0.0
    dnabdz  = 0.0
     
    sound  = 0.0
    dcsdd  = 0.0
    dcsdt  = 0.0
    dcsda  = 0.0
    dcsdz  = 0.0
     
    chit_gas     = 0.0
    dchit_gasdd  = 0.0
    dchit_gasdt  = 0.0
    dchit_gasda  = 0.0
    dchit_gasdz  = 0.0
    
    chid_gas     = 0.0
    dchid_gasdd  = 0.0
    dchid_gasdt  = 0.0
    dchid_gasda  = 0.0
    dchid_gasdz  = 0.0
     
    cv_gas     = 0.0
    dcv_gasdd  = 0.0
    dcv_gasdt  = 0.0
    dcv_gasda  = 0.0
    dcv_gasdz  = 0.0
    
    cp_gas     = 0.0
    dcp_gasdd  = 0.0
    dcp_gasdt  = 0.0
    dcp_gasda  = 0.0
    dcp_gasdz  = 0.0
     
    gam1_gas     = 0.0
    dgam1_gasdd  = 0.0
    dgam1_gasdt  = 0.0
    dgam1_gasda  = 0.0
    dgam1_gasdz  = 0.0
     
    gam2_gas     = 0.0
    dgam2_gasdd  = 0.0
    dgam2_gasdt  = 0.0
    dgam2_gasda  = 0.0
    dgam2_gasdz  = 0.0
     
    gam3_gas     = 0.0
    dgam3_gasdd  = 0.0
    dgam3_gasdt  = 0.0
    dgam3_gasda  = 0.0
    dgam3_gasdz  = 0.0
     
    nabad_gas   = 0.0
    dnab_gasdd  = 0.0
    dnab_gasdt  = 0.0
    dnab_gasda  = 0.0
    dnab_gasdz  = 0.0
     
    sound_gas  = 0.0
    dcs_gasdd  = 0.0
    dcs_gasdt  = 0.0
    dcs_gasda  = 0.0
    dcs_gasdz  = 0.0
     
    dse = 0.0
    dpe = 0.0
    dsp = 0.0
    
     
    
    
    ## Other Variables
    temp_var = beta_var*mecc/kerg
    
    
    
    #################### Radiation Terms #######################################
    ############################################################################
    
    if radmult == 1.0:
        
        ytot1   = 1.0/abar_var
        deni    = 1.0/den_var
        kt      = beta_var*mecc
        temp_var       = kt/kerg
        tempi       = 1/temp_var
        kti     = 1.0/kt
        beta_var32 = beta_var**1.5
        beta_var52 = beta_var*beta_var32
        dbetadt    = kerg/mecc
        
        # pressure in erg/cm**3
        prad    = asol * temp_var * temp_var * temp_var * temp_var/3.0
        dpraddd = 0.0
        dpraddt = 4.0 * prad/temp_var
        dpradda = 0.0
        dpraddz = 0.0
        dprad_dbeta = 4.0*asol*(beta_var**3.0)*(mecc)**4.0/kerg**4.0
        
        dpradddd = 0.0
        dpradddt = 0.0
        dpraddda = 0.0
        dpradddz = 0.0
        dpraddtt = 4.0*asol*temp_var
        dpraddta = 0.0
        dpraddtz = 0.0
        dpraddaa = 0.0
        dpraddaz = 0.0
        dpraddzz = 0.0
        
        # energy in erg/gr
        erad    = 3.0 * prad * deni
        
        deraddd = -erad * deni
        deraddt = 3.0 * dpraddt * deni
        deradda = 0.0
        deraddz = 0.0
        derad_deta = 0.0
        derad_dbeta = 3.0*dprad_dbeta/den_var
        
        deradddd = (3.0*dpradddd - 2.0*deraddd)*deni
        deradddt = (3.0*dpradddt - deraddt)*deni
        deraddda = 0.0
        deradddz = 0.0
        deraddtt = 3.0*dpraddtt*deni
        deraddta = 0.0
        deraddtz = 0.0
        deraddaa = 0.0
        deraddaz = 0.0
        deraddzz = 0.0
     
        # entropy in erg/g/kelvin
        srad    = (prad*deni + erad) * tempi
        dsraddd = (dpraddd*deni - prad*deni**2 + deraddd) * tempi
        dsraddt = (dpraddt*deni + deraddt - srad)  * tempi
        dsradda = 0.0
        dsraddz = 0.0
        
        dsraddtt = -dsraddt*tempi + (dpraddtt*deni + deraddtt - dsraddt)  * tempi
        
        dsradddd = (2.0*prad*deni**3.0 - 2.0*dpraddd*deni**2.0 + dpradddd*deni + deradddd)*tempi
        dsradddt = (-dpraddt*deni**2.0 + dpradddt*deni + deradddt - dsraddd)*tempi
        dsraddda = 0.0
        dsradddz = 0.0
        dsraddtt = (dpraddtt*deni + deraddtt - 2.0*dsraddt)*tempi
        dsraddta = 0.0
        dsraddtz = 0.0
        dsraddaa = 0.0
        dsraddaz = 0.0
        dsraddzz = 0.0
    
    ########################### Ion Terms ######################################
    ############################################################################
        
    ytot1   = 1.0/abar_var
    temp_var = beta_var * mecc / kerg
    tempinv  = 1/temp_var
    kt       = kerg * temp_var
    sifac    = (avo*(h**2.0)/2.0/pi)**1.5
    deninv   = 1.0/den_var
    
    # number density in 1/cm**3,
    xni     = avo * den_var / abar_var
    #print(xni,avo,den_var,abar_var)
    dxnidd  = avo / abar_var
    dxnidt  = 0.0
    dxnida  = -xni / abar_var
    dxnidz  = 0.0
    
    dxniddd = 0.0
    dxniddt = 0.0
    dxnidda = -avo/abar_var/abar_var
    dxniddz = 0.0
    dxnidtt = 0.0
    dxnidta = 0.0
    dxnidtz = 0.0
    dxnidaa = -2.0*dxnida/abar_var
    dxnidaz = 0.0
    dxnidzz = 0.0
    
    # set returned elements to zero
    pion     = 0.0
    
    dpiondd  = 0.0
    dpiondt  = 0.0
    dpionda  = 0.0
    dpiondz  = 0.0
    
    dpionddd = 0.0
    dpionddt = 0.0
    dpiondda = 0.0
    dpionddz = 0.0
    dpiondtt = 0.0
    dpiondta = 0.0
    dpiondtz = 0.0
    dpiondaa = 0.0
    dpiondaz = 0.0
    dpiondzz = 0.0
    
    eion        = 0.0
    deion_deta  = 0.0
    deion_deta2  = 0.0
    deion_dbeta = 0.0 
    sion        = 0.0 
    dsiondd     = 0.0 
    dsiondt     = 0.0 
    dsionddd    = 0.0
    
    if (ionmult != 0.0):
        # pressure in erg/cm**3
        pion    = xni * kt
        #print('kt',kt)
        dpiondd = dxnidd * kt
        dpiondt = xni * kerg
        dpionda = -pion / abar_var
        dpiondz = 0.0
        
        dpionddd = dxniddd*kt
        dpionddt = dxniddt*kt + dxnidd*kerg
        dpiondda = dxnidda*kt
        dpionddz = 0.0
        dpiondtt = dxnidtt*kt + 2.0*dxnidt*kerg
        dpiondta = dxnidta*kt + 2.0*dxnida*kerg
        dpiondtz = 0.0
        dpiondaa = dxnidaa*kt
        dpiondaz = 0.0
        dpiondzz = 0.0
    
       # energy in erg/gr
        eion    = 1.5 * pion / den_var
        deiondd = (1.5 * dpiondd - eion)/den_var
        deiondt = 1.5 * dpiondt / den_var
        deionda = 1.5 * dpionda / den_var
        deiondz = 0.0
        deion_deta  = 0.0
        deion_deta2  = 0.0
        deion_dbeta = 1.5 * xni * mecc / den_var
        
        deionddd = (1.5*dpionddd - 2.0*deiondd)/den_var
        deionddt = -deiondt/den_var + 1.5*dpionddt/den_var
        deiondda = -deionda/den_var + 1.5*dpiondda/den_var
        deionddz = 0.0
        deiondtt = 1.5*dpiondtt/den_var
        deiondta = 1.5*dpiondta/den_var
        deiondtz = 0.0
        deiondaa = 1.5*dpiondaa/den_var
        deiondaz = 0.0
        deiondzz = 0.0
        
       # ! ion degeneracy parameter (c&g 9.60)
        y       = 1.0/(abar_var*kt)
        #print(abar_var,kt)
        y32     = y * sqrt(y)
        
        
        dy32dd = 0.0
        dy32dt = -1.5*y*y32*abar_var*kerg
        dy32da = -1.5*y*y32*kt
        dy32dz = 0.0
        
        dy32ddd = 0.0
        dy32ddt = 0.0
        dy32dda = 0.0
        dy32ddz = 0.0
        dy32dtt = 1.5*2.5*y*y*y32*(abar_var*kerg)**2.0
        dy32dta = 1.5*2.5*y*y*y32*abar_var*temp_var*kerg**2.0 - 1.5*y*y32*kerg
        dy32dtz = 0.0
        dy32daa = 1.5*2.5*y*y*y32*kt**2.0
        dy32daz = 0.0
        dy32dzz = 0.0
        
        z       = xni * sifac * y32
        #print(xni,sifac,y32)
        dzdd = sifac*(dxnidd*y32 + xni*dy32dd)
        dzdt = sifac*(dxnidt*y32 + xni*dy32dt)
        dzda = sifac*(dxnida*y32 + xni*dy32da)
        dzdz = sifac*(dxnidz*y32 + xni*dy32dz)
        
        dzddd = sifac*(dxniddd*y32 + 2.0*dxnidd*dy32dd + xni*dy32ddd)
        dzddt = sifac*(dxniddt*y32 + dxnidd*dy32dt + dxnidt*dy32dd + xni*dy32ddt)
        dzdda = sifac*(dxnidda*y32 + dxnidd*dy32da + dxnida*dy32dd + xni*dy32dda)
        dzddz = sifac*(dxniddz*y32 + dxnidd*dy32dz + dxnidz*dy32dd + xni*dy32ddz)
        dzdtt = sifac*(dxnidtt*y32 + 2.0*dxnidt*dy32dt + xni*dy32dtt)
        dzdta = sifac*(dxnidta*y32 + dxnidt*dy32da + dxnida*dy32dt + xni*dy32dta)
        dzdtz = sifac*(dxnidtz*y32 + dxnidt*dy32dz + dxnidz*dy32dt + xni*dy32dtz)
        dzdaa = sifac*(dxnidaa*y32 + 2.0*dxnida*dy32da + xni*dy32daa)
        dzdaz = sifac*(dxnidaz*y32 + dxnidz*dy32da + dxnida*dy32dz + xni*dy32daz)
        dzdzz = sifac*(dxnidzz*y32 + 2.0*dxnidz*dy32dz + xni*dy32dzz)
        
        etaion  = log(z)
    
        zinv     = 1.0/z
        detaidd = dzdd*zinv
        detaidt = dzdt*zinv
        detaida = dzda*zinv
        detaidz = dzdz*zinv
        
        detaiddd = -detaidd**2.0 + zinv*dzddd
        detaiddt = -dzdd*dzdt*zinv**2.0 + dzddt*zinv
        detaidda = -dzdd*dzda*zinv**2.0 + dzdda*zinv
        detaiddz = -dzdd*dzdz*zinv**2.0 + dzddz*zinv
        detaidtt = -detaidt**2.0 + zinv*dzdtt
        detaidta = -dzdt*dzda*zinv**2.0 + dzdta*zinv
        detaidtz = -dzdt*dzdz*zinv**2.0 + dzdtz*zinv
        detaidaa = -detaida**2.0 + zinv*dzdaa
        detaidaz = -dzda*dzdz*zinv**2.0 + dzdaz*zinv
        detaidzz = -detaidz**2.0 + zinv*dzdzz
     
       # entropy in erg/gr/kelvin
       # the last term is the usual  etaion * kerg * xni/den
       # sometimes called the sacker-tetrode equation
    
        sion    = (eion + pion*deninv)*tempinv - etaion * kerg*avo*ytot1
     
        dsiondd = (deiondd + dpiondd*deninv - pion*deninv**2)*tempinv \
                  - detaidd * kerg * avo*ytot1
     
        dsiondt = (deiondt + dpiondt*deninv)*tempinv \
                   - (eion + pion*deninv)*tempinv**2 \
                   - detaidt * kerg * avo*ytot1
     
        dsionda = (deionda + dpionda*deninv)*tempinv \
                   - detaida * kerg * avo*ytot1 \
                   + etaion * kerg * avo * ytot1**2
    
        dsiondz = 0.0
        
        ## Second derivatives
        
    #    dsionddd = (deionddd + 2.0*pion*deninv**3 - 2.0*dpiondd*deninv**2  \
    #               + dpionddd*deninv)*tempinv  - detaiddd * kerg * avo*ytot1
        ## Note: I set all z derivatives to zero because of the lack of dependence       
        dsionddd = (2*pion/den_var**3.0 - 2.0*dpiondd/den_var**2.0 + dpionddd/den_var + deionddd)/temp_var \
                    - kerg*(xni*etaion/den_var - etaion*dxnidd - xni*detaidd)/den_var**2.0 \
                    + kerg*(-xni*etaion/den_var**2.0 + dxnidd*etaion/den_var + xni*detaidd/den_var \
                           -(detaidd*dxnidd + etaion*dxnidd) -(dxnidd*detaidd + xni*detaiddd))/den_var
        dsionddt = -(-pion/den_var**2.0 + dpiondd/den_var + deiondd)/temp_var**2.0 \
                    + (-dpiondt/den_var**2.0 + dpionddt/den_var + deionddt)/temp_var \
                    + kerg*((dxnidt*etaion + xni*detaidt)/den_var -(dxnidd*detaidt + etaion*dxniddt)
                    - (dxnidt*detaidd + xni*detaiddt))/den_var
        dsiondda = (-dpionda/den_var**2.0 + dpiondda/den_var + deiondda)/temp_var \
                    +kerg*((dxnida*etaion + xni*detaida)/den_var - (detaida*dxnidd + etaion*dxnidda)
                    - (dxnida*detaidd + xni*detaidda))/den_var
        dsionddz = 0.0
        dsiondtt = 2.0*(pion/den_var + eion)/temp_var**3.0 - 2.0*(dpiondt/den_var + deiondt)/temp_var**2.0 \
                    + (dpiondtt/den_var + deiondtt)/temp_var \
                    - kerg*(dxnidtt*etaion + dxnidt*detaidt + dxnidt*detaidt + xni*detaidtt)/den_var
        dsiondta = -(dpionda/den_var + deionda)/temp_var**2.0 + (dpiondta/den_var + deiondta)/temp_var \
                   - kerg*(dxnidta*etaion + dxnidt*detaida + dxnida*detaidt + xni*detaidta)/den_var
        dsiondtz = 0.0
        dsiondaa = (dpiondaa/den_var + deiondaa)/temp_var  \
                    - kerg*(dxnidaa*etaion + 2.0*dxnida*detaida + xni*detaidaa)/den_var
        dsiondaz = 0.0
        dsiondzz = 0.0
    
    a1 = -0.898004
    b1 =  0.96786
    c1 =  0.220703
    d1 = -0.86097
    e1 =  2.5269
    a2 =  0.29561
    b2 =  1.9885
    c2 =  0.288675
    
    ########################### Coulomb Corrections ############################
    ############################################################################
        
    # coulomb corrections section:
    if (coulmult != 0):
#        print('coulsection')
        
        # this code fragment implments coulomb corrections
        # see yakovlev & shalybkov 1989, uniform background corrections
     
        # input:
        # den  = density g/cc
        # temp = temperature k
        # abar = avarage atomic weight
        # zbar = avarge charge
        # pion and all its derivatives = ion pressure through common block
        # xne and all its derivatives = electron number density through common block
        
        # output:
        # plasg = ion coupling parameter
        # pcoul and all its derivatives = coulomb pressure through common block
        # ecoul and all its derivatives = coulomb energy through common block
        # scoul and all its derivatives = coulomb entropy through common block
         
        # yakovlev & shalybkov eqs 5, 9 and 10
        # use the ion number density instead of the free electron number desnity
        # to avoid issues with positrons being included or not (eosfxt vs helmeos)
        temp_var = beta_var * mecc / kerg
        tempinv = 1/temp_var
        kt       = kerg * temp_var
        ktinv    = 1.0/kt
        fiveth  = 5.0/3.0
        teninth = 10.0/9.0
        
        # number density in 1/cm**3,
        xni     = avo * den_var / abar_var
        dxnidd  = avo / abar_var
        dxnidt  = 0.0
        dxnida  = -xni / abar_var
        dxnidz  = 0.0
        dxniddd = 0.0
        dxniddt = 0.0
        dxnidda = -dxnidd*ytot1
        dxniddz = 0.0
        dxnidtt = 0.0
        dxnidta = 0.0
        dxnidtz = 0.0
        dxnidaa = -2.0 * dxnida * ytot1
        dxnidaz = 0.0
        dxnidzz = 0.0
        
        forth    = 4.0/3.0
        forthpi  = forth*pi
        third    = 1.0/3.0
     
        y        = forthpi * zbar_var
        s        = y * xni
        sinv     = 1.0/s
     
        # first derivatives
        dsdd     = y * dxnidd
        dsdt     = y * dxnidt
        dsda     = y * dxnida
        dsdz     = y * dxnidz + forthpi*xni
     
        # second derivatives
        dsddd   = y * dxniddd
        dsddt   = y * dxniddt
        dsdda   = y * dxnidda
        dsddz   = y * dxniddz + forthpi*dxnidd
        dsdtt   = y * dxnidtt
        dsdta   = y * dxnidta
        dsdtz   = y * dxnidtz + forthpi*dxnidt
        dsdaa   = y * dxnidaa
        dsdaz   = y * dxnidaz + forthpi*dxnida
        dsdzz   = y * dxnidzz + 2.0*forthpi*dxnidz
     
     
        #electron-sphere radius aele
        aele     = sinv**third
        aeleinv  = 1.0/aele
        z        = -third * aele * sinv
        y        = -forth * z * sinv
         
        # first derivatives
        daeledd  = z * dsdd
        daeledt  = z * dsdt
        daeleda  = z * dsda
        daeledz  = z * dsdz
    
        # second derivatives
        daeleddd = y*dsdd*dsdd + z*dsddd
        daeleddt = y*dsdt*dsdd + z*dsddt
        daeledda = y*dsda*dsdd + z*dsdda
        daeleddz = y*dsdz*dsdd + z*dsddz
        daeledtt = y*dsdt*dsdt + z*dsdtt
        daeledta = y*dsda*dsdt + z*dsdta
        daeledtz = y*dsdz*dsdt + z*dsdtz
        daeledaa = y*dsda*dsda + z*dsdaa
        daeledaz = y*dsdz*dsda + z*dsdaz
        daeledzz = y*dsdz*dsdz + z*dsdzz
     
     
        # electron coupling parameter eplasg
        eplasg   = esqu * ktinv * aeleinv
        z        = -eplasg * aeleinv
        y        = -2.0 * z * aeleinv
     
        # first derivatives
        deplasgdd = z * daeledd
        deplasgdt = z * daeledt - eplasg*tempinv
        deplasgda = z * daeleda
        deplasgdz = z * daeledz
     
        # second derivatives
        deplasgddd = y*daeledd*daeledd + z*daeleddd
        deplasgddt = y*daeledt*daeledd - deplasgdd*tempinv + z*daeleddt
        deplasgdda = y*daeleda*daeledd + z*daeledda
        deplasgddz = y*daeledz*daeledd + z*daeleddz
        deplasgdtt = y*daeledt*daeledt + z*daeledtt \
                     + (2.0*z*daeledt + 2.0*eplasg*tempinv)*tempinv
        deplasgdta = y*daeleda*daeledt + z*daeledta \
                     - deplasgda*tempinv
        deplasgdtz = y*daeledz*daeledt + z*daeledtz \
                     - deplasgdz*tempinv
        deplasgdaa = y*daeleda*daeleda + z*daeledaa
        deplasgdaz = y*daeledz*daeleda + z*daeledaz
        deplasgdzz = y*daeledz*daeledz + z*daeledzz
    
        # ion-sphere radius aion
        x          = zbar_var**third
        z          = x*x*x*x*x
        ww         = fiveth * x * x
        ion_radius = x * aele
        #print(x,z,ww,ion_radius,aele)
     
        # ion coupling parameter plasg
        plasg    = z * eplasg
        plasginv = 1.0/plasg
     
     
        # first derivatives
        dplasgdd  = z * deplasgdd
        dplasgdt  = z * deplasgdt
        dplasgda  = z * deplasgda
        dplasgdz  = z * deplasgdz + ww*eplasg
     
        # second derivatives
        dplasgddd  = z * deplasgddd
        dplasgddt  = z * deplasgddt
        dplasgdda  = z * deplasgdda
        dplasgddz  = z * deplasgddz + ww*deplasgdd
        dplasgdtt  = z * deplasgdtt
        dplasgdta  = z * deplasgdta
        dplasgdtz  = z * deplasgdtz + ww*deplasgdt
        dplasgdaa  = z * deplasgdaa
        dplasgdaz  = z * deplasgdaz + ww*deplasgda
        dplasgdzz  = z * deplasgdzz + 2.0*ww*deplasgdz +teninth/x*eplasg
     
     
        # yakovlev & shalybkov 1989 equations 82, 85, 86, 87
        if (plasg >= 1.0):
            
            x        = sqrt(sqrt(plasg))
            p1       = x
            p2       = 1.0/x
            p3       = p1*plasginv
            p4       = p2*plasginv
            p5       = p3*plasginv
            p6       = p4*plasginv
     
            u0       = a1*plasg + b1*p1 + c1*p2 + d1
            du0      = a1 + 0.25*b1*p3 - 0.25*c1*p4
            ddu0     = -0.1875*b1*p5 + 0.3125*c1*p6
     
     
            # energy in erg/gr
            z        = pion * deninv
            ecoul    = z * u0
     
            x        = deninv*u0
            y        = deninv*du0
            ww       = z*du0
            dfk       = z*ddu0
     
            # first derivatives
            decouldd = dpiondd*x - z*x + ww*dplasgdd
            decouldt = dpiondt*x + ww*dplasgdt
            decoulda = dpionda*x + ww*dplasgda
            decouldz = dpiondz*x + ww*dplasgdz
            decoul_deta = 0.0
            decoul_dbeta = decouldt*mecc/kerg
    
             # second derivatives
            decoulddd = dpionddd*x + 2.0*dpiondd*(y*dplasgdd - deninv*x) \
                         + z*(2.0*deninv*x \
                         + (ddu0*dplasgdd - 2.0*y)*dplasgdd + du0*dplasgddd)
            decoulddt = dpionddt*x + dpiondd*y*dplasgdt \
                         - dpiondt*deninv*x  - z*y*dplasgdt \
                         + (dpiondt*y + dfk*dplasgdt)*dplasgdd + ww*dplasgddt
            decouldda = dpiondda*x + dpiondd*y*dplasgda \
                         - dpionda*deninv*x - z*y*dplasgda \
                         + (dpionda*y + dfk*dplasgda)*dplasgdd + ww*dplasgdda
            decoulddz = dpionddz*x + dpiondd*y*dplasgdz \
                         - dpiondz*deninv*x - z*y*dplasgdz \
                         + (dpiondz*y + dfk*dplasgdz)*dplasgdd + ww*dplasgddz
            decouldtt = dpiondtt*x + (2.0*dpiondt*y \
                         + dfk*dplasgdt)*dplasgdt + ww*dplasgdtt
            decouldta = dpiondta*x + dpiondt*y*dplasgda \
                         + (dpionda*y + dfk*dplasgda)*dplasgdt + ww*dplasgdta
            decouldtz = dpiondtz*x + dpiondt*y*dplasgdz \
                         + (dpiondz*y + dfk*dplasgdz)*dplasgdt + ww*dplasgdtz
            decouldaa = dpiondaa*x + (2.0*dpionda*y \
                         + dfk*dplasgda)*dplasgda + ww*dplasgdaa
            decouldaz = dpiondaz*x + dpionda*y*dplasgdz \
                         + (dpiondz*y + dfk*dplasgdz)*dplasgda + ww*dplasgdaz
            decouldzz = dpiondzz*x + (2.0*dpiondz*y \
                         + dfk*dplasgdz)*dplasgdz + ww*dplasgdzz
    
            # pressure in erg/cc
            y        = third * den_var
            pcoul    = y * ecoul
     
            # first derivatives
            dpcouldd = third*ecoul + y*decouldd
            dpcouldt = y * decouldt
            dpcoulda = y * decoulda
            dpcouldz = y * decouldz
     
            # second derivatives
            dpcoulddd = 2.0*third*decouldd + y*decoulddd
            dpcoulddt = third*decouldt + y*decoulddt
            dpcouldda = third*decoulda + y*decouldda
            dpcoulddz = third*decouldz + y*decoulddz
            dpcouldtt = y * decouldtt
            dpcouldta = y * decouldta
            dpcouldtz = y * decouldtz
            dpcouldaa = y * decouldaa
            dpcouldaz = y * decouldaz
            dpcouldzz = y * decouldzz
    
    
            # entropy in erg/g/kelvin
            u0   = 3.0*b1*p1 - 5.0*c1*p2 + d1*(log(plasg) - 1.0) - e1
            du0  = 0.75*b1*p3 + 1.25*c1*p4 + d1*plasginv
            ddu0 = -0.5625*b1*p5 - 1.5625*c1*p6 - d1*plasginv*plasginv
     
            z    = -avo*ytot1*kerg
     
            scoul = z*u0
            ww    = z*du0
            x     = z*ddu0
     
            # first derivatives
            dscouldd = ww*dplasgdd
            dscouldt = ww*dplasgdt
            dscoulda = ww*dplasgda - scoul*ytot1
            dscouldz = ww*dplasgdz
     
            # second derivatives
            dscoulddd = x*dplasgdd*dplasgdd + ww*dplasgddd
            dscoulddt = x*dplasgdt*dplasgdd + ww*dplasgddt
            dscouldda = x*dplasgda*dplasgdd + ww*dplasgdda - x*ytot1*dplasgdd
            dscoulddz = x*dplasgdz*dplasgdd + ww*dplasgddz
            dscouldtt = x*dplasgdt*dplasgdt + ww*dplasgdtt
            dscouldta = x*dplasgda*dplasgdt + ww*dplasgdta - x*ytot1*dplasgdt
            dscouldtz = x*dplasgdz*dplasgdt + ww*dplasgdtz
            dscouldaa = x*dplasgda*dplasgda + ww*dplasgdaa - x*ytot1*dplasgda \
                         - ww*dplasgda*ytot1 + 2.0*scoul*ytot1*ytot1
            dscouldaz = x*dplasgdz*dplasgda + ww*dplasgdaz \
                         - ww*dplasgdz*ytot1
            dscouldzz = x*dplasgdz*dplasgdz + ww*dplasgdzz
             
        elif (plasg < 1.0):
            
            x        = sqrt(plasg)
            p1       = plasg*x
            p2       = plasg**b2
            p3       = x
            p4       = p2*plasginv
            p5       = p3*plasginv
            p6       = p4*plasginv
     
            u0   = c2*p1 - third*a2*p2
            du0  = 1.5*c2*p3 - third*a2*b2*p4
            ddu0 = 0.75*c2*p5 - third*a2*b2*(b2-1.0)*p6
      
            # pressure
            pcoul    = -pion * u0
     
            x        = pion*du0
            y        = pion*ddu0
     
            # first derivatives
            dpcouldd = -dpiondd*u0 - x*dplasgdd
            dpcouldt = -dpiondt*u0 - x*dplasgdt
            dpcoulda = -dpionda*u0 - x*dplasgda
            dpcouldz = -dpiondz*u0 - x*dplasgdz
    
     
            # second derivatives
            dpcoulddd = -dpionddd*u0 - (2.0*dpiondd*du0 \
                         + y*dplasgdd)*dplasgdd - x*dplasgddd
            dpcoulddt = -dpionddt*u0 - dpiondd*du0*dplasgdt \
                         - (dpiondt*du0 + y*dplasgdt)*dplasgdd - x*dplasgddt
            dpcouldda = -dpiondda*u0 - dpiondd*du0*dplasgda \
                         - (dpionda*du0 + y*dplasgda)*dplasgdd - x*dplasgdda
            dpcoulddz = -dpionddz*u0 - dpiondd*du0*dplasgdz \
                         - (dpiondz*du0 + y*dplasgdz)*dplasgdd - x*dplasgddz
            dpcouldtt = -dpiondtt*u0 - (2.0*dpiondt*du0 \
                         + y*dplasgdt)*dplasgdt - x*dplasgdtt
            dpcouldta = -dpiondta*u0 - dpiondt*du0*dplasgda \
                         - (dpionda*du0 + y*dplasgda)*dplasgdt - x*dplasgdta
            dpcouldtz = -dpiondtz*u0 - dpiondt*du0*dplasgdz \
                         - (dpiondz*du0 + y*dplasgdz)*dplasgdt - x*dplasgdtz
            dpcouldaa = -dpiondaa*u0 - (2.0*dpionda*du0 \
                         + y*dplasgda)*dplasgda - x*dplasgdaa
            dpcouldaz = -dpiondaz*u0 - dpionda*du0*dplasgdz \
                         - (dpiondz*du0 +  y*dplasgdz)*dplasgda - x*dplasgdaz
            dpcouldzz = -dpiondzz*u0 - (2.0*dpiondz*du0 \
                         + y*dplasgdz)*dplasgdz - x*dplasgdzz
     
     
            # energy in erg/gr
            z        = 3.0*deninv
            y        = -z*deninv
     
            ecoul    = z*pcoul
     
            x = deninv*deninv
    
            # first derivatives
            decouldd = z*dpcouldd - ecoul*deninv
            decouldt = z*dpcouldt
            decoulda = z*dpcoulda
            decouldz = z*dpcouldz
            decoul_deta = 0.0
            decoul_dbeta = decouldt*mecc/kerg
    
            # second derivatives
            decoulddd = z*dpcoulddd + y*dpcouldd \
                         + (2.0*ecoul - 3.0*dpcouldd)*x
            decoulddt = z*dpcoulddt - 3.0*dpcouldt*x
            decouldda = z*dpcouldda - 3.0*dpcoulda*x
            decoulddz = z*dpcoulddz - 3.0*dpcouldz*x
            decouldtt = z*dpcouldtt
            decouldta = z*dpcouldta
            decouldtz = z*dpcouldtz
            decouldaa = z*dpcouldaa
            decouldaz = z*dpcouldaz
            decouldzz = z*dpcouldzz
     
     
     
             # entropy in erg/g/kelvin
            u0    = c2*p1 - a2/b2*(b2-1.0)*p2
            du0   = 1.5*c2*p3 - a2*(b2-1.0)*p4
            ddu0  = 0.75*c2*p5 - a2*(b2-1.0)*(b2-1.0)*p6
            z     = -avo*ytot1*kerg
            y     = -z*ytot1
     
            scoul = z*u0
            x     = z*du0
            y     = z*ddu0
    
            # first derivatives
            dscouldd = x*dplasgdd
            dscouldt = x*dplasgdt
            dscoulda = x*dplasgda - scoul*ytot1
            dscouldz = x*dplasgdz
     
            # second derivatives
            dscoulddd = y*dplasgdd*dplasgdd + x*dplasgddd
            dscoulddt = y*dplasgdt*dplasgdd + x*dplasgddt
            dscouldda = y*dplasgda*dplasgdd + x*dplasgdda - x*ytot1*dplasgdd
            dscoulddz = y*dplasgdz*dplasgdd + x*dplasgddz
            dscouldtt = y*dplasgdt*dplasgdt + x*dplasgdtt
            dscouldta = y*dplasgda*dplasgdt + x*dplasgdta - x*ytot1*dplasgdt
            dscouldtz = y*dplasgdz*dplasgdt + x*dplasgdtz
            dscouldaa = y*dplasgda*dplasgda + x*dplasgdaa - x*ytot1*dplasgda \
                         - x*dplasgda*ytot1 + 2.0*scoul*ytot1*ytot1
            dscouldaz = y*dplasgdz*dplasgda + x*dplasgdaz - x*dplasgdz*ytot1
            dscouldzz = y*dplasgdz*dplasgdz + x*dplasgdzz
            
    ########################### Electron-Positron Terms  #######################
    ############################################################################
            
    ytot1   = 1.0/abar_var
    deni    = 1.0/den_var
    kt      = beta_var*mecc
    kti     = 1.0/kt
    t       = kt/kerg
    #print('temp=',t)
    tempinv = 1/t
    beta_var32 = beta_var**1.5
    beta_var52 = beta_var*beta_var32
    dbetadt    = kerg/mecc
    
    # ion number density in 1/cm**3
    xni     = avo * ytot1 * den_var
    
    # first derivatve of ion number density
    dxnidd  = avo * ytot1
    dxnidt  = 0.0
    dxnida  = -xni * ytot1
    dxnidz  = 0.0
    
    # second derivative of ion number density
    dxniddd = 0.0
    dxniddt = 0.0
    dxnidda = -dxnidd*ytot1
    dxniddz = 0.0
    dxnidtt = 0.0
    dxnidta = 0.0
    dxnidtz = 0.0
    dxnidaa = -2.0 * dxnida * ytot1
    dxnidaz = 0.0
    dxnidzz = 0.0
        
    
      # assume fully ionized
    chi        = 0.0
    chifac     = 0.0
    dchifacdt  = 0.0
    dchifacdz  = 0.0
    dchifacdtt = 0.0
    dchifacdtz = 0.0
    dchifacdzz = 0.0
    saha       = 0.0
    dsaha_dd   = 0.0
    dsaha_dt   = 0.0
    dsaha_da   = 0.0
    dsaha_dz   = 0.0
    dsaha_deta = 0.0
    dsaha_dbeta = 0.0
    dsaha_ddd  = 0.0
    dsaha_ddt  = 0.0
    dsaha_dda  = 0.0
    dsaha_ddz  = 0.0
    dsaha_dtt  = 0.0
    dsaha_dta  = 0.0
    dsaha_dtz  = 0.0
    dsaha_daa  = 0.0
    dsaha_daz  = 0.0
    dsaha_dzz  = 0.0
    dsaha_deta_dd = 0.0
    dsaha_deta_dt = 0.0
    dsaha_deta_da = 0.0
    dsaha_deta_dz = 0.0
    dsaha_deta2   = 0.0
        
        #do a simple saha approach
        
    if (ionized == 0.0):
        #print('ionized')
    
        denion     = 0.1
        chi        = hion * ev2erg * zbar_var
        chifac     = chi*kti
        chifacbeta = chi/mecc/beta_var
        dchifacbeta_dbeta = -chifacbeta/beta_var
        dchifacdt  = -chifac/t
        dchifacdz  = chifac/zbar_var
        dchifacdtt = -2.0*dchifacdt/t
        dchifacdtz = -dchifacdz/t
        dchifacdzz = 0.0
        
        yy         = chifac - den_var/denion
        #print('yy',yy,chifac,chi)
        # completely neutral; set it for a fake convergence
        if (yy > 200.0):
            saha       = 1.0e90
            f          = 0.0
            df         = 1.0
            etaele     = -100.0
         #  if (mode == 0):
          #      return
         # ionization possible
        elif (yy > -200.0):
            ww   = min(200.0,chifac + eta_var - den_var/denion)
            saha = 2.0 * exp(ww)
            if (ww != 200.0):
                
                # first derivative
                dsaha_dd   = -saha/denion
                dsaha_dt   = saha * dchifacdt
                dsaha_da   = 0.0
                dsaha_dz   = saha * dchifacdz
                dsaha_deta  = saha
                dsaha_dbeta = saha * dchifacbeta_dbeta
    # saha factor, the fraction ionized
    sfac       = 1.0/(1.0 + saha)
    
    # first derivatives
    y          = -sfac*sfac
    dsfac_dd   = y*dsaha_dd
    dsfac_dt   = y*dsaha_dt
    dsfac_da   = y*dsaha_da
    dsfac_dz   = y*dsaha_dz
    dsfac_deta = y*dsaha_deta
    dsfac_dbeta= y*dsaha_dbeta
        
    # second derivatives
    ww            = -2.0*sfac
    dsfac_ddd     = ww*dsfac_dd*dsaha_dd + y*dsaha_ddd
    dsfac_ddt     = ww*dsfac_dt*dsaha_dd + y*dsaha_ddt
    dsfac_dda     = ww*dsfac_da*dsaha_dd + y*dsaha_dda
    dsfac_ddz     = ww*dsfac_dz*dsaha_dd + y*dsaha_ddz
    dsfac_dtt     = ww*dsfac_dt*dsaha_dt + y*dsaha_dtt
    dsfac_dta     = ww*dsfac_da*dsaha_dt + y*dsaha_dta
    dsfac_dtz     = ww*dsfac_dz*dsaha_dt + y*dsaha_dtz
    dsfac_daa     = ww*dsfac_da*dsaha_da + y*dsaha_daa
    dsfac_daz     = ww*dsfac_dz*dsaha_da + y*dsaha_daz
    dsfac_dzz     = ww*dsfac_dz*dsaha_dz + y*dsaha_dzz
    dsfac_deta_dd = ww*dsfac_dd*dsaha_deta + y*dsaha_deta_dd
    dsfac_deta_dt = ww*dsfac_dt*dsaha_deta + y*dsaha_deta_dt
    dsfac_deta_da = ww*dsfac_da*dsaha_deta + y*dsaha_deta_da
    dsfac_deta_dz = ww*dsfac_dz*dsaha_deta + y*dsaha_deta_dz
    dsfac_deta2   = ww*dsfac_deta*dsaha_deta + y*dsaha_deta2
    
    # effective charge
    zeff       = zbar_var * sfac
    #print(zeff)
    # first derivatives
    dzeff_dd   = zbar_var * dsfac_dd
    #print('dzdd',dzeff_dd)
    dzeff_dt   = zbar_var * dsfac_dt
    #print('dzdt',dzeff_dt)
    dzeff_da   = zbar_var * dsfac_da
    dzeff_dz   = sfac
    dzeff_deta  = zbar_var * dsfac_deta
    dzeff_dbeta = zbar_var * dsfac_dbeta
    #print(dzeff_dd,dzeff_dt,dzeff_da,dzeff_dz,dzeff_deta,dzeff_dbeta)
    # second derivatives
    dzeff_ddd  = zbar_var*dsfac_ddd
    dzeff_ddt  = zbar_var*dsfac_ddt
    dzeff_dda  = zbar_var*dsfac_dda
    dzeff_ddz  = zbar_var*dsfac_ddz
    dzeff_dtt  = zbar_var*dsfac_dtt
    dzeff_dta  = zbar_var*dsfac_dta
    dzeff_dtz  = zbar_var*dsfac_dtz
    dzeff_daa  = zbar_var*dsfac_daa
    dzeff_daz  = zbar_var*dsfac_daz
    dzeff_dzz  = dsfac_dz
    dzeff_deta_dd = zbar_var*dsfac_deta_dd
    dzeff_deta_dt = zbar_var*dsfac_deta_dt
    dzeff_deta_da = zbar_var*dsfac_deta_da
    dzeff_deta_dz = dsfac_deta + zbar_var*dsfac_deta_dz
    dzeff_deta2   = zbar_var * dsfac_deta2
    
        
    # number density of free electrons
    xne        = xni * zeff
    
    # first derivatives
    dxne_dd    = dxnidd * zeff + xni*dzeff_dd
    #print('dxnidd',dxnidd)
    dxne_dt    = dxnidt * zeff + xni*dzeff_dt
    #print('dxnidt',dxnidt)
    dxne_da    = dxnida * zeff + xni*dzeff_da
    dxne_dz    = dxnidz * zeff + xni*dzeff_dz
    dxne_deta  = xni*dzeff_deta
    dxne_dbeta = xni*dzeff_dbeta
    # second derivatives
    dxne_ddd = dxniddd*zeff + 2.0*dxnidd*dzeff_dd + xni*dzeff_ddd
    dxne_ddt = dxniddt*zeff + dxnidd*dzeff_dt \
               + dxnidt*dzeff_dt + xni*dzeff_ddt
    dxne_dda = dxnidda*zeff + dxnidd*dzeff_da \
               + dxnida*dzeff_dt + xni*dzeff_dda
    dxne_ddz = dxniddz*zeff + dxnidd*dzeff_dz \
               + dxnidz*dzeff_dt + xni*dzeff_ddz
    dxne_dtt = dxnidtt*zeff + 2.0*dxnidt*dzeff_dt + xni*dzeff_dtt
    dxne_dta = dxnidta*zeff + dxnidt*dzeff_da \
               + dxnida*dzeff_dt + xni*dzeff_dta
    dxne_dtz = dxnidtz*zeff + dxnidt*dzeff_dz \
               + dxnidz*dzeff_dt + xni*dzeff_dtz
    dxne_daa = dxnidaa*zeff + 2.0*dxnida*dzeff_da + xni*dzeff_daa
    dxne_daz = dxnidaz*zeff + dxnida*dzeff_dz \
               + dxnidz*dzeff_da + xni*dzeff_daz
    dxne_dzz = dxnidzz*zeff + 2.0*dxnidz*dzeff_dz + xni*dzeff_dzz
    dxne_deta_dd = dxnidd*dzeff_deta + xni*dzeff_deta_dd
    dxne_deta_dt = dxnidt*dzeff_deta + xni*dzeff_deta_dt
    dxne_deta_dz = dxnidz*dzeff_deta + xni*dzeff_deta_dz
    dxne_deta_da = dxnida*dzeff_deta + xni*dzeff_deta_da
    dxne_deta2   = xni*dzeff_deta2
    
    #print (xne,dsfac_deta,dsfac_dbeta,y)
    #################################################################
    # get the fermi-dirac integral electron contribution
    [f12, f12eta, f12beta, f12eta2, f12beta2, f12etabeta]= FermiDiracIntegralCalculator(0.5, eta_var, beta_var)
    [f32, f32eta, f32beta, f32eta2, f32beta2, f32etabeta]= FermiDiracIntegralCalculator(1.5, eta_var, beta_var)
    beta_var12 = sqrt(beta_var)
    beta_var32 = beta_var*beta_var12
     
    zz   = xconst * beta_var32
    ww   = xconst * 1.5 * beta_var12
    yy   = f12 + beta_var * f32
    dum1 = f12eta + beta_var*f32eta
    dum2 = f12beta + f32 + beta_var*f32beta
    #print(f12,f32,f12beta,f32,beta_var,f32beta)
    xnefer             = zz * yy
    #print(xconst,beta_var32,f12,beta_var,f32)
    dxnefer_deta       = zz * dum1
    dxnefer_dbeta      = ww*yy + zz*dum2
    dxnefer_deta2      = zz * (f12eta2 + beta_var * f32eta2)
    dxnefer_dbeta2     = 0.5*ww/beta_var*yy + 2.0*ww*dum2 + zz*(f12beta2 + 2.0*f32beta + beta_var*f32beta2)
    dxnefer_deta_dbeta = ww*dum1 + zz*(f12etabeta + f32eta + beta_var*f32etabeta)
    #################################################################
    # if the temperature is not too low, get the positron contributions
    # chemical equilibrium means etaele + etapos = eta_photon = 0.
    etapos             = 0.0
    detap_deta         = 0.0
    detap_dbeta        = 0.0
    detap_deta2        = 0.0
    detap_dbeta2       = 0.0
    detap_deta_dbeta   = 0.0
    xnpfer             = 0.0
    dxnpfer_deta       = 0.0
    dxnpfer_dbeta      = 0.0
    dxnpfer_deta2      = 0.0
    dxnpfer_dbeta2     = 0.0
    dxnpfer_deta_dbeta = 0.0
    positron_start     = 0.02
    aa                 = eta_var
    ppos              = 0.0
    dppos_detap       = 0.0
    dppos_dbeta       = 0.0
    dppos_detap2      = 0.0
    dppos_dbeta2      = 0.0
    dppos_detap_dbeta = 0.0
    epos              = 0.0
    depos_detap       = 0.0
    depos_dbeta       = 0.0
    depos_deta        = 0.0
    depos_detap2      = 0.0
    depos_dbeta2      = 0.0
    depos_detap_dbeta = 0.0
         
    if (beta_var > positron_start):
        etapos           = -aa - 2.0/beta_var
        detap_deta       = -1.0
        detap_deta2      = 0.0
        detap_dbeta      = 2.0/beta_var**2
        detap_dbeta2     = -4.0/beta_var**3
        detap_deta_dbeta = 0.0
     
        [f12, f12eta, f12beta, f12eta2, f12beta2, f12etabeta]= FermiDiracIntegralCalculator(0.5, etapos, beta_var)
        [f32, f32eta, f32beta, f32eta2, f32beta2, f32etabeta]= FermiDiracIntegralCalculator(1.5, etapos, beta_var)
     
        zz   = xconst * beta_var32
        ww   = xconst * 1.5 * beta_var12
        yy   = f12 + beta_var * f32
        dum1 = f12eta + beta_var*f32eta
        dum2 = f12beta + f32 + beta_var*f32beta
     
        xnpfer              = zz * yy
        dxnpfer_detap       = zz * dum1
        dxnpfer_dbeta       = ww*yy + zz*dum2
        dxnpfer_detap2      = zz * (f12eta2 + beta_var * f32eta2)
        dxnpfer_dbeta2      = 0.5*ww/beta_var*yy + 2.0*ww*dum2 \
                             + zz*(f12beta2 + 2.0*f32beta + beta_var*f32beta2)
        dxnpfer_detap_dbeta = ww*dum1 \
                             + zz*(f12etabeta + f32eta + beta_var*f32etabeta)
       
        # convert the etap derivatives to eta derivatives
        # all derived from the operator dxp = dxp/detap detap + dxp/dbeta dbeta
     
        dxnpfer_deta  = dxnpfer_detap * detap_deta
        dxnpfer_dbeta = dxnpfer_dbeta + dxnpfer_detap * detap_dbeta
        dxnpfer_deta2 = dxnpfer_detap2 * detap_deta**2 \
                        + dxnpfer_detap * detap_deta2
        dxnpfer_dbeta2 = dxnpfer_dbeta2 \
                        + 2.0 * dxnpfer_detap_dbeta * detap_dbeta \
                        + dxnpfer_detap2 * detap_dbeta**2 \
                        + dxnpfer_detap * detap_dbeta2
        dxnpfer_deta_dbeta = dxnpfer_detap2 * detap_dbeta * detap_deta \
                        + dxnpfer_detap_dbeta * detap_deta \
                        + dxnpfer_detap * detap_deta_dbeta
    
       
    
    
    # all the derivatives are in terms of eta and beta.
    # we want to convert to temperature, density, abar and zbar derivatives.
    # so, after the root find above on eta we have xne = xnefer - xnpfer.
    # taking the derivative of this for property p
    # dxne/deta deta + dxne/dp dp = dxnefer/deta deta + dxnefer/dbeta dbeta
    # solving for the unknown eta derivative yields
    # deta/dp = (dxne_dp - dxnefer/dbeta dbeta/dp) / (dxnefer/deta - dxne/deta)
     
     
    dxep_deta       = dxnefer_deta  - dxnpfer_deta
    dxep_dbeta      = dxnefer_dbeta - dxnpfer_dbeta
    dxep_deta2      = dxnefer_deta2  - dxnpfer_deta2
    dxep_dbeta2     = dxnefer_dbeta2 - dxnpfer_dbeta2
    dxep_deta_dbeta = dxnefer_deta_dbeta  - dxnpfer_deta_dbeta
    #print('dxep_deta', dxep_deta,'dxep_dbeta',dxep_dbeta )
     
    y          = 1.0/(dxep_deta - dxne_deta)
    #print('dxep',dxep_deta,'dxne',dxne_deta)
    # the all important first derivatives of eta
    detadd = dxne_dd * y
    #print ('detadd',detadd)
    detadt = (dxne_dt - dxep_dbeta*dbetadt) * y
    #print ('detadt',detadt)
    detada = dxne_da * y
    #print ('detada',detada)
    detadz = dxne_dz * y
    #print ('detadz',detadz)
    # second derivatives
    detaddd = dxne_ddd*y - detadd*y*(dxep_deta2 - dxne_deta2)*detadd
    detaddt = dxne_ddt*y - detadd*y*(dxep_deta2*detadt \
              + dxep_deta_dbeta*dbetadt - dxne_deta2*detadt)
    detadda = dxne_dda*y - detadd*y*(dxep_deta2 - dxne_deta2)*detada
    detaddz = dxne_ddz*y - detadd*y*(dxep_deta2 - dxne_deta2)*detadz
     
    detadtt = (dxne_dtt - (dxep_deta_dbeta*detadt \
              + dxep_dbeta2*dbetadt)*dbetadt)*y \
              - detadt*y*(dxep_deta2*detadt \
              + dxep_deta_dbeta*dbetadt - dxne_deta2*detadt)
     
    detadta = (dxne_dta - dxep_deta_dbeta*detada*dbetadt)*y \
              - detadt*y*(dxep_deta2*detada - dxne_deta2*detada)
     
    detadtz = (dxne_dtz - dxep_deta_dbeta*detadz*dbetadt)*y \
              - detadt*y*(dxep_deta2*detadz - dxne_deta2*detadz)
     
    detadaa = dxne_daa*y - detada*y*(dxep_deta2 - dxne_deta2)*detada
    detadaz = dxne_daz*y - detada*y*(dxep_deta2 - dxne_deta2)*detadz
    detadzz = dxne_dzz*y - detadz*y*(dxep_deta2 - dxne_deta2)*detadz
    
     
    # first derivatives of the effective charge
    dzeffdd = dzeff_deta*detadd + dzeff_dd
    dzeffdt = dzeff_deta*detadt + dzeff_dt
    dzeffda = dzeff_deta*detada + dzeff_da
    dzeffdz = dzeff_deta*detadz + dzeff_dz
     
    # second derivatives
    dzeffddd = dzeff_deta_dd*detadd + dzeff_deta*detaddd + dzeff_ddd
    dzeffddt = dzeff_deta_dt*detadd + dzeff_deta*detaddt + dzeff_ddt
    dzeffdda = dzeff_deta_da*detadd + dzeff_deta*detadda + dzeff_dda
    dzeffddz = dzeff_deta_dz*detadd + dzeff_deta*detaddz + dzeff_ddz
    dzeffdtt = dzeff_deta_dt*detadt + dzeff_deta*detadtt + dzeff_dtt
    dzeffdta = dzeff_deta_da*detadt + dzeff_deta*detadta + dzeff_dta
    dzeffdtz = dzeff_deta_dz*detadt + dzeff_deta*detadtz + dzeff_dtz
    dzeffdaa = dzeff_deta_da*detada + dzeff_deta*detadaa + dzeff_daa
    dzeffdaz = dzeff_deta_dz*detada + dzeff_deta*detadaz + dzeff_daz
    dzeffdzz = dzeff_deta_dz*detadz + dzeff_deta*detadzz + dzeff_dzz
     
     
       # first derivatives of the electron number density
    dxnedd = dxnidd * zeff + xni * dzeffdd
    #print('zeff',zeff)
    dxnedt = dxnidt * zeff + xni * dzeffdt
    #print('dxnedt',dxnedt)
    dxneda = dxnida * zeff + xni * dzeffda
    dxnedz = dxnidz * zeff + xni * dzeffdz
    #print(dxnidt,zeff,xni,dzeffdt)
     
       # second derivatives
    dxneddd = dxniddd*zeff+dxnidd*dzeffdd+dxnidd*dzeffdd+xni*dzeffddd
    dxneddt = dxniddt*zeff+dxnidd*dzeffdt+dxnidt*dzeffdd+xni*dzeffddt
    dxnedda = dxnidda*zeff+dxnidd*dzeffda+dxnida*dzeffdd+xni*dzeffdda
    dxneddz = dxniddz*zeff+dxnidd*dzeffdz+dxnidz*dzeffdd+xni*dzeffddz
    dxnedtt = dxnidtt*zeff+dxnidt*dzeffdt+dxnidt*dzeffdt+xni*dzeffdtt
    dxnedta = dxnidta*zeff+dxnidt*dzeffda+dxnida*dzeffdt+xni*dzeffdta
    dxnedtz = dxnidtz*zeff+dxnidt*dzeffdz+dxnidz*dzeffdt+xni*dzeffdtz
    dxnedaa = dxnidaa*zeff+dxnida*dzeffda+dxnida*dzeffda+xni*dzeffdaa
    dxnedaz = dxnidaz*zeff+dxnida*dzeffdz+dxnidz*dzeffda+xni*dzeffdaz
    dxnedzz = dxnidzz*zeff+dxnidz*dzeffdz+dxnidz*dzeffdz+xni*dzeffdzz
      
       # first derivatives of the fermi integral electron number densities
    dxneferdd = dxnefer_deta * detadd
    dxneferdt = dxnefer_deta * detadt + dxnefer_dbeta * dbetadt
    dxneferda = dxnefer_deta * detada
    dxneferdz = dxnefer_deta * detadz
     
      # second derivatives
    dxneferddd = dxnefer_deta2*detadd*detadd + dxnefer_deta*detaddd
    dxneferddt = dxnefer_deta2*detadt*detadd \
                 + dxnefer_deta_dbeta * dbetadt * detadd \
                 + dxnefer_deta *detaddt
    dxneferdda = dxnefer_deta2*detada*detadd + dxnefer_deta*detadda
    dxneferddz = dxnefer_deta2*detadz*detadd + dxnefer_deta*detaddz
    dxneferdtt = dxnefer_deta2*detadt*detadt + dxnefer_deta*detadtt \
                 + 2.0*dxnefer_deta_dbeta * detadt * dbetadt \
                 + dxnefer_dbeta2 * dbetadt * dbetadt 
    dxneferdta = dxnefer_deta2 * detada * detadt \
                 + dxnefer_deta_dbeta * dbetadt * detada \
                 + dxnefer_deta * detadta 
    dxneferdtz = dxnefer_deta2 * detadz * detadt \
                 + dxnefer_deta_dbeta * dbetadt * detadz \
                 + dxnefer_deta * detadtz
    dxneferdaa = dxnefer_deta2*detada*detada + dxnefer_deta*detadaa
    dxneferdaz = dxnefer_deta2*detadz*detada + dxnefer_deta*detadaz
    dxneferdzz = dxnefer_deta2*detadz*detadz + dxnefer_deta*detadzz
    
    # first derivatives of the fermi integral positron number densities
    dxnpferdd = dxnpfer_deta * detadd
    dxnpferdt = dxnpfer_deta * detadt + dxnpfer_dbeta * dbetadt
    dxnpferda = dxnpfer_deta * detada
    dxnpferdz = dxnpfer_deta * detadz
     
    
    # second derivatives
    dxnpferddd = dxnpfer_deta2*detadd*detadd + dxnpfer_deta*detaddd
    dxnpferddt = dxnpfer_deta2*detadt*detadd \
                 + dxnpfer_deta_dbeta * dbetadt * detadd \
                 + dxnpfer_deta *detaddt
    dxnpferdda = dxnpfer_deta2*detada*detadd + dxnpfer_deta*detadda
    dxnpferddz = dxnpfer_deta2*detadz*detadd + dxnpfer_deta*detaddz
    dxnpferdtt = dxnpfer_deta2*detadt*detadt + dxnpfer_deta*detadtt \
                 + 2.0*dxnpfer_deta_dbeta * detadt * dbetadt \
                 + dxnpfer_dbeta2 * dbetadt * dbetadt 
    dxnpferdta = dxnpfer_deta2 * detada * detadt \
                 + dxnpfer_deta_dbeta * dbetadt * detada \
                 + dxnpfer_deta * detadta
    dxnpferdtz = dxnpfer_deta2 * detadz * detadt \
                 + dxnpfer_deta_dbeta * dbetadt * detadz \
                 + dxnpfer_deta * detadtz
    dxnpferdaa = dxnpfer_deta2*detada*detada + dxnpfer_deta*detadaa
    dxnpferdaz = dxnpfer_deta2*detadz*detada + dxnpfer_deta*detadaz
    dxnpferdzz = dxnpfer_deta2*detadz*detadz + dxnpfer_deta*detadzz
    
    ## Need Entropy derivative dS2dde=-t^2(drhodt)^-1
    zeffsqi = 1/(zeff*zeff)
    drhodt = abar_var*((dxneferdt-dxnpferdt)/zeff - zeffsqi*dzeffdt*xne)/avo
    #print('drhodt',drhodt,dxneferdt,dxnpferdt,dzeffdt)
    tempinvsq = tempinv * tempinv
    #dsded  = -tempinvsq/drhodt
    #print(dsded,drhodt,dxneferdt)
    
       # ion section:
     
       # number density in 1/cm**3,
    xni     = avo * ytot1 * den_var
    dxnidd  = avo * ytot1
    dxnidt  = 0.0
    dxnida  = -xni * ytot1
    dxnidz  = 0.0
    
       
       # ! for the electrons
    [f32, f32eta, f32beta, f32eta2, f32beta2, f32etabeta]= FermiDiracIntegralCalculator(1.5, eta_var, beta_var)
    [f52, f52eta, f52beta, f52eta2, f52beta2, f52etabeta]= FermiDiracIntegralCalculator(2.5, eta_var, beta_var)
    
       # pressure in erg/cm**3
    yy   = pconst * beta_var52
    ww   = pconst * 2.5 * beta_var32
    dum1 = f32 + 0.5*beta_var*f52
    dum2 = f32beta + 0.5*f52 + 0.5*beta_var*f52beta
    dum3 = f32eta + 0.5*beta_var*f52eta
     
    pele         = yy * dum1
       # print(pele)
    dpele_deta   = yy * dum3
    dpele_dbeta  = ww * dum1 + yy * dum2
    dpele_deta2  = yy * (f32eta2 + 0.5*beta_var*f52eta2)
    dpele_dbeta2 = pconst*beta_var12*3.75*dum1 + 2.0*ww*dum2 \
                   + yy * (f32beta2 + f52beta + 0.5*beta_var*f52beta2)
    dpele_deta_dbeta = ww*dum3 + yy*(f32etabeta +0.5*f52eta \
                       + 0.5*beta_var*f52etabeta)
     
       # first derivatives of the electron pressure
    dpeledd = dpele_deta * detadd
    dpeledt = dpele_deta * detadt + dpele_dbeta * dbetadt
    dpeleda = dpele_deta * detada
    dpeledz = dpele_deta * detadz
     
       # second derivatives
    dpeleddd = dpele_deta2*detadd*detadd + dpele_deta*detaddd
    dpeleddt = dpele_deta2*detadt*detadd \
               + dpele_deta_dbeta * dbetadt * detadd \
               + dpele_deta *detaddt
    dpeledda = dpele_deta2*detada*detadd + dpele_deta*detadda
    dpeleddz = dpele_deta2*detadz*detadd + dpele_deta*detaddz
    dpeledtt = dpele_deta2*detadt*detadt + dpele_deta*detadtt \
               + 2.0*dpele_deta_dbeta * detadt * dbetadt \
               + dpele_dbeta2 * dbetadt * dbetadt
    dpeledta = dpele_deta2 * detada * detadt \
               + dpele_deta_dbeta * dbetadt * detada \
               + dpele_deta * detadta
    dpeledtz = dpele_deta2 * detadz * detadt \
               + dpele_deta_dbeta * dbetadt * detadz \
               + dpele_deta * detadtz
    dpeledaa = dpele_deta2*detada*detada + dpele_deta*detadaa
    dpeledaz = dpele_deta2*detadz*detada + dpele_deta*detadaz
    dpeledzz = dpele_deta2*detadz*detadz + dpele_deta*detadzz
     
       #! NOTE!!: energy in erg/cm**3 and will be converted to erg/gr later
    yy   = econst * beta_var52
    ww   = econst * 2.5 * beta_var32
    dum1 = f32 + beta_var*f52
    dum2 = f32beta + f52 + beta_var*f52beta
    dum3 = f32eta + beta_var*f52eta
    
    eele        = yy * dum1
    deele_deta  = yy * dum3
    #print(yy,dum3)
    deele_dbeta = ww * dum1 + yy * dum2
    deele_deta2  = yy * (f32eta2 + beta_var*f52eta2)
    deele_dbeta2 = econst*beta_var12*3.75*dum1 + 2.0*ww*dum2 \
                   + yy * (f32beta2 + 2.0*f52beta + beta_var*f52beta2)
    deele_deta_dbeta = ww*dum3 \
                       + yy*(f32etabeta + f52eta + beta_var*f52etabeta)
    
    # first derivatives of the electron energy
    deeledd = deele_deta * detadd
    deeledt = deele_deta * detadt + deele_dbeta * dbetadt
    #print(deele_deta,detadt,deele_dbeta,dbetadt)
    deeleda = deele_deta * detada
    deeledz = deele_deta * detadz
     
    # second derivatives
    deeleddd = deele_deta2*detadd*detadd + deele_deta*detaddd
    deeleddt = deele_deta2*detadt*detadd + deele_deta_dbeta * dbetadt * detadd + deele_deta *detaddt
    deeledda = deele_deta2*detada*detadd + deele_deta*detadda
    deeleddz = deele_deta2*detadz*detadd + deele_deta*detaddz
    deeledtt = deele_deta2*detadt*detadt + deele_deta*detadtt  \
              + 2.0*deele_deta_dbeta * detadt * dbetadt + deele_dbeta2 * dbetadt * dbetadt
    deeledta = deele_deta2 * detada * detadt + deele_deta_dbeta * dbetadt * detada + deele_deta * detadta
    deeledtz = deele_deta2 * detadz * detadt + deele_deta_dbeta * dbetadt * detadz + deele_deta * detadtz
    deeledaa = deele_deta2*detada*detada + deele_deta*detadaa
    deeledaz = deele_deta2*detadz*detada + deele_deta*detadaz
    deeledzz = deele_deta2*detadz*detadz + deele_deta*detadzz
     
    # ! for the positrons
    ppos              = 0.0
    dppos_detap       = 0.0
    dppos_dbeta       = 0.0
    dppos_detap2      = 0.0
    dppos_dbeta2      = 0.0
    dppos_detap_dbeta = 0.0
    epos              = 0.0
    depos_detap       = 0.0
    depos_dbeta       = 0.0
    depos_detap2      = 0.0
    depos_dbeta2      = 0.0
    depos_detap_dbeta = 0.0
    depos_deta        = 0.0
    depos_dbeta       = 0.0
        
    if (beta_var > positron_start):
        
        [f32, f32eta, f32beta, f32eta2, f32beta2, f32etabeta]= FermiDiracIntegralCalculator(1.5, etapos, beta_var)
        [f52, f52eta, f52beta, f52eta2, f52beta2, f52etabeta]= FermiDiracIntegralCalculator(2.5, etapos, beta_var)
    
    # pressure
        yy   = pconst * beta_var52
        ww   = pconst * 2.5 * beta_var32
        dum1 = f32 + 0.5*beta_var*f52
        dum2 = f32beta + 0.5*f52 + 0.5*beta_var*f52beta
        dum3 = f32eta + 0.5*beta_var*f52eta
     
        ppos         = yy * dum1
        dppos_detap  = yy * dum3
        dppos_dbeta  = ww * dum1 + yy * dum2
        dppos_detap2 = yy * (f32eta2 + 0.5*beta_var*f52eta2)
        dppos_dbeta2 = pconst*beta_var12*3.75*dum1 + 2.0*ww*dum2 \
                       + yy * (f32beta2 + f52beta + 0.5*beta_var*f52beta2)
        dppos_detap_dbeta = ww*dum3 + yy*(f32etabeta +0.5*f52eta \
                            + 0.5*beta_var*f52etabeta)
     
    #! NOTE!!: energy in erg/cm**3 and will be converted to erg/gr later
        yy   = econst * beta_var52 
        ww   = econst * 2.5 * beta_var32
        dum1 = f32 + beta_var*f52
        dum2 = f32beta + f52 + beta_var*f52beta
        dum3 = f32eta + beta_var*f52eta
     
        epos          = yy * dum1
        depos_detap   = yy * dum3
        depos_dbeta   = ww * dum1 + yy * dum2
        depos_detap2  = yy * (f32eta2 + beta_var*f52eta2)
        depos_dbeta2  = econst*beta_var12*3.75*dum1 + 2.0*ww*dum2 \
                        + yy * (f32beta2 + 2.0*f52beta + beta_var*f52beta2)
        depos_detap_dbeta = ww*dum3 \
                            + yy*(f32etabeta + f52eta + beta_var*f52etabeta)
                            
    # convert the etap derivatives to eta derivatives
    # all derived from the operator dxp = dxp/detap detap + dxp/dbeta dbeta
     
    dppos_deta  = dppos_detap * detap_deta
    dppos_dbeta = dppos_dbeta + dppos_detap * detap_dbeta
    dppos_deta2 = dppos_detap2 * detap_deta**2 \
                  + dppos_detap * detap_deta2
    dppos_dbeta2 = dppos_dbeta2 \
                   + 2.0 * dppos_detap_dbeta * detap_dbeta \
                   + dppos_detap2 * detap_dbeta**2 \
                   + dppos_detap * detap_dbeta2
    dppos_deta_dbeta = dppos_detap2 * detap_dbeta * detap_deta \
                       + dppos_detap_dbeta * detap_deta \
                       + dppos_detap * detap_deta_dbeta
    
    # first derivatives of the positron pressure
    dpposdd     = dppos_deta * detadd
    dpposdt     = dppos_deta * detadt + dppos_dbeta * dbetadt
    dpposda     = dppos_deta * detada
    dpposdz     = dppos_deta * detadz
     
    # second derivatives
    dpposddd = dppos_deta2*detadd*detadd + dppos_deta*detaddd
    dpposddt = dppos_deta2*detadt*detadd + dppos_deta_dbeta * dbetadt * detadd + dppos_deta *detaddt
    dpposdda = dppos_deta2*detada*detadd + dppos_deta*detadda
    dpposddz = dppos_deta2*detadz*detadd + dppos_deta*detaddz
    dpposdtt = dppos_deta2*detadt*detadt + dppos_deta*detadtt + 2.0*dppos_deta_dbeta * detadt * dbetadt + dppos_dbeta2 * dbetadt * dbetadt
    dpposdta = dppos_deta2 * detada * detadt + dppos_deta_dbeta * dbetadt * detada + dppos_deta * detadta
    dpposdtz = dppos_deta2 * detadz * detadt + dppos_deta_dbeta * dbetadt * detadz + dppos_deta * detadtz
    dpposdaa = dppos_deta2*detada*detada + dppos_deta*detadaa
    dpposdaz = dppos_deta2*detadz*detada + dppos_deta*detadaz
    dpposdzz = dppos_deta2*detadz*detadz + dppos_deta*detadzz
    
    # convert the etap derivatives to eta derivatives
    # all derived from the operator dxp = dxp/detap detap + dxp/dbeta dbeta
     
    depos_deta  = depos_detap * detap_deta
    depos_dbeta = depos_dbeta + depos_detap * detap_dbeta
    depos_deta2 = depos_detap2 * detap_deta**2 + depos_detap * detap_deta2
    depos_dbeta2 = depos_dbeta2 \
                   + 2.0 * depos_detap_dbeta * detap_dbeta \
                   + depos_detap2 * detap_dbeta**2  \
                   + depos_detap * detap_dbeta2
    depos_deta_dbeta = depos_detap2 * detap_dbeta * detap_deta \
                       + depos_detap_dbeta * detap_deta \
                       + depos_detap * detap_deta_dbeta

     
     # first derivatives of the positron energy
    deposdd     = depos_deta * detadd
    deposdt     = depos_deta * detadt + depos_dbeta * dbetadt
    deposda     = depos_deta * detada
    deposdz     = depos_deta * detadz
     
     # second derivatives
    deposddd = depos_deta2*detadd*detadd + depos_deta*detaddd
    deposddt = depos_deta2*detadt*detadd \
               + depos_deta_dbeta * dbetadt * detadd \
               + depos_deta *detaddt
    deposdda = depos_deta2*detada*detadd + depos_deta*detadda
    deposddz = depos_deta2*detadz*detadd + depos_deta*detaddz
    deposdtt = depos_deta2*detadt*detadt + depos_deta*detadtt \
               + 2.0*depos_deta_dbeta * detadt * dbetadt \
               + depos_dbeta2 * dbetadt * dbetadt
    deposdta = depos_deta2 * detada * detadt \
               + depos_deta_dbeta * dbetadt * detada \
               + depos_deta * detadta
    deposdtz = depos_deta2 * detadz * detadt \
               + depos_deta_dbeta * dbetadt * detadz \
               + depos_deta * detadtz
    deposdaa = depos_deta2*detada*detada + depos_deta*detadaa
    deposdaz = depos_deta2*detadz*detada + depos_deta*detadaz
    deposdzz = depos_deta2*detadz*detadz + depos_deta*detadzz
    
    # electron+positron pressure and its derivatives
    # note: at high temperatures and low densities, dpepdd is very small
    # and can go negative, so limit it to be positive definite
     
    pep     = pele    + ppos
    dpepdd  = max(dpeledd + dpposdd, 1.0e-30)
    dpepdt  = dpeledt + dpposdt
    dpepda  = dpeleda + dpposda
    dpepdz  = dpeledz + dpposdz
    dpepddd = dpeleddd + dpposddd
    dpepddt = dpeleddt + dpposddt
    dpepdda = dpeledda + dpposdda
    dpepddz = dpeleddz + dpposddz
    dpepdtt = dpeledtt + dpposdtt
    dpepdta = dpeledta + dpposdta
    dpepdtz = dpeledtz + dpposdtz
    dpepdaa = dpeledaa + dpposdaa
    dpepdaz = dpeledaz + dpposdaz
    dpepdzz = dpeledzz + dpposdzz
     
    # electron+positron thermal energy and its derivatives
    eep     = eele    + epos
    deepdd  = deeledd + deposdd
    deepdt  = deeledt + deposdt
    deepda  = deeleda + deposda
    deepdz  = deeledz + deposdz
    deepddd = deeleddd + deposddd
    deepddt = deeleddt + deposddt
    deepdda = deeledda + deposdda
    deepddz = deeleddz + deposddz
    deepdtt = deeledtt + deposdtt
    deepdta = deeledta + deposdta
    deepdtz = deeledtz + deposdtz
    deepdaa = deeledaa + deposdaa
    deepdaz = deeledaz + deposdaz
    deepdzz = deeledzz + deposdzz
    
    # electron entropy in erg/gr/kelvin and its derivatives
    y       = kerg*deni
    sele    = ((pele + eele)*kti - eta_var*xnefer) * y
     
    # first derivatives
    dseledd = ((dpeledd + deeledd)*kti - detadd*xnefer \
              - eta_var*dxneferdd)*y - sele*deni
    dseledt = ((dpeledt + deeledt)*kti - (pele + eele)/(kt*t) \
                   - detadt*xnefer - eta_var*dxneferdt)*y
    dseleda = ((dpeleda + deeleda)*kti - detada*xnefer \
                    - eta_var*dxneferda)*y
    dseledz = ((dpeledz + deeledz)*kti - detadz*xnefer \
                   - eta_var*dxneferdz)*y
     
    # second derivatives
    dseleddd = ((dpeleddd + deeleddd)*kti - detaddd*xnefer \
               - 2.0*detadd*dxneferdd - eta_var*dxneferddd)*y \
               - 2.0*((dpeledd + deeledd)*kti - detadd*xnefer \
               - eta_var*dxneferdd)*y*deni \
               + 2.0*sele*deni**2
    dseleddt = ((dpeleddt + deeleddt)*kti \
               - (dpeledd + deeledd)*kti/t \
               - detaddt*xnefer - detadd*dxneferdt \
               - detadt*dxneferdd - eta_var*dxneferddt)*y \
               - dseledt*deni
    dseledda = ((dpeledda + deeledda)*kti \
               - detadda*xnefer - detadd*dxneferda \
               - detada*dxneferdd - eta_var*dxneferdda)*y \
               - dseleda*deni
    dseleddz = ((dpeleddz + deeleddz)*kti \
               - detaddz*xnefer - detadd*dxneferdz \
               - detadz*dxneferdd - eta_var*dxneferddz)*y \
               - dseledz*deni
    dseledtt = ((dpeledtt + deeledtt)*kti \
               - 2.0*(dpeledt + deeledt)*kti/t \
               + 2.0*(pele + eele)*kti/t**2 \
               - detadtt*xnefer - 2.0*detadt*dxneferdt \
               - eta_var*dxneferdtt)*y
    dseledta = ((dpeledta + deeledta)*kti \
                - (dpeleda + deeleda)*kti/t \
                - detadta*xnefer - detadt*dxneferda \
                - detada*dxneferdt - eta_var*dxneferdta)*y
    dseledtz = ((dpeledtz + deeledtz)*kti \
                - (dpeledz + deeledz)*kti/t \
                - detadtz*xnefer - detadt*dxneferdz \
                - detadz*dxneferdt - eta_var*dxneferdtz)*y
    dseledaa = ((dpeledaa + deeledaa)*kti - detadaa*xnefer \
               - 2.0*detada*dxneferda - eta_var*dxneferdaa)*y
    dseledaz = ((dpeledaz + deeledaz)*kti - detadaz*xnefer \
               - detada*dxneferdz \
               - detadz*dxneferda  - eta_var*dxneferdaz)*y
    dseledzz = ((dpeledzz + deeledzz)*kti - detadzz*xnefer \
               - 2.0*detadz*dxneferdz - eta_var*dxneferdzz)*y
     
    # positron entropy in erg/gr/kelvin and its derivatives
    spos    = ((ppos + epos)/kt - etapos*xnpfer) * y
     
    # first derivatives
    dsposdd = ((dpposdd + deposdd)*kti \
              - detap_deta*detadd*xnpfer \
              - etapos*dxnpferdd)*y - spos*deni
    #print(etapos,dxnpferdd,y,spos,deni)
    dsposdt = ((dpposdt + deposdt)*kti - (ppos + epos)/(kt*t) \
              - (detap_deta*detadt + detap_dbeta*dbetadt)*xnpfer \
              - etapos*dxnpferdt)*y
    dsposda = ((dpposda + deposda)*kti \
              - detap_deta*detada*xnpfer \
              - etapos*dxnpferda)*y
    dsposdz = ((dpposdz + deposdz)*kti - detap_deta*detadz*xnpfer \
              - etapos*dxnpferdz)*y
    
    # second derivatives
    dsposddd = ((dpposddd + deposddd)*kti - detap_deta*detaddd*xnpfer \
               - 2.0*detap_deta*detadd*dxnpferdd \
               - etapos*dxnpferddd)*y \
               - 2.0*((dpposdd + deposdd)*kti \
               - detap_deta*detadd*xnpfer \
               - etapos*dxnpferdd)*y*deni \
               + 2.0*spos*deni**2
    dsposddt = ((dpposddt + deposddt)*kti \
               - (dpposdd + deposdd)*kti/t \
               - detap_deta*detaddt*xnpfer \
               - detap_deta*detadd*dxnpferdt \
               - (detap_deta*detadt + detap_dbeta*dbetadt)*dxnpferdd \
               - etapos*dxnpferddt)*y \
               - dsposdt*deni
    dsposdda = ((dpposdda + deposdda)*kti \
               - detap_deta*detadda*xnpfer \
               - detap_deta*detadd*dxnpferda \
               - detap_deta*detada*dxnpferdd \
               - etapos*dxnpferdda)*y \
               - dsposda*deni
    dsposddz = ((dpposddz + deposddz)*kti \
               - detap_deta*detaddz*xnpfer \
               - detap_deta*detadd*dxnpferdz \
               - detap_deta*detadz*dxnpferdd \
               - etapos*dxnpferddz)*y \
               - dsposdz*deni
    dsposdtt = ((dpposdtt + deposdtt)*kti \
                - 2.0*(dpposdt + deposdt)*kti/t \
                + 2.0*(ppos + epos)*kti/t**2 \
                - (detap_deta2*detadt**2 \
                + 2.0*detap_deta_dbeta*detadt*dbetadt \
                + detap_deta*detadtt \
                + detap_dbeta2*dbetadt**2)*xnpfer \
                - 2.0*(detap_deta*detadt \
                + detap_dbeta*dbetadt)*dxnpferdt \
                - etapos*dxnpferdtt)*y
    #      dsposdt = ((dpposdt + deposdt)*kti - (ppos + epos)/(kt*temp)
    #     1           - (detap_deta*detadt + detap_dbeta*dbetadt)*xnpfer
    #     2           - etapos*dxnpferdt)*y
     
    dsposdta = ((dpposdta + deposdta)*kti \
               - (dpposda + deposda)*kti/t \
               - (detap_deta2*detadt*detada \
               + detap_deta_dbeta*dbetadt*detada \
               + detap_deta*detadta)*xnpfer \
               - (detap_deta*detadt \
               + detap_dbeta*dbetadt)*dxnpferda \
               - detap_deta*detada*dxnpferdt \
               - etapos*dxnpferdta)*y
    dsposdtz = ((dpposdtz + deposdtz)*kti \
               - (dpposdz + deposdz)*kti/t \
               - (detap_deta2*detadt*detadz \
               + detap_deta_dbeta*dbetadt*detadz \
               + detap_deta*detadtz)*xnpfer \
               - (detap_deta*detadt \
               + detap_dbeta*dbetadt)*dxnpferdz \
               - detap_deta*detadz*dxnpferdt \
               - etapos*dxnpferdtz)*y
    dsposdaa = ((dpposdaa + deposdaa)*kti \
               - detap_deta*detadaa*xnpfer \
               - 2.0*detap_deta*detada*dxnpferda \
               - etapos*dxnpferdaa)*y
    dsposdaz = ((dpposdaz + deposdaz)*kti \
               - detap_deta*detadaz*xnpfer \
               - detap_deta*detada*dxnpferdz \
               - detap_deta*detadz*dxnpferda \
               - etapos*dxnpferdaz)*y
    dsposdzz = ((dpposdzz + deposdzz)*kti \
               - detap_deta*detadzz*xnpfer \
               - 2.0*detap_deta*detadz*dxnpferdz \
               - etapos*dxnpferdzz)*y
     
    # and their sum
    sep      = sele + spos
    dsepdd   = dseledd + dsposdd
    dsepdt   = dseledt + dsposdt
    dsepda   = dseleda + dsposda
    dsepdz   = dseledz + dsposdz
    dsepddd  = dseleddd + dsposddd
    dsepddt  = dseleddt + dsposddt
    dsepdda  = dseledda + dsposdda
    dsepddz  = dseleddz + dsposddz
    dsepdtt  = dseledtt + dsposdtt
    dsepdta  = dseledta + dsposdta
    dsepdtz  = dseledtz + dsposdtz
    dsepdaa  = dseledaa + dsposdaa
    dsepdaz  = dseledaz + dsposdaz
    dsepdzz  = dseledzz + dsposdzz
     
    # adjust for the rest mass energy of the positrons
    y                = 2.0 * mecc
    epos             = epos + y * xnpfer
    depos_deta       = depos_deta  + y * dxnpfer_deta
    depos_dbeta      = depos_dbeta + y * dxnpfer_dbeta
    #depos_deta2      = depos_deta2 + y * dxnpferdeta2
    #depos_dbeta2     = 
    #depos_deta_dbeta = 
    deposdd  = deposdd  + y * dxnpferdd
    deposdt  = deposdt  + y * dxnpferdt
    deposda  = deposda  + y * dxnpferda
    deposdz  = deposdz  + y * dxnpferdz
    deposddd = deposddd + y * dxnpferddd
    deposddt = deposddt + y * dxnpferddt
    deposdda = deposdda + y * dxnpferdda
    deposddz = deposddz + y * dxnpferddz
    deposdtt = deposdtt + y * dxnpferdtt
    deposdta = deposdta + y * dxnpferdta
    deposdtz = deposdtz + y * dxnpferdtz
    deposdaa = deposdaa + y * dxnpferdaa
    deposdaz = deposdaz + y * dxnpferdaz
    deposdzz = deposdzz + y * dxnpferdzz
     
    # and resum
    deepdd  = deeledd + deposdd
    deepdt  = deeledt + deposdt
    deepda  = deeleda + deposda
    deepdz  = deeledz + deposdz
    deepddd = deeleddd + deposddd
    deepddt = deeleddt + deposddt
    deepdda = deeledda + deposdda
    deepddz = deeleddz + deposddz
    deepdtt = deeledtt + deposdtt
    deepdta = deeledta + deposdta
    deepdtz = deeledtz + deposdtz
    deepdaa = deeledaa + deposdaa
    deepdaz = deeledaz + deposdaz
    deepdzz = deeledzz + deposdzz
    
    # convert the electron-positron thermal energy in erg/cm**3
    # to a specific thermal energy in erg/gr
    
    eele             = eele*deni
    deele_deta       = deele_deta*deni
    deele_dbeta      = deele_dbeta*deni
    deele_deta2      = deele_deta2*deni
    deele_dbeta2     = deele_dbeta2*deni
    deele_deta_dbeta = deele_deta_dbeta*deni
    deeledd  = (deeledd - eele)*deni
    deeledt  = deeledt*deni
    deeleda  = deeleda*deni
    deeledz  = deeledz*deni
    deeleddd = (deeleddd - 2.0*deeledd)*deni
    deeleddt = (deeleddt - deeledt)*deni
    deeledda = (deeledda - deeleda)*deni
    deeleddz = (deeleddz - deeledz)*deni
    deeledtt = deeledtt*deni
    deeledta = deeledta*deni
    deeledtz = deeledtz*deni
    deeledaa = deeledaa*deni
    deeledaz = deeledaz*deni
    deeledzz = deeledzz*deni
     
    epos         = epos*deni
    depos_deta   = depos_deta*deni
    depos_dbeta  = depos_dbeta*deni
    deposdd  = (deposdd - epos)*deni
    deposdt  = deposdt*deni
    deposda  = deposda*deni
    deposdz  = deposdz*deni
    deposddd = (deposddd - 2.0*deposdd)*deni
    deposddt = (deposddt - deposdt)*deni
    deposdda = (deposdda - deposda)*deni
    deposddz = (deposddz - deposdz)*deni
    deposdtt = deposdtt*deni
    deposdta = deposdta*deni
    deposdtz = deposdtz*deni
    deposdaa = deposdaa*deni
    deposdaz = deposdaz*deni
    deposdzz = deposdzz*deni
     
       # and resum
    deepdd  = deeledd + deposdd
    deepdt  = deeledt + deposdt
    deepda  = deeleda + deposda
    deepdz  = deeledz + deposdz
    deepddd = deeleddd + deposddd
    deepddt = deeleddt + deposddt
    deepdda = deeledda + deposdda
    deepddz = deeleddz + deposddz
    deepdtt = deeledtt + deposdtt
    deepdta = deeledta + deposdta
    deepdtz = deeledtz + deposdtz
    deepdaa = deeledaa + deposdaa
    deepdaz = deeledaz + deposdaz
    deepdzz = deeledzz + deposdzz
    
    ####################### Ionization Potential Terms #########################
    ############################################################################
    
    if (potmult == 0):
        
        eip     = 0.0
        deipdd  = 0.0
        deipdt  = 0.0
        deipda  = 0.0
        deipdz  = 0.0
        deipddd = 0.0
        deipddt = 0.0
        deipdda = 0.0
        deipddz = 0.0
        deipdtt = 0.0
        deipdta = 0.0
        deipdtz = 0.0
        deipdaa = 0.0
        deipdaz = 0.0
        deipdzz = 0.0
     
        sip     = 0.0
        dsipdd  = 0.0
        dsipdt  = 0.0
        dsipda  = 0.0
        dsipdz  = 0.0
        dsipddd = 0.0
        dsipddt = 0.0
        dsipdda = 0.0
        dsipddz = 0.0
        dsipdtt = 0.0
        dsipdta = 0.0
        dsipdtz = 0.0
        dsipdaa = 0.0
        dsipdaz = 0.0
        dsipdzz = 0.0
     
    else:
        eip     = chi * xne
        deipdd  = chi * dxnedd
        deipdt  = chi * dxnedt
        deipda  = chi * dxneda
        deipdz  = chi * dxnedz + hion*ev2erg*xne
        deipddd = chi * dxneddd
        deipddt = chi * dxneddt
        deipdda = chi * dxnedda
        deipddz = chi * dxneddz
        deipdtt = chi * dxnedtt
        deipdta = chi * dxnedta
        deipdtz = chi * dxnedtz
        deipdaa = chi * dxnedaa
        deipdaz = chi * dxnedaz
        deipdzz = chi * dxnedzz + 2.0*hion*ev2erg*dxnedz
     
        # the ionization entropy in erg/gr/kelvin and its derivatives
        y       = kerg*deni
        sip     = eip*kti*y
        dsipdd  = deipdd*kti*y - sip*deni
        dsipdt  = (deipdt*kti - eip*kti/t)*y
        dsipda  = deipda*kti*y
        dsipdz  = deipdz*kti*y
        dsipddd = deipddd*kti*y - dsipdd*deni + sip*deni*deni
        dsipddt = deipddt*kti*y - dsipdt*deni
        dsipddt = deipdda*kti*y - dsipda*deni
        dsipddt = deipddz*kti*y - dsipdz*deni
        dsipdtt = (deipdtt*kti - 2.0*deipdt*kti/t \
                   + 2.0*eip*kti/t**2)*y
        dsipdta = (deipdta*kti - deipda*kti/t)*y
        dsipdtz = (deipdtz*kti - deipdz*kti/t)*y
        dsipdaa = deipdaa*kti*y
        dsipdaz = deipdaz*kti*y
        dsipdzz = deipdzz*kti*y
    
    #      sip    = (eip/kt - etaele*xne) * y
    #       dsipdd = (deipdd/kt
    #    1            - detadd*xne)*y
    #    2            - etaele*dxnedd*y
    #    3            - sip*deni
    #       dsipdt = (deipdt/kt
    #    1             - detadt*xne
    #     2             - etaele*dxnedt
    #     3             - eip/(kt*temp))*y
    
      # convert the ionization energy from erg/cm**3 to  erg/gr
        eip    = eip*deni
        deipdd = (deipdd - eip)*deni
        deipdt = deipdt*deni
        deipda = deipda*deni
        deipdz = deipdz*deni
        deipddd = (deipddd - 2.0*deipdd)*deni
        deipddt = (deipddt - deipdt)*deni
        deipdda = (deipdda - deipda)*deni
        deipddz = (deipddz - deipdz)*deni
        deipdtt = deipdtt*deni
        deipdta = deipdta*deni
        deipdtz = deipdtz*deni
        deipdaa = deipdaa*deni
        deipdaz = deipdaz*deni
        deipdzz = deipdzz*deni
        
    ########################### Creating Totals Terms ##########################
    ############################################################################
    
     
    # bomb proof
    # wish we didn't have to do this if statement
    x   = prad + pion + pele + pcoul
    y   = erad + eion + eele + ecoul
    z   = srad + sion + sele + scoul
    
    if (x < 0.0 or y < 0.0 or z < 0.0):
        
    #         write(6,*)
    #         write(6,*) 'coulomb corrections are causing a negative pressure'
    #         write(6,*) 'setting all coulomb corrections to zero'
    #        write(6,*)
     
        pcoul    = 0.0
        dpcouldd = 0.0
        dpcouldt = 0.0
        dpcoulda = 0.0
        dpcouldz = 0.0
        ecoul    = 0.0
        decouldd = 0.0
        decouldt = 0.0
        decoulda = 0.0
        decouldz = 0.0
        scoul    = 0.0
        dscouldd = 0.0
        dscouldt = 0.0
        dscoulda = 0.0
        dscouldz = 0.0
        decoul_deta = 0.0
        decoul_dbeta = decouldt*mecc/kerg
    
    
    # bomb proof the coulomb corrections
    x   = prad + pion + pele + ppos + pcoul
    y   = erad + eion + eele + epos + ecoul
    z   = srad + sion + sele + spos + scoul
    if (x < 0.0 or y < 0.0 or z < 0.0):
     
    #         write(6,*)
    #         write(6,*) 'coulomb corrections are causing a negative pressure'
    #         write(6,*) 'setting all coulomb corrections to zero'
    #         write(6,*)
     
        pcoul    = 0.0
        dpcouldd = 0.0
        dpcouldt = 0.0
        dpcoulda = 0.0
        dpcouldz = 0.0
        ecoul    = 0.0
        decouldd = 0.0
        decouldt = 0.0
        decoulda = 0.0
        decouldz = 0.0
        scoul    = 0.0
        dscouldd = 0.0
        dscouldt = 0.0
        dscoulda = 0.0
        dscouldz = 0.0
        decoul_deta = 0.0
        decoul_dbeta = decouldt*mecc/kerg
        
        
#    if (NewtonItmode == 0):
#        #print('inside NR')
#        
#        ## Newton Rapshon Guesses
#        ## 
#        xguess = eta_var
#        yguess =  beta_var
#        print('ges0',xguess,yguess)
#        ## Newton-Raphson
#        eostol = 1e-13
#        fpmin  = 1e-14
#        itermax= 1000
#        #print ('%%')
#            
#        for i in range(itermax):
#            print(i)
            # charge neutrality means ne_ionization = ne_electrons - ne_positrons
    f1       = xnefer        - xnpfer        - xne
    # derivative of f with eta and beta for newton-like root finders
    df1deta  = dxnefer_deta  - dxnpfer_deta  - dxne_deta
    df1dbeta = dxnefer_dbeta - dxnpfer_dbeta - dxne_dbeta
    #print('dxne_deta',dxne_deta,dxne_dbeta)
    #print('xne',xne)
     # Second fitting equation for energy with eta and beta derivatives
            
    f2       = (erad + eion + epos + eele + ecoul) - (einput_var)
    #print(den_var)
    df2deta  = (derad_deta + deion_deta + depos_deta + deele_deta + decoul_deta)
    df2dbeta = (derad_dbeta + deion_dbeta + depos_dbeta + deele_dbeta + decoul_dbeta)
            #print(f2)
            ##
            
#            [new_xguess,new_yguess,diff_x,diff_y,x_ratio,y_ratio]= \
#            TwoDNewtonRaphson(xguess,yguess,f1,df1deta,df1dbeta,f2,df2deta,df2dbeta)
#            #print(new_xguess)
#            xguess = new_xguess
#            yguess = new_yguess
#            print('ges',xguess,yguess)
#            
#            if (max(diff_x,diff_y) < eostol) or ((max(abs(x_ratio),abs(y_ratio))) < fpmin):
#                eta_var  = xguess
#                beta_var = yguess
#                #print('xges',xguess,'yges',yguess)
#                break
#            elif (i==itermax-1):
#                print('The Newton-Raphson Method reached the maximum number of iterations')
#            else:
#                continue
#                print('continue') 
    
       
       # sum all the gas components
    pgas    = pion + pele + ppos + pcoul
    egas    = eion + eele + epos + ecoul
    sgas    = sion + sele + spos + scoul
    #First derivs of gas components
    dpgasdd = dpiondd + dpepdd + dpcouldd
    dpgasdt = dpiondt + dpepdt + dpcouldt
    dpgasda = dpionda + dpepda + dpcoulda
    dpgasdz = dpiondz + dpepdz + dpcouldz
     
    degasdd = deiondd + deepdd + decouldd
    degasdt = deiondt + deepdt + decouldt
    degasda = deionda + deepda + decoulda
    degasdz = deiondz + deepdz + decouldz
    
    dsgasdd = dsiondd + dsepdd + dscouldd
    dsgasdt = dsiondt + dsepdt + dscouldt
    dsgasda = dsionda + dsepda + dscoulda
    dsgasdz = dsiondz + dsepdz + dscouldz
    
    # Second derivs of gas components
    dpgasddd = dpionddd + dpepddd + dpcoulddd
    dpgasddt = dpionddt + dpepddt + dpcoulddt
    dpgasdda = dpiondda + dpepdda + dpcouldda
    dpgasddz = dpionddz + dpepddz + dpcoulddz
    dpgasdtt = dpiondtt + dpepdtt + dpcouldtt
    dpgasdta = dpiondta + dpepdta + dpcouldta
    dpgasdtz = dpiondtz + dpepdtz + dpcouldtz
    dpgasdaa = dpiondaa + dpepdaa + dpcouldaa
    dpgasdaz = dpiondaz + dpepdaz + dpcouldaz
    dpgasdzz = dpiondzz + dpepdzz + dpcouldzz
    #
    degasddd = deionddd + deepddd + decoulddd
    degasddt = deionddt + deepddt + decoulddt
    degasdda = deiondda + deepdda + decouldda
    degasddz = deionddz + deepddz + decoulddz
    degasdtt = deiondtt + deepdtt + decouldtt
    degasdta = deiondta + deepdta + decouldta
    degasdtz = deiondtz + deepdtz + decouldtz
    degasdaa = deiondaa + deepdaa + decouldaa
    degasdaz = deiondaz + deepdaz + decouldaz
    degasdzz = deiondzz + deepdzz + decouldzz
    #
    dsgasddd = dsionddd + dsepddd + dscoulddd
    dsgasddt = dsionddt + dsepddt + dscoulddt
    dsgasdda = dsiondda + dsepdda + dscouldda
    dsgasddz = dsionddz + dsepddz + dscoulddz
    dsgasdtt = dsiondtt + dsepdtt + dscouldtt
    dsgasdta = dsiondta + dsepdta + dscouldta
    dsgasdtz = dsiondtz + dsepdtz + dscouldtz
    dsgasdaa = dsiondaa + dsepdaa + dscouldaa
    dsgasdaz = dsiondaz + dsepdaz + dscouldaz
    dsgasdzz = dsiondzz + dsepdzz + dscouldzz
     
     
    # add in radiation to get the total
    pres    = prad + pgas
    ener    = erad + egas
    entr    = srad + sgas
     
    dpresdd = dpraddd + dpgasdd
    dpresdt = dpraddt + dpgasdt
    dpresda = dpradda + dpgasda
    dpresdz = dpraddz + dpgasdz
     
    denerdd = deraddd + degasdd
    denerdt = deraddt + degasdt
    denerda = deradda + degasda
    denerdz = deraddz + degasdz
     
    dentrdd = dsraddd + dsgasdd
    dentrdt = dsraddt + dsgasdt
    dentrda = dsradda + dsgasda
    dentrdz = dsraddz + dsgasdz
    
    #Second Derivatives
    
    dpresddd = dpradddd + dpgasddd
    dpresddt = dpradddt + dpgasddt
    dpresdda = dpraddda + dpgasdda
    dpresddz = dpradddz + dpgasddz 
    dpresdtt = dpraddtt + dpgasdtt
    dpresdta = dpraddta + dpgasdta 
    dpresdtz = dpraddtz + dpgasdtz 
    dpresdaa = dpraddaa + dpgasdaa 
    dpresdaz = dpraddaz + dpgasdaz 
    dpresdzz = dpraddzz + dpgasdzz 
    #
    denerddd = deradddd + degasddd 
    denerddt = deradddt + degasddt 
    denerdda = deraddda + degasdda 
    denerddz = deradddz + degasddz 
    denerdtt = deraddtt + degasdtt 
    denerdta = deraddta + degasdta 
    denerdtz = deraddtz + degasdtz 
    denerdaa = deraddaa + degasdaa 
    denerdaz = deraddaz + degasdaz 
    denerdzz = deraddzz + degasdzz 
    #
    dentrddd = dsradddd + dsgasddd 
    dentrddt = dsradddt + dsgasddt
    dentrdda = dsraddda + dsgasdda 
    dentrddz = dsradddz + dsgasddz 
    dentrdtt = dsraddtt + dsgasdtt 
    dentrdta = dsraddta + dsgasdta 
    dentrdtz = dsraddtz + dsgasdtz 
    dentrdaa = dsraddaa + dsgasdaa 
    dentrdaz = dsraddaz + dsgasdaz 
    dentrdzz = dsraddzz + dsgasdzz
    
    
#    print("%2s %10.6e %2s %10.6e %2s %10.6e %2s %10.6e"% ('temp =',temp_var,'den = ',den_var,'abar=',abar_var,'zbar=',zbar_var))
#    print("%2s %10.6e %2s %10.6e"% ('eta =',eta_var,'beta = ',beta_var))
#    print()
#    print("%2s %16s %16s %16s %16s %16s"% ('teos','value','d/dd','d/dt','d/da','d/dz'))
#    print("%2s %16.8e %16.8e %16.8e %16.8e %16.8e"% ('ptot=',pres,dpresdd,dpresdt,dpresda,dpresdz))
#    print("%2s %16.8e %16.8e %16.8e %16.8e %16.8e"% ('pgas=',pgas,dpgasdd,dpgasdt,dpgasda,dpgasdz))
#    print()
#    print("%2s %16.8e %16.8e %16.8e %16.8e %16.8e"% ('etot=',ener,denerdd,denerdt,denerda,denerdz))
#    print()
#    print("%2s %16.8e %16.8e %16.8e %16.8e %16.8e"% ('stot=',entr,dentrdd,dentrdt,dentrda,dentrdz))
#    print("%2s %16.8e %16.8e %16.8e %16.8e %16.8e"% ('sgas=',sgas,dsgasdd,dsgasdt,dsgasda,dsgasdz))
#    print("%2s %16.8e %16.8e %16.8e %16.8e %16.8e"% ('srad=',srad,dsraddd,dsraddt,dsradda,dsraddz))
#    print("%2s %16.8e %16.8e %16.8e %16.8e %16.8e"% ('sion=',sion,dsiondd,dsiondt,dsionda,dsiondz))
#    print("%2s %16.8e %16.8e %16.8e %16.8e %16.8e"% ('s e-=',sele,dseledd,dseledt,dseleda,dseledz))
#    print("%2s %16.8e %16.8e %16.8e %16.8e %16.8e"% ('s e+=',spos,dsposdd,dsposdt,dsposda,dsposdz))
#    print("%2s %16.8e %16.8e %16.8e %16.8e %16.8e"% ('scoul',scoul,dscouldd,dscouldt,dscoulda,dscouldz))
#    print()
    
    ####################### Entropy Based EOS ###################################
    #############################################################################
    s = entr
    # First Derivatives
    dsde = dentrdt/denerdt
    dsdd = dentrdd - dentrdt*denerdd/denerdt
    dsda = dentrda - dentrdt*denerda/denerdt
    dsdz = dentrdz - dentrdt*denerdz/denerdt
    
    #Second Derivatives
    
    #dsdee = dentrdtt/denerdtt
    dsdee = -dentrdt*denerdtt/denerdt**3.0 + dentrdtt/denerdt**2.0
    #dsdee_check = -1.0/denerdt/temp_var**2.0
    dsded = -dentrdt*(denerddt - denerdtt*denerdd/denerdt)/denerdt**2.0 + (dentrddt - dentrdtt*denerdd/denerdt)/denerdt
    dsdea = -dentrdt*(denerdta - denerdtt*denerda/denerdt)/denerdt**2.0 + (dentrdta - dentrdtt*denerda/denerdt)/denerdt
    dsdez = -dentrdt*(denerdtz - denerdtt*denerdz/denerdt)/denerdt**2.0 + (dentrdtz - dentrdtt*denerdz/denerdt)/denerdt
    dsddd = dentrddd - dentrddt*denerdd/denerdt
    dsdda = dentrdda - dentrdta*denerdd/denerdt
    dsddz = dentrddz - dentrdtz*denerdd/denerdt
    dsdaa = dentrdaa - dentrdta*denerda/denerdt
    dsdaz = dentrdaz - dentrdtz*denerda/denerdt
    dsdzz = dentrdzz - dentrdtz*denerdz/denerdt
    #print('First Energy deriv',dsde - (1/temp_var))
    #print('First density deriv',dsdd + (pres/temp_var/den_var**2.0))
    #print(dsdd,pres,temp_var,den_var)
    
    #print('Second Energy deriv',dsdee_test - dsdee_check)
    
    
    ## Thermodynamic Consistency (tc)
    maxw1 = dsdd + dsde*pres/den_var**2.0
    maxw2 = (dpresdt/denerdt/temp_var - pres/denerdt/temp_var**2.0)/den_var**2.0 - denerdd/denerdt/temp_var**2.0
    #print('maxw1=',maxw1,'maxw2=',maxw2)
    
    return(s,dsde,dsdd,dsda,dsdz,dsdee,dsded,dsdea,dsdez,dsddd,dsdda,dsddz,\
           dsdaa,dsdaz,dsdzz,eta_var,beta_var,f1,df1deta,df1dbeta,f2,df2deta,\
           df2dbeta,maxw1,maxw2)
        
#def EntropicEOS_singlefile_test(beta_var,den_var,eta_var,abar_var,zbar_var,\
#                                einput_var,radmult,ionmult,ionized,elemult,\
#                                coulmult,potmult,1):
#
#    ## inputs: Temp, den, abar, zbar, all modes
#    beta_var = kerg*temp_var/mecc
#    
#    if (NewtonItmode==0):
#    
#        [s,dsde,dsdd,dsda,dsdz,dsdee,dsded,dsdea,dsdez,dsddd,dsdda,dsddz,dsdaa,\
#         dsdaz,dsdzz,eta_out,beta_out]=EntropicEOS_singlefile(beta_var,den_var,\
#                                         eta_var,abar_var,zbar_var,einput_var,\
#                                         radmult,ionmult,ionized,elemult,coulmult,\
#                                         potmult,0.0)
#    
#    else:
#        [s,dsde,dsdd,dsda,dsdz,dsdee,dsded,dsdea,dsdez,dsddd,dsdda,dsddz,dsdaa,\
#         dsdaz,dsdzzeta_out,beta_out]=EntropicEOS_singlefile(beta_var,den_var,\
#                                         eta_var,abar_var,zbar_var,einput_var,\
#                                         radmult,ionmult,ionized,elemult,coulmult,\
#                                         potmult,1)
#    print(s,dsdd,dsda,dsdz)    
  



EntropicEOS_singlefile(beta_var,den_var,eta_var,abar_var,zbar_var,\
                                einput_var,radmult,ionmult,ionized,elemult,\
                                coulmult,potmult)