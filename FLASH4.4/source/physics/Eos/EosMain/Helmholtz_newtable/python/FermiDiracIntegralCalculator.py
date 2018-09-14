#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 15:39:59 2018

@author: alexgrannan
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt

def FermiDiracIntegralCalculator(dkvar,etavar,thetavar):

    # routine dfermi gets the fermi-dirac functions and their derivaties
    # routine fdfunc1 forms the integrand of the fermi-dirac functions
    # routine fdfunc2 same as fdfunc but with the change of variable z**2=x
    # routine dqleg010 does 10 point gauss-legendre integration  9 fig accuracy
    # routine dqleg020 does 20 point gauss-legendre integration 14 fig accuracy
    # routine dqleg040 does 40 point gauss-legendre integration 18 fig accuracy
    # routine dqleg080 does 80 point gauss-legendre integration 32 fig accuracy
    # routine dqlag010 does 10 point gauss-laguerre integration  9 fig accuracy
    # routine dqlag020 does 20 point gauss-laguerre integration 14 fig accuracy
    # routine dqlag040 does 40 point gauss-laguerre integration 18 fig accuracy
    # routine dqlag080 does 80 point gauss-laguerre integration 32 fig accuracy
    
    #   parameters defining the location of the breakpoints for the
    #   subintervals of integration:
    d=  3.3609e0
    sg= 9.1186e-2 
    a1= 6.7774e0 
    b1= 1.1418e0 
    c1= 2.9826e0 
    a2= 3.7601e0 
    b2= 9.3719e-2 
    c2= 2.1063e-2 
    d2= 3.1084e1 
    e2= 1.0056e0 
    a3= 7.5669e0 
    b3= 1.1695e0 
    c3= 7.5416e-1 
    d3= 6.6558e0 
    e3=-1.2819e-1
    
    # Integrand Parameters
    partemp = [dkvar, etavar, thetavar]
    
    # Definition of the x_i:
    eta1=sg*(etavar-d)
    if (eta1 <= 50):
        xi=log(1.0 + exp(eta1))/sg
    else:
        xi=etavar-d
        
    xi2=xi*xi
    
    # definition of the x_i:
    x1=(a1 + b1*xi + c1*xi2)/(1.0 + c1*xi)
    x2=(a2 + b2*xi + c2*d2*xi2)/(1.0 + e2*xi + c2*xi2)
    x3=(a3 + b3*xi + c3*d3*xi2)/(1.0 + e3*xi + c3*xi2)
    
    # breakpoints:
    s1=x1-x2
    s2=x1
    s3=x1+x3
    s12=sqrt(s1)
    #print (dkvar,etavar,thetavar,s1,s2,s3)
    
    def fdfunc1(xtemp,partemp):
        dkvar    = partemp[0]
        etavar   = partemp[1]
        thetavar = partemp[2]
        xdk   = xtemp**dkvar
        xst   = 1.0 + 0.5*xtemp*thetavar
        dxst  = sqrt(xst)
        
        #   avoid overflow in the exponentials at large xtemp
        if ((xtemp-etavar) < 100.0):
           factor      = exp(xtemp-etavar)
           denom       = factor + 1.0
           denomi      = 1.0/denom
           fd          = xdk*dxst*denomi
           fdeta       = fd * factor * denomi
           fdeta2      = (2.0 * factor * denomi - 1.0)*fdeta
           denom2      = 1.0/(4.0 * xst)
           fdtheta     = fd * xtemp * denom2
           fdtheta2    = -fdtheta * xtemp * denom2
           fdetadtheta = fdtheta * factor * denomi
        else:
           factor      = exp(etavar-xtemp)
           fd          = xdk*dxst*factor
           fdeta       = fd
           fdeta2      = fd
           denom2      = 1.0/(4.0 * xst)
           fdtheta     = fd * xtemp * denom2
           fdtheta2    = -fdtheta * xtemp * denom2
           fdetadtheta = fdtheta
        #print(1,xtemp,etavar,thetavar,fd,fdeta,fdtheta)
        #print(1,xtemp,etavar,xtemp-etavar)
        return (fd,fdeta,fdtheta,fdeta2,fdtheta2,fdetadtheta)
    
    def fdfunc2(xtemp,partemp):
        dkvar    = partemp[0]
        etavar   = partemp[1]
        thetavar = partemp[2]
        xsq   = xtemp*xtemp
        xdk   = xtemp**(2.0*dkvar + 1.0)
        xst   = 1.0 + 0.5*xsq*thetavar
        dxst  = sqrt(xst)
        
        #   avoid overflow in the exponentials at large x
        if ((xsq-etavar) < 100.0):
           factor      = exp(xsq - etavar)
           denom       = factor + 1.0
           denomi      = 1.0/denom
           fd          = 2.0*xdk*dxst*denomi
           fdeta       = fd * factor * denomi
           fdeta2      = (2.0 * factor * denomi - 1.0)*fdeta
           denom2      = 1.0/(4.0 * xst)
           fdtheta     = fd * xsq * denom2
           fdtheta2    = -fdtheta * xsq * denom2
           fdetadtheta = fdtheta * factor * denomi
        else:
           factor      = exp(etavar-xsq)
           fd          = 2.0*xdk*dxst*factor
           fdeta       = fd
           fdeta2      = fd
           denom2      = 1.0/(4.0 * xst)
           fdtheta     = fd * xsq * denom2
           fdtheta2    = -fdtheta * xsq * denom2
           fdetadtheta = fdtheta
        #print(2,xtemp,etavar,thetavar,fd,fdeta,fdtheta)
        #print(dkvar,etavar,thetavar,factor,fd)
        return (fd,fdeta,fdtheta,fdeta2,fdtheta2,fdetadtheta)
        
    def dqleg020(f, avar, bvar, partemp):
        xg=np.zeros(10)
        xg[0]  = 7.65265211334973337546404093988382110e-2
        xg[1]  = 2.27785851141645078080496195368574624e-1 
        xg[2]  = 3.73706088715419560672548177024927237e-1
        xg[3]  = 5.10867001950827098004364050955250998e-1 
        xg[4]  = 6.36053680726515025452836696226285936e-1
        xg[5]  = 7.46331906460150792614305070355641590e-1
        xg[6]  = 8.39116971822218823394529061701520685e-1
        xg[7]  = 9.12234428251325905867752441203298113e-1
        xg[8]  = 9.63971927277913791267666131197277221e-1
        xg[9]  = 9.93128599185094924786122388471320278e-1
        
        wg=np.zeros(10)
        wg[0]  = 1.52753387130725850698084331955097593e-1
        wg[1]  = 1.49172986472603746787828737001969436e-1
        wg[2]  = 1.42096109318382051329298325067164933e-1
        wg[3]  = 1.31688638449176626898494499748163134e-1
        wg[4]  = 1.18194531961518417312377377711382287e-1
        wg[5]  = 1.01930119817240435036750135480349876e-1
        wg[6]  = 8.32767415767047487247581432220462061e-2
        wg[7]  = 6.26720483341090635695065351870416063e-2
        wg[8]  = 4.06014298003869413310399522749321098e-2
        wg[9]  = 1.76140071391521183118619623518528163e-2
        
        #           absc   - abscissa
       #       fval*  - function value
       #        result - result of the 20-point gauss formula
        center       = 0.5*(avar+bvar)
        hlfrun       = 0.5*(bvar-avar)
        result       = 0.0
        drdeta       = 0.0
        drdtheta     = 0.0
        drdeta2      = 0.0
        drdtheta2    = 0.0
        drdetadtheta = 0.0
        
        for j in range(0,len(wg)):
            absc1 = center + hlfrun*xg[j]
            absc2 = center - hlfrun*xg[j]
            
            [fval_1,deta_1,dtheta_1,deta2_1,dtheta2_1,detadtheta_1]=f(absc1, partemp)
            [fval_2,deta_2,dtheta_2,deta2_2,dtheta2_2,detadtheta_2]=f(absc2, partemp)
            #print ('dqleg',fval_1,deta_1,dtheta_1)
            result       = result    + (fval_1    + fval_2)            *wg[j]
            drdeta       = drdeta    + (deta_1    + deta_2)            *wg[j]
            drdtheta     = drdtheta  + (dtheta_1  + dtheta_2)          *wg[j]
            #print ('dqleg',dtheta_1,dtheta_2)
            drdeta2      = drdeta2   + (deta2_1   + deta2_2)           *wg[j]
            drdtheta2    = drdtheta2 + (dtheta2_1 + dtheta2_2)         *wg[j]
            drdetadtheta = drdetadtheta + (detadtheta_1 + detadtheta_2)*wg[j]
            #print (j,avar,bvar,absc1,absc2,fval1,fval2)
    
        result       = result * hlfrun
        drdeta       = drdeta * hlfrun
        drdtheta     = drdtheta * hlfrun
        drdeta2      = drdeta2 * hlfrun
        drdtheta2    = drdtheta2 * hlfrun
        drdetadtheta = drdetadtheta * hlfrun
        #print ('dqleg',result,drdeta,drdtheta)
        return (result,drdeta,drdtheta,drdeta2,drdtheta2,drdetadtheta)
    
    def dqlag020(f,avar,bvar, partemp):
        
        xg=np.zeros(20)
        xg[0]  = 7.05398896919887533666890045842150958e-2
        xg[1]  = 3.72126818001611443794241388761146636e-1
        xg[2]  = 9.16582102483273564667716277074183187e-1 
        xg[3]  = 1.70730653102834388068768966741305070e0
        xg[4]  = 2.74919925530943212964503046049481338e0
        xg[5]  = 4.04892531385088692237495336913333219e0
        xg[6]  = 5.61517497086161651410453988565189234e0 
        xg[7]  = 7.45901745367106330976886021837181759e0 
        xg[8]  = 9.59439286958109677247367273428279837e0 
        xg[9]  = 1.20388025469643163096234092988655158e1 
        xg[10] = 1.48142934426307399785126797100479756e1 
        xg[11] = 1.79488955205193760173657909926125096e1 
        xg[12] = 2.14787882402850109757351703695946692e1 
        xg[13] = 2.54517027931869055035186774846415418e1 
        xg[14] = 2.99325546317006120067136561351658232e1 
        xg[15] = 3.50134342404790000062849359066881395e1 
        xg[16] = 4.08330570567285710620295677078075526e1 
        xg[17] = 4.76199940473465021399416271528511211e1 
        xg[18] = 5.58107957500638988907507734444972356e1 
        xg[19] = 6.65244165256157538186403187914606659e1 
        wg=np.zeros(20)
        wg[0]  = 1.81080062418989255451675405913110644e-1 
        wg[1]  = 4.22556767878563974520344172566458197e-1 
        wg[2]  = 6.66909546701848150373482114992515927e-1 
        wg[3]  = 9.15352372783073672670604684771868067e-1 
        wg[4]  = 1.16953970719554597380147822239577476e0 
        wg[5]  = 1.43135498592820598636844994891514331e0 
        wg[6]  = 1.70298113798502272402533261633206720e0
        wg[7]  = 1.98701589079274721410921839275129020e0 
        wg[8]  = 2.28663578125343078546222854681495651e0 
        wg[9]  = 2.60583472755383333269498950954033323e0 
        wg[10] = 2.94978373421395086600235416827285951e0 
        wg[11] = 3.32539578200931955236951937421751118e0 
        wg[12] = 3.74225547058981092111707293265377811e0 
        wg[13] = 4.21423671025188041986808063782478746e0 
        wg[14] = 4.76251846149020929695292197839096371e0 
        wg[15] = 5.42172604424557430380308297989981779e0 
        wg[16] = 6.25401235693242129289518490300707542e0 
        wg[17] = 7.38731438905443455194030019196464791e0 
        wg[18] = 9.15132873098747960794348242552950528e0 
        wg[19] = 1.28933886459399966710262871287485278e1 
        
        #           absc   - abscissa
       #       fval*  - function value
       #        result - result of the 20-point gauss formula
        result       = 0.0
        drdeta       = 0.0
        drdtheta     = 0.0
        drdeta2      = 0.0
        drdtheta2    = 0.0
        drdetadtheta = 0.0
        
        for j in range(0,len(wg)):
            absc = avar + bvar*xg[j]
    
            [fval,deta,dtheta,deta2,dtheta2,detadtheta]=f(absc, partemp)
            #print ('dqlag',fval,deta,dtheta)
            result       = result    + fval*wg[j]
            drdeta       = drdeta    + deta*wg[j]
            drdtheta     = drdtheta  + dtheta*wg[j]
            #print ('dqlag',dtheta)
            drdeta2      = drdeta2   + deta2*wg[j]
            drdtheta2    = drdtheta2 + dtheta2*wg[j]
            drdetadtheta = drdetadtheta + detadtheta*wg[j]
            #print (j,avar,bvar,absc)
            #print (j,avar,bvar,absc,fval)
        result       = result*bvar
        drdeta       = drdeta      *bvar
        drdtheta     = drdtheta    *bvar
        drdeta2      = drdeta2     *bvar
        drdtheta2    = drdtheta2   *bvar
        drdetadtheta = drdetadtheta*bvar
        #print ('dqlag',result,drdeta,drdtheta)
        return (result,drdeta,drdtheta,drdeta2,drdtheta2,drdetadtheta)
    
    [res_1,drde_1,drdt_1,drde2_1,drdt2_1,drdet_1] = dqleg020(fdfunc2,  0.0,  s12, partemp)
    
    [res_2,drde_2,drdt_2,drde2_2,drdt2_2,drdet_2] = dqleg020(fdfunc1,   s1,   s2, partemp)
    
    [res_3,drde_3,drdt_3,drde2_3,drdt2_3,drdet_3] = dqleg020(fdfunc1,   s2,   s3, partemp)
    
    [res_4,drde_4,drdt_4,drde2_4,drdt2_4,drdet_4] = dqlag020(fdfunc1,   s3,  1.0, partemp)
    
    # sum the contributions
    fd          = res_1 + res_2 + res_3 + res_4
    #print (res_1,res_2,res_3,res_4)
    fdeta       = drde_1  + drde_2  + drde_3  + drde_4
    fdtheta     = drdt_1  + drdt_2  + drdt_3  + drdt_4
    fdeta2      = drde2_1 + drde2_2 + drde2_3 + drde2_4
    fdtheta2    = drdt2_1 + drdt2_2 + drdt2_3 + drdt2_4
    fdetadtheta = drdet_1 + drdet_2 + drdet_3 + drdet_4
    #print (dkvar,etavar,thetavar,fd,fdeta,fdtheta)
    return (fd,fdeta,fdtheta,fdeta2,fdtheta2,fdetadtheta)