#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 11:24:19 2018

@author: alexgrannan
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt

## Constants
avo     = 6.0221417930e23
kerg    = 1.3806504240000000e-16
me      = 9.1093821545e-28
clight  = 2.99792458e10
mecc    = me * clight * clight
h       = 6.6260689633e-27
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

def Entropic_EOS(beta_var,den_var,eta_var,abar_var,zbar_var,einput_var,\
                           radmult,ionmult,ionized,elemult,coulmult,potmult,\
                           mode):
    
    if (mode==0):
        
        ## Determines the best eta and beta
       
    
        ## Newton Rapshon Guesses
        
        xguess = eta_var
        yguess = beta_var
#        print('ges0',xguess,yguess)
        ## Newton-Raphson
        eostol = 1e-13
        fpmin  = 1e-14
        itermax=100
        #print ('%%')
#        
        for i in range(itermax):
#            print(i)
            
            [s,dsde,dsdd,dsda,dsdz,dsdee,dsded,dsdea,dsdez,dsddd,dsdda,dsddz,\
             dsdaa,dsdaz,dsdzz,eta_out,beta_out,f1,df1deta,df1dbeta,f2,df2deta,\
             df2dbeta,maxw1,maxw2]=EntropicEOS_singlefile(yguess,den_var,xguess,\
                            abar_var,zbar_var,einput_var,\
                            radmult,ionmult,ionized,elemult,coulmult,potmult)
            
            [new_xguess,new_yguess,diff_x,diff_y,x_ratio,y_ratio]= \
            TwoDNewtonRaphson(xguess,yguess,f1,df1deta,df1dbeta,f2,df2deta,df2dbeta)
            
            xguess = new_xguess
            yguess = new_yguess
#            print('ges',xguess,yguess)
            
            if (max(diff_x,diff_y) < eostol) or ((max(abs(x_ratio),abs(y_ratio))) < fpmin):
                mode=1.0
                eta_var  = xguess
                beta_var = yguess
#                print('gesfin',eta_var,beta_var)
                break
            elif (i==itermax-1):
                print('The Newton-Raphson Method reached the maximum number of iterations')
            else:
                continue
                
                
                
    [s,dsde,dsdd,dsda,dsdz,dsdee,dsded,dsdea,dsdez,dsddd,dsdda,dsddz,\
         dsdaa,dsdaz,dsdzz,eta_out,beta_out,f1,df1deta,df1dbeta,f2,df2deta,\
         df2dbeta,maxw1,maxw2]=EntropicEOS_singlefile(beta_var,den_var,eta_var,\
                            abar_var,zbar_var,einput_var,\
                            radmult,ionmult,ionized,elemult,coulmult,potmult)
#            
#
#            #print (i)
#            #print ('A',xguess,yguess)
#            
#            ## Make the root finding functions
#        
#            # Electron-Positron terms
#            [xnefer,dxnefer_deta,dxnefer_dbeta,\
#            xnpfer,dxnpfer_deta,dxnpfer_dbeta,\
#            xne,dxne_deta,dxne_dbeta,\
#            epos,depos_deta,depos_dbeta,spos,dsposdd,dsposddd,deposdt,depos_deta2,\
#            eele,deele_deta,deele_dbeta,sele,dseledd,dseleddd,deeledt,deele_deta2,\
#            dsded,eip,sip,detadd,detaddd]=\
#                    EOS_ele_pos_ionizepotential(yguess,xguess,abar_var,zbar_var,den_var,ionized,potmult)
#                    
#            # Ion and Coulomb Correction terms
#            [xni,eion,deion_deta,deion_deta2,deion_dbeta, sion, dsiondd, deiondt, dsionddd,pcoul,ecoul,scoul] = EOS_ion_coul_corr(yguess,xguess,abar_var,zbar_var,den_var,ionmult,coulmult)
#        
#            # First fitting equation for number densities with eta and beta derivatives
#              
#            # charge neutrality means ne_ionization = ne_electrons - ne_positrons
#            f1       = xnefer        - xnpfer        - xne
#            
#            df1deta  = dxnefer_deta  - dxnpfer_deta  - dxne_deta
#            df1dbeta = dxnefer_dbeta - dxnpfer_dbeta - dxne_dbeta
#            
#            # Second fitting equation for energy with eta and beta derivatives
#            
#            f2       = (eion + epos + eele) - einput_var
#            #print(den_var)
#            df2deta  = (deion_deta + depos_deta + deele_deta)
#            df2dbeta = (deion_dbeta + depos_dbeta + deele_dbeta)
#            
#            ##
#        
#            [new_xguess,new_yguess,diff_x,diff_y,x_ratio,y_ratio]= \
#            TwoDNewtonRaphson(xguess,yguess,f1,df1deta,df1dbeta,f2,df2deta,df2dbeta)
#            
#            xguess = new_xguess
#            yguess = new_yguess
#            
#            if (max(diff_x,diff_y) < eostol) or ((max(abs(x_ratio),abs(y_ratio))) < fpmin):
#                mode=1.0
#                eta_var  = xguess
#                beta_var = yguess
#                break
#            elif (i==itermax):
#                print('The Newton-Raphson Method reached the maximum number of iterations')
#            else:
#                continue
#                print('continue') 
    #print (xguess,yguess)   
    # Electron-Positron terms
#    [xnefer,dxnefer_deta,dxnefer_dbeta,\
#            xnpfer,dxnpfer_deta,dxnpfer_dbeta,\
#            xne,dxne_deta,dxne_dbeta,\
#            epos,depos_deta,depos_dbeta,spos,dsposdd,dsposddd,deposdt,depos_deta2,\
#            eele,deele_deta,deele_dbeta,sele,dseledd,dseleddd,deeledt,deele_deta2,\
#            dsded,eip,sip,detadd,detaddd]=\
#            EOS_ele_pos_ionizepotential(beta_var,eta_var,abar_var,zbar_var,den_var,ionized,potmult)
#          
#    # Ion and Coulomb Correction terms
#    [xni,eion,deion_deta,deion_deta2,deion_dbeta, sion, dsiondd, deiondt, dsionddd,\
#           pcoul,ecoul,scoul] = EOS_ion_coul_corr(beta_var,eta_var,abar_var,zbar_var,den_var,ionmult,coulmult)
#    #print('e_in',eion+epos+eele)
#    ## Radiation Term  
#    [prad,erad,srad] = EOS_rad(beta_var,eta_var,abar_var,zbar_var,den_var,1.0)
#    #print ('B',eta_var,beta_var)
##    print('drhodbeta',(dxnefer_dbeta - dxnpfer_dbeta)/avo)
##    print('drhodeta',(dxnefer_deta - dxnpfer_deta)/avo)
##    print('dedeta',df2deta)
##    print('dedbeta',df2dbeta)
#    #print(xne,xnefer,xnpfer)
    return (s,dsde,dsdd,dsda,dsdz,dsdee,dsded,dsdea,dsdez,dsddd,dsdda,dsddz,\
            dsdaa,dsdaz,dsdzz,eta_var,beta_var,maxw1,maxw2)

#    # Total Entropy
#    s3     = sion + sele + spos
#    # First deriviative of total entropy with density
#    d3dd  = dsiondd + dseledd + dsposdd
#    # Second deriviative of total entropy with density
#    d3ddd = dsionddd + dseleddd + dsposddd
#    # First derivative of entropy with energy
#    d3de  = temp
#    # Second drivative of entropy with energy
#    d3dee = -1/((deiondt + deeledt + deposdt)*temp**2)
#    # Mixed derivative of entropy with temperature and energy
#    d3ded = dsded
#    
#    # Chemical Potential
#    eta_totable  = etaele
#    # First Derivative of Chemical Potential with density
#    detadd_totable = detadd
#    # Second Derivative of Chemical Potential with density
#    detaddd_totable = detaddd
#    # First Derivative of Chemical Potential with energy
#    detade_totable = 1/(deion_deta + deele_deta + depos_deta)
#    # Second Derivative of Chemical Potential with energy
#    detadee_totable = 1/(deion_deta2 + deele_deta2 + depos_deta2)
#    # Mixed Derivative of Chemical Potential
#    detaded_totable = 0.0
#    
#    return(s3,d3dd,d3ddd,d3de,d3dee,d3ded,eta_totable,detadd_totable,\
#           detaddd_totable,detade_totable,detadee_totable,detaded_totable)


def Entropic_EOS_test():
      
    ## Inputs
    temp_0    = 1.0000E3
    beta_var      = kerg * temp_0 / mecc
    eta_var  = -21.653832766028515  
    abar_var     = 1.0
    zbar_var      = 0.1*1.5
    den_var     = 1.0000E-12
    einput_var  = 144754322608.70731
    # set the on/off switches
    radmult  = 1
    ionmult  = 1
    ionized  = 1
    elemult  = 1
    coulmult = 1
    potmult  = 0
    mode = 0.0
    # set the on/off switches
#    radmult  = 1 # Turn radiation on (1) or off (0)
#    ionmult  = 1 # Turn ion contributionon (1) or off (0)
#    ionized  = 1 # Turn ionization with Saha on (0)** or off (1)
#    elemult  = 1 # Turn electron/positron contribution on (1) or off (0)
#    coulmult = 1 # Turn coulomb corrections on (1) or off (0)
#    potmult  = 0 # Turn potential ionization on (1) or off (0)
#    mode     = 0 # Turn Newton Raphson on (0) or off (1)

    
    [s,dsde,dsdd,dsda,dsdz,dsdee,dsded,dsdea,dsdez,dsddd,dsdda,dsddz,\
            dsdaa,dsdaz,dsdzz,eta_out,beta_out,maxw1,maxw2]=\
                Entropic_EOS(beta_var,den_var,eta_var,abar_var,zbar_var,einput_var,\
                           radmult,ionmult,ionized,elemult,coulmult,potmult,mode)
        
        
    print ('beta',beta_out,'etaele',eta_out)  

Entropic_EOS_test()
