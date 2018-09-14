#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 13:19:09 2018

@author: alexgrannan
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import colors, ticker, cm
from numpy.linalg import inv
from scipy import misc, optimize
from scipy.optimize import minimize, rosen
import scipy.linalg
from mpl_toolkits.mplot3d import Axes3D

## Construct 4D log liklihood functions

################################### test function ############################
##############################################################################
    


def TestFuncandDerivs(w_data,x_data,y_data,z_data):
        
        f = np.zeros(len(x_data))
        dfdw = np.zeros(len(x_data))
        dfdx = np.zeros(len(x_data))
        dfdy = np.zeros(len(x_data))
        dfdz = np.zeros(len(x_data))
        dfdww = np.zeros(len(x_data))
        dfdwx = np.zeros(len(x_data))
        dfdwy = np.zeros(len(x_data))
        dfdwz = np.zeros(len(x_data))
        dfdxx = np.zeros(len(x_data))
        dfdxy = np.zeros(len(x_data))
        dfdxz = np.zeros(len(x_data))
        dfdyy = np.zeros(len(x_data))
        dfdyz = np.zeros(len(x_data))
        dfdzz = np.zeros(len(x_data))
        
        for i in range(0,len(x_data)):
            f_coef=w_data[i]**2.0  + x_data[i]**2.0 + y_data[i]**2.0  + z_data[i]**2.0
            f[i] = sqrt(f_coef)
            dfdw[i] = w_data[i]/sqrt(f_coef)
            dfdx[i] = x_data[i]/sqrt(f_coef)
            dfdy[i] = y_data[i]/sqrt(f_coef)
            dfdz[i] = z_data[i]/sqrt(f_coef)
            dfdww[i] = -w_data[i]**2.0/f_coef**1.5 + 1.0/sqrt(f_coef)
            dfdwx[i] = -w_data[i]*x_data[i]/f_coef**1.5
            dfdwy[i] = -w_data[i]*y_data[i]/f_coef**1.5
            dfdwz[i] = -w_data[i]*z_data[i]/f_coef**1.5
            dfdxx[i] = -x_data[i]**2.0/f_coef**1.5 + 1.0/sqrt(f_coef)
            dfdxy[i] = -x_data[i]*y_data[i]/f_coef**1.5
            dfdxz[i] = -x_data[i]*z_data[i]/f_coef**1.5
            dfdyy[i] = -y_data[i]**2.0/f_coef**1.5 + 1.0/sqrt(f_coef)
            dfdyz[i] = -y_data[i]*z_data[i]/f_coef**1.5
            dfdzz[i] = -z_data[i]**2.0/f_coef**1.5 + 1.0/sqrt(f_coef)
        
        f_df_data = np.concatenate((f,dfdw,dfdx,dfdy,dfdz,dfdww,dfdwx,dfdwy,dfdwz,\
                                    dfdxx,dfdxy,dfdxz,dfdyy,dfdyz,dfdzz))
        #print('f_df_data',f_df_data)
        return(f_df_data)
############################################################################
############################################################################
        


############################ Function #######################################


def EOS_MLL_Solver(w_data,x_data,y_data,z_data,w_star,x_star,y_star,z_star,f_df_data,gph_sigma_f,sigma_n,hyperparams_in):
    
    fold = f_df_data

    def BuildCovarMatrix(hyperparams):
      
   
        gph_l1 = hyperparams[0]
        gph_l2 = hyperparams[1]
        gph_l3 = hyperparams[2]
        gph_l4 = hyperparams[3]
        
        k1_1 = np.zeros((dim_x_data,dim_x_data))
        k1_2 = np.zeros((dim_x_data,dim_x_data))
        k1_3 = np.zeros((dim_x_data,dim_x_data))
        k1_4 = np.zeros((dim_x_data,dim_x_data))
        k1_5 = np.zeros((dim_x_data,dim_x_data))
        k1_6 = np.zeros((dim_x_data,dim_x_data))
        k1_7 = np.zeros((dim_x_data,dim_x_data))
        k1_8 = np.zeros((dim_x_data,dim_x_data))
        k1_9 = np.zeros((dim_x_data,dim_x_data))
        k1_10 = np.zeros((dim_x_data,dim_x_data))
        k1_11 = np.zeros((dim_x_data,dim_x_data))
        k1_12 = np.zeros((dim_x_data,dim_x_data))
        k1_13 = np.zeros((dim_x_data,dim_x_data))
        k1_14 = np.zeros((dim_x_data,dim_x_data))
        k1_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k2_1 = np.zeros((dim_x_data,dim_x_data))
        k2_2 = np.zeros((dim_x_data,dim_x_data))
        k2_3 = np.zeros((dim_x_data,dim_x_data))
        k2_4 = np.zeros((dim_x_data,dim_x_data))
        k2_5 = np.zeros((dim_x_data,dim_x_data))
        k2_6 = np.zeros((dim_x_data,dim_x_data))
        k2_7 = np.zeros((dim_x_data,dim_x_data))
        k2_8 = np.zeros((dim_x_data,dim_x_data))
        k2_9 = np.zeros((dim_x_data,dim_x_data))
        k2_10 = np.zeros((dim_x_data,dim_x_data))
        k2_11 = np.zeros((dim_x_data,dim_x_data))
        k2_12 = np.zeros((dim_x_data,dim_x_data))
        k2_13 = np.zeros((dim_x_data,dim_x_data))
        k2_14 = np.zeros((dim_x_data,dim_x_data))
        k2_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k3_1 = np.zeros((dim_x_data,dim_x_data))
        k3_2 = np.zeros((dim_x_data,dim_x_data))
        k3_3 = np.zeros((dim_x_data,dim_x_data))
        k3_4 = np.zeros((dim_x_data,dim_x_data))
        k3_5 = np.zeros((dim_x_data,dim_x_data))
        k3_6 = np.zeros((dim_x_data,dim_x_data))
        k3_7 = np.zeros((dim_x_data,dim_x_data))
        k3_8 = np.zeros((dim_x_data,dim_x_data))
        k3_9 = np.zeros((dim_x_data,dim_x_data))
        k3_10 = np.zeros((dim_x_data,dim_x_data))
        k3_11 = np.zeros((dim_x_data,dim_x_data))
        k3_12 = np.zeros((dim_x_data,dim_x_data))
        k3_13 = np.zeros((dim_x_data,dim_x_data))
        k3_14 = np.zeros((dim_x_data,dim_x_data))
        k3_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k4_1 = np.zeros((dim_x_data,dim_x_data))
        k4_2 = np.zeros((dim_x_data,dim_x_data))
        k4_3 = np.zeros((dim_x_data,dim_x_data))
        k4_4 = np.zeros((dim_x_data,dim_x_data))
        k4_5 = np.zeros((dim_x_data,dim_x_data))
        k4_6 = np.zeros((dim_x_data,dim_x_data))
        k4_7 = np.zeros((dim_x_data,dim_x_data))
        k4_8 = np.zeros((dim_x_data,dim_x_data))
        k4_9 = np.zeros((dim_x_data,dim_x_data))
        k4_10 = np.zeros((dim_x_data,dim_x_data))
        k4_11 = np.zeros((dim_x_data,dim_x_data))
        k4_12 = np.zeros((dim_x_data,dim_x_data))
        k4_13 = np.zeros((dim_x_data,dim_x_data))
        k4_14 = np.zeros((dim_x_data,dim_x_data))
        k4_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k5_1 = np.zeros((dim_x_data,dim_x_data))
        k5_2 = np.zeros((dim_x_data,dim_x_data))
        k5_3 = np.zeros((dim_x_data,dim_x_data))
        k5_4 = np.zeros((dim_x_data,dim_x_data))
        k5_5 = np.zeros((dim_x_data,dim_x_data))
        k5_6 = np.zeros((dim_x_data,dim_x_data))
        k5_7 = np.zeros((dim_x_data,dim_x_data))
        k5_8 = np.zeros((dim_x_data,dim_x_data))
        k5_9 = np.zeros((dim_x_data,dim_x_data))
        k5_10 = np.zeros((dim_x_data,dim_x_data))
        k5_11 = np.zeros((dim_x_data,dim_x_data))
        k5_12 = np.zeros((dim_x_data,dim_x_data))
        k5_13 = np.zeros((dim_x_data,dim_x_data))
        k5_14 = np.zeros((dim_x_data,dim_x_data))
        k5_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k6_1 = np.zeros((dim_x_data,dim_x_data))
        k6_2 = np.zeros((dim_x_data,dim_x_data))
        k6_3 = np.zeros((dim_x_data,dim_x_data))
        k6_4 = np.zeros((dim_x_data,dim_x_data))
        k6_5 = np.zeros((dim_x_data,dim_x_data))
        k6_6 = np.zeros((dim_x_data,dim_x_data))
        k6_7 = np.zeros((dim_x_data,dim_x_data))
        k6_8 = np.zeros((dim_x_data,dim_x_data))
        k6_9 = np.zeros((dim_x_data,dim_x_data))
        k6_10 = np.zeros((dim_x_data,dim_x_data))
        k6_11 = np.zeros((dim_x_data,dim_x_data))
        k6_12 = np.zeros((dim_x_data,dim_x_data))
        k6_13 = np.zeros((dim_x_data,dim_x_data))
        k6_14 = np.zeros((dim_x_data,dim_x_data))
        k6_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k7_1 = np.zeros((dim_x_data,dim_x_data))
        k7_2 = np.zeros((dim_x_data,dim_x_data))
        k7_3 = np.zeros((dim_x_data,dim_x_data))
        k7_4 = np.zeros((dim_x_data,dim_x_data))
        k7_5 = np.zeros((dim_x_data,dim_x_data))
        k7_6 = np.zeros((dim_x_data,dim_x_data))
        k7_7 = np.zeros((dim_x_data,dim_x_data))
        k7_8 = np.zeros((dim_x_data,dim_x_data))
        k7_9 = np.zeros((dim_x_data,dim_x_data))
        k7_10 = np.zeros((dim_x_data,dim_x_data))
        k7_11 = np.zeros((dim_x_data,dim_x_data))
        k7_12 = np.zeros((dim_x_data,dim_x_data))
        k7_13 = np.zeros((dim_x_data,dim_x_data))
        k7_14 = np.zeros((dim_x_data,dim_x_data))
        k7_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k8_1 = np.zeros((dim_x_data,dim_x_data))
        k8_2 = np.zeros((dim_x_data,dim_x_data))
        k8_3 = np.zeros((dim_x_data,dim_x_data))
        k8_4 = np.zeros((dim_x_data,dim_x_data))
        k8_5 = np.zeros((dim_x_data,dim_x_data))
        k8_6 = np.zeros((dim_x_data,dim_x_data))
        k8_7 = np.zeros((dim_x_data,dim_x_data))
        k8_8 = np.zeros((dim_x_data,dim_x_data))
        k8_9 = np.zeros((dim_x_data,dim_x_data))
        k8_10 = np.zeros((dim_x_data,dim_x_data))
        k8_11 = np.zeros((dim_x_data,dim_x_data))
        k8_12 = np.zeros((dim_x_data,dim_x_data))
        k8_13 = np.zeros((dim_x_data,dim_x_data))
        k8_14 = np.zeros((dim_x_data,dim_x_data))
        k8_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k9_1 = np.zeros((dim_x_data,dim_x_data))
        k9_2 = np.zeros((dim_x_data,dim_x_data))
        k9_3 = np.zeros((dim_x_data,dim_x_data))
        k9_4 = np.zeros((dim_x_data,dim_x_data))
        k9_5 = np.zeros((dim_x_data,dim_x_data))
        k9_6 = np.zeros((dim_x_data,dim_x_data))
        k9_7 = np.zeros((dim_x_data,dim_x_data))
        k9_8 = np.zeros((dim_x_data,dim_x_data))
        k9_9 = np.zeros((dim_x_data,dim_x_data))
        k9_10 = np.zeros((dim_x_data,dim_x_data))
        k9_11 = np.zeros((dim_x_data,dim_x_data))
        k9_12 = np.zeros((dim_x_data,dim_x_data))
        k9_13 = np.zeros((dim_x_data,dim_x_data))
        k9_14 = np.zeros((dim_x_data,dim_x_data))
        k9_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k10_1 = np.zeros((dim_x_data,dim_x_data))
        k10_2 = np.zeros((dim_x_data,dim_x_data))
        k10_3 = np.zeros((dim_x_data,dim_x_data))
        k10_4 = np.zeros((dim_x_data,dim_x_data))
        k10_5 = np.zeros((dim_x_data,dim_x_data))
        k10_6 = np.zeros((dim_x_data,dim_x_data))
        k10_7 = np.zeros((dim_x_data,dim_x_data))
        k10_8 = np.zeros((dim_x_data,dim_x_data))
        k10_9 = np.zeros((dim_x_data,dim_x_data))
        k10_10 = np.zeros((dim_x_data,dim_x_data))
        k10_11 = np.zeros((dim_x_data,dim_x_data))
        k10_12 = np.zeros((dim_x_data,dim_x_data))
        k10_13 = np.zeros((dim_x_data,dim_x_data))
        k10_14 = np.zeros((dim_x_data,dim_x_data))
        k10_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k11_1 = np.zeros((dim_x_data,dim_x_data))
        k11_2 = np.zeros((dim_x_data,dim_x_data))
        k11_3 = np.zeros((dim_x_data,dim_x_data))
        k11_4 = np.zeros((dim_x_data,dim_x_data))
        k11_5 = np.zeros((dim_x_data,dim_x_data))
        k11_6 = np.zeros((dim_x_data,dim_x_data))
        k11_7 = np.zeros((dim_x_data,dim_x_data))
        k11_8 = np.zeros((dim_x_data,dim_x_data))
        k11_9 = np.zeros((dim_x_data,dim_x_data))
        k11_10 = np.zeros((dim_x_data,dim_x_data))
        k11_11 = np.zeros((dim_x_data,dim_x_data))
        k11_12 = np.zeros((dim_x_data,dim_x_data))
        k11_13 = np.zeros((dim_x_data,dim_x_data))
        k11_14 = np.zeros((dim_x_data,dim_x_data))
        k11_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k12_1 = np.zeros((dim_x_data,dim_x_data))
        k12_2 = np.zeros((dim_x_data,dim_x_data))
        k12_3 = np.zeros((dim_x_data,dim_x_data))
        k12_4 = np.zeros((dim_x_data,dim_x_data))
        k12_5 = np.zeros((dim_x_data,dim_x_data))
        k12_6 = np.zeros((dim_x_data,dim_x_data))
        k12_7 = np.zeros((dim_x_data,dim_x_data))
        k12_8 = np.zeros((dim_x_data,dim_x_data))
        k12_9 = np.zeros((dim_x_data,dim_x_data))
        k12_10 = np.zeros((dim_x_data,dim_x_data))
        k12_11 = np.zeros((dim_x_data,dim_x_data))
        k12_12 = np.zeros((dim_x_data,dim_x_data))
        k12_13 = np.zeros((dim_x_data,dim_x_data))
        k12_14 = np.zeros((dim_x_data,dim_x_data))
        k12_15 = np.zeros((dim_x_data,dim_x_data))
         #
        k13_1 = np.zeros((dim_x_data,dim_x_data))
        k13_2 = np.zeros((dim_x_data,dim_x_data))
        k13_3 = np.zeros((dim_x_data,dim_x_data))
        k13_4 = np.zeros((dim_x_data,dim_x_data))
        k13_5 = np.zeros((dim_x_data,dim_x_data))
        k13_6 = np.zeros((dim_x_data,dim_x_data))
        k13_7 = np.zeros((dim_x_data,dim_x_data))
        k13_8 = np.zeros((dim_x_data,dim_x_data))
        k13_9 = np.zeros((dim_x_data,dim_x_data))
        k13_10 = np.zeros((dim_x_data,dim_x_data))
        k13_11 = np.zeros((dim_x_data,dim_x_data))
        k13_12 = np.zeros((dim_x_data,dim_x_data))
        k13_13 = np.zeros((dim_x_data,dim_x_data))
        k13_14 = np.zeros((dim_x_data,dim_x_data))
        k13_15 = np.zeros((dim_x_data,dim_x_data))
         #
        k14_1 = np.zeros((dim_x_data,dim_x_data))
        k14_2 = np.zeros((dim_x_data,dim_x_data))
        k14_3 = np.zeros((dim_x_data,dim_x_data))
        k14_4 = np.zeros((dim_x_data,dim_x_data))
        k14_5 = np.zeros((dim_x_data,dim_x_data))
        k14_6 = np.zeros((dim_x_data,dim_x_data))
        k14_7 = np.zeros((dim_x_data,dim_x_data))
        k14_8 = np.zeros((dim_x_data,dim_x_data))
        k14_9 = np.zeros((dim_x_data,dim_x_data))
        k14_10 = np.zeros((dim_x_data,dim_x_data))
        k14_11 = np.zeros((dim_x_data,dim_x_data))
        k14_12 = np.zeros((dim_x_data,dim_x_data))
        k14_13 = np.zeros((dim_x_data,dim_x_data))
        k14_14 = np.zeros((dim_x_data,dim_x_data))
        k14_15 = np.zeros((dim_x_data,dim_x_data))
         #
        k15_1 = np.zeros((dim_x_data,dim_x_data))
        k15_2 = np.zeros((dim_x_data,dim_x_data))
        k15_3 = np.zeros((dim_x_data,dim_x_data))
        k15_4 = np.zeros((dim_x_data,dim_x_data))
        k15_5 = np.zeros((dim_x_data,dim_x_data))
        k15_6 = np.zeros((dim_x_data,dim_x_data))
        k15_7 = np.zeros((dim_x_data,dim_x_data))
        k15_8 = np.zeros((dim_x_data,dim_x_data))
        k15_9 = np.zeros((dim_x_data,dim_x_data))
        k15_10 = np.zeros((dim_x_data,dim_x_data))
        k15_11 = np.zeros((dim_x_data,dim_x_data))
        k15_12 = np.zeros((dim_x_data,dim_x_data))
        k15_13 = np.zeros((dim_x_data,dim_x_data))
        k15_14 = np.zeros((dim_x_data,dim_x_data))
        k15_15 = np.zeros((dim_x_data,dim_x_data))
        
        for i in range(dim_x_data):
            for j in range(dim_x_data):
                if (i==j):
                    sigma_tmp = 0.0
                else:
                    sigma_tmp = 0.0
                
                dw = w_data[i] - w_data[j]
                dx = x_data[i] - x_data[j]
                dy = y_data[i] - y_data[j]
                dz = z_data[i] - z_data[j]
                coef = gph_sigma_f * gph_sigma_f * exp( - 0.5 * ( dw**2.0/gph_l1**2.0 + dx**2.0/gph_l2**2.0 + dy**2.0/gph_l3**2.0 + dz**2.0/gph_l4**2.0 ) )
                c1 = dw / gph_l1**2.0 
                c2 = dx / gph_l2**2.0
                c3 = dy / gph_l3**2.0
                c4 = dz / gph_l4**2.0
                #print('c1=',c1,'c2=',c2,'c3=',c3,'c4=',c4)
                # Effectively creates an upper triangular matrix
                k1_1[i,j] = coef + sigma_tmp
                k1_2[i,j] = c1*k1_1[i,j] + sigma_tmp
                k1_3[i,j] = c2*k1_1[i,j] + sigma_tmp
                k1_4[i,j] = c3*k1_1[i,j] + sigma_tmp
                k1_5[i,j] = c4*k1_1[i,j] + sigma_tmp
                k1_6[i,j] = (c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k1_7[i,j] = c1*c2*k1_1[i,j] + sigma_tmp
                k1_8[i,j] = c1*c3*k1_1[i,j] + sigma_tmp
                k1_9[i,j] = c1*c4*k1_1[i,j] + sigma_tmp
                k1_10[i,j] = (c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k1_11[i,j] = c2*c3*k1_1[i,j] + sigma_tmp
                k1_12[i,j] = c2*c4*k1_1[i,j] + sigma_tmp
                k1_13[i,j] = (c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k1_14[i,j] = c3*c4*k1_1[i,j] + sigma_tmp
                k1_15[i,j] = (c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k2_1[i,j] = -c1*k1_1[i,j] + sigma_tmp
                k2_2[i,j] = (1.0/gph_l1**2.0 - c1**2.0)*k1_1[i,j] + sigma_tmp
                k2_3[i,j] = -c1*c2*k1_1[i,j] + sigma_tmp
                k2_4[i,j] = -c1*c3*k1_1[i,j] + sigma_tmp
                k2_5[i,j] = -c1*c4*k1_1[i,j] + sigma_tmp
                k2_6[i,j] = -c1*(c1**2.0 - 3.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k2_7[i,j] = -c2*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k2_8[i,j] = -c3*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k2_9[i,j] = -c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k2_10[i,j]= -c1*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k2_11[i,j]= -c1*c2*c3*k1_1[i,j] + sigma_tmp
                k2_12[i,j]= -c1*c2*c4*k1_1[i,j] + sigma_tmp
                k2_13[i,j]= -c1*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k2_14[i,j]= -c1*c3*c4*k1_1[i,j] + sigma_tmp
                k2_15[i,j]= -c1*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k3_1[i,j] = -c2*k1_1[i,j] + sigma_tmp
                k3_2[i,j] = -c1*c2*k1_1[i,j] + sigma_tmp
                k3_3[i,j] = (1.0/gph_l2**2.0 - c2**2.0)*k1_1[i,j] + sigma_tmp
                k3_4[i,j] = -c2*c3*k1_1[i,j] + sigma_tmp
                k3_5[i,j] = -c2*c4*k1_1[i,j] + sigma_tmp
                k3_6[i,j] = -c2*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k3_7[i,j] = -c1*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k3_8[i,j] = -c1*c2*c3*k1_1[i,j] + sigma_tmp
                k3_9[i,j] = -c1*c2*c4*k1_1[i,j] + sigma_tmp
                k3_10[i,j] =-c2*(c2**2.0 - 3.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k3_11[i,j] =-c3*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k3_12[i,j] =-c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k3_13[i,j] =-c2*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k3_14[i,j] =-c2*c3*c4*k1_1[i,j] + sigma_tmp
                k3_15[i,j] =-c2*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k4_1[i,j] = -c3*k1_1[i,j] + sigma_tmp
                k4_2[i,j] = -c1*c3*k1_1[i,j] + sigma_tmp
                k4_3[i,j] = -c2*c3*k1_1[i,j] + sigma_tmp
                k4_4[i,j] = (1.0/gph_l3**2.0 - c3**2.0)*k1_1[i,j] + sigma_tmp
                k4_5[i,j] = -c3*c4*k1_1[i,j] + sigma_tmp
                k4_6[i,j] = -c3*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k4_7[i,j] = -c1*c2*c3*k1_1[i,j] + sigma_tmp
                k4_8[i,j] = -c1*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k4_9[i,j] = -c1*c3*c4*k1_1[i,j] + sigma_tmp
                k4_10[i,j]= -c3*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k4_11[i,j]= -c2*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k4_12[i,j]= -c2*c3*c4*k1_1[i,j] + sigma_tmp
                k4_13[i,j]= -c3*(c3**2.0 - 3.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k4_14[i,j]= -c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k4_15[i,j]= -c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k5_1[i,j] = -c4*k1_1[i,j] + sigma_tmp
                k5_2[i,j] = -c1*c4*k1_1[i,j] + sigma_tmp
                k5_3[i,j] = -c2*c4*k1_1[i,j] + sigma_tmp
                k5_4[i,j] = -c3*c4*k1_1[i,j] + sigma_tmp
                k5_5[i,j] = (1.0/gph_l4**2.0 - c4**2.0)*k1_1[i,j] + sigma_tmp
                k5_6[i,j] = -c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k5_7[i,j]= -c1*c2*c4*k1_1[i,j] + sigma_tmp
                k5_8[i,j] = -c1*c3*c4*k1_1[i,j] + sigma_tmp
                k5_9[i,j] = -c1*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k5_10[i,j]= -c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k5_11[i,j]= -c2*c3*c4*k1_1[i,j] + sigma_tmp
                k5_12[i,j]= -c2*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k5_13[i,j]= -c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k5_14[i,j]= -c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k5_15[i,j]= -c4*(c4**2.0 - 3.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k6_1[i,j] =       (c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_2[i,j] =    c1*(c1**2.0 - 3.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_3[i,j] =    c2*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_4[i,j] =    c3*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_5[i,j] =    c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_6[i,j] = (c1**4.0 - 6.0*c1**2.0/gph_l1**2.0 + 3.0/gph_l1**4.0)*k1_1[i,j] + sigma_tmp
                k6_7[i,j] = c1*c2*(c1**2.0 - 3.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_8[i,j] = c1*c3*(c1**2.0 - 3.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_9[i,j] = c1*c4*(c1**2.0 - 3.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_10[i,j]=       (c1**2.0 - 1.0/gph_l1**2.0)*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k6_11[i,j]= c2*c3*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_12[i,j]= c2*c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_13[i,j]=       (c1**2.0 - 1.0/gph_l1**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k6_14[i,j]= c3*c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_15[i,j]=       (c1**2.0 - 1.0/gph_l1**2.0)*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k7_1[i,j] = c1*c2*k1_1[i,j] + sigma_tmp
                k7_2[i,j] = c2*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k7_3[i,j] = c1*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k7_4[i,j] = c1*c2*c3*k1_1[i,j] + sigma_tmp
                k7_5[i,j] = c1*c2*c4*k1_1[i,j] + sigma_tmp
                k7_6[i,j] = c1*c2*(c1**2.0 - 3.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k7_7[i,j] =       (c1**2.0 - 1.0/gph_l1**2.0)*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k7_8[i,j] = c2*c3*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k7_9[i,j] = c2*c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k7_10[i,j]= c1*c2*(c2**2.0 - 3.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k7_11[i,j]= c1*c3*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k7_12[i,j]= c1*c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k7_13[i,j]= c1*c2*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k7_14[i,j]= c1*c2*c3*c4*k1_1[i,j] + sigma_tmp
                k7_15[i,j]= c1*c2*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k8_1[i,j] = c1*c3*k1_1[i,j] + sigma_tmp
                k8_2[i,j] = c3*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k8_3[i,j] = c1*c2*c3*k1_1[i,j] + sigma_tmp
                k8_4[i,j] = c1*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k8_5[i,j] = c1*c3*c4*k1_1[i,j] + sigma_tmp
                k8_6[i,j] = c1*c3*(c1**2.0 - 3.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k8_7[i,j] = c2*c3*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k8_8[i,j] =       (c1**2.0 - 1.0/gph_l1**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k8_9[i,j] = c3*c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k8_10[i,j]= c1*c3*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k8_11[i,j]= c1*c2*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k8_12[i,j]= c1*c2*c3*c4*k1_1[i,j] + sigma_tmp
                k8_13[i,j]= c1*c3*(c3**2.0 - 3.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k8_14[i,j]= c1*c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k8_15[i,j]= c1*c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k9_1[i,j] = c1*c4*k1_1[i,j] + sigma_tmp
                k9_2[i,j] = c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k9_3[i,j] = c1*c2*c4*k1_1[i,j] + sigma_tmp
                k9_4[i,j] = c1*c3*c4*k1_1[i,j] + sigma_tmp
                k9_5[i,j] = c1*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k9_6[i,j] = c1*c4*(c1**2.0 - 3.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k9_7[i,j] = c2*c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k9_8[i,j] = c3*c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k9_9[i,j] =       (c4**2.0 - 1.0/gph_l4**2.0)*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k9_10[i,j]= c1*c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k9_11[i,j]= c1*c2*c3*c4*k1_1[i,j] + sigma_tmp
                k9_12[i,j]= c1*c2*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k9_13[i,j]= c1*c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k9_14[i,j]= c1*c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k9_15[i,j]= c1*c4*(c4**2.0 - 3.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k10_1[i,j] =       (c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_2[i,j] =    c1*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_3[i,j] =    c2*(c2**2.0 - 3.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_4[i,j] =    c3*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_5[i,j] =    c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_6[i,j] =       (c2**2.0 - 1.0/gph_l2**2.0)*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k10_7[i,j] = c2*c1*(c2**2.0 - 3.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_8[i,j] = c3*c1*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_9[i,j] = c4*c1*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_10[i,j]= (c2**4.0 - 6.0*c2**2.0/gph_l2**2.0 + 3.0/gph_l2**4.0)*k1_1[i,j] + sigma_tmp
                k10_11[i,j]= c2*c3*(c2**2.0 - 3.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_12[i,j]= c2*c4*(c2**2.0 - 3.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_13[i,j]=       (c2**2.0 - 1.0/gph_l2**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k10_14[i,j]= c3*c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_15[i,j]=       (c2**2.0 - 1.0/gph_l2**2.0)*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k11_1[i,j] = c2*c3*k1_1[i,j] + sigma_tmp
                k11_2[i,j] = c1*c2*c3*k1_1[i,j] + sigma_tmp
                k11_3[i,j] = c3*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k11_4[i,j] = c2*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k11_5[i,j] = c2*c3*c4*k1_1[i,j] + sigma_tmp
                k11_6[i,j] = c2*c3*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k11_7[i,j] = c1*c3*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k11_8[i,j] = c1*c2*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k11_9[i,j] = c1*c2*c3*c4*k1_1[i,j] + sigma_tmp
                k11_10[i,j]= c2*c3*(c2**2.0 - 3.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k11_11[i,j]=       (c2**2.0 - 1.0/gph_l2**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k11_12[i,j]= c3*c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k11_13[i,j]= c2*c3*(c3**2.0 - 3.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k11_14[i,j]= c2*c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k11_15[i,j]= c2*c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k12_1[i,j] = c2*c4*k1_1[i,j] + sigma_tmp
                k12_2[i,j] = c1*c2*c4*k1_1[i,j] + sigma_tmp
                k12_3[i,j] = c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k12_4[i,j] = c2*c3*c4*k1_1[i,j] + sigma_tmp
                k12_5[i,j] = c2*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k12_6[i,j] = c2*c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k12_7[i,j] = c1*c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k12_8[i,j] = c1*c2*c3*c4*k1_1[i,j] + sigma_tmp
                k12_9[i,j] =  c1*c2*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k12_10[i,j] = c2*c4*(c2**2.0 - 3.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k12_11[i,j] = c3*c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k12_12[i,j] =       (c4**2.0 - 1.0/gph_l4**2.0)*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k12_13[i,j] = c2*c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k12_14[i,j] = c2*c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k12_15[i,j] = c2*c4*(c4**2.0 - 3.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                 #
                k13_1[i,j] =        (c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_2[i,j] =     c1*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_3[i,j] =     c2*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_4[i,j] =     c3*(c3**2.0 - 3.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_5[i,j] =     c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_6[i,j] =        (c3**2.0 - 1.0/gph_l3**2.0)*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k13_7[i,j] =  c1*c2*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_8[i,j] =  c1*c3*(c3**2.0 - 3.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_9[i,j] =  c1*c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_10[i,j] =       (c3**2.0 - 1.0/gph_l3**2.0)*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k13_11[i,j] = c2*c3*(c3**2.0 - 3.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_12[i,j] = c2*c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_13[i,j] = (c3**4.0 - 6.0*c3**2.0/gph_l3**2.0 + 3.0/gph_l3**4.0)*k1_1[i,j] + sigma_tmp
                k13_14[i,j] = c3*c4*(c3**2.0 - 3.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_15[i,j] =       (c3**2.0 - 1.0/gph_l3**2.0)*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                 #
                k14_1[i,j] =    c3*c4*k1_1[i,j] + sigma_tmp
                k14_2[i,j] = c1*c3*c4*k1_1[i,j] + sigma_tmp
                k14_3[i,j] = c2*c3*c4*k1_1[i,j] + sigma_tmp
                k14_4[i,j] =    c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k14_5[i,j] =    c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k14_6[i,j] = c3*c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k14_7[i,j] = c1*c2*c3*c4*k1_1[i,j] + sigma_tmp
                k14_8[i,j] = c1*c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k14_9[i,j] = c1*c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k14_10[i,j]= c3*c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k14_11[i,j]= c2*c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k14_12[i,j]= c2*c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k14_13[i,j]= c3*c4*(c3**2.0 - 3.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k14_14[i,j]=       (c4**2.0 - 1.0/gph_l4**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k14_15[i,j]= c3*c4*(c4**2.0 - 3.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                 #
                k15_1[i,j] =       (c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_2[i,j] =    c1*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_3[i,j] =    c2*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_4[i,j] =    c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_5[i,j] =    c4*(c4**2.0 - 3.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_6[i,j] =       (c4**2.0 - 1.0/gph_l4**2.0)*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k15_7[i,j] = c1*c2*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_8[i,j] = c1*c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_9[i,j] = c1*c4*(c4**2.0 - 3.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_10[i,j] =      (c4**2.0 - 1.0/gph_l4**2.0)*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k15_11[i,j] =c2*c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_12[i,j] =c2*c4*(c4**2.0 - 3.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_13[i,j] =      (c4**2.0 - 1.0/gph_l4**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k15_14[i,j] =c3*c4*(c4**2.0 - 3.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_15[i,j] = (c4**4.0 - 6.0*c4**2.0/gph_l4**2.0 + 3.0/gph_l4**4.0)*k1_1[i,j] + sigma_tmp
        
        
        Kgrad =    np.concatenate((np.concatenate((k1_1,k1_2,k1_3,k1_4,k1_5),axis=1),\
                                   np.concatenate((k2_1,k2_2,k2_3,k2_4,k2_5),axis=1),\
                                   np.concatenate((k3_1,k3_2,k3_3,k3_4,k3_5),axis=1),\
                                   np.concatenate((k4_1,k4_2,k4_3,k4_4,k4_5),axis=1),\
                                   np.concatenate((k5_1,k5_2,k5_3,k5_4,k5_5),axis=1)), axis=0)
        
       
        
        
        Kgradlap = np.concatenate((np.concatenate((k1_1, k1_2, k1_3, k1_4, k1_5, k1_6, k1_7, k1_8, k1_9, k1_10, k1_11, k1_12, k1_13, k1_14, k1_15),axis=1),\
                                   np.concatenate((k2_1, k2_2, k2_3, k2_4, k2_5, k2_6, k2_7, k2_8, k2_9, k2_10, k2_11, k2_12, k2_13, k2_14, k2_15),axis=1),\
                                   np.concatenate((k3_1, k3_2, k3_3, k3_4, k3_5, k3_6, k3_7, k3_8, k3_9, k3_10, k3_11, k3_12, k3_13, k3_14, k3_15),axis=1),\
                                   np.concatenate((k4_1, k4_2, k4_3, k4_4, k4_5, k4_6, k4_7, k4_8, k4_9, k4_10, k4_11, k4_12, k4_13, k4_14, k4_15),axis=1),\
                                   np.concatenate((k5_1, k5_2, k5_3, k5_4, k5_5, k5_6, k5_7, k5_8, k5_9, k5_10, k5_11, k5_12, k5_13, k5_14, k5_15),axis=1),\
                                   np.concatenate((k6_1, k6_2, k6_3, k6_4, k6_5, k6_6, k6_7, k6_8, k6_9, k6_10, k6_11, k6_12, k6_13, k6_14, k6_15),axis=1),\
                                   np.concatenate((k7_1, k7_2, k7_3, k7_4, k7_5, k7_6, k7_7, k7_8, k7_9, k7_10, k7_11, k7_12, k7_13, k7_14, k7_15),axis=1),\
                                   np.concatenate((k8_1, k8_2, k8_3, k8_4, k8_5, k8_6, k8_7, k8_8, k8_9, k8_10, k8_11, k8_12, k8_13, k8_14, k8_15),axis=1),\
                                   np.concatenate((k9_1, k9_2, k9_3, k9_4, k9_5, k9_6, k9_7, k9_8, k9_9, k9_10, k9_11, k9_12, k9_13, k9_14, k9_15),axis=1),\
                                   np.concatenate((k10_1,k10_2,k10_3,k10_4,k10_5,k10_6,k10_7,k10_8,k10_9,k10_10,k10_11,k10_12,k10_13,k10_14,k10_15),axis=1),\
                                   np.concatenate((k11_1,k11_2,k11_3,k11_4,k11_5,k11_6,k11_7,k11_8,k11_9,k11_10,k11_11,k11_12,k11_13,k11_14,k11_15),axis=1),\
                                   np.concatenate((k12_1,k12_2,k12_3,k12_4,k12_5,k12_6,k12_7,k12_8,k12_9,k12_10,k12_11,k12_12,k12_13,k12_14,k12_15),axis=1),\
                                   np.concatenate((k13_1,k13_2,k13_3,k13_4,k13_5,k13_6,k13_7,k13_8,k13_9,k13_10,k13_11,k13_12,k13_13,k13_14,k13_15),axis=1),\
                                   np.concatenate((k14_1,k14_2,k14_3,k14_4,k14_5,k14_6,k14_7,k14_8,k14_9,k14_10,k14_11,k14_12,k14_13,k14_14,k14_15),axis=1),\
                                   np.concatenate((k15_1,k15_2,k15_3,k15_4,k15_5,k15_6,k15_7,k15_8,k15_9,k15_10,k15_11,k15_12,k15_13,k15_14,k15_15),axis=1)), axis=0)
       
        K0 = k1_1 + sigma_n*np.identity(len(k1_1))
        K1 = Kgrad + sigma_n*np.identity(len(Kgrad))
        K2 = Kgradlap + sigma_n*np.identity(len(Kgradlap))
    #    fold = f_df_data
    #    #print('postsym',K[15,16])
    #    n = len(K)
    #    L=np.linalg.cholesky(K)
    #    cho = scipy.linalg.cho_factor(K)
    #    alpha_lin = scipy.linalg.cho_solve(cho, fold)
    ##    #K_inv = inv(K)
    ##    #m1 = np.dot(K_inv,fold)
    ##    #part1 = np.dot(y_data.T,m1)
    ##    #part2 = math.log(np.linalg.det(K))
    #    part1 = np.dot(fold.T,alpha_lin)
    #    part2 = 2.0*sum(np.log(np.diag(L)))
    #    return(0.5 * part1 + 0.5 * part2 + n*math.log(2.0*math.pi)/2.0)
        return(K0,K1,K2)
        
    def GP_posterior(hyperparams,K0,K1,K2):
    
     
        gph_l1 = hyperparams[0]
        gph_l2 = hyperparams[1]
        gph_l3 = hyperparams[2]
        gph_l4 = hyperparams[3]
        
        kstar1_1 = np.zeros((dim_x_data))
        kstar1_2 = np.zeros((dim_x_data))
        kstar1_3 = np.zeros((dim_x_data))
        kstar1_4 = np.zeros((dim_x_data))
        kstar1_5 = np.zeros((dim_x_data))
        kstar1_6 = np.zeros((dim_x_data))
        kstar1_7 = np.zeros((dim_x_data))
        kstar1_8 = np.zeros((dim_x_data))
        kstar1_9 = np.zeros((dim_x_data))
        kstar1_10 = np.zeros((dim_x_data))
        kstar1_11 = np.zeros((dim_x_data))
        kstar1_12 = np.zeros((dim_x_data))
        kstar1_13 = np.zeros((dim_x_data))
        kstar1_14 = np.zeros((dim_x_data))
        kstar1_15 = np.zeros((dim_x_data))
        #
        kstar2_1 = np.zeros((dim_x_data))
        kstar2_2 = np.zeros((dim_x_data))
        kstar2_3 = np.zeros((dim_x_data))
        kstar2_4 = np.zeros((dim_x_data))
        kstar2_5 = np.zeros((dim_x_data))
        kstar2_6 = np.zeros((dim_x_data))
        kstar2_7 = np.zeros((dim_x_data))
        kstar2_8 = np.zeros((dim_x_data))
        kstar2_9 = np.zeros((dim_x_data))
        kstar2_10 = np.zeros((dim_x_data))
        kstar2_11 = np.zeros((dim_x_data))
        kstar2_12 = np.zeros((dim_x_data))
        kstar2_13 = np.zeros((dim_x_data))
        kstar2_14 = np.zeros((dim_x_data))
        kstar2_15 = np.zeros((dim_x_data))
        #
        kstar3_1 = np.zeros((dim_x_data))
        kstar3_2 = np.zeros((dim_x_data))
        kstar3_3 = np.zeros((dim_x_data))
        kstar3_4 = np.zeros((dim_x_data))
        kstar3_5 = np.zeros((dim_x_data))
        kstar3_6 = np.zeros((dim_x_data))
        kstar3_7 = np.zeros((dim_x_data))
        kstar3_8 = np.zeros((dim_x_data))
        kstar3_9 = np.zeros((dim_x_data))
        kstar3_10 = np.zeros((dim_x_data))
        kstar3_11 = np.zeros((dim_x_data))
        kstar3_12 = np.zeros((dim_x_data))
        kstar3_13 = np.zeros((dim_x_data))
        kstar3_14 = np.zeros((dim_x_data))
        kstar3_15 = np.zeros((dim_x_data))
        #
        kstar4_1 = np.zeros((dim_x_data))
        kstar4_2 = np.zeros((dim_x_data))
        kstar4_3 = np.zeros((dim_x_data))
        kstar4_4 = np.zeros((dim_x_data))
        kstar4_5 = np.zeros((dim_x_data))
        kstar4_6 = np.zeros((dim_x_data))
        kstar4_7 = np.zeros((dim_x_data))
        kstar4_8 = np.zeros((dim_x_data))
        kstar4_9 = np.zeros((dim_x_data))
        kstar4_10 = np.zeros((dim_x_data))
        kstar4_11 = np.zeros((dim_x_data))
        kstar4_12 = np.zeros((dim_x_data))
        kstar4_13 = np.zeros((dim_x_data))
        kstar4_14 = np.zeros((dim_x_data))
        kstar4_15 = np.zeros((dim_x_data))
        #
        kstar5_1 = np.zeros((dim_x_data))
        kstar5_2 = np.zeros((dim_x_data))
        kstar5_3 = np.zeros((dim_x_data))
        kstar5_4 = np.zeros((dim_x_data))
        kstar5_5 = np.zeros((dim_x_data))
        kstar5_6 = np.zeros((dim_x_data))
        kstar5_7 = np.zeros((dim_x_data))
        kstar5_8 = np.zeros((dim_x_data))
        kstar5_9 = np.zeros((dim_x_data))
        kstar5_10 = np.zeros((dim_x_data))
        kstar5_11 = np.zeros((dim_x_data))
        kstar5_12 = np.zeros((dim_x_data))
        kstar5_13 = np.zeros((dim_x_data))
        kstar5_14 = np.zeros((dim_x_data))
        kstar5_15 = np.zeros((dim_x_data))
        #
        kstar6_1 = np.zeros((dim_x_data))
        kstar6_2 = np.zeros((dim_x_data))
        kstar6_3 = np.zeros((dim_x_data))
        kstar6_4 = np.zeros((dim_x_data))
        kstar6_5 = np.zeros((dim_x_data))
        kstar6_6 = np.zeros((dim_x_data))
        kstar6_7 = np.zeros((dim_x_data))
        kstar6_8 = np.zeros((dim_x_data))
        kstar6_9 = np.zeros((dim_x_data))
        kstar6_10 = np.zeros((dim_x_data))
        kstar6_11 = np.zeros((dim_x_data))
        kstar6_12 = np.zeros((dim_x_data))
        kstar6_13 = np.zeros((dim_x_data))
        kstar6_14 = np.zeros((dim_x_data))
        kstar6_15 = np.zeros((dim_x_data))
        #
        kstar7_1 = np.zeros((dim_x_data))
        kstar7_2 = np.zeros((dim_x_data))
        kstar7_3 = np.zeros((dim_x_data))
        kstar7_4 = np.zeros((dim_x_data))
        kstar7_5 = np.zeros((dim_x_data))
        kstar7_6 = np.zeros((dim_x_data))
        kstar7_7 = np.zeros((dim_x_data))
        kstar7_8 = np.zeros((dim_x_data))
        kstar7_9 = np.zeros((dim_x_data))
        kstar7_10 = np.zeros((dim_x_data))
        kstar7_11 = np.zeros((dim_x_data))
        kstar7_12 = np.zeros((dim_x_data))
        kstar7_13 = np.zeros((dim_x_data))
        kstar7_14 = np.zeros((dim_x_data))
        kstar7_15 = np.zeros((dim_x_data))
        #
        kstar8_1 = np.zeros((dim_x_data))
        kstar8_2 = np.zeros((dim_x_data))
        kstar8_3 = np.zeros((dim_x_data))
        kstar8_4 = np.zeros((dim_x_data))
        kstar8_5 = np.zeros((dim_x_data))
        kstar8_6 = np.zeros((dim_x_data))
        kstar8_7 = np.zeros((dim_x_data))
        kstar8_8 = np.zeros((dim_x_data))
        kstar8_9 = np.zeros((dim_x_data))
        kstar8_10 = np.zeros((dim_x_data))
        kstar8_11 = np.zeros((dim_x_data))
        kstar8_12 = np.zeros((dim_x_data))
        kstar8_13 = np.zeros((dim_x_data))
        kstar8_14 = np.zeros((dim_x_data))
        kstar8_15 = np.zeros((dim_x_data))
        #
        kstar9_1 = np.zeros((dim_x_data))
        kstar9_2 = np.zeros((dim_x_data))
        kstar9_3 = np.zeros((dim_x_data))
        kstar9_4 = np.zeros((dim_x_data))
        kstar9_5 = np.zeros((dim_x_data))
        kstar9_6 = np.zeros((dim_x_data))
        kstar9_7 = np.zeros((dim_x_data))
        kstar9_8 = np.zeros((dim_x_data))
        kstar9_9 = np.zeros((dim_x_data))
        kstar9_10 = np.zeros((dim_x_data))
        kstar9_11 = np.zeros((dim_x_data))
        kstar9_12 = np.zeros((dim_x_data))
        kstar9_13 = np.zeros((dim_x_data))
        kstar9_14 = np.zeros((dim_x_data))
        kstar9_15 = np.zeros((dim_x_data))
        #
        kstar10_1 = np.zeros((dim_x_data))
        kstar10_2 = np.zeros((dim_x_data))
        kstar10_3 = np.zeros((dim_x_data))
        kstar10_4 = np.zeros((dim_x_data))
        kstar10_5 = np.zeros((dim_x_data))
        kstar10_6 = np.zeros((dim_x_data))
        kstar10_7 = np.zeros((dim_x_data))
        kstar10_8 = np.zeros((dim_x_data))
        kstar10_9 = np.zeros((dim_x_data))
        kstar10_10 = np.zeros((dim_x_data))
        kstar10_11 = np.zeros((dim_x_data))
        kstar10_12 = np.zeros((dim_x_data))
        kstar10_13 = np.zeros((dim_x_data))
        kstar10_14 = np.zeros((dim_x_data))
        kstar10_15 = np.zeros((dim_x_data))
        #
        kstar11_1 = np.zeros((dim_x_data))
        kstar11_2 = np.zeros((dim_x_data))
        kstar11_3 = np.zeros((dim_x_data))
        kstar11_4 = np.zeros((dim_x_data))
        kstar11_5 = np.zeros((dim_x_data))
        kstar11_6 = np.zeros((dim_x_data))
        kstar11_7 = np.zeros((dim_x_data))
        kstar11_8 = np.zeros((dim_x_data))
        kstar11_9 = np.zeros((dim_x_data))
        kstar11_10 = np.zeros((dim_x_data))
        kstar11_11 = np.zeros((dim_x_data))
        kstar11_12 = np.zeros((dim_x_data))
        kstar11_13 = np.zeros((dim_x_data))
        kstar11_14 = np.zeros((dim_x_data))
        kstar11_15 = np.zeros((dim_x_data))
        #
        kstar12_1 = np.zeros((dim_x_data))
        kstar12_2 = np.zeros((dim_x_data))
        kstar12_3 = np.zeros((dim_x_data))
        kstar12_4 = np.zeros((dim_x_data))
        kstar12_5 = np.zeros((dim_x_data))
        kstar12_6 = np.zeros((dim_x_data))
        kstar12_7 = np.zeros((dim_x_data))
        kstar12_8 = np.zeros((dim_x_data))
        kstar12_9 = np.zeros((dim_x_data))
        kstar12_10 = np.zeros((dim_x_data))
        kstar12_11 = np.zeros((dim_x_data))
        kstar12_12 = np.zeros((dim_x_data))
        kstar12_13 = np.zeros((dim_x_data))
        kstar12_14 = np.zeros((dim_x_data))
        kstar12_15 = np.zeros((dim_x_data))
         #
        kstar13_1 = np.zeros((dim_x_data))
        kstar13_2 = np.zeros((dim_x_data))
        kstar13_3 = np.zeros((dim_x_data))
        kstar13_4 = np.zeros((dim_x_data))
        kstar13_5 = np.zeros((dim_x_data))
        kstar13_6 = np.zeros((dim_x_data))
        kstar13_7 = np.zeros((dim_x_data))
        kstar13_8 = np.zeros((dim_x_data))
        kstar13_9 = np.zeros((dim_x_data))
        kstar13_10 = np.zeros((dim_x_data))
        kstar13_11 = np.zeros((dim_x_data))
        kstar13_12 = np.zeros((dim_x_data))
        kstar13_13 = np.zeros((dim_x_data))
        kstar13_14 = np.zeros((dim_x_data))
        kstar13_15 = np.zeros((dim_x_data))
         #
        kstar14_1 = np.zeros((dim_x_data))
        kstar14_2 = np.zeros((dim_x_data))
        kstar14_3 = np.zeros((dim_x_data))
        kstar14_4 = np.zeros((dim_x_data))
        kstar14_5 = np.zeros((dim_x_data))
        kstar14_6 = np.zeros((dim_x_data))
        kstar14_7 = np.zeros((dim_x_data))
        kstar14_8 = np.zeros((dim_x_data))
        kstar14_9 = np.zeros((dim_x_data))
        kstar14_10 = np.zeros((dim_x_data))
        kstar14_11 = np.zeros((dim_x_data))
        kstar14_12 = np.zeros((dim_x_data))
        kstar14_13 = np.zeros((dim_x_data))
        kstar14_14 = np.zeros((dim_x_data))
        kstar14_15 = np.zeros((dim_x_data))
         #
        kstar15_1 = np.zeros((dim_x_data))
        kstar15_2 = np.zeros((dim_x_data))
        kstar15_3 = np.zeros((dim_x_data))
        kstar15_4 = np.zeros((dim_x_data))
        kstar15_5 = np.zeros((dim_x_data))
        kstar15_6 = np.zeros((dim_x_data))
        kstar15_7 = np.zeros((dim_x_data))
        kstar15_8 = np.zeros((dim_x_data))
        kstar15_9 = np.zeros((dim_x_data))
        kstar15_10 = np.zeros((dim_x_data))
        kstar15_11 = np.zeros((dim_x_data))
        kstar15_12 = np.zeros((dim_x_data))
        kstar15_13 = np.zeros((dim_x_data))
        kstar15_14 = np.zeros((dim_x_data))
        kstar15_15 = np.zeros((dim_x_data))
        
        for i in range(dim_x_data):
                
                dw = w_star - w_data[i]
                dx = x_star - x_data[i]
                dy = y_star - y_data[i]
                dz = z_star - z_data[i]
                coef = gph_sigma_f * gph_sigma_f * exp( - 0.5 * ( dw**2.0/gph_l1**2.0 + dx**2.0/gph_l2**2.0 + dy**2.0/gph_l3**2.0 + dz**2.0/gph_l4**2.0 ) )
                c1 = dw / gph_l1**2.0 
                c2 = dx / gph_l2**2.0
                c3 = dy / gph_l3**2.0
                c4 = dz / gph_l4**2.0
                #print('c1=',c1,'c2=',c2,'c3=',c3,'c4=',c4)
                # Effectively creates an upper triangular matrix
                kstar1_1[i] = coef 
                kstar1_2[i] = c1*kstar1_1[i]
                kstar1_3[i] = c2*kstar1_1[i] 
                kstar1_4[i] = c3*kstar1_1[i] 
                kstar1_5[i] = c4*kstar1_1[i] 
                kstar1_6[i] = (c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar1_7[i] = c1*c2*kstar1_1[i]
                kstar1_8[i] = c1*c3*kstar1_1[i] 
                kstar1_9[i] = c1*c4*kstar1_1[i] 
                kstar1_10[i] = (c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar1_11[i] = c2*c3*kstar1_1[i] 
                kstar1_12[i] = c2*c4*kstar1_1[i] 
                kstar1_13[i] = (c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar1_14[i] = c3*c4*kstar1_1[i] 
                kstar1_15[i] = (c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                #
                kstar2_1[i] = -c1*kstar1_1[i] 
                kstar2_2[i] = (1.0/gph_l1**2.0 - c1**2.0)*kstar1_1[i] 
                kstar2_3[i] = -c1*c2*kstar1_1[i]
                kstar2_4[i] = -c1*c3*kstar1_1[i] 
                kstar2_5[i] = -c1*c4*kstar1_1[i] 
                kstar2_6[i] = -c1*(c1**2.0 - 3.0/gph_l1**2.0)*kstar1_1[i] 
                kstar2_7[i] = -c2*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar2_8[i] = -c3*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar2_9[i] = -c4*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar2_10[i]= -c1*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar2_11[i]= -c1*c2*c3*kstar1_1[i] 
                kstar2_12[i]= -c1*c2*c4*kstar1_1[i] 
                kstar2_13[i]= -c1*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar2_14[i]= -c1*c3*c4*kstar1_1[i] 
                kstar2_15[i]= -c1*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                #
                kstar3_1[i] = -c2*kstar1_1[i] 
                kstar3_2[i] = -c1*c2*kstar1_1[i]
                kstar3_3[i] = (1.0/gph_l2**2.0 - c2**2.0)*kstar1_1[i] 
                kstar3_4[i] = -c2*c3*kstar1_1[i] 
                kstar3_5[i] = -c2*c4*kstar1_1[i] 
                kstar3_6[i] = -c2*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar3_7[i] = -c1*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar3_8[i] = -c1*c2*c3*kstar1_1[i] 
                kstar3_9[i] = -c1*c2*c4*kstar1_1[i] 
                kstar3_10[i] =-c2*(c2**2.0 - 3.0/gph_l2**2.0)*kstar1_1[i] 
                kstar3_11[i] =-c3*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar3_12[i] =-c4*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar3_13[i] =-c2*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar3_14[i] =-c2*c3*c4*kstar1_1[i] 
                kstar3_15[i] =-c2*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                #
                kstar4_1[i] = -c3*kstar1_1[i] 
                kstar4_2[i] = -c1*c3*kstar1_1[i] 
                kstar4_3[i] = -c2*c3*kstar1_1[i] 
                kstar4_4[i] = (1.0/gph_l3**2.0 - c3**2.0)*kstar1_1[i] 
                kstar4_5[i] = -c3*c4*kstar1_1[i] 
                kstar4_6[i] = -c3*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar4_7[i] = -c1*c2*c3*kstar1_1[i] 
                kstar4_8[i] = -c1*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar4_9[i] = -c1*c3*c4*kstar1_1[i] 
                kstar4_10[i]= -c3*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar4_11[i]= -c2*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar4_12[i]= -c2*c3*c4*kstar1_1[i] 
                kstar4_13[i]= -c3*(c3**2.0 - 3.0/gph_l3**2.0)*kstar1_1[i] 
                kstar4_14[i]= -c4*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar4_15[i]= -c3*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                #
                kstar5_1[i] = -c4*kstar1_1[i]
                kstar5_2[i] = -c1*c4*kstar1_1[i] 
                kstar5_3[i] = -c2*c4*kstar1_1[i] 
                kstar5_4[i] = -c3*c4*kstar1_1[i] 
                kstar5_5[i] = (1.0/gph_l4**2.0 - c4**2.0)*kstar1_1[i] 
                kstar5_6[i] = -c4*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar5_7[i]= -c1*c2*c4*kstar1_1[i]
                kstar5_8[i] = -c1*c3*c4*kstar1_1[i] 
                kstar5_9[i] = -c1*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                kstar5_10[i]= -c4*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar5_11[i]= -c2*c3*c4*kstar1_1[i] 
                kstar5_12[i]= -c2*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                kstar5_13[i]= -c4*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar5_14[i]= -c3*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                kstar5_15[i]= -c4*(c4**2.0 - 3.0/gph_l4**2.0)*kstar1_1[i]
                #
                kstar6_1[i] =       (c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i]
                kstar6_2[i] =    c1*(c1**2.0 - 3.0/gph_l1**2.0)*kstar1_1[i] 
                kstar6_3[i] =    c2*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar6_4[i] =    c3*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar6_5[i] =    c4*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar6_6[i] = (c1**4.0 - 6.0*c1**2.0/gph_l1**2.0 + 3.0/gph_l1**4.0)*kstar1_1[i] 
                kstar6_7[i] = c1*c2*(c1**2.0 - 3.0/gph_l1**2.0)*kstar1_1[i] 
                kstar6_8[i] = c1*c3*(c1**2.0 - 3.0/gph_l1**2.0)*kstar1_1[i] 
                kstar6_9[i] = c1*c4*(c1**2.0 - 3.0/gph_l1**2.0)*kstar1_1[i] 
                kstar6_10[i]=       (c1**2.0 - 1.0/gph_l1**2.0)*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar6_11[i]= c2*c3*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar6_12[i]= c2*c4*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar6_13[i]=       (c1**2.0 - 1.0/gph_l1**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar6_14[i]= c3*c4*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar6_15[i]=       (c1**2.0 - 1.0/gph_l1**2.0)*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                #
                kstar7_1[i] = c1*c2*kstar1_1[i] 
                kstar7_2[i] = c2*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar7_3[i] = c1*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar7_4[i] = c1*c2*c3*kstar1_1[i] 
                kstar7_5[i] = c1*c2*c4*kstar1_1[i] 
                kstar7_6[i] = c1*c2*(c1**2.0 - 3.0/gph_l1**2.0)*kstar1_1[i] 
                kstar7_7[i] =       (c1**2.0 - 1.0/gph_l1**2.0)*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i]
                kstar7_8[i] = c2*c3*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar7_9[i] = c2*c4*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar7_10[i]= c1*c2*(c2**2.0 - 3.0/gph_l2**2.0)*kstar1_1[i] 
                kstar7_11[i]= c1*c3*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i]
                kstar7_12[i]= c1*c4*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar7_13[i]= c1*c2*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar7_14[i]= c1*c2*c3*c4*kstar1_1[i] 
                kstar7_15[i]= c1*c2*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                #
                kstar8_1[i] = c1*c3*kstar1_1[i]
                kstar8_2[i] = c3*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar8_3[i] = c1*c2*c3*kstar1_1[i]
                kstar8_4[i] = c1*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar8_5[i] = c1*c3*c4*kstar1_1[i]
                kstar8_6[i] = c1*c3*(c1**2.0 - 3.0/gph_l1**2.0)*kstar1_1[i]
                kstar8_7[i] = c2*c3*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar8_8[i] =       (c1**2.0 - 1.0/gph_l1**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar8_9[i] = c3*c4*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar8_10[i]= c1*c3*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar8_11[i]= c1*c2*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar8_12[i]= c1*c2*c3*c4*kstar1_1[i] 
                kstar8_13[i]= c1*c3*(c3**2.0 - 3.0/gph_l3**2.0)*kstar1_1[i] 
                kstar8_14[i]= c1*c4*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar8_15[i]= c1*c3*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                #
                kstar9_1[i] = c1*c4*kstar1_1[i] 
                kstar9_2[i] = c4*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar9_3[i] = c1*c2*c4*kstar1_1[i]
                kstar9_4[i] = c1*c3*c4*kstar1_1[i] 
                kstar9_5[i] = c1*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i]
                kstar9_6[i] = c1*c4*(c1**2.0 - 3.0/gph_l1**2.0)*kstar1_1[i] 
                kstar9_7[i] = c2*c4*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar9_8[i] = c3*c4*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar9_9[i] =       (c4**2.0 - 1.0/gph_l4**2.0)*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar9_10[i]= c1*c4*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i]
                kstar9_11[i]= c1*c2*c3*c4*kstar1_1[i] 
                kstar9_12[i]= c1*c2*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                kstar9_13[i]= c1*c4*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar9_14[i]= c1*c3*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                kstar9_15[i]= c1*c4*(c4**2.0 - 3.0/gph_l4**2.0)*kstar1_1[i] 
                #
                kstar10_1[i] =       (c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i]
                kstar10_2[i] =    c1*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar10_3[i] =    c2*(c2**2.0 - 3.0/gph_l2**2.0)*kstar1_1[i]
                kstar10_4[i] =    c3*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar10_5[i] =    c4*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar10_6[i] =       (c2**2.0 - 1.0/gph_l2**2.0)*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar10_7[i] = c2*c1*(c2**2.0 - 3.0/gph_l2**2.0)*kstar1_1[i] 
                kstar10_8[i] = c3*c1*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar10_9[i] = c4*c1*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar10_10[i]= (c2**4.0 - 6.0*c2**2.0/gph_l2**2.0 + 3.0/gph_l2**4.0)*kstar1_1[i] 
                kstar10_11[i]= c2*c3*(c2**2.0 - 3.0/gph_l2**2.0)*kstar1_1[i] 
                kstar10_12[i]= c2*c4*(c2**2.0 - 3.0/gph_l2**2.0)*kstar1_1[i] 
                kstar10_13[i]=       (c2**2.0 - 1.0/gph_l2**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar10_14[i]= c3*c4*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar10_15[i]=       (c2**2.0 - 1.0/gph_l2**2.0)*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i]
                #
                kstar11_1[i] = c2*c3*kstar1_1[i]
                kstar11_2[i] = c1*c2*c3*kstar1_1[i] 
                kstar11_3[i] = c3*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i]
                kstar11_4[i] = c2*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i]
                kstar11_5[i] = c2*c3*c4*kstar1_1[i] 
                kstar11_6[i] = c2*c3*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar11_7[i] = c1*c3*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i]
                kstar11_8[i] = c1*c2*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i]
                kstar11_9[i] = c1*c2*c3*c4*kstar1_1[i] 
                kstar11_10[i]= c2*c3*(c2**2.0 - 3.0/gph_l2**2.0)*kstar1_1[i] 
                kstar11_11[i]=       (c2**2.0 - 1.0/gph_l2**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i]
                kstar11_12[i]= c3*c4*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar11_13[i]= c2*c3*(c3**2.0 - 3.0/gph_l3**2.0)*kstar1_1[i]
                kstar11_14[i]= c2*c4*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar11_15[i]= c2*c3*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                #
                kstar12_1[i] = c2*c4*kstar1_1[i] 
                kstar12_2[i] = c1*c2*c4*kstar1_1[i] 
                kstar12_3[i] = c4*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i]
                kstar12_4[i] = c2*c3*c4*kstar1_1[i] 
                kstar12_5[i] = c2*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                kstar12_6[i] = c2*c4*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar12_7[i] = c1*c4*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar12_8[i] = c1*c2*c3*c4*kstar1_1[i] 
                kstar12_9[i] =  c1*c2*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                kstar12_10[i] = c2*c4*(c2**2.0 - 3.0/gph_l2**2.0)*kstar1_1[i]
                kstar12_11[i] = c3*c4*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i]
                kstar12_12[i] =       (c4**2.0 - 1.0/gph_l4**2.0)*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i]
                kstar12_13[i] = c2*c4*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar12_14[i] = c2*c3*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                kstar12_15[i] = c2*c4*(c4**2.0 - 3.0/gph_l4**2.0)*kstar1_1[i] 
                 #
                kstar13_1[i] =        (c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar13_2[i] =     c1*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i]
                kstar13_3[i] =     c2*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar13_4[i] =     c3*(c3**2.0 - 3.0/gph_l3**2.0)*kstar1_1[i] 
                kstar13_5[i] =     c4*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar13_6[i] =        (c3**2.0 - 1.0/gph_l3**2.0)*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar13_7[i] =  c1*c2*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar13_8[i] =  c1*c3*(c3**2.0 - 3.0/gph_l3**2.0)*kstar1_1[i] 
                kstar13_9[i] =  c1*c4*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar13_10[i] =       (c3**2.0 - 1.0/gph_l3**2.0)*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar13_11[i] = c2*c3*(c3**2.0 - 3.0/gph_l3**2.0)*kstar1_1[i] 
                kstar13_12[i] = c2*c4*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar13_13[i] = (c3**4.0 - 6.0*c3**2.0/gph_l3**2.0 + 3.0/gph_l3**4.0)*kstar1_1[i]
                kstar13_14[i] = c3*c4*(c3**2.0 - 3.0/gph_l3**2.0)*kstar1_1[i] 
                kstar13_15[i] =       (c3**2.0 - 1.0/gph_l3**2.0)*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i]
                 #
                kstar14_1[i] =    c3*c4*kstar1_1[i] 
                kstar14_2[i] = c1*c3*c4*kstar1_1[i] 
                kstar14_3[i] = c2*c3*c4*kstar1_1[i] 
                kstar14_4[i] =    c4*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar14_5[i] =    c3*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                kstar14_6[i] = c3*c4*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i] 
                kstar14_7[i] = c1*c2*c3*c4*kstar1_1[i] 
                kstar14_8[i] = c1*c4*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar14_9[i] = c1*c3*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                kstar14_10[i]= c3*c4*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar14_11[i]= c2*c4*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar14_12[i]= c2*c3*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                kstar14_13[i]= c3*c4*(c3**2.0 - 3.0/gph_l3**2.0)*kstar1_1[i] 
                kstar14_14[i]=       (c4**2.0 - 1.0/gph_l4**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i]
                kstar14_15[i]= c3*c4*(c4**2.0 - 3.0/gph_l4**2.0)*kstar1_1[i] 
                 #
                kstar15_1[i] =       (c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i]
                kstar15_2[i] =    c1*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                kstar15_3[i] =    c2*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i]
                kstar15_4[i] =    c3*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                kstar15_5[i] =    c4*(c4**2.0 - 3.0/gph_l4**2.0)*kstar1_1[i] 
                kstar15_6[i] =       (c4**2.0 - 1.0/gph_l4**2.0)*(c1**2.0 - 1.0/gph_l1**2.0)*kstar1_1[i]
                kstar15_7[i] = c1*c2*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                kstar15_8[i] = c1*c3*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                kstar15_9[i] = c1*c4*(c4**2.0 - 3.0/gph_l4**2.0)*kstar1_1[i] 
                kstar15_10[i] =      (c4**2.0 - 1.0/gph_l4**2.0)*(c2**2.0 - 1.0/gph_l2**2.0)*kstar1_1[i] 
                kstar15_11[i] =c2*c3*(c4**2.0 - 1.0/gph_l4**2.0)*kstar1_1[i] 
                kstar15_12[i] =c2*c4*(c4**2.0 - 3.0/gph_l4**2.0)*kstar1_1[i] 
                kstar15_13[i] =      (c4**2.0 - 1.0/gph_l4**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*kstar1_1[i] 
                kstar15_14[i] =c3*c4*(c4**2.0 - 3.0/gph_l4**2.0)*kstar1_1[i] 
                kstar15_15[i] = (c4**4.0 - 6.0*c4**2.0/gph_l4**2.0 + 3.0/gph_l4**4.0)*kstar1_1[i] 
        
        
        Kstargrad =    np.stack((np.concatenate((kstar1_1,kstar1_2,kstar1_3,kstar1_4,kstar1_5)),\
                                   np.concatenate((kstar2_1,kstar2_2,kstar2_3,kstar2_4,kstar2_5)),\
                                   np.concatenate((kstar3_1,kstar3_2,kstar3_3,kstar3_4,kstar3_5)),\
                                   np.concatenate((kstar4_1,kstar4_2,kstar4_3,kstar4_4,kstar4_5)),\
                                   np.concatenate((kstar5_1,kstar5_2,kstar5_3,kstar5_4,kstar5_5))))
        
       
        
        
        Kstargradlap = np.stack((np.concatenate((kstar1_1, kstar1_2, kstar1_3, kstar1_4, kstar1_5, kstar1_6, kstar1_7, kstar1_8, kstar1_9, kstar1_10, kstar1_11, kstar1_12, kstar1_13, kstar1_14, kstar1_15)),\
                                   np.concatenate((kstar2_1, kstar2_2, kstar2_3, kstar2_4, kstar2_5, kstar2_6, kstar2_7, kstar2_8, kstar2_9, kstar2_10, kstar2_11, kstar2_12, kstar2_13, kstar2_14, kstar2_15)),\
                                   np.concatenate((kstar3_1, kstar3_2, kstar3_3, kstar3_4, kstar3_5, kstar3_6, kstar3_7, kstar3_8, kstar3_9, kstar3_10, kstar3_11, kstar3_12, kstar3_13, kstar3_14, kstar3_15)),\
                                   np.concatenate((kstar4_1, kstar4_2, kstar4_3, kstar4_4, kstar4_5, kstar4_6, kstar4_7, kstar4_8, kstar4_9, kstar4_10, kstar4_11, kstar4_12, kstar4_13, kstar4_14, kstar4_15)),\
                                   np.concatenate((kstar5_1, kstar5_2, kstar5_3, kstar5_4, kstar5_5, kstar5_6, kstar5_7, kstar5_8, kstar5_9, kstar5_10, kstar5_11, kstar5_12, kstar5_13, kstar5_14, kstar5_15)),\
                                   np.concatenate((kstar6_1, kstar6_2, kstar6_3, kstar6_4, kstar6_5, kstar6_6, kstar6_7, kstar6_8, kstar6_9, kstar6_10, kstar6_11, kstar6_12, kstar6_13, kstar6_14, kstar6_15)),\
                                   np.concatenate((kstar7_1, kstar7_2, kstar7_3, kstar7_4, kstar7_5, kstar7_6, kstar7_7, kstar7_8, kstar7_9, kstar7_10, kstar7_11, kstar7_12, kstar7_13, kstar7_14, kstar7_15)),\
                                   np.concatenate((kstar8_1, kstar8_2, kstar8_3, kstar8_4, kstar8_5, kstar8_6, kstar8_7, kstar8_8, kstar8_9, kstar8_10, kstar8_11, kstar8_12, kstar8_13, kstar8_14, kstar8_15)),\
                                   np.concatenate((kstar9_1, kstar9_2, kstar9_3, kstar9_4, kstar9_5, kstar9_6, kstar9_7, kstar9_8, kstar9_9, kstar9_10, kstar9_11, kstar9_12, kstar9_13, kstar9_14, kstar9_15)),\
                                   np.concatenate((kstar10_1,kstar10_2,kstar10_3,kstar10_4,kstar10_5,kstar10_6,kstar10_7,kstar10_8,kstar10_9,kstar10_10,kstar10_11,kstar10_12,kstar10_13,kstar10_14,kstar10_15)),\
                                   np.concatenate((kstar11_1,kstar11_2,kstar11_3,kstar11_4,kstar11_5,kstar11_6,kstar11_7,kstar11_8,kstar11_9,kstar11_10,kstar11_11,kstar11_12,kstar11_13,kstar11_14,kstar11_15)),\
                                   np.concatenate((kstar12_1,kstar12_2,kstar12_3,kstar12_4,kstar12_5,kstar12_6,kstar12_7,kstar12_8,kstar12_9,kstar12_10,kstar12_11,kstar12_12,kstar12_13,kstar12_14,kstar12_15)),\
                                   np.concatenate((kstar13_1,kstar13_2,kstar13_3,kstar13_4,kstar13_5,kstar13_6,kstar13_7,kstar13_8,kstar13_9,kstar13_10,kstar13_11,kstar13_12,kstar13_13,kstar13_14,kstar13_15)),\
                                   np.concatenate((kstar14_1,kstar14_2,kstar14_3,kstar14_4,kstar14_5,kstar14_6,kstar14_7,kstar14_8,kstar14_9,kstar14_10,kstar14_11,kstar14_12,kstar14_13,kstar14_14,kstar14_15)),\
                                   np.concatenate((kstar15_1,kstar15_2,kstar15_3,kstar15_4,kstar15_5,kstar15_6,kstar15_7,kstar15_8,kstar15_9,kstar15_10,kstar15_11,kstar15_12,kstar15_13,kstar15_14,kstar15_15))))
       
        Kstar0 = kstar1_1 
        Kstar1 = Kstargrad 
        Kstar2 = Kstargradlap
        
        
        #m1 = np.dot(Kstar,K_inv)
        fnew0 = np.dot(Kstar0,scipy.linalg.cho_solve(scipy.linalg.cho_factor(K0), fold[0:1*dim_x_data]))
        fnew1 = np.dot(Kstar1,scipy.linalg.cho_solve(scipy.linalg.cho_factor(K1), fold[0:5*dim_x_data]))
        fnew2 = np.dot(Kstar2,scipy.linalg.cho_solve(scipy.linalg.cho_factor(K2), fold))
    #    fold = f_df_data
    #    #print('postsym',K[15,16])
    #    n = len(K)
    #    L=np.linalg.cholesky(K)
    #    cho = scipy.linalg.cho_factor(K)
    #    alpha_lin = scipy.linalg.cho_solve(cho, fold)
    ##    #K_inv = inv(K)
    ##    #m1 = np.dot(K_inv,fold)
    ##    #part1 = np.dot(y_data.T,m1)
    ##    #part2 = math.log(np.linalg.det(K))
    #    part1 = np.dot(fold.T,alpha_lin)
    #    part2 = 2.0*sum(np.log(np.diag(L)))
    #    return(0.5 * part1 + 0.5 * part2 + n*math.log(2.0*math.pi)/2.0)
        return(fnew0,fnew1,fnew2)
    
    def log_likelihood_function_0(hyperparams):
 
        gph_l1 = hyperparams[0]
        gph_l2 = hyperparams[1]
        gph_l3 = hyperparams[2]
        gph_l4 = hyperparams[3]
        
        k1_1 = np.zeros((dim_x_data,dim_x_data))
        
        for i in range(dim_x_data):
            for j in range(dim_x_data):
                
                dw = w_data[i] - w_data[j]
                dx = x_data[i] - x_data[j]
                dy = y_data[i] - y_data[j]
                dz = z_data[i] - z_data[j]
                coef = - 0.5 * ( dw**2.0/gph_l1**2.0 + dx**2.0/gph_l2**2.0 + dy**2.0/gph_l3**2.0 + dz**2.0/gph_l4**2.0 ) 
                c1 = dw / gph_l1**2.0 
                c2 = dx / gph_l2**2.0
                c3 = dy / gph_l3**2.0
                c4 = dz / gph_l4**2.0
                k1_1[i,j] = gph_sigma_f * gph_sigma_f * exp(coef)
          
        K = k1_1 + sigma_n*np.identity(len(k1_1))
        fold = f_df_data[0:1*dim_x_data]
        
        n = len(K)
        L=np.linalg.cholesky(K)
        cho = scipy.linalg.cho_factor(K)
        alpha_lin = scipy.linalg.cho_solve(cho, fold)
        #K_inv = inv(K)
        #m1 = np.dot(K_inv,fold)
        #part1 = np.dot(y_data.T,m1)
        #part2 = math.log(np.linalg.det(K))
        part1 = np.dot(fold.T,alpha_lin)
        part2 = 2.0*sum(np.log(np.diag(L)))
        return(0.5 * part1 + 0.5 * part2 + n*log(2.0*pi)/2.0)   
        
    def log_likelihood_function_1(hyperparams):
        
        gph_l1 = hyperparams[0]
        gph_l2 = hyperparams[1]
        gph_l3 = hyperparams[2]
        gph_l4 = hyperparams[3]
        
        k1_1 = np.zeros((dim_x_data,dim_x_data))
        k1_2 = np.zeros((dim_x_data,dim_x_data))
        k1_3 = np.zeros((dim_x_data,dim_x_data))
        k1_4 = np.zeros((dim_x_data,dim_x_data))
        k1_5 = np.zeros((dim_x_data,dim_x_data))
        k2_1 = np.zeros((dim_x_data,dim_x_data))
        k2_2 = np.zeros((dim_x_data,dim_x_data))
        k2_3 = np.zeros((dim_x_data,dim_x_data))
        k2_4 = np.zeros((dim_x_data,dim_x_data))
        k2_5 = np.zeros((dim_x_data,dim_x_data))
        k3_1 = np.zeros((dim_x_data,dim_x_data))
        k3_2 = np.zeros((dim_x_data,dim_x_data))
        k3_3 = np.zeros((dim_x_data,dim_x_data))
        k3_4 = np.zeros((dim_x_data,dim_x_data))
        k3_5 = np.zeros((dim_x_data,dim_x_data))
        k4_1 = np.zeros((dim_x_data,dim_x_data))
        k4_2 = np.zeros((dim_x_data,dim_x_data))
        k4_3 = np.zeros((dim_x_data,dim_x_data))
        k4_4 = np.zeros((dim_x_data,dim_x_data))
        k4_5 = np.zeros((dim_x_data,dim_x_data))
        k5_1 = np.zeros((dim_x_data,dim_x_data))
        k5_2 = np.zeros((dim_x_data,dim_x_data))
        k5_3 = np.zeros((dim_x_data,dim_x_data))
        k5_4 = np.zeros((dim_x_data,dim_x_data))
        k5_5 = np.zeros((dim_x_data,dim_x_data))
        
        for i in range(dim_x_data):
            for j in range(dim_x_data):
                
                dw = w_data[i] - w_data[j]
                dx = x_data[i] - x_data[j]
                dy = y_data[i] - y_data[j]
                dz = z_data[i] - z_data[j]
                coef = - 0.5 * ( dw**2.0/gph_l1**2.0 + dx**2.0/gph_l2**2.0 + dy**2.0/gph_l3**2.0 + dz**2.0/gph_l4**2.0 ) 
                c1 = dw / gph_l1**2.0 
                c2 = dx / gph_l2**2.0
                c3 = dy / gph_l3**2.0
                c4 = dz / gph_l4**2.0
                k1_1[i,j] = gph_sigma_f * gph_sigma_f * exp(coef)
                k1_2[i,j] = c1*k1_1[i,j]
                k1_3[i,j] = c2*k1_1[i,j]
                k1_4[i,j] = c3*k1_1[i,j]
                k1_5[i,j] = c4*k1_1[i,j]
                k2_1[i,j] = -c1*k1_1[i,j]
                k2_2[i,j] = (1.0/gph_l1**2.0 - c1**2.0)*k1_1[i,j]
                k2_3[i,j] = -c1*c2*k1_1[i,j]
                k2_4[i,j] = -c1*c3*k1_1[i,j]
                k2_5[i,j] = -c1*c4*k1_1[i,j]
                k3_1[i,j] = -c2*k1_1[i,j]
                k3_2[i,j] = -c1*c2*k1_1[i,j]
                k3_3[i,j] = (1.0/gph_l2**2.0 - c2**2.0)*k1_1[i,j]
                k3_4[i,j] = -c2*c3*k1_1[i,j]
                k3_5[i,j] = -c2*c4*k1_1[i,j]
                k4_1[i,j] = -c3*k1_1[i,j]
                k4_2[i,j] = -c1*c3*k1_1[i,j]
                k4_3[i,j] = -c2*c3*k1_1[i,j]
                k4_4[i,j] = (1.0/gph_l3**2.0 - c3**2.0)*k1_1[i,j]
                k4_5[i,j] = -c3*c4*k1_1[i,j]
                k5_1[i,j] = -c4*k1_1[i,j]
                k5_2[i,j] = -c1*c4*k1_1[i,j]
                k5_3[i,j] = -c2*c4*k1_1[i,j]
                k5_4[i,j] = -c3*c4*k1_1[i,j]
                k5_5[i,j] = (1.0/gph_l4**2.0 - c4**2.0)*k1_1[i,j]
        
        
        Kgrad =    np.concatenate((np.concatenate((k1_1,k1_2,k1_3,k1_4,k1_5),axis=1),\
                                   np.concatenate((k2_1,k2_2,k2_3,k2_4,k2_5),axis=1),\
                                   np.concatenate((k3_1,k3_2,k3_3,k3_4,k3_5),axis=1),\
                                   np.concatenate((k4_1,k4_2,k4_3,k4_4,k4_5),axis=1),\
                                   np.concatenate((k5_1,k5_2,k5_3,k5_4,k5_5),axis=1)), axis=0)
        #print(Kgrad[15,16])
        # Joint covariance with function, gradient, and Laplacian
        K = Kgrad + sigma_n*np.identity(len(Kgrad))
        fold = f_df_data[0:5*dim_x_data]
        #print(fold)
        
        n = len(K)
        L=np.linalg.cholesky(K)
        cho = scipy.linalg.cho_factor(K)
        alpha_lin = scipy.linalg.cho_solve(cho, fold)
    #    K_inv = inv(K)
    #    m1 = np.dot(K_inv,fold)
    #    part1 = np.dot(y_data.T,m1)
    #    part2 = math.log(np.linalg.det(K))
        part1 = np.dot(fold.T,alpha_lin)
        part2 = 2.0*sum(np.log(np.diag(L)))
        return(0.5 * part1 + 0.5 * part2 + n*log(2.0*pi)/2.0)
    
        
    def log_likelihood_function_2(hyperparams):
        
        gph_l1 = hyperparams[0]
        gph_l2 = hyperparams[1]
        gph_l3 = hyperparams[2]
        gph_l4 = hyperparams[3]
        
        k1_1 = np.zeros((dim_x_data,dim_x_data))
        k1_2 = np.zeros((dim_x_data,dim_x_data))
        k1_3 = np.zeros((dim_x_data,dim_x_data))
        k1_4 = np.zeros((dim_x_data,dim_x_data))
        k1_5 = np.zeros((dim_x_data,dim_x_data))
        k1_6 = np.zeros((dim_x_data,dim_x_data))
        k1_7 = np.zeros((dim_x_data,dim_x_data))
        k1_8 = np.zeros((dim_x_data,dim_x_data))
        k1_9 = np.zeros((dim_x_data,dim_x_data))
        k1_10 = np.zeros((dim_x_data,dim_x_data))
        k1_11 = np.zeros((dim_x_data,dim_x_data))
        k1_12 = np.zeros((dim_x_data,dim_x_data))
        k1_13 = np.zeros((dim_x_data,dim_x_data))
        k1_14 = np.zeros((dim_x_data,dim_x_data))
        k1_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k2_1 = np.zeros((dim_x_data,dim_x_data))
        k2_2 = np.zeros((dim_x_data,dim_x_data))
        k2_3 = np.zeros((dim_x_data,dim_x_data))
        k2_4 = np.zeros((dim_x_data,dim_x_data))
        k2_5 = np.zeros((dim_x_data,dim_x_data))
        k2_6 = np.zeros((dim_x_data,dim_x_data))
        k2_7 = np.zeros((dim_x_data,dim_x_data))
        k2_8 = np.zeros((dim_x_data,dim_x_data))
        k2_9 = np.zeros((dim_x_data,dim_x_data))
        k2_10 = np.zeros((dim_x_data,dim_x_data))
        k2_11 = np.zeros((dim_x_data,dim_x_data))
        k2_12 = np.zeros((dim_x_data,dim_x_data))
        k2_13 = np.zeros((dim_x_data,dim_x_data))
        k2_14 = np.zeros((dim_x_data,dim_x_data))
        k2_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k3_1 = np.zeros((dim_x_data,dim_x_data))
        k3_2 = np.zeros((dim_x_data,dim_x_data))
        k3_3 = np.zeros((dim_x_data,dim_x_data))
        k3_4 = np.zeros((dim_x_data,dim_x_data))
        k3_5 = np.zeros((dim_x_data,dim_x_data))
        k3_6 = np.zeros((dim_x_data,dim_x_data))
        k3_7 = np.zeros((dim_x_data,dim_x_data))
        k3_8 = np.zeros((dim_x_data,dim_x_data))
        k3_9 = np.zeros((dim_x_data,dim_x_data))
        k3_10 = np.zeros((dim_x_data,dim_x_data))
        k3_11 = np.zeros((dim_x_data,dim_x_data))
        k3_12 = np.zeros((dim_x_data,dim_x_data))
        k3_13 = np.zeros((dim_x_data,dim_x_data))
        k3_14 = np.zeros((dim_x_data,dim_x_data))
        k3_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k4_1 = np.zeros((dim_x_data,dim_x_data))
        k4_2 = np.zeros((dim_x_data,dim_x_data))
        k4_3 = np.zeros((dim_x_data,dim_x_data))
        k4_4 = np.zeros((dim_x_data,dim_x_data))
        k4_5 = np.zeros((dim_x_data,dim_x_data))
        k4_6 = np.zeros((dim_x_data,dim_x_data))
        k4_7 = np.zeros((dim_x_data,dim_x_data))
        k4_8 = np.zeros((dim_x_data,dim_x_data))
        k4_9 = np.zeros((dim_x_data,dim_x_data))
        k4_10 = np.zeros((dim_x_data,dim_x_data))
        k4_11 = np.zeros((dim_x_data,dim_x_data))
        k4_12 = np.zeros((dim_x_data,dim_x_data))
        k4_13 = np.zeros((dim_x_data,dim_x_data))
        k4_14 = np.zeros((dim_x_data,dim_x_data))
        k4_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k5_1 = np.zeros((dim_x_data,dim_x_data))
        k5_2 = np.zeros((dim_x_data,dim_x_data))
        k5_3 = np.zeros((dim_x_data,dim_x_data))
        k5_4 = np.zeros((dim_x_data,dim_x_data))
        k5_5 = np.zeros((dim_x_data,dim_x_data))
        k5_6 = np.zeros((dim_x_data,dim_x_data))
        k5_7 = np.zeros((dim_x_data,dim_x_data))
        k5_8 = np.zeros((dim_x_data,dim_x_data))
        k5_9 = np.zeros((dim_x_data,dim_x_data))
        k5_10 = np.zeros((dim_x_data,dim_x_data))
        k5_11 = np.zeros((dim_x_data,dim_x_data))
        k5_12 = np.zeros((dim_x_data,dim_x_data))
        k5_13 = np.zeros((dim_x_data,dim_x_data))
        k5_14 = np.zeros((dim_x_data,dim_x_data))
        k5_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k6_1 = np.zeros((dim_x_data,dim_x_data))
        k6_2 = np.zeros((dim_x_data,dim_x_data))
        k6_3 = np.zeros((dim_x_data,dim_x_data))
        k6_4 = np.zeros((dim_x_data,dim_x_data))
        k6_5 = np.zeros((dim_x_data,dim_x_data))
        k6_6 = np.zeros((dim_x_data,dim_x_data))
        k6_7 = np.zeros((dim_x_data,dim_x_data))
        k6_8 = np.zeros((dim_x_data,dim_x_data))
        k6_9 = np.zeros((dim_x_data,dim_x_data))
        k6_10 = np.zeros((dim_x_data,dim_x_data))
        k6_11 = np.zeros((dim_x_data,dim_x_data))
        k6_12 = np.zeros((dim_x_data,dim_x_data))
        k6_13 = np.zeros((dim_x_data,dim_x_data))
        k6_14 = np.zeros((dim_x_data,dim_x_data))
        k6_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k7_1 = np.zeros((dim_x_data,dim_x_data))
        k7_2 = np.zeros((dim_x_data,dim_x_data))
        k7_3 = np.zeros((dim_x_data,dim_x_data))
        k7_4 = np.zeros((dim_x_data,dim_x_data))
        k7_5 = np.zeros((dim_x_data,dim_x_data))
        k7_6 = np.zeros((dim_x_data,dim_x_data))
        k7_7 = np.zeros((dim_x_data,dim_x_data))
        k7_8 = np.zeros((dim_x_data,dim_x_data))
        k7_9 = np.zeros((dim_x_data,dim_x_data))
        k7_10 = np.zeros((dim_x_data,dim_x_data))
        k7_11 = np.zeros((dim_x_data,dim_x_data))
        k7_12 = np.zeros((dim_x_data,dim_x_data))
        k7_13 = np.zeros((dim_x_data,dim_x_data))
        k7_14 = np.zeros((dim_x_data,dim_x_data))
        k7_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k8_1 = np.zeros((dim_x_data,dim_x_data))
        k8_2 = np.zeros((dim_x_data,dim_x_data))
        k8_3 = np.zeros((dim_x_data,dim_x_data))
        k8_4 = np.zeros((dim_x_data,dim_x_data))
        k8_5 = np.zeros((dim_x_data,dim_x_data))
        k8_6 = np.zeros((dim_x_data,dim_x_data))
        k8_7 = np.zeros((dim_x_data,dim_x_data))
        k8_8 = np.zeros((dim_x_data,dim_x_data))
        k8_9 = np.zeros((dim_x_data,dim_x_data))
        k8_10 = np.zeros((dim_x_data,dim_x_data))
        k8_11 = np.zeros((dim_x_data,dim_x_data))
        k8_12 = np.zeros((dim_x_data,dim_x_data))
        k8_13 = np.zeros((dim_x_data,dim_x_data))
        k8_14 = np.zeros((dim_x_data,dim_x_data))
        k8_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k9_1 = np.zeros((dim_x_data,dim_x_data))
        k9_2 = np.zeros((dim_x_data,dim_x_data))
        k9_3 = np.zeros((dim_x_data,dim_x_data))
        k9_4 = np.zeros((dim_x_data,dim_x_data))
        k9_5 = np.zeros((dim_x_data,dim_x_data))
        k9_6 = np.zeros((dim_x_data,dim_x_data))
        k9_7 = np.zeros((dim_x_data,dim_x_data))
        k9_8 = np.zeros((dim_x_data,dim_x_data))
        k9_9 = np.zeros((dim_x_data,dim_x_data))
        k9_10 = np.zeros((dim_x_data,dim_x_data))
        k9_11 = np.zeros((dim_x_data,dim_x_data))
        k9_12 = np.zeros((dim_x_data,dim_x_data))
        k9_13 = np.zeros((dim_x_data,dim_x_data))
        k9_14 = np.zeros((dim_x_data,dim_x_data))
        k9_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k10_1 = np.zeros((dim_x_data,dim_x_data))
        k10_2 = np.zeros((dim_x_data,dim_x_data))
        k10_3 = np.zeros((dim_x_data,dim_x_data))
        k10_4 = np.zeros((dim_x_data,dim_x_data))
        k10_5 = np.zeros((dim_x_data,dim_x_data))
        k10_6 = np.zeros((dim_x_data,dim_x_data))
        k10_7 = np.zeros((dim_x_data,dim_x_data))
        k10_8 = np.zeros((dim_x_data,dim_x_data))
        k10_9 = np.zeros((dim_x_data,dim_x_data))
        k10_10 = np.zeros((dim_x_data,dim_x_data))
        k10_11 = np.zeros((dim_x_data,dim_x_data))
        k10_12 = np.zeros((dim_x_data,dim_x_data))
        k10_13 = np.zeros((dim_x_data,dim_x_data))
        k10_14 = np.zeros((dim_x_data,dim_x_data))
        k10_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k11_1 = np.zeros((dim_x_data,dim_x_data))
        k11_2 = np.zeros((dim_x_data,dim_x_data))
        k11_3 = np.zeros((dim_x_data,dim_x_data))
        k11_4 = np.zeros((dim_x_data,dim_x_data))
        k11_5 = np.zeros((dim_x_data,dim_x_data))
        k11_6 = np.zeros((dim_x_data,dim_x_data))
        k11_7 = np.zeros((dim_x_data,dim_x_data))
        k11_8 = np.zeros((dim_x_data,dim_x_data))
        k11_9 = np.zeros((dim_x_data,dim_x_data))
        k11_10 = np.zeros((dim_x_data,dim_x_data))
        k11_11 = np.zeros((dim_x_data,dim_x_data))
        k11_12 = np.zeros((dim_x_data,dim_x_data))
        k11_13 = np.zeros((dim_x_data,dim_x_data))
        k11_14 = np.zeros((dim_x_data,dim_x_data))
        k11_15 = np.zeros((dim_x_data,dim_x_data))
        #
        k12_1 = np.zeros((dim_x_data,dim_x_data))
        k12_2 = np.zeros((dim_x_data,dim_x_data))
        k12_3 = np.zeros((dim_x_data,dim_x_data))
        k12_4 = np.zeros((dim_x_data,dim_x_data))
        k12_5 = np.zeros((dim_x_data,dim_x_data))
        k12_6 = np.zeros((dim_x_data,dim_x_data))
        k12_7 = np.zeros((dim_x_data,dim_x_data))
        k12_8 = np.zeros((dim_x_data,dim_x_data))
        k12_9 = np.zeros((dim_x_data,dim_x_data))
        k12_10 = np.zeros((dim_x_data,dim_x_data))
        k12_11 = np.zeros((dim_x_data,dim_x_data))
        k12_12 = np.zeros((dim_x_data,dim_x_data))
        k12_13 = np.zeros((dim_x_data,dim_x_data))
        k12_14 = np.zeros((dim_x_data,dim_x_data))
        k12_15 = np.zeros((dim_x_data,dim_x_data))
         #
        k13_1 = np.zeros((dim_x_data,dim_x_data))
        k13_2 = np.zeros((dim_x_data,dim_x_data))
        k13_3 = np.zeros((dim_x_data,dim_x_data))
        k13_4 = np.zeros((dim_x_data,dim_x_data))
        k13_5 = np.zeros((dim_x_data,dim_x_data))
        k13_6 = np.zeros((dim_x_data,dim_x_data))
        k13_7 = np.zeros((dim_x_data,dim_x_data))
        k13_8 = np.zeros((dim_x_data,dim_x_data))
        k13_9 = np.zeros((dim_x_data,dim_x_data))
        k13_10 = np.zeros((dim_x_data,dim_x_data))
        k13_11 = np.zeros((dim_x_data,dim_x_data))
        k13_12 = np.zeros((dim_x_data,dim_x_data))
        k13_13 = np.zeros((dim_x_data,dim_x_data))
        k13_14 = np.zeros((dim_x_data,dim_x_data))
        k13_15 = np.zeros((dim_x_data,dim_x_data))
         #
        k14_1 = np.zeros((dim_x_data,dim_x_data))
        k14_2 = np.zeros((dim_x_data,dim_x_data))
        k14_3 = np.zeros((dim_x_data,dim_x_data))
        k14_4 = np.zeros((dim_x_data,dim_x_data))
        k14_5 = np.zeros((dim_x_data,dim_x_data))
        k14_6 = np.zeros((dim_x_data,dim_x_data))
        k14_7 = np.zeros((dim_x_data,dim_x_data))
        k14_8 = np.zeros((dim_x_data,dim_x_data))
        k14_9 = np.zeros((dim_x_data,dim_x_data))
        k14_10 = np.zeros((dim_x_data,dim_x_data))
        k14_11 = np.zeros((dim_x_data,dim_x_data))
        k14_12 = np.zeros((dim_x_data,dim_x_data))
        k14_13 = np.zeros((dim_x_data,dim_x_data))
        k14_14 = np.zeros((dim_x_data,dim_x_data))
        k14_15 = np.zeros((dim_x_data,dim_x_data))
         #
        k15_1 = np.zeros((dim_x_data,dim_x_data))
        k15_2 = np.zeros((dim_x_data,dim_x_data))
        k15_3 = np.zeros((dim_x_data,dim_x_data))
        k15_4 = np.zeros((dim_x_data,dim_x_data))
        k15_5 = np.zeros((dim_x_data,dim_x_data))
        k15_6 = np.zeros((dim_x_data,dim_x_data))
        k15_7 = np.zeros((dim_x_data,dim_x_data))
        k15_8 = np.zeros((dim_x_data,dim_x_data))
        k15_9 = np.zeros((dim_x_data,dim_x_data))
        k15_10 = np.zeros((dim_x_data,dim_x_data))
        k15_11 = np.zeros((dim_x_data,dim_x_data))
        k15_12 = np.zeros((dim_x_data,dim_x_data))
        k15_13 = np.zeros((dim_x_data,dim_x_data))
        k15_14 = np.zeros((dim_x_data,dim_x_data))
        k15_15 = np.zeros((dim_x_data,dim_x_data))
        
        for i in range(dim_x_data):
            for j in range(dim_x_data):
                if (i==j):
                    sigma_tmp = 0.0
                else:
                    sigma_tmp = 0.0
                
                dw = w_data[i] - w_data[j]
                dx = x_data[i] - x_data[j]
                dy = y_data[i] - y_data[j]
                dz = z_data[i] - z_data[j]
                coef = gph_sigma_f * gph_sigma_f * exp( - 0.5 * ( dw**2.0/gph_l1**2.0 + dx**2.0/gph_l2**2.0 + dy**2.0/gph_l3**2.0 + dz**2.0/gph_l4**2.0 ) )
                c1 = dw / gph_l1**2.0 
                c2 = dx / gph_l2**2.0
                c3 = dy / gph_l3**2.0
                c4 = dz / gph_l4**2.0
                #print('c1=',c1,'c2=',c2,'c3=',c3,'c4=',c4)
                # Effectively creates an upper triangular matrix
                k1_1[i,j] = coef + sigma_tmp
                k1_2[i,j] = c1*k1_1[i,j] + sigma_tmp
                k1_3[i,j] = c2*k1_1[i,j] + sigma_tmp
                k1_4[i,j] = c3*k1_1[i,j] + sigma_tmp
                k1_5[i,j] = c4*k1_1[i,j] + sigma_tmp
                k1_6[i,j] = (c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k1_7[i,j] = c1*c2*k1_1[i,j] + sigma_tmp
                k1_8[i,j] = c1*c3*k1_1[i,j] + sigma_tmp
                k1_9[i,j] = c1*c4*k1_1[i,j] + sigma_tmp
                k1_10[i,j] = (c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k1_11[i,j] = c2*c3*k1_1[i,j] + sigma_tmp
                k1_12[i,j] = c2*c4*k1_1[i,j] + sigma_tmp
                k1_13[i,j] = (c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k1_14[i,j] = c3*c4*k1_1[i,j] + sigma_tmp
                k1_15[i,j] = (c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k2_1[i,j] = -c1*k1_1[i,j] + sigma_tmp
                k2_2[i,j] = (1.0/gph_l1**2.0 - c1**2.0)*k1_1[i,j] + sigma_tmp
                k2_3[i,j] = -c1*c2*k1_1[i,j] + sigma_tmp
                k2_4[i,j] = -c1*c3*k1_1[i,j] + sigma_tmp
                k2_5[i,j] = -c1*c4*k1_1[i,j] + sigma_tmp
                k2_6[i,j] = -c1*(c1**2.0 - 3.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k2_7[i,j] = -c2*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k2_8[i,j] = -c3*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k2_9[i,j] = -c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k2_10[i,j]= -c1*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k2_11[i,j]= -c1*c2*c3*k1_1[i,j] + sigma_tmp
                k2_12[i,j]= -c1*c2*c4*k1_1[i,j] + sigma_tmp
                k2_13[i,j]= -c1*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k2_14[i,j]= -c1*c3*c4*k1_1[i,j] + sigma_tmp
                k2_15[i,j]= -c1*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k3_1[i,j] = -c2*k1_1[i,j] + sigma_tmp
                k3_2[i,j] = -c1*c2*k1_1[i,j] + sigma_tmp
                k3_3[i,j] = (1.0/gph_l2**2.0 - c2**2.0)*k1_1[i,j] + sigma_tmp
                k3_4[i,j] = -c2*c3*k1_1[i,j] + sigma_tmp
                k3_5[i,j] = -c2*c4*k1_1[i,j] + sigma_tmp
                k3_6[i,j] = -c2*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k3_7[i,j] = -c1*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k3_8[i,j] = -c1*c2*c3*k1_1[i,j] + sigma_tmp
                k3_9[i,j] = -c1*c2*c4*k1_1[i,j] + sigma_tmp
                k3_10[i,j] =-c2*(c2**2.0 - 3.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k3_11[i,j] =-c3*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k3_12[i,j] =-c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k3_13[i,j] =-c2*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k3_14[i,j] =-c2*c3*c4*k1_1[i,j] + sigma_tmp
                k3_15[i,j] =-c2*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k4_1[i,j] = -c3*k1_1[i,j] + sigma_tmp
                k4_2[i,j] = -c1*c3*k1_1[i,j] + sigma_tmp
                k4_3[i,j] = -c2*c3*k1_1[i,j] + sigma_tmp
                k4_4[i,j] = (1.0/gph_l3**2.0 - c3**2.0)*k1_1[i,j] + sigma_tmp
                k4_5[i,j] = -c3*c4*k1_1[i,j] + sigma_tmp
                k4_6[i,j] = -c3*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k4_7[i,j] = -c1*c2*c3*k1_1[i,j] + sigma_tmp
                k4_8[i,j] = -c1*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k4_9[i,j] = -c1*c3*c4*k1_1[i,j] + sigma_tmp
                k4_10[i,j]= -c3*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k4_11[i,j]= -c2*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k4_12[i,j]= -c2*c3*c4*k1_1[i,j] + sigma_tmp
                k4_13[i,j]= -c3*(c3**2.0 - 3.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k4_14[i,j]= -c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k4_15[i,j]= -c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k5_1[i,j] = -c4*k1_1[i,j] + sigma_tmp
                k5_2[i,j] = -c1*c4*k1_1[i,j] + sigma_tmp
                k5_3[i,j] = -c2*c4*k1_1[i,j] + sigma_tmp
                k5_4[i,j] = -c3*c4*k1_1[i,j] + sigma_tmp
                k5_5[i,j] = (1.0/gph_l4**2.0 - c4**2.0)*k1_1[i,j] + sigma_tmp
                k5_6[i,j] = -c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k5_7[i,j]= -c1*c2*c4*k1_1[i,j] + sigma_tmp
                k5_8[i,j] = -c1*c3*c4*k1_1[i,j] + sigma_tmp
                k5_9[i,j] = -c1*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k5_10[i,j]= -c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k5_11[i,j]= -c2*c3*c4*k1_1[i,j] + sigma_tmp
                k5_12[i,j]= -c2*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k5_13[i,j]= -c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k5_14[i,j]= -c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k5_15[i,j]= -c4*(c4**2.0 - 3.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k6_1[i,j] =       (c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_2[i,j] =    c1*(c1**2.0 - 3.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_3[i,j] =    c2*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_4[i,j] =    c3*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_5[i,j] =    c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_6[i,j] = (c1**4.0 - 6.0*c1**2.0/gph_l1**2.0 + 3.0/gph_l1**4.0)*k1_1[i,j] + sigma_tmp
                k6_7[i,j] = c1*c2*(c1**2.0 - 3.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_8[i,j] = c1*c3*(c1**2.0 - 3.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_9[i,j] = c1*c4*(c1**2.0 - 3.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_10[i,j]=       (c1**2.0 - 1.0/gph_l1**2.0)*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k6_11[i,j]= c2*c3*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_12[i,j]= c2*c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_13[i,j]=       (c1**2.0 - 1.0/gph_l1**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k6_14[i,j]= c3*c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k6_15[i,j]=       (c1**2.0 - 1.0/gph_l1**2.0)*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k7_1[i,j] = c1*c2*k1_1[i,j] + sigma_tmp
                k7_2[i,j] = c2*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k7_3[i,j] = c1*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k7_4[i,j] = c1*c2*c3*k1_1[i,j] + sigma_tmp
                k7_5[i,j] = c1*c2*c4*k1_1[i,j] + sigma_tmp
                k7_6[i,j] = c1*c2*(c1**2.0 - 3.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k7_7[i,j] =       (c1**2.0 - 1.0/gph_l1**2.0)*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k7_8[i,j] = c2*c3*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k7_9[i,j] = c2*c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k7_10[i,j]= c1*c2*(c2**2.0 - 3.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k7_11[i,j]= c1*c3*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k7_12[i,j]= c1*c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k7_13[i,j]= c1*c2*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k7_14[i,j]= c1*c2*c3*c4*k1_1[i,j] + sigma_tmp
                k7_15[i,j]= c1*c2*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k8_1[i,j] = c1*c3*k1_1[i,j] + sigma_tmp
                k8_2[i,j] = c3*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k8_3[i,j] = c1*c2*c3*k1_1[i,j] + sigma_tmp
                k8_4[i,j] = c1*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k8_5[i,j] = c1*c3*c4*k1_1[i,j] + sigma_tmp
                k8_6[i,j] = c1*c3*(c1**2.0 - 3.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k8_7[i,j] = c2*c3*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k8_8[i,j] =       (c1**2.0 - 1.0/gph_l1**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k8_9[i,j] = c3*c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k8_10[i,j]= c1*c3*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k8_11[i,j]= c1*c2*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k8_12[i,j]= c1*c2*c3*c4*k1_1[i,j] + sigma_tmp
                k8_13[i,j]= c1*c3*(c3**2.0 - 3.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k8_14[i,j]= c1*c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k8_15[i,j]= c1*c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k9_1[i,j] = c1*c4*k1_1[i,j] + sigma_tmp
                k9_2[i,j] = c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k9_3[i,j] = c1*c2*c4*k1_1[i,j] + sigma_tmp
                k9_4[i,j] = c1*c3*c4*k1_1[i,j] + sigma_tmp
                k9_5[i,j] = c1*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k9_6[i,j] = c1*c4*(c1**2.0 - 3.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k9_7[i,j] = c2*c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k9_8[i,j] = c3*c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k9_9[i,j] =       (c4**2.0 - 1.0/gph_l4**2.0)*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k9_10[i,j]= c1*c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k9_11[i,j]= c1*c2*c3*c4*k1_1[i,j] + sigma_tmp
                k9_12[i,j]= c1*c2*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k9_13[i,j]= c1*c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k9_14[i,j]= c1*c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k9_15[i,j]= c1*c4*(c4**2.0 - 3.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k10_1[i,j] =       (c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_2[i,j] =    c1*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_3[i,j] =    c2*(c2**2.0 - 3.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_4[i,j] =    c3*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_5[i,j] =    c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_6[i,j] =       (c2**2.0 - 1.0/gph_l2**2.0)*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k10_7[i,j] = c2*c1*(c2**2.0 - 3.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_8[i,j] = c3*c1*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_9[i,j] = c4*c1*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_10[i,j]= (c2**4.0 - 6.0*c2**2.0/gph_l2**2.0 + 3.0/gph_l2**4.0)*k1_1[i,j] + sigma_tmp
                k10_11[i,j]= c2*c3*(c2**2.0 - 3.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_12[i,j]= c2*c4*(c2**2.0 - 3.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_13[i,j]=       (c2**2.0 - 1.0/gph_l2**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k10_14[i,j]= c3*c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k10_15[i,j]=       (c2**2.0 - 1.0/gph_l2**2.0)*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k11_1[i,j] = c2*c3*k1_1[i,j] + sigma_tmp
                k11_2[i,j] = c1*c2*c3*k1_1[i,j] + sigma_tmp
                k11_3[i,j] = c3*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k11_4[i,j] = c2*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k11_5[i,j] = c2*c3*c4*k1_1[i,j] + sigma_tmp
                k11_6[i,j] = c2*c3*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k11_7[i,j] = c1*c3*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k11_8[i,j] = c1*c2*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k11_9[i,j] = c1*c2*c3*c4*k1_1[i,j] + sigma_tmp
                k11_10[i,j]= c2*c3*(c2**2.0 - 3.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k11_11[i,j]=       (c2**2.0 - 1.0/gph_l2**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k11_12[i,j]= c3*c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k11_13[i,j]= c2*c3*(c3**2.0 - 3.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k11_14[i,j]= c2*c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k11_15[i,j]= c2*c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                #
                k12_1[i,j] = c2*c4*k1_1[i,j] + sigma_tmp
                k12_2[i,j] = c1*c2*c4*k1_1[i,j] + sigma_tmp
                k12_3[i,j] = c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k12_4[i,j] = c2*c3*c4*k1_1[i,j] + sigma_tmp
                k12_5[i,j] = c2*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k12_6[i,j] = c2*c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k12_7[i,j] = c1*c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k12_8[i,j] = c1*c2*c3*c4*k1_1[i,j] + sigma_tmp
                k12_9[i,j] =  c1*c2*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k12_10[i,j] = c2*c4*(c2**2.0 - 3.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k12_11[i,j] = c3*c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k12_12[i,j] =       (c4**2.0 - 1.0/gph_l4**2.0)*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k12_13[i,j] = c2*c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k12_14[i,j] = c2*c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k12_15[i,j] = c2*c4*(c4**2.0 - 3.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                 #
                k13_1[i,j] =        (c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_2[i,j] =     c1*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_3[i,j] =     c2*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_4[i,j] =     c3*(c3**2.0 - 3.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_5[i,j] =     c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_6[i,j] =        (c3**2.0 - 1.0/gph_l3**2.0)*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k13_7[i,j] =  c1*c2*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_8[i,j] =  c1*c3*(c3**2.0 - 3.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_9[i,j] =  c1*c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_10[i,j] =       (c3**2.0 - 1.0/gph_l3**2.0)*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k13_11[i,j] = c2*c3*(c3**2.0 - 3.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_12[i,j] = c2*c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_13[i,j] = (c3**4.0 - 6.0*c3**2.0/gph_l3**2.0 + 3.0/gph_l3**4.0)*k1_1[i,j] + sigma_tmp
                k13_14[i,j] = c3*c4*(c3**2.0 - 3.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k13_15[i,j] =       (c3**2.0 - 1.0/gph_l3**2.0)*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                 #
                k14_1[i,j] =    c3*c4*k1_1[i,j] + sigma_tmp
                k14_2[i,j] = c1*c3*c4*k1_1[i,j] + sigma_tmp
                k14_3[i,j] = c2*c3*c4*k1_1[i,j] + sigma_tmp
                k14_4[i,j] =    c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k14_5[i,j] =    c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k14_6[i,j] = c3*c4*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k14_7[i,j] = c1*c2*c3*c4*k1_1[i,j] + sigma_tmp
                k14_8[i,j] = c1*c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k14_9[i,j] = c1*c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k14_10[i,j]= c3*c4*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k14_11[i,j]= c2*c4*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k14_12[i,j]= c2*c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k14_13[i,j]= c3*c4*(c3**2.0 - 3.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k14_14[i,j]=       (c4**2.0 - 1.0/gph_l4**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k14_15[i,j]= c3*c4*(c4**2.0 - 3.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                 #
                k15_1[i,j] =       (c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_2[i,j] =    c1*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_3[i,j] =    c2*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_4[i,j] =    c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_5[i,j] =    c4*(c4**2.0 - 3.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_6[i,j] =       (c4**2.0 - 1.0/gph_l4**2.0)*(c1**2.0 - 1.0/gph_l1**2.0)*k1_1[i,j] + sigma_tmp
                k15_7[i,j] = c1*c2*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_8[i,j] = c1*c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_9[i,j] = c1*c4*(c4**2.0 - 3.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_10[i,j] =      (c4**2.0 - 1.0/gph_l4**2.0)*(c2**2.0 - 1.0/gph_l2**2.0)*k1_1[i,j] + sigma_tmp
                k15_11[i,j] =c2*c3*(c4**2.0 - 1.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_12[i,j] =c2*c4*(c4**2.0 - 3.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_13[i,j] =      (c4**2.0 - 1.0/gph_l4**2.0)*(c3**2.0 - 1.0/gph_l3**2.0)*k1_1[i,j] + sigma_tmp
                k15_14[i,j] =c3*c4*(c4**2.0 - 3.0/gph_l4**2.0)*k1_1[i,j] + sigma_tmp
                k15_15[i,j] = (c4**4.0 - 6.0*c4**2.0/gph_l4**2.0 + 3.0/gph_l4**4.0)*k1_1[i,j] + sigma_tmp
        
        # By transposing what we have already done we create a symmetric matrix
    #    plt.figure(1)
    #    #plt.rc('text', usetex=True)
    #    #plt.rc('font', family='serif')
    #    plt.contourf(k13_13, 100, cmap='rainbow')
    #    plt.axis('equal')
    #    plt.colorbar()
    #    plt.grid()
    #    print(k10_13-k13_10.T)
    #    print(k1_1[0,1])
    #    
    #    plt.figure(2)
    #    #plt.rc('text', usetex=True)
    #    #plt.rc('font', family='serif')
    #    plt.contourf(k13_15, 100, cmap='rainbow')
    #    plt.colorbar()
    #    plt.grid()
            
        
        Kgradlap = np.concatenate((np.concatenate((k1_1, k1_2, k1_3, k1_4, k1_5, k1_6, k1_7, k1_8, k1_9, k1_10, k1_11, k1_12, k1_13, k1_14, k1_15),axis=1),\
                                   np.concatenate((k2_1, k2_2, k2_3, k2_4, k2_5, k2_6, k2_7, k2_8, k2_9, k2_10, k2_11, k2_12, k2_13, k2_14, k2_15),axis=1),\
                                   np.concatenate((k3_1, k3_2, k3_3, k3_4, k3_5, k3_6, k3_7, k3_8, k3_9, k3_10, k3_11, k3_12, k3_13, k3_14, k3_15),axis=1),\
                                   np.concatenate((k4_1, k4_2, k4_3, k4_4, k4_5, k4_6, k4_7, k4_8, k4_9, k4_10, k4_11, k4_12, k4_13, k4_14, k4_15),axis=1),\
                                   np.concatenate((k5_1, k5_2, k5_3, k5_4, k5_5, k5_6, k5_7, k5_8, k5_9, k5_10, k5_11, k5_12, k5_13, k5_14, k5_15),axis=1),\
                                   np.concatenate((k6_1, k6_2, k6_3, k6_4, k6_5, k6_6, k6_7, k6_8, k6_9, k6_10, k6_11, k6_12, k6_13, k6_14, k6_15),axis=1),\
                                   np.concatenate((k7_1, k7_2, k7_3, k7_4, k7_5, k7_6, k7_7, k7_8, k7_9, k7_10, k7_11, k7_12, k7_13, k7_14, k7_15),axis=1),\
                                   np.concatenate((k8_1, k8_2, k8_3, k8_4, k8_5, k8_6, k8_7, k8_8, k8_9, k8_10, k8_11, k8_12, k8_13, k8_14, k8_15),axis=1),\
                                   np.concatenate((k9_1, k9_2, k9_3, k9_4, k9_5, k9_6, k9_7, k9_8, k9_9, k9_10, k9_11, k9_12, k9_13, k9_14, k9_15),axis=1),\
                                   np.concatenate((k10_1,k10_2,k10_3,k10_4,k10_5,k10_6,k10_7,k10_8,k10_9,k10_10,k10_11,k10_12,k10_13,k10_14,k10_15),axis=1),\
                                   np.concatenate((k11_1,k11_2,k11_3,k11_4,k11_5,k11_6,k11_7,k11_8,k11_9,k11_10,k11_11,k11_12,k11_13,k11_14,k11_15),axis=1),\
                                   np.concatenate((k12_1,k12_2,k12_3,k12_4,k12_5,k12_6,k12_7,k12_8,k12_9,k12_10,k12_11,k12_12,k12_13,k12_14,k12_15),axis=1),\
                                   np.concatenate((k13_1,k13_2,k13_3,k13_4,k13_5,k13_6,k13_7,k13_8,k13_9,k13_10,k13_11,k13_12,k13_13,k13_14,k13_15),axis=1),\
                                   np.concatenate((k14_1,k14_2,k14_3,k14_4,k14_5,k14_6,k14_7,k14_8,k14_9,k14_10,k14_11,k14_12,k14_13,k14_14,k14_15),axis=1),\
                                   np.concatenate((k15_1,k15_2,k15_3,k15_4,k15_5,k15_6,k15_7,k15_8,k15_9,k15_10,k15_11,k15_12,k15_13,k15_14,k15_15),axis=1)), axis=0)
       
    
        K = Kgradlap + sigma_n*np.identity(len(Kgradlap))
        fold = f_df_data
        #print('postsym',K[15,16])
        n = len(K)
        L=np.linalg.cholesky(K)
        cho = scipy.linalg.cho_factor(K)
        alpha_lin = scipy.linalg.cho_solve(cho, fold)
    #    #K_inv = inv(K)
    #    #m1 = np.dot(K_inv,fold)
    #    #part1 = np.dot(y_data.T,m1)
    #    #part2 = math.log(np.linalg.det(K))
        part1 = np.dot(fold.T,alpha_lin)
        part2 = 2.0*sum(np.log(np.diag(L)))
        return(0.5 * part1 + 0.5 * part2 + n*log(2.0*pi)/2.0)
        #return(K)
    
    dim_x_data = len(x_data)
    
    res0 = minimize(log_likelihood_function_0, hyperparams_in, method='Nelder-Mead', tol=1e-6)
    hyperparams_out_0 = res0.x
    [K0tmp,K1tmp,K2tmp] = BuildCovarMatrix(hyperparams_out_0)
    K0 = K0tmp
    [fnew0_tmp,fnew1_tmp,fnew2_tmp]  = GP_posterior(hyperparams_out_0,K0,K1tmp,K2tmp)
    fnew0=fnew0_tmp
       
    res1 = minimize(log_likelihood_function_1, hyperparams_in, method='Nelder-Mead', tol=1e-6)
    hyperparams_out_1 = res1.x
    [K0tmp,K1tmp,K2tmp] = BuildCovarMatrix(hyperparams_out_1)
    K1 = K1tmp
    [fnew0_tmp,fnew1_tmp,fnew2_tmp]  = GP_posterior(hyperparams_out_1,K0tmp,K1,K2tmp)
    fnew1=fnew1_tmp
    
    res2 = minimize(log_likelihood_function_2, hyperparams_in, method='Nelder-Mead', tol=1e-6)
    hyperparams_out_2 = res2.x
    [K0tmp,K1tmp,K2tmp] = BuildCovarMatrix(hyperparams_out_2)
    K2 = K2tmp
    [fnew0_tmp,fnew1_tmp,fnew2_tmp]  = GP_posterior(hyperparams_out_2,K0tmp,K1tmp,K2)
    fnew2=fnew2_tmp
    
    return(hyperparams_out_0,K0,fnew0,hyperparams_out_1,K1,fnew1,hyperparams_out_2,K2,fnew2)

#def EOS_MLL_Solver_test():

w_data = (1.0e2*(1.0 + 0.3*np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0])))
x_data = (2.0 + 0.4*np.array([1.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0]))
y_data = (0.75 + 0.5*np.array([1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0]))
z_data = (0.3 + 0.2*np.array([1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0]))


#w_star = ([np.mean(w_data)])
#x_star = ([np.mean(x_data)])
#y_star = ([np.mean(y_data)])
#z_star = ([np.mean(z_data)])
w_star = [w_data[0]]
x_star = [x_data[0]]
y_star = [y_data[0]]
z_star = [z_data[0]]

f_df_data = TestFuncandDerivs(w_data,x_data,y_data,z_data)
f_df_star = TestFuncandDerivs(w_star,x_star,y_star,z_star)


gph_sigma_f = 1.0
sigma_n = 1.0e-10
hyperparams_in=[1.0*(np.max(w_data)-np.min(w_data)),1.0*(np.max(x_data)-np.min(x_data)),1.0*(np.max(y_data)-np.min(y_data)),1.0*(np.max(z_data)-np.min(z_data))]
print(hyperparams_in)

[hyperparams_out_0,K0,fnew0,hyperparams_out_1,K1,fnew1,\
 hyperparams_out_2,K2,fnew2]=EOS_MLL_Solver((w_data),(x_data),(y_data),(z_data),\
 (w_star),(x_star),(y_star),(z_star),(f_df_data),gph_sigma_f,sigma_n,hyperparams_in)


print('Test Function Results')
print('Mode 0:')
print('Hyperparameters',hyperparams_out_0)
print('Error0:',sqrt((f_df_star[0]-fnew0)**2.0)/fnew0)
print('#####')
print('Answer = ',f_df_star[0])
print('Interpolant =',fnew0)
print('#########################################################')
print('Mode 1:')
print('Hyperparameters',hyperparams_out_1)
print('Error1:',np.mean(np.sqrt((f_df_star[0:5]-fnew1)**2.0)/fnew1))
print('#####')
print('Answer = ',f_df_star[0:5])
print('Interpolant =',fnew1)
print('#########################################################')
print('Mode 2:')
print('Hyperparameters',hyperparams_out_2)
print('Error2:',np.mean(np.sqrt((f_df_star-fnew2)**2.0)/fnew2))
print('#####')
print('Answer = ',f_df_star)
print('Interpolant =',fnew2)

plt.figure(1)
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
plt.contourf(K0, 100, cmap='rainbow')
plt.axis('equal')
plt.colorbar()
plt.grid()
#    plt.xticks(np.arange(0, 240, step=15))
#    plt.yticks(np.arange(0, 240, step=15))
#    L=np.linalg.cholesky(K)
#    print('Is K symmetric?:',np.allclose(K,K.T))
#    print('Are all eigenvalues of K greater than zero?:',np.all(np.linalg.eigvals(K) >= 0))
#    print('Does Cholesky-Decomposition work?:',np.allclose(K,np.dot(L,L.T)))

plt.figure(2)
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
plt.contourf(K1, 100, cmap='rainbow')
plt.axis('equal')
plt.colorbar()
plt.grid()
#    plt.xticks(np.arange(0, 240, step=15))
#    plt.yticks(np.arange(0, 240, step=15))

plt.figure(3)
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
plt.contourf(K2, 100, cmap='rainbow')
plt.axis('equal')
plt.colorbar()
plt.grid()
plt.xticks(np.arange(0, 240, step=15))
plt.yticks(np.arange(0, 240, step=15))
#print(len(K2))
#print(K1)
#    return()

    
#EOS_MLL_Solver_test()

 
#gridsize=16
#l1_range=np.linspace(1.0, 3.0, gridsize)
#l2_range=np.linspace(1.0, 3.0, gridsize)
#l3_range=np.linspace(1.0, 3.0, gridsize)
#l4_range=np.linspace(1.0, 3.0, gridsize)
#
#
#l1_mesh, l2_mesh, l3_mesh,l4_mesh = np.meshgrid(l1_range,l2_range,l3_range,l4_range)
#
#MLL_0 = np.zeros_like(l1_mesh)
#MLL_1 = np.zeros_like(l1_mesh)
#MLL_2 = np.zeros_like(l1_mesh)
#
#
#for i in range(0,len(l1_range)):
#    #print('i=',i)
#    for j in range(0,len(l2_range)):
#        for k in range(0,len(l3_range)):
#            for l in range(0,len(l4_range)):
#                print('i=',i,'j=',j,'k=',k,'l=',l)
#                hyperparams=[l1_mesh[i,j,k,l],l2_mesh[i,j,k,l],l3_mesh[i,j,k,l],l4_mesh[i,j,k,l]]
#                print(hyperparams)
#                MLL_0[i,j,k,l] = log_likelihood_function_0(hyperparams)
#                #print(0)
#                MLL_1[i,j,k,l] = log_likelihood_function_1(hyperparams)
#                #print(1)
#                MLL_2[i,j,k,l] = log_likelihood_function_2(hyperparams)
#                #print(2)
#                
#
#ind_0 = np.unravel_index(np.argmin(MLL_0, axis=None), MLL_0.shape)
#ind_1 = np.unravel_index(np.argmin(MLL_1, axis=None), MLL_1.shape)
#ind_2 = np.unravel_index(np.argmax(MLL_2, axis=None), MLL_2.shape)
#print('For mode0: l_1=',l1_mesh[ind_0],'l_2=',l2_mesh[ind_0],'l_3=',l3_mesh[ind_0],'l_4=',l4_mesh[ind_0])
#print('For mode1: l_1=',l1_mesh[ind_1],'l_2=',l2_mesh[ind_1],'l_3=',l3_mesh[ind_1],'l_4=',l4_mesh[ind_1])
#print('For mode2: l_1=',l1_mesh[ind_2],'l_2=',l2_mesh[ind_2],'l_3=',l3_mesh[ind_2],'l_4=',l4_mesh[ind_2])
#
#
#hyperparams=[1.0,1.0,1.0,1.0]
#
#res0 = minimize(log_likelihood_function_0, hyperparams, method='Nelder-Mead', tol=1e-6)
#res1 = minimize(log_likelihood_function_1, hyperparams, method='Nelder-Mead', tol=1e-6)
#res2 = minimize(log_likelihood_function_2, hyperparams, method='Nelder-Mead', tol=1e-6)
#
#print('From Nelder-Mead_0',res0.x)
#print('From Nelder-Mead_1',res1.x)
#print('From Nelder-Mead_2',res2.x)


#hyperparams=[1.0,1.0,1.0,1.0]
#
#Kcheck= log_likelihood_function_2(hyperparams)
#
#matsize=17
#K=Kcheck[0:(16*matsize - 1),0:(16*matsize - 1)]

#plt.figure(3)
##plt.rc('text', usetex=True)
##plt.rc('font', family='serif')
#plt.contourf(K, 100, cmap='rainbow')
#plt.axis('equal')
#plt.colorbar()
#plt.grid()
#plt.xticks(np.arange(0, 240, step=15))
#plt.yticks(np.arange(0, 240, step=15))
#L=np.linalg.cholesky(K)
#print('Is K symmetric?:',np.allclose(K,K.T))
#print('Are all eigenvalues of K greater than zero?:',np.all(np.linalg.eigvals(K) >= 0))
#print('Does Cholesky-Decomposition work?:',np.allclose(K,np.dot(L,L.T)))

#plt.figure(3)
##plt.rc('text', usetex=True)
##plt.rc('font', family='serif')
#plt.contourf(Ktmp.T , 100, cmap='rainbow')
#plt.axis('equal')
#plt.colorbar()
#plt.grid()
#L=np.linalg.cholesky(K)
#print('Is K symmetric?:',np.allclose(K,K.T))
#print('Are all eigenvalues of K greater than zero?:',np.all(np.linalg.eigvals(K) > 0))
#print('Does Cholesky-Decomposition work?:',np.allclose(K,np.dot(L,L.T)))

#plt.figure(4)
##plt.rc('text', usetex=True)
##plt.rc('font', family='serif')
#plt.contourf(Ktmp + Ktmp.T , 100, cmap='rainbow')
#plt.axis('equal')
#plt.colorbar()
#plt.grid()
##L=np.linalg.cholesky(K)
#print('Is K symmetric?:',np.allclose(K,K.T))
#print('Are all eigenvalues of K greater than zero?:',np.all(np.linalg.eigvals(K) > 0))
#print('Does Cholesky-Decomposition work?:',np.allclose(K,np.dot(L,L.T)))

#for i in range(0,14,1):
#    for j in range(0,14,1):
#        print('i=',i*16,'j=',j*16)
#        L=np.linalg.cholesky(K[i*16:i*16 + 15,j*16:j*16 + 15])