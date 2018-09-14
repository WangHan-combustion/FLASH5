#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 15:31:26 2018

@author: alexgrannan
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt



def TwoDNewtonRaphson(x0_var,y0_var,f1,df1dx,df1dy,f2,df2dx,df2dy):
  
    xguess = x0_var
    yguess = y0_var
    
    Jacob_det = df1dx*df2dy - df1dy*df2dx

    x_ratio = ( df2dy*f1 - df1dy*f2)/Jacob_det
    y_ratio = (-df2dx*f1 + df1dx*f2)/Jacob_det
   
    new_x  = xguess - x_ratio
    new_y  = yguess - y_ratio
    
    diff_x  = abs(new_x - xguess)/xguess
    diff_y  = abs(new_y - yguess)/yguess

    return (new_x,new_y,diff_x,diff_y,x_ratio,y_ratio)



def TwoDNewtonRaphson_test():
    ## For this test the answer is 
    ## x=0.06177012633860736, y=0.7244905153472167
        
    def func1(x_var,y_var):
        
        func1 = 1.0 - 4.0*x_var + 2.0*x_var**2.0 - 2.0*y_var**3.0
        dfunc1dx = -4.0 + 4.0*x_var
        dfunc1dy = -6.0*y_var**2.0
        
        return (func1,dfunc1dx,dfunc1dy)
    
    def func2(x_var,y_var):
        
        func2 = -4.0 + x_var**4.0 + 4.0*y_var + 4.0*y_var**4.0
        dfunc2dx = 4.0*x_var**3.0
        dfunc2dy = 4.0 + 16.0*y_var**3.0
        
        return (func2,dfunc2dx,dfunc2dy)
    
    ## Initial Guesses 
    x0 = 0.1
    y0 = 0.7
    
    ## Limits
    eostol = 1e-13
    fpmin  = 1e-14
    itermax=100
    
    xguess = x0
    yguess = y0
    
    for i in range(1,itermax):
        
        ## Insert specific functions
        
        [f1,df1dx,df1dy]=func1(xguess,yguess)
        [f2,df2dx,df2dy]=func2(xguess,yguess)
        
        ##
    
        [new_xguess,new_yguess,diff_x,diff_y,x_ratio,y_ratio]=TwoDNewtonRaphson(xguess,yguess,f1,df1dx,df1dy,f2,df2dx,df2dy)
        
        xguess = new_xguess
        yguess = new_yguess
        
        if (max(diff_x,diff_y) < eostol) or ((max(abs(x_ratio),abs(y_ratio))) < fpmin):
            print('x',xguess,'y',yguess)
            print('break')
            break
        elif (i==itermax):
            print('The Newton-Raphson Method reached the maximum number of iterations')
        else:
            continue
            print('continue') 
        