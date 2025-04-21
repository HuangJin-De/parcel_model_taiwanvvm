#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt 
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
import csv
from scipy import interpolate
import matplotlib.gridspec as gridspec
#
#import metpy.calc as mpcalc
#from metpy.plots import Hodograph, SkewT
#from metpy.units import units
#
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys
sys.path.append('/Users/jfp5948/Desktop/PYTHON_PROJECTS/')
from ECAPE_FUNCTIONS import *
import os
import os.path
from os import path
import subprocess
import warnings
from scipy.optimize import fsolve

def create_profile(Z_TR,Z_top,P_SFC,RH_SFC,RH_TROP,CAPE,T_SFC,T1,T2,prate_global,heavy,dzz=2.5):

    
    # Read the text file
    with open('CAPE500_LOW.txt', 'r') as file:
        lines = file.readlines()
    
    # Initialize an empty matrix
    matrix = np.zeros((len(lines), 2))
    
    # Parse the coordinates and populate the matrix
    for i, line in enumerate(lines):
        x, y = line.strip().split(',')
        matrix[i] = [float(x), float(y)]
        
    BUOY_295 = matrix[:,0]
    
    
    # Read the text file
    with open('CAPE500_HIGH.txt', 'r') as file:
        lines = file.readlines()
    
    # Initialize an empty matrix
    matrix = np.zeros((len(lines), 2))
    
    # Parse the coordinates and populate the matrix
    for i, line in enumerate(lines):
        x, y = line.strip().split(',')
        matrix[i] = [float(x), float(y)]
        
    BUOY_305 = matrix[:,0]
    dz=100
    dbuoy = (BUOY_305 - BUOY_295)/10
    B_prof_500 = BUOY_295 + (T_SFC - 295)*dbuoy



    # Read the text file
    with open('CAPE1000_LOW.txt', 'r') as file:
        lines = file.readlines()
    
    # Initialize an empty matrix
    matrix = np.zeros((len(lines), 2))
    
    # Parse the coordinates and populate the matrix
    for i, line in enumerate(lines):
        x, y = line.strip().split(',')
        matrix[i] = [float(x), float(y)]
        
    BUOY_295 = matrix[:,0]
    
    
    # Read the text file
    with open('CAPE1000_HIGH.txt', 'r') as file:
        lines = file.readlines()
    
    # Initialize an empty matrix
    matrix = np.zeros((len(lines), 2))
    
    # Parse the coordinates and populate the matrix
    for i, line in enumerate(lines):
        x, y = line.strip().split(',')
        matrix[i] = [float(x), float(y)]
        
    BUOY_305 = matrix[:,0]
    dz=100
    dbuoy = (BUOY_305 - BUOY_295)/10
    B_prof_1000 = BUOY_295 + (T_SFC - 295)*dbuoy





    
    # Read the text file
    with open('CAPE2000_LOW.txt', 'r') as file:
        lines = file.readlines()
    
    # Initialize an empty matrix
    matrix = np.zeros((len(lines), 2))
    
    # Parse the coordinates and populate the matrix
    for i, line in enumerate(lines):
        x, y = line.strip().split(',')
        matrix[i] = [float(x), float(y)]
        
    BUOY_295 = matrix[:,0]
    
    
    # Read the text file
    with open('CAPE2000_HIGH.txt', 'r') as file:
        lines = file.readlines()
    
    # Initialize an empty matrix
    matrix = np.zeros((len(lines), 2))
    
    # Parse the coordinates and populate the matrix
    for i, line in enumerate(lines):
        x, y = line.strip().split(',')
        matrix[i] = [float(x), float(y)]
        
    BUOY_305 = matrix[:,0]
    dz=100
    dbuoy = (BUOY_305 - BUOY_295)/10
    B_prof_2000 = BUOY_295 + (T_SFC - 295)*dbuoy
    
    
    
    # Read the text file
    with open('CAPE3000_LOW.txt', 'r') as file:
        lines = file.readlines()
    
    # Initialize an empty matrix
    matrix = np.zeros((len(lines), 2))
    
    # Parse the coordinates and populate the matrix
    for i, line in enumerate(lines):
        x, y = line.strip().split(',')
        matrix[i] = [float(x), float(y)]
        
    BUOY_295 = matrix[:,0]
    
    
    # Read the text file
    with open('CAPE3000_HIGH.txt', 'r') as file:
        lines = file.readlines()
    
    # Initialize an empty matrix
    matrix = np.zeros((len(lines), 2))
    
    # Parse the coordinates and populate the matrix
    for i, line in enumerate(lines):
        x, y = line.strip().split(',')
        matrix[i] = [float(x), float(y)]
        
    BUOY_305 = matrix[:,0]
    
    dbuoy = (BUOY_305 - BUOY_295)/10
    B_prof_3000 = BUOY_295 + (T_SFC - 295)*dbuoy
    
    if CAPE >= 2000:
        dCAPE = (B_prof_3000 - B_prof_2000)/1000
        B_prof = B_prof_2000 + (CAPE - 2000)*dCAPE
    else:
        dCAPE = (B_prof_2000)/2000
        B_prof = (CAPE)*dCAPE
    
    aninds = np.where(np.isfinite(B_prof))[0]
    ln = len(aninds)
    mxB = np.max(np.where(~np.isnan(B_prof))[0])
    maxval = np.nanmax(B_prof)
    if B_prof[aninds[ln-1]]>maxval*.6:
        
        if B_prof[mxB]>0:      
            rlen = 10
            mxbp = mxB+rlen
            for iz in np.arange(mxB+1,mxbp+1,1):
                m = -B_prof[mxB]/rlen
                B_prof[iz] = B_prof[mxB] + m*(iz - mxB)
            
    naninds = np.where(np.isnan(B_prof))[0]
    aninds = np.where(np.isfinite(B_prof))[0]
    
    z_old = matrix[:,1]
    
    f = interpolate.interp1d(z_old[aninds],B_prof[aninds],fill_value="extrapolate",kind="linear"); B_prof[naninds] = f(z_old[naninds])
    
    B_prof[np.where(B_prof<0)[0]]=np.NaN
        
    B_prof = CAPE*B_prof/(np.nansum(np.maximum(B_prof,0))*dz)
    
    dz=dzz
    
    
    

    #CONSTANTS
    Rd=287.04#%dry gas constant
    Rv=461.5 #water vapor gas constant
    epsilon=Rd/Rv
    cp=1005 #specific heat of dry air at constant pressure
    g=9.81 #gravitational acceleration
    xlv=2501000 #reference latent heat of vaporization at the triple point temperature
    xls=2834000 #reference latent heat of sublimation at the triple point temperature
    cpv=1870 #specific heat of water vapor at constant pressure
    cpl=4190 #specific heat of liquid water
    cpi=2106 #specific heat of ice
    ttrip=273.15; #triple point temperature
    eref=611.2 #reference pressure at the triple point temperature
    pi = 3.1415926535
    
    #Z_BL = 1000 #HEIGHT OF THE BOUNDARY LAYER TOP
    
    #
    iz=0
    r = compute_rsat(T_SFC,P_SFC,1,T1,T2)
    qsfc=RH_SFC*r/(1 + r)
    #Z_BL=compute_LCL(T_SFC,qsfc,P_SFC)
    
    Z_BL=compute_LCL_NUMERICAL(T_SFC,qsfc,P_SFC,dz)
    
    low_ind = np.where(np.isfinite(B_prof))[0][0]
    z_old = z_old - z_old[low_ind] + Z_BL
    
    #
    Z = np.arange(0,Z_top,dz)
    f = interpolate.interp1d(z_old,B_prof,fill_value="extrapolate",kind="linear"); B_prof = f(Z)
    B_prof = np.maximum(B_prof,0)
    B_prof[np.where(np.isnan(B_prof))[0]]=0
    #Z_BL = Z[np.min(np.where(B_prof>0)[0])]
    bbb = np.where(B_prof>0)[0]
    if len(bbb)>0:
        Z_TR = Z[np.max(bbb)]
    else:
        Z_TR = 12*1000
    z0 = Z
    inds_BL = np.where(Z<=Z_BL)
    inds_FT = np.where( np.logical_and(Z>Z_BL,Z<=Z_TR) )
    inds_ST = np.where(Z>Z_TR )
    #
    

    
    T0 = np.zeros(Z.shape)
    q0 = np.zeros(Z.shape)
    p0 = np.zeros(Z.shape)
    T_lif = np.zeros(Z.shape)
    Qv_lif = np.zeros(Z.shape)
    Qt_lif = np.zeros(Z.shape)
    
    for iz in inds_BL[0]:
        if iz==0:
            T0[iz]=T_SFC
            p0[iz]=P_SFC
            r = compute_rsat(T0[iz],p0[iz],1,T1,T2)
            q0[iz]=RH_SFC*r/(1 + r)
            Qv_lif[iz] = q0[iz]
            Qt_lif[iz] = q0[iz]   
            T_lif[iz] = T0[iz]            
        else:
            T0[iz]= T0[iz-1] + dz*drylift(T0[iz-1],q0[iz-1],T0[iz-1],q0[iz-1],0)
            p0[iz]= p0[iz-1] - dz*(p0[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*q0[iz-1] )*T0[iz-1] )
            q0[iz]=q0[iz-1]
            Qv_lif[iz] = q0[0]
            Qt_lif[iz] = q0[0]  
            T_lif[iz] = T_lif[iz-1] + dz*drylift(T_lif[iz-1],Qv_lif[iz-1],T0[iz-1],q0[iz-1],0)             
                   
    for iz in inds_FT[0]:
        qvv=(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],0,T1,T2)
        qvi=(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],2,T1,T2)
        T_lif[iz] = T_lif[iz-1] + dz*moislif(T_lif[iz-1],Qv_lif[iz-1],qvv,qvi,p0[iz-1],T0[iz-1],q0[iz-1],Qt_lif[iz-1],0,0,T1,T2)
        p0[iz]= p0[iz-1] - dz*(p0[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*q0[iz-1] )*T0[iz-1] )
        r = compute_rsat(T_lif[iz],p0[iz],1,T1,T2)
        #Qv_lif[iz]=r/(1 + r)
        Qv_lif[iz]=(1-Qt_lif[iz-1])*r
        Qt_lif[iz]=Qt_lif[iz-1]
        def equations(x):
            y = np.zeros_like(x)
            #y[0] = np.abs( B_prof[iz] - ( g*(T_lif[iz]-x[0])/x[0] + g*(Rv/Rd - 1)*(Qv_lif[iz]-RH_TROP*compute_rsat(x[0],p0[iz],1,T1,T2)/(1+RH_TROP*compute_rsat(x[0],p0[iz],1,T1,T2))) - g*(Qt_lif[iz]-Qv_lif[iz]) ) )
            
            T_rho_par =  T_lif[iz]*(1 + (Rv/Rd)*Qv_lif[iz] - Qt_lif[iz] )
            T_rho_env = x[0]*(1 + (Rv/Rd - 1)*RH_TROP*compute_rsat(x[0],p0[iz],1,T1,T2)/(1+RH_TROP*compute_rsat(x[0],p0[iz],1,T1,T2)) )
            
            y[0] = np.abs( B_prof[iz] - g*(T_rho_par - T_rho_env)/T_rho_env )
            return y
        
        # Initial guess for the solution
        initial_guess = T0[iz-1]#,q0[iz-1]
        solution = fsolve(equations, initial_guess)
        T0[iz]=solution[0]
        #q0[iz]=solution[1]
        q0[iz] = RH_TROP*compute_rsat(T0[iz],p0[iz],1,T1,T2)/(1+RH_TROP*compute_rsat(T0[iz],p0[iz],1,T1,T2))
        
    T0[inds_ST[0]]=T0[np.max(inds_FT)]
    
    for iz in inds_ST[0]:
        p0[iz]= p0[iz-1] - dz*(p0[iz-1]*g)/(Rd*T0[iz-1] )
    
    return z0,p0,T0,q0


