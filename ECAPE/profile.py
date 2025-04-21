import numpy as np
import matplotlib.pyplot as plt
from ECAPE_FUNCTIONS import *
from scipy import interpolate
from metpy.plots import SkewT
from metpy.calc import dewpoint_from_relative_humidity
from metpy.units import units
import sys

#===CONSTANTS
T1 = 273.15
T2 = 233.15
dz = 2.5

#===READ SOUNDING FROM VVM fort.98========================================================
file_path = "/data/cysu/taiwanVVM/tpe20050712nor/fort.98"
header_line = " K, RHO(K),THBAR(K),PBAR(K),PIBAR(K),QVBAR(K)"
header_line1= " K, ZZ(K),ZT(K),FNZ(K),FNT(K)"
stop_indicator = "="

with open(file_path, "r") as file:
    lines = file.readlines()

#Find the header line
start_index = None
for i, line in enumerate(lines):
    if header_line1 in line:
        start_index = i + 2
        break

#Parse and collect the data
parts = []
for line in lines[start_index:]:
    if line.strip().startswith(stop_indicator):
        break
    parts.append(line.strip().split())
#True data exists from level 2 to level N-1 
Z_orig=np.array(parts).T[:,1:-1].astype(float)[2,:]

#Find the header line
start_index = None
for i, line in enumerate(lines):
    if header_line in line:
        start_index = i + 2
        break

#Parse and collect the data
parts = []
for line in lines[start_index:]:
    if line.strip().startswith(stop_indicator):
        break
    parts.append(line.strip().split())
#True data exists from level 2 to level N-1 
sounding=np.array(parts).T[:,1:-1].astype(float)

T_orig=sounding[2,:]*sounding[4,:]
P_orig=sounding[3,:]
RS_orig = compute_rsat(T_orig,P_orig,1,T1,T2)
RW_orig = np.minimum(sounding[5,:],0.999*compute_rsat(T_orig,P_orig,1,T1,T2))#ensure no supersaturation
RH_orig = RW_orig/RS_orig

#===INTERPOLATE ORIGINAL SOUNDING TO FINE DZ==============================================
Z_fine=np.arange(Z_orig[0],Z_orig[-1]+2.5,2.5)

f = interpolate.interp1d(Z_orig,T_orig,fill_value="extrapolate",kind="linear")
T_fine = f(Z_fine)

f = interpolate.interp1d(Z_orig,RW_orig,fill_value="extrapolate",kind="linear")
RW_fine = f(Z_fine)#mixing ratio

P_fine = np.zeros(Z_fine.shape)
RS_fine = np.zeros(Z_fine.shape)
Q_fine = np.zeros(Z_fine.shape)

P_fine[0] = P_orig[0]
Q_fine[0] = RW_fine[0]/(1+RW_fine[0])#mass fraction
RS_fine[0] = compute_rsat(T_fine[0],P_fine[0],1,T1,T2)

for iz in range(1,len(Z_fine)):
    P_fine[iz]= P_fine[iz-1] - dz*(P_fine[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_fine[iz-1] )*T_fine[iz-1] )
    RS_fine[iz] = compute_rsat(T_fine[iz],P_fine[iz],1,T1,T2)
    RW_fine[iz] = np.minimum(RW_fine[iz],0.999*RS_fine[iz])
    Q_fine[iz] = RW_fine[iz]/(1+RW_fine[iz])

RH_fine = RW_fine/RS_fine

#===LIFT FROM SURFACE=====================================================================
#===B_lif&RH WILL BE KEPT THROUGHOUT THE CALCULATION======================================
print("T_fine:",compute_CAPE_AND_CIN(T_fine,P_fine,Q_fine,0,0,0,Z_fine,T1,T2))
T_lif,Qv_lif,Qt_lif,B_lif=lift_parcel_adiabatic(T_fine,P_fine,Q_fine,0,0,0,Z_fine,T1,T2)

#===FIND EL HEIGHT AS THE HIGHEST LEVEL OF B>0============================================
bbb = np.where(B_lif>0)[0]
Z_fine_TR = Z_fine[np.max(bbb)]
inds_ST = np.where(Z_fine>Z_fine_TR )

#===GET WELL-MIXED PBL(DSE&MIXING RATIO) BY ENERGY METHOD=================================
dse=cp*T_fine+g*Z_fine

nz=Z_fine.shape[0]

dse_pbl=np.zeros(nz)
qv_pbl=np.zeros(nz)

dse_pbl[0]=dse[0]
qv_pbl[0]=RW_fine[0]

T_fine_mix=np.copy(T_fine)
RW_fine_mix=np.copy(RW_fine)
Q_fine_mix=np.copy(Q_fine)

for k in np.arange(1,nz):
    dse_pbl[k]=(np.sum(0.5*(dse[0:k]+dse[1:k+1]))*2.5+dse_pbl[0]*Z_fine[0])/Z_fine[k]
    qv_pbl[k]=(np.sum(0.5*(RW_fine[0:k]+ RW_fine[1:k+1]))*2.5+qv_pbl[0]*Z_fine[0])/Z_fine[k]

    top_T=(dse_pbl[k]-g*Z_fine[k])/cp
    top_qs = compute_rsat(top_T,P_fine[k],1,T1,T2)
    
    if top_qs<=qv_pbl[k]:
        T_fine_mix[0:k]=(dse_pbl[k-1]-g*Z_fine[0:k])/cp
        RW_fine_mix[0:k]=qv_pbl[k-1]
        Q_fine_mix[0:k]=RW_fine_mix[0:k]/(1+RW_fine_mix[0:k])
        Z_mix_BL=Z_fine[k-1]
        T_mix_dev=T_fine_mix[0:k]-T_fine[0:k]
        Q_mix_dev=Q_fine_mix[0:k]-Q_fine[0:k]
        break

#===FIND PBL AND FREE TROPO. HEIGHT ACOORDING TO MIXING==================================
inds_BL = np.where(Z_fine<=Z_mix_BL)
inds_FT = np.where( np.logical_and(Z_fine>Z_mix_BL,Z_fine<=Z_fine_TR) )

#===RECALCULATE VAR. AFTER MIXING========================================================
P_fine_mix=np.copy(P_fine)
RH_fine_mix=np.copy(RH_fine)
RH_fine_mix[0]=RW_fine_mix[0]/compute_rsat(T_fine_mix[0],P_fine_mix[0],1,T1,T2)

for iz in range(1,len(Z_fine)):
    if iz in inds_BL[0]:
        P_fine_mix[iz]= P_fine_mix[iz-1] - dz*(P_fine_mix[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_fine_mix[iz-1] )*T_fine_mix[iz-1] )
        RH_fine_mix[iz]=RW_fine_mix[iz]/compute_rsat(T_fine_mix[iz],P_fine_mix[iz],1,T1,T2)
    else:
        P_fine_mix[iz]= P_fine_mix[iz-1] - dz*(P_fine_mix[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_fine_mix[iz-1] )*T_fine_mix[iz-1] )
        RW_fine_mix[iz]=RH_fine_mix[iz]*compute_rsat(T_fine_mix[iz],P_fine_mix[iz],1,T1,T2)#in FT RH_fine = RH_fine_mix
        Q_fine_mix[iz]=RW_fine_mix[iz]/(1+RW_fine_mix[iz])

#===RECONSTRUCT PROFILE in FT==================================================================
T_fine_re = np.copy(T_fine_mix)
Q_fine_re = np.copy(Q_fine_mix)
RW_fine_re = np.copy(RW_fine_mix)
P_fine_re = np.copy(P_fine_mix)
RH_fine_re = np.copy(RH_fine_mix)
T_lif_re = np.copy(T_fine_mix)
Qv_lif_re = np.copy(Q_fine_mix)
Qt_lif_re = np.copy(Q_fine_mix)

for iz in inds_FT[0]:
    qvv=(1-Qt_lif_re[iz-1])*compute_rsat(T_lif_re[iz-1],P_fine_re[iz-1],0,T1,T2)
    qvi=(1-Qt_lif_re[iz-1])*compute_rsat(T_lif_re[iz-1],P_fine_re[iz-1],2,T1,T2)
    T_lif_re[iz] = T_lif_re[iz-1] + dz*moislif(T_lif_re[iz-1],Qv_lif_re[iz-1],qvv,qvi,P_fine_re[iz-1],T_fine_re[iz-1],Q_fine_re[iz-1],Qt_lif_re[iz-1],0,0,T1,T2)
    P_fine_re[iz]= P_fine_re[iz-1] - dz*(P_fine_re[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_fine_re[iz-1] )*T_fine_re[iz-1] )
    r = compute_rsat(T_lif_re[iz],P_fine_re[iz],1,T1,T2)
    Qv_lif_re[iz]=(1-Qt_lif_re[iz-1])*r
    Qt_lif_re[iz]=Qt_lif_re[iz-1]
    def equations(x):
        y = np.zeros_like(x)
        
        T_rho_par =  T_lif_re[iz]*(1 + (Rv/Rd)*Qv_lif_re[iz] - Qt_lif_re[iz] )
        T_rho_env = x[0]*(1 + (Rv/Rd - 1)*RH_fine_re[iz]*compute_rsat(x[0],P_fine_re[iz],1,T1,T2)/(1+RH_fine_re[iz]*compute_rsat(x[0],P_fine_re[iz],1,T1,T2)) )
        
        y[0] = np.abs( B_lif[iz] - g*(T_rho_par - T_rho_env)/T_rho_env )
        return y
    
    # Initial guess for the solution
    initial_guess = T_fine_re[iz-1]#,q0[iz-1]
    solution = fsolve(equations, initial_guess)
    T_fine_re[iz]=solution[0]
    RW_fine_re[iz] = RH_fine_re[iz]*compute_rsat(T_fine_re[iz],P_fine_re[iz],1,T1,T2)
    Q_fine_re[iz] = RW_fine_re[iz]/(1+RW_fine_re[iz])
    
T_fine_re[inds_ST[0]]=T_fine_re[np.max(inds_FT)]

for iz in inds_ST[0]:
    P_fine_re[iz]= P_fine_re[iz-1] - dz*(P_fine_re[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_fine_re[iz-1] )*T_fine_re[iz-1] )
    RW_fine_re[iz]=RH_fine_re[iz]*compute_rsat(T_fine_re[iz],P_fine_re[iz],1,T1,T2)#in ST RH_fine_re = RH_fine_mix
    Q_fine_re[iz]=RW_fine_re[iz]/(1+RW_fine_re[iz])
    qvv=(1-Qt_lif_re[iz-1])*compute_rsat(T_lif_re[iz-1],P_fine_re[iz-1],0,T1,T2)
    qvi=(1-Qt_lif_re[iz-1])*compute_rsat(T_lif_re[iz-1],P_fine_re[iz-1],2,T1,T2)
    T_lif_re[iz] = T_lif_re[iz-1] + dz*moislif(T_lif_re[iz-1],Qv_lif_re[iz-1],qvv,qvi,P_fine_re[iz-1],T_fine_re[iz-1],Q_fine_re[iz-1],Qt_lif_re[iz-1],0,0,T1,T2)
    r = compute_rsat(T_lif_re[iz],P_fine_re[iz],1,T1,T2)
    Qv_lif_re[iz]=(1-Qt_lif_re[iz-1])*r
    Qt_lif_re[iz]=Qt_lif_re[iz-1]

#RW_fine_re=Q_fine_re/(1-Q_fine_re)

print("T_fine_re:",compute_CAPE_AND_CIN(T_fine_re,P_fine_re,Q_fine_re,0,0,0,Z_fine,T1,T2))

Td_fine_re=dewpoint_from_relative_humidity((T_fine_re-273.15)*units.degC,(RH_fine_re*100)*units.percent)

#=== warming scenario
T_fine_warm = np.copy(T_fine_re)
RW_fine_warm= np.copy(RW_fine_re)
RH_fine_warm= np.copy(RH_fine_re)
P_fine_warm = np.copy(P_fine_re)
Q_fine_warm = np.copy(Q_fine_re)
B_lif_warm  = np.copy(B_lif)

T_fine_warm[0]=T_fine_warm[0]+3
RW_fine_warm[0] = compute_rsat(T_fine_warm[0],P_fine_warm[0],1,T1,T2)*RH_fine_warm[0]
Q_fine_warm[0]=RW_fine_warm[0]/(1+RW_fine_warm[0])
#Z_warm_BL=compute_LCL_NUMERICAL(T_fine_warm[0],Q_fine_warm[0],P_fine_warm[0],2.5)
#print(Z_mix_BL,Z_warm_BL)

#inds_warm_BL = np.where(Z_fine<=Z_mix_BL)
#inds_warm_FT = np.where( np.logical_and(Z_fine>Z_mix_BL,Z_fine<=Z_fine_TR) )
#inds_warm_ST = np.where(Z_fine>Z_fine_TR )
#print(len(B_lif_warm[inds_warm_BL[0][-1]:inds_FT[0][0]+1]))
#B_lif_warm[inds_warm_BL[0][-1]:inds_FT[0][0]+1]=np.linspace(0,B_lif[inds_FT[0][0]],len(B_lif_warm[inds_warm_BL[0][-1]:inds_FT[0][0]+1]))

T_lif_warm   = np.zeros(Z_fine.shape)
Qv_lif_warm  = np.zeros(Z_fine.shape)
Qt_lif_warm  = np.zeros(Z_fine.shape)

for iz in inds_BL[0]:
    if iz==0:
        Qv_lif_warm[iz] = Q_fine_warm[iz]
        Qt_lif_warm[iz] = Q_fine_warm[iz]   
        T_lif_warm[iz]  = T_fine_warm[iz]            
    else:
        T_fine_warm[iz]= T_fine_warm[iz-1] + dz*drylift(T_fine_warm[iz-1],Q_fine_warm[iz-1],T_fine_warm[iz-1],Q_fine_warm[iz-1],0)
        P_fine_warm[iz]= P_fine_warm[iz-1] - dz*(P_fine_warm[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_fine_warm[iz-1] )*T_fine_warm[iz-1] )
        Q_fine_warm[iz]=Q_fine_warm[iz-1]
        Qv_lif_warm[iz] = Q_fine_warm[0]
        Qt_lif_warm[iz] = Q_fine_warm[0]  
        T_lif_warm[iz] = T_lif_warm[iz-1] + dz*drylift(T_lif_warm[iz-1],Qv_lif_warm[iz-1],T_fine_warm[iz-1],Q_fine_warm[iz-1],0)             
    RH_fine_warm[iz]=(Q_fine_warm[iz]/(1-Q_fine_warm[iz]))/compute_rsat(T_fine_warm[iz],P_fine_warm[iz],1,T1,T2)

#RH_fine_warm[inds_warm_BL[0][-1]:inds_FT[0][0]+1]=np.linspace(RH_fine_warm[inds_warm_BL[0][-1]],RH_fine[inds_FT[0][0]],len(B_lif_warm[inds_warm_BL[0][-1]:inds_FT[0][0]+1]))

for iz in inds_FT[0]:
    qvv=(1-Qt_lif_warm[iz-1])*compute_rsat(T_lif_warm[iz-1],P_fine_warm[iz-1],0,T1,T2)
    qvi=(1-Qt_lif_warm[iz-1])*compute_rsat(T_lif_warm[iz-1],P_fine_warm[iz-1],2,T1,T2)
    T_lif_warm[iz] = T_lif_warm[iz-1] + dz*moislif(T_lif_warm[iz-1],Qv_lif_warm[iz-1],qvv,qvi,P_fine_warm[iz-1],T_fine_warm[iz-1],Q_fine_warm[iz-1],Qt_lif_warm[iz-1],0,0,T1,T2)
    P_fine_warm[iz]= P_fine_warm[iz-1] - dz*(P_fine_warm[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_fine_warm[iz-1] )*T_fine_warm[iz-1] )
    r = compute_rsat(T_lif_warm[iz],P_fine_warm[iz],1,T1,T2)
    Qv_lif_warm[iz]=(1-Qt_lif_warm[iz-1])*r
    Qt_lif_warm[iz]=Qt_lif_warm[iz-1]
    #print(iz,B_lif_warm[iz])
    def equations(x):
        y = np.zeros_like(x)
        #y[0] = np.abs( B_prof[iz] - ( g*(T_lif[iz]-x[0])/x[0] + g*(Rv/Rd - 1)*(Qv_lif[iz]-RH_TROP*compute_rsat(x[0],p0[iz],1,T1,T2)/(1+RH_TROP*compute_rsat(x[0],p0[iz],1,T1,T2))) - g*(Qt_lif[iz]-Qv_lif[iz]) ) )
        
        T_rho_par =  T_lif_warm[iz]*(1 + (Rv/Rd)*Qv_lif_warm[iz] - Qt_lif_warm[iz] )
        T_rho_env = x[0]*(1 + (Rv/Rd - 1)*RH_fine_warm[iz]*compute_rsat(x[0],P_fine_warm[iz],1,T1,T2)/(1+RH_fine_warm[iz]*compute_rsat(x[0],P_fine_warm[iz],1,T1,T2)) )
        
        y[0] = np.abs( B_lif_warm[iz] - g*(T_rho_par - T_rho_env)/T_rho_env )
#        y[0] = B[iz] - g*(T_rho_par - T_rho_env)/T_rho_env
        return y
    
    # Initial guess for the solution
    initial_guess = T_fine_warm[iz-1]#,q0[iz-1]
    solution = fsolve(equations, initial_guess)
    T_fine_warm[iz]=solution[0]
    #q0[iz]=solution[1]
    Q_fine_warm[iz] = RH_fine_warm[iz]*compute_rsat(T_fine_warm[iz],P_fine_warm[iz],1,T1,T2)/(1+RH_fine_warm[iz]*compute_rsat(T_fine_warm[iz],P_fine_warm[iz],1,T1,T2))
    
T_fine_warm[inds_ST[0]]=T_fine_warm[np.max(inds_FT)]
#T_fine_warm[inds_ST[0][1:]]=np.nan

for iz in inds_ST[0]:
    P_fine_warm[iz]= P_fine_warm[iz-1] - dz*(P_fine_warm[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_fine_warm[iz-1] )*T_fine_warm[iz-1] )
    RW_fine_warm[iz]=RH_fine_warm[iz]*compute_rsat(T_fine_warm[iz],P_fine_warm[iz],1,T1,T2)#in ST RH_fine_warm = RH_fine_mix
    Q_fine_warm[iz]=RW_fine_warm[iz]/(1+RW_fine_warm[iz])
#    P_fine_warm[iz]= P_fine_warm[iz-1] - dz*(P_fine_warm[iz-1]*g)/(Rd*T_fine_warm[iz-1] )
    qvv=(1-Qt_lif_warm[iz-1])*compute_rsat(T_lif_warm[iz-1],P_fine_warm[iz-1],0,T1,T2)
    qvi=(1-Qt_lif_warm[iz-1])*compute_rsat(T_lif_warm[iz-1],P_fine_warm[iz-1],2,T1,T2)
    T_lif_warm[iz] = T_lif_warm[iz-1] + dz*moislif(T_lif_warm[iz-1],Qv_lif_warm[iz-1],qvv,qvi,P_fine_warm[iz-1],T_fine_warm[iz-1],Q_fine_warm[iz-1],Qt_lif_warm[iz-1],0,0,T1,T2)
    r = compute_rsat(T_lif_warm[iz],P_fine_warm[iz],1,T1,T2)
    Qv_lif_warm[iz]=(1-Qt_lif_warm[iz-1])*r
    Qt_lif_warm[iz]=Qt_lif_warm[iz-1]


#for iz in inds_ST[0]:
#    P_fine_warm[iz]= P_fine_warm[iz-1] - dz*(P_fine_warm[iz-1]*g)/(Rd*T_fine_warm[iz-1] )
#    qvv=(1-Qt_lif_warm[iz-1])*compute_rsat(T_lif_warm[iz-1],P_fine_warm[iz-1],0,T1,T2)
#    qvi=(1-Qt_lif_warm[iz-1])*compute_rsat(T_lif_warm[iz-1],P_fine_warm[iz-1],2,T1,T2)
#    T_lif_warm[iz] = T_lif_warm[iz-1] + dz*moislif(T_lif_warm[iz-1],Qv_lif_warm[iz-1],qvv,qvi,P_fine_warm[iz-1],T_fine_warm[iz-1],Q_fine_warm[iz-1],Qt_lif_warm[iz-1],0,0,T1,T2)

print("T_fine_warm:",compute_CAPE_AND_CIN(T_fine_warm,P_fine_warm,Q_fine_warm,0,0,0,Z_fine,T1,T2))

Td_fine_warm=dewpoint_from_relative_humidity((T_fine_warm-273.15)*units.degC,(RH_fine_warm*100)*units.percent)

T_warm=np.copy(T_fine_warm)
Q_warm=np.copy(Q_fine_warm)
P_warm=np.copy(P_fine_warm)
for iz in inds_BL[0]:
    T_warm[iz]=T_fine_warm[iz]-T_mix_dev[iz]
    Q_warm[iz]=Q_fine_warm[iz]-Q_mix_dev[iz]
for iz in range(1,len(Z_fine)):
    P_warm[iz]= P_warm[iz-1] - dz*(P_warm[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_warm[iz-1] )*T_warm[iz-1] )

print("T_warm:",compute_CAPE_AND_CIN(T_warm,P_warm,Q_warm,0,0,0,Z_fine,T1,T2))
T_lif_warm_re,Qv_lif_warm_re,Qt_lif_warm_re,_=lift_parcel_adiabatic(T_warm,P_warm,Q_warm,0,0,0,Z_fine,T1,T2)

T_warm_re = np.copy(T_warm)
Q_warm_re = np.copy(Q_warm)
RW_warm_re=Q_warm_re/(1-Q_warm_re)
P_warm_re = np.copy(P_warm)

for iz in inds_FT[0]:
    qvv=(1-Qt_lif_warm_re[iz-1])*compute_rsat(T_lif_warm_re[iz-1],P_warm_re[iz-1],0,T1,T2)
    qvi=(1-Qt_lif_warm_re[iz-1])*compute_rsat(T_lif_warm_re[iz-1],P_warm_re[iz-1],2,T1,T2)
    T_lif_warm_re[iz] = T_lif_warm_re[iz-1] + dz*moislif(T_lif_warm_re[iz-1],Qv_lif_warm_re[iz-1],qvv,qvi,P_warm_re[iz-1],T_warm_re[iz-1],Q_warm_re[iz-1],Qt_lif_warm_re[iz-1],0,0,T1,T2)
    P_warm_re[iz]= P_warm_re[iz-1] - dz*(P_warm_re[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_warm_re[iz-1] )*T_warm_re[iz-1] )
    r = compute_rsat(T_lif_warm_re[iz],P_warm_re[iz],1,T1,T2)
    Qv_lif_warm_re[iz]=(1-Qt_lif_warm_re[iz-1])*r
    Qt_lif_warm_re[iz]=Qt_lif_warm_re[iz-1]
    #print(iz,B_lif_warm[iz])
    def equations(x):
        y = np.zeros_like(x)
        #y[0] = np.abs( B_prof[iz] - ( g*(T_lif[iz]-x[0])/x[0] + g*(Rv/Rd - 1)*(Qv_lif[iz]-RH_TROP*compute_rsat(x[0],p0[iz],1,T1,T2)/(1+RH_TROP*compute_rsat(x[0],p0[iz],1,T1,T2))) - g*(Qt_lif[iz]-Qv_lif[iz]) ) )
        
        T_rho_par =  T_lif_warm_re[iz]*(1 + (Rv/Rd)*Qv_lif_warm_re[iz] - Qt_lif_warm_re[iz] )
        T_rho_env = x[0]*(1 + (Rv/Rd - 1)*RH_fine_warm[iz]*compute_rsat(x[0],P_warm_re[iz],1,T1,T2)/(1+RH_fine_warm[iz]*compute_rsat(x[0],P_warm_re[iz],1,T1,T2)) )
        
        y[0] = np.abs( B_lif_warm[iz] - g*(T_rho_par - T_rho_env)/T_rho_env )
#        y[0] = B[iz] - g*(T_rho_par - T_rho_env)/T_rho_env
        return y
    
    # Initial guess for the solution
    initial_guess = T_warm_re[iz-1]#,q0[iz-1]
    solution = fsolve(equations, initial_guess)
    T_warm_re[iz]=solution[0]
    #q0[iz]=solution[1]
    RW_warm_re[iz] = RH_fine_warm[iz]*compute_rsat(T_warm_re[iz],P_warm_re[iz],1,T1,T2)
    Q_warm_re[iz] = RW_warm_re[iz]/(1+RW_warm_re[iz])
    
T_warm_re[inds_ST[0]]=T_warm_re[np.max(inds_FT)]
#T_fine_warm[inds_ST[0][1:]]=np.nan

for iz in inds_ST[0]:
    P_warm_re[iz]= P_warm_re[iz-1] - dz*(P_warm_re[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_warm_re[iz-1] )*T_warm_re[iz-1] )
    RW_warm_re[iz]=RH_fine_warm[iz]*compute_rsat(T_warm_re[iz],P_warm_re[iz],1,T1,T2)#in ST RH_fine_warm = RH_fine_mix
    Q_warm_re[iz]=RW_warm_re[iz]/(1+RW_warm_re[iz])
#    P_warm_re[iz]= P_warm_re[iz-1] - dz*(P_warm_re[iz-1]*g)/(Rd*T_warm_re[iz-1] )
    qvv=(1-Qt_lif_warm_re[iz-1])*compute_rsat(T_lif_warm_re[iz-1],P_warm_re[iz-1],0,T1,T2)
    qvi=(1-Qt_lif_warm_re[iz-1])*compute_rsat(T_lif_warm_re[iz-1],P_warm_re[iz-1],2,T1,T2)
    T_lif_warm_re[iz] = T_lif_warm_re[iz-1] + dz*moislif(T_lif_warm_re[iz-1],Qv_lif_warm_re[iz-1],qvv,qvi,P_warm_re[iz-1],T_warm_re[iz-1],Q_warm_re[iz-1],Qt_lif_warm_re[iz-1],0,0,T1,T2)
    r = compute_rsat(T_lif_warm_re[iz],P_warm_re[iz],1,T1,T2)
    Qv_lif_warm_re[iz]=(1-Qt_lif_warm_re[iz-1])*r
    Qt_lif_warm_re[iz]=Qt_lif_warm_re[iz-1]

print("T_warm_re:",compute_CAPE_AND_CIN(T_warm_re,P_warm_re,Q_warm_re,0,0,0,Z_fine,T1,T2))
#for iz in inds_BL[0]:
#    print("BL: ",iz,RH_fine[iz],RH_fine_re[iz],RH_fine_warm[iz])
#for iz in inds_FT[0]:
#    print("FT: ",iz,RH_fine[iz],RH_fine_re[iz],RH_fine_warm[iz])
#for iz in inds_ST[0]:
#    print("ST: ",iz,RH_fine[iz],RH_fine_re[iz],RH_fine_warm[iz])
#sys.exit()
#print(T0_lif[inds_BL[0]])
#print(T_lif[inds_BL[0]])
#sys.exit()
#T1_lif,Qv1_lif,Qt1_lif,B_lif=lift_parcel_adiabatic(T0,p0,q0,0,0,0,Z,T1,T2)
#for iz in inds_FT[0]:
    #print(iz,B_lif[iz],B_lif_warm[iz])
T_lif_2,_,_,B_lif_2=lift_parcel_adiabatic(T_fine_mix,P_fine_mix,Q_fine_mix,0,0,0,Z_fine,T1,T2)

Td_orig=dewpoint_from_relative_humidity((T_orig-273.15)*units.degC,(RH_orig*100)*units.percent)
Td_fine=dewpoint_from_relative_humidity((T_fine-273.15)*units.degC,(RH_fine*100)*units.percent)
Td_fine_mix=dewpoint_from_relative_humidity((T_fine_mix-273.15)*units.degC,(RH_fine_mix*100)*units.percent)

#print(compute_CAPE_AND_CIN(T_fine_re,P_fine_re,Q_fine_re,0,0,0,Z_fine,T1,T2))
#print(compute_CAPE_AND_CIN(T_fine_mix,P_fine_mix,Q_fine_mix,0,0,0,Z_fine,T1,T2))
#print(compute_CAPE_AND_CIN(T_fine,P_fine,Q_fine,0,0,0,Z_fine,T1,T2))
#print(compute_CAPE_AND_CIN(T_fine_warm,P_fine_warm,Q_fine_warm,0,0,0,Z_fine,T1,T2))
#print(B_lif[inds_BL[0]])
#print(B_lif_warm[inds_BL[0]])

#print(P_fine[inds_FT[0][-1]],P_fine_re[inds_FT[0][-1]],P_fine_warm[inds_FT[0][-1]])
skew = SkewT()

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot
skew.plot(P_fine/100, T_lif-273.15, 'r')
#skew.plot(P_fine_mix/100, T_lif_2-273.15, 'r')
skew.plot(P_fine/100, T_fine-273.15, 'k')
#skew.plot(P_fine_mix/100, T_fine_mix-273.15, 'k')
#skew.plot(P_fine_re/100, T_fine_re-273.15, 'g')
#skew.plot(P_fine_warm/100, T_fine_warm-273.15, 'b')
#skew.plot(P_fine_warm/100, T_lif_warm-273.15, 'r')
#skew.plot(P_warm/100, T_warm-273.15, 'k')
skew.plot(P_warm_re/100, T_warm_re-273.15, 'k')
skew.plot(P_warm/100, T_lif_warm_re-273.15, 'r')
#skew.plot(P_fine_re/100, T_fine_re-273.15, 'r')
#skew.plot(P_fine_re/100, T_lif_re-273.15, 'b')
#skew.plot(P_fine_warm/100, T_lif_warm-273.15, 'b')
#skew.plot(P_fine_warm/100, T_fine_warm-273.15, 'k')

#skew.plot(P_fine/100, Td_fine, 'g',ls='--')
#skew.plot(P_fine_mix/100, Td_fine_mix, 'b',ls='--')
#skew.plot(P_fine_re/100, Td_fine_re, 'k',ls='--')
#skew.plot(P_fine_warm/100, Td_fine_warm, 'k',ls='--')

#skew.plot(P_orig/100, Td_orig, 'r', ls='--')
#skew.plot(P_fine/100, Td_fine, 'g', ls='--')
#skew.plot(P_fine_mix/100, Td_fine_mix, 'b', ls='--')
#skew.plot(P_fine_re/100, Td_fine_re, 'r', ls='--')
#skew.plot(p0/100, T0-273.15, 'r')
#skew.plot(p0/100, Td0, 'g')
#skew.plot(p0/100, T0-273.15, 'r')
#skew.plot(P_orig_hydrostatic/100, Td, 'g', ls='--')
#skew.plot(P_orig_hydrostatic/100, T-273.15, 'r', ls='--')
#skew.plot(p0/100, T_lif-273.15, 'b',ls='-')
#skew.plot(p0/100, T0_lif-273.15, 'b',ls='--')
#skew.plot(p0/100, Td, 'g')
#skew.plot(P_orig_hydrostatic/100, T-273.15, 'r',ls='--')
#skew.plot(P_orig/100, Td_orig, 'g',ls='--')

# Set some better labels than the default
skew.ax.set_xlabel('Temperature (\N{DEGREE CELSIUS})')
skew.ax.set_ylabel('Pressure (mb)')

# Add the relevant special lines
#skew.plot_dry_adiabats()
#skew.plot_moist_adiabats()
#skew.plot_mixing_lines()
skew.ax.set_ylim(1000, P_fine_warm[-1]/100)
skew.ax.set_xlim(-30, 35)
plt.show()

#r0=sounding[5,:]
#rs0 = compute_rsat(T0,p0,1,T1,T2)
#print(r0,r0/rs0)



###r=mixing ratio, q=mass fraction or specific humidity
#T_SFC=sounding[2,0]*sounding[4,0]
#P_SFC=sounding[3,0]
#r_SFC=sounding[5,0]
##rs = compute_rsat(T_SFC,P_SFC,1,T1,T2)
##RH_SFC=r_SFC/rs
##
#q_sfc=r_SFC/(1+r_SFC)
##
##T_lif  = np.zeros(z0.shape)
##Qv_lif = np.zeros(z0.shape)
##Qt_lif = np.zeros(z0.shape)
##B_lif  = np.zeros(z0.shape)
#
#for k in range(len(z0)):
#    print(k,z0[k],p0[k],T_lif[k],Qv_lif[k],B_lif[k])
