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
file_path = "/data/cysu/taiwanVVM/tpe20140525nor/fort.98"
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
Td_orig=dewpoint_from_relative_humidity((T_orig-273.15)*units.degC,(RH_orig*100)*units.percent)

#===INTERPOLATE ORIGINAL SOUNDING TO FINE DZ==============================================
Z_fine=np.arange(Z_orig[0],Z_orig[-1]+2.5,2.5)

f = interpolate.interp1d(Z_orig,T_orig,fill_value="extrapolate",kind="linear")
T_fine = f(Z_fine)

f = interpolate.interp1d(Z_orig,RW_orig,fill_value="extrapolate",kind="linear")
RW_fine = f(Z_fine)#mixing ratio

P_fine = np.zeros(Z_fine.shape)
Q_fine = np.zeros(Z_fine.shape)
RS_fine = np.zeros(Z_fine.shape)

P_fine[0] = P_orig[0]
Q_fine[0] = RW_fine[0]/(1+RW_fine[0])#mass fraction
RS_fine[0] = compute_rsat(T_fine[0],P_fine[0],1,T1,T2)

for iz in range(1,len(Z_fine)):
    P_fine[iz]= P_fine[iz-1] - dz*(P_fine[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_fine[iz-1] )*T_fine[iz-1] )
    RS_fine[iz] = compute_rsat(T_fine[iz],P_fine[iz],1,T1,T2)
    RW_fine[iz] = np.minimum(RW_fine[iz],0.999*RS_fine[iz])#ensure no supersaturation
    Q_fine[iz] = RW_fine[iz]/(1+RW_fine[iz])

RH_fine = RW_fine/RS_fine
Td_fine=dewpoint_from_relative_humidity((T_fine-273.15)*units.degC,(RH_fine*100)*units.percent)

#===LIFT FROM SURFACE=====================================================================
print("T_fine:",compute_CAPE_AND_CIN(T_fine,P_fine,Q_fine,0,0,0,Z_fine,T1,T2))
T_lif_fine,Qv_lif_fine,Qt_lif_fine,B_lif_fine=lift_parcel_adiabatic(T_fine,P_fine,Q_fine,0,0,0,Z_fine,T1,T2)

##===FIND EL HEIGHT AS THE HIGHEST LEVEL OF B>0============================================
#bbb = np.where(B_lif_fine>0)[0]
#Z_fine_TR = Z_fine[np.max(bbb)]
#inds_ST = np.where(Z_fine>Z_fine_TR )
#
#dT_ST=np.zeros(Z_fine.shape)
#for iz in inds_ST[0]:
#    dT_ST[iz]=T_fine[iz]-T_fine[iz-1]
#
##===GET WELL-MIXED TOP-SATURATED BL BY ENERGY METHOD=====================================
#dse=cp*T_fine+g*Z_fine
#
#nz=Z_fine.shape[0]
#
#dse_pbl=np.zeros(nz)
#qv_pbl=np.zeros(nz)
#
#dse_pbl[0]=dse[0]
#qv_pbl[0]=RW_fine[0]
#
#T_fine_mix=np.copy(T_fine)
#RW_fine_mix=np.copy(RW_fine)
#Q_fine_mix=np.copy(Q_fine)
#
#for k in np.arange(1,nz):
#    dse_pbl[k]=(np.sum(0.5*(dse[0:k]+dse[1:k+1]))*dz+dse_pbl[0]*Z_fine[0])/Z_fine[k]
#    qv_pbl[k]=(np.sum(0.5*(RW_fine[0:k]+ RW_fine[1:k+1]))*dz+qv_pbl[0]*Z_fine[0])/Z_fine[k]
#
#    top_T=(dse_pbl[k]-g*Z_fine[k])/cp
#    top_qs = compute_rsat(top_T,P_fine[k],1,T1,T2)
#    
##mixing from the surface to the level just below saturation
#    if top_qs<=qv_pbl[k]:
#        Z_mix_BL=Z_fine[k-1]
#        T_fine_mix[0:k]=(dse_pbl[k-1]-g*Z_fine[0:k])/cp
#        RW_fine_mix[0:k]=qv_pbl[k-1]
#        Q_fine_mix[0:k]=RW_fine_mix[0:k]/(1+RW_fine_mix[0:k])
#        T_mix_dev=T_fine_mix[0:k]-T_fine[0:k]
#        RW_mix_dev=RW_fine_mix[0:k]-RW_fine[0:k]
#        break
#
##===FIND PBL AND FREE TROPO. HEIGHT ACOORDING TO MIXING==================================
#inds_BL = np.where(Z_fine<=Z_mix_BL)
#inds_FT = np.where( np.logical_and(Z_fine>Z_mix_BL,Z_fine<=Z_fine_TR) )
#
##===CALCULATE P, RH, Td AFTER MIXING=====================================================
##===Note: RH_fine above BL is fixed and copied throughout the code=======================
#P_fine_mix=np.copy(P_fine)
#RH_fine_mix=np.copy(RH_fine)
#RH_fine_mix[0]=RW_fine_mix[0]/compute_rsat(T_fine_mix[0],P_fine_mix[0],1,T1,T2)#surface RW and RH will decrease after mixing
#
#for iz in range(1,len(Z_fine)):
#    if iz in inds_BL[0]:
#        P_fine_mix[iz]= P_fine_mix[iz-1] - dz*(P_fine_mix[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_fine_mix[iz-1] )*T_fine_mix[iz-1] )
#        RH_fine_mix[iz]=RW_fine_mix[iz]/compute_rsat(T_fine_mix[iz],P_fine_mix[iz],1,T1,T2)
#    else:
#        P_fine_mix[iz]= P_fine_mix[iz-1] - dz*(P_fine_mix[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_fine_mix[iz-1] )*T_fine_mix[iz-1] )
#        RW_fine_mix[iz]=RH_fine_mix[iz]*compute_rsat(T_fine_mix[iz],P_fine_mix[iz],1,T1,T2)#in FT RH_fine = RH_fine_mix
#        Q_fine_mix[iz]=RW_fine_mix[iz]/(1+RW_fine_mix[iz])
#
#Td_fine_mix=dewpoint_from_relative_humidity((T_fine_mix-273.15)*units.degC,(RH_fine_mix*100)*units.percent)
#
##===PBL WARMING========================================================
#T_fine_warm = np.copy(T_fine_mix)
#RW_fine_warm= np.copy(RW_fine_mix)
#P_fine_warm = np.copy(P_fine_mix)
#Q_fine_warm = np.copy(Q_fine_mix)
#
#for iz in inds_BL[0]:
#    if iz==0:
#        T_fine_warm[iz]=T_fine_warm[iz]+3#+3k
##       Get surface RW based on the surface RH after mixing 
#        RW_fine_warm[iz] = compute_rsat(T_fine_warm[iz],P_fine_warm[iz],1,T1,T2)*RH_fine_mix[iz]
#        Q_fine_warm[iz]=RW_fine_warm[iz]/(1+RW_fine_warm[iz])
#    else:
#        P_fine_warm[iz]= P_fine_warm[iz-1] - dz*(P_fine_warm[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_fine_warm[iz-1] )*T_fine_warm[iz-1] )
#        T_fine_warm[iz]= T_fine_warm[iz-1] + dz*drylift(T_fine_warm[iz-1],Q_fine_warm[iz-1],T_fine_warm[iz-1],Q_fine_warm[iz-1],0)
#        RW_fine_warm[iz]=RW_fine_warm[iz-1]
#        Q_fine_warm[iz]=RW_fine_warm[iz]/(1+RW_fine_warm[iz])
#
##===PUT BACK DEVIATION OF T&RW FROM WELL-MIXED PROFILE===================================
#RH_fine_warm = np.copy(RH_fine_mix)
#
#for iz in inds_BL[0]:
#    T_fine_warm[iz]=T_fine_warm[iz]-T_mix_dev[iz]
#    RW_fine_warm[iz]=RW_fine_warm[iz]-RW_mix_dev[iz]
#    Q_fine_warm[iz]=RW_fine_warm[iz]/(1+RW_fine_warm[iz])
#    
#for iz in range(1,len(Z_fine)):
#    P_fine_warm[iz]= P_fine_warm[iz-1] - dz*(P_fine_warm[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_fine_warm[iz-1] )*T_fine_warm[iz-1] )
#
#for iz in inds_BL[0]:
#    RH_fine_warm[iz]=RW_fine_warm[iz]/compute_rsat(T_fine_warm[iz],P_fine_warm[iz],1,T1,T2)
#    
#
##===TAKE PARCEL PROPERTIES IN BL THEN RECONSTRUCT FT&ST PROFILE===============================
#T_lif_warm,Qv_lif_warm,Qt_lif_warm,_=lift_parcel_adiabatic(T_fine_warm,P_fine_warm,Q_fine_warm,0,0,0,Z_fine,T1,T2)
#
#B_lif_warm=np.copy(B_lif_fine)
#for iz in inds_FT[0]:
#
#    if B_lif_fine[iz]>0:
#        B_lif_warm[iz]=B_lif_fine[iz]*np.exp(0.18)
#
#    qvv=(1-Qt_lif_warm[iz-1])*compute_rsat(T_lif_warm[iz-1],P_fine_warm[iz-1],0,T1,T2)
#    qvi=(1-Qt_lif_warm[iz-1])*compute_rsat(T_lif_warm[iz-1],P_fine_warm[iz-1],2,T1,T2)
#    T_lif_warm[iz] = T_lif_warm[iz-1] + dz*moislif(T_lif_warm[iz-1],Qv_lif_warm[iz-1],qvv,qvi,P_fine_warm[iz-1],T_fine_warm[iz-1],Q_fine_warm[iz-1],Qt_lif_warm[iz-1],0,0,T1,T2)
#    P_fine_warm[iz]= P_fine_warm[iz-1] - dz*(P_fine_warm[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_fine_warm[iz-1] )*T_fine_warm[iz-1] )
#    r = compute_rsat(T_lif_warm[iz],P_fine_warm[iz],1,T1,T2)
#    Qv_lif_warm[iz]=(1-Qt_lif_warm[iz-1])*r
#    Qt_lif_warm[iz]=Qt_lif_warm[iz-1]
#    def equations(x):
#        y = np.zeros_like(x)
#        
#        T_rho_par =  T_lif_warm[iz]*(1 + (Rv/Rd)*Qv_lif_warm[iz] - Qt_lif_warm[iz] )
##       RH_fine_warm is equilavent to RH_fine_mix above BL
#        T_rho_env = x[0]*(1 + (Rv/Rd - 1)*RH_fine_warm[iz]*compute_rsat(x[0],P_fine_warm[iz],1,T1,T2)/(1+RH_fine_warm[iz]*compute_rsat(x[0],P_fine_warm[iz],1,T1,T2)) )
##       fixed B_lif_fine         
#        y[0] = np.abs( B_lif_warm[iz] - g*(T_rho_par - T_rho_env)/T_rho_env )
#        return y
#    
#    # Initial guess for the solution
#    initial_guess = T_fine_warm[iz-1]#,q0[iz-1]
#    solution = fsolve(equations, initial_guess)
#    T_fine_warm[iz]=solution[0]
#    RW_fine_warm[iz] = RH_fine_warm[iz]*compute_rsat(T_fine_warm[iz],P_fine_warm[iz],1,T1,T2)
#    Q_fine_warm[iz] = RW_fine_warm[iz]/(1+RW_fine_warm[iz])
#
#for iz in inds_ST[0]:
#    P_fine_warm[iz]= P_fine_warm[iz-1] - dz*(P_fine_warm[iz-1]*g)/(Rd*(1 + (Rv/Rd - 1)*Q_fine_warm[iz-1] )*T_fine_warm[iz-1] )
#    T_fine_warm[iz]=T_fine_warm[iz-1]+dT_ST[iz]
##   RH_fine_warm is equilavent to RH_fine_mix above BL
#    RW_fine_warm[iz]=RH_fine_warm[iz]*compute_rsat(T_fine_warm[iz],P_fine_warm[iz],1,T1,T2)
#    Q_fine_warm[iz]=RW_fine_warm[iz]/(1+RW_fine_warm[iz])
#
#    qvv=(1-Qt_lif_warm[iz-1])*compute_rsat(T_lif_warm[iz-1],P_fine_warm[iz-1],0,T1,T2)
#    qvi=(1-Qt_lif_warm[iz-1])*compute_rsat(T_lif_warm[iz-1],P_fine_warm[iz-1],2,T1,T2)
#    T_lif_warm[iz] = T_lif_warm[iz-1] + dz*moislif(T_lif_warm[iz-1],Qv_lif_warm[iz-1],qvv,qvi,P_fine_warm[iz-1],T_fine_warm[iz-1],Q_fine_warm[iz-1],Qt_lif_warm[iz-1],0,0,T1,T2)
#    r = compute_rsat(T_lif_warm[iz],P_fine_warm[iz],1,T1,T2)
#    Qv_lif_warm[iz]=(1-Qt_lif_warm[iz-1])*r
#    Qt_lif_warm[iz]=Qt_lif_warm[iz-1]
#
#print("T_fine_warm:",compute_CAPE_AND_CIN(T_fine_warm,P_fine_warm,Q_fine_warm,0,0,0,Z_fine,T1,T2))
#
#Td_fine_warm=dewpoint_from_relative_humidity((T_fine_warm-273.15)*units.degC,(RH_fine_warm*100)*units.percent)
#
#Td_fine[inds_ST[0]]=np.nan
#Td_fine_warm[inds_ST[0]]=np.nan


#fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(4,4),dpi=300)
#ax.plot(B_lif_fine,Z_fine,np.zeros(10),np.arange(0,100000,10000))
#ax.plot(B_lif_warm,Z_fine,np.zeros(10),np.arange(0,100000,10000))
##ax.plot(np.arange(-0,300,50),np.arange(-0,300,50),'k-',lw=1)
##im=ax.scatter(vvm_w,w_cape_m*np.ones(vvm_w.shape),s=0.1,alpha=0.6,c=vvm_w_prec,vmax=50,vmin=0,cmap='jet',edgecolors=None,linewidths=0.)
##im=ax.scatter(w_vvm_m,w_cape_m,s=8,c=dcape_max,vmax=0.6,vmin=0,cmap='turbo',edgecolors='gray',linewidths=0.3)
#
##cbar=fig.colorbar(im)
##cbar.set_ticks(np.arange(0,1+.1,.1))
#
##[bar.set_alpha(0.2) for bar in bars]
##[cap.set_alpha(0.2) for cap in caps]
#
#ax.set_xlim([-0.1,1])
#ax.set_ylim([0,15000])
##ax.set_xticks(np.arange(-0,80.1,20))
##ax.set_yticks(np.arange(-0,80.1,20))
##ax.set_xlabel('W in VVM [m/s]',fontsize=8)
##ax.set_ylabel('W from CAPE [m/s]',fontsize=8)
##ax.tick_params(labelsize=5)
#
#plt.show()
##plt.savefig('./figure/w_cape_max_obj_prec.png')
#===PLOT

skew = SkewT()

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot
skew.plot(P_fine/100, T_lif_fine-273.15, 'k')
skew.plot(P_fine/100, T_fine-273.15, 'b')
skew.plot(P_fine/100, Td_fine, 'b',ls='--')

#skew.plot(P_fine_warm/100, T_lif_warm-273.15, 'k')
#skew.plot(P_fine_warm/100, T_fine_warm-273.15, 'r')
#skew.plot(P_fine_warm/100, Td_fine_warm, 'r',ls='--')

# Set some better labels than the default
skew.ax.set_xlabel('Temperature (\N{DEGREE CELSIUS})')
skew.ax.set_ylabel('Pressure (mb)')

# Add the relevant special lines
#skew.plot_dry_adiabats()
#skew.plot_moist_adiabats()
#skew.plot_mixing_lines()
skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-30, 35)
plt.show()
