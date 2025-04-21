import numpy as np
from scipy.special import lambertw
from scipy.special import sici, exp1
from CONSTANTS import *
from scipy.optimize import fsolve



#==============================================================================
#==============================================================================
#==============================================================================
#FUNCTION DEFINITIONS==========================================================
#==============================================================================
#==============================================================================
#==============================================================================


#==============================================================================
def comp_cdwave(F):
   gamma_em = np.euler_gamma
   return 4*(F*np.sin(2/F) - (F**2)*(np.sin(1/F)**2) - cosint(2/F) + np.log(2/F) - 1 + gamma_em )
#==============================================================================

#==============================================================================
def cosint(x):
   si, ci = sici(x)
   return ci
#==============================================================================


#==============================================================================
 #descriminator function between liquid and ice (i.e., omega defined in the
 #beginning of section 2e in Peters et al. 2022)
def omega(T,T1,T2):
    return ((T - T1)/(T2-T1))*np.heaviside((T - T1)/(T2-T1),1)*np.heaviside((1 - (T - T1)/(T2-T1)),1) + np.heaviside(-(1 - (T - T1)/(T2-T1)),1);
def domega(T,T1,T2):
    return (np.heaviside(T1-T,1) - np.heaviside(T2-T,1))/(T2-T1)
#==============================================================================

#==============================================================================
#FUNCTION THAT CALCULATES THE SATURATION MIXING RATIO
def compute_rsat(T,p,iceflag,T1,T2):
    
    #THIS FUNCTION COMPUTES THE SATURATION MIXING RATIO, USING THE INTEGRATED
    #CLAUSIUS CLAPEYRON EQUATION (eq. 7-12 in Peters et al. 2022).
    #https://doi-org.ezaccess.libraries.psu.edu/10.1175/JAS-D-21-0118.1 

    #input arguments
    #T temperature (in K)
    #p pressure (in Pa)
    #iceflag (give mixing ratio with respect to liquid (0), combo liquid and
    #ice (2), or ice (3)
    #T1 warmest mixed-phase temperature
    #T2 coldest mixed-phase temperature
    
    #NOTE: most of my scripts and functions that use this function need
    #saturation mass fraction qs, not saturation mixing ratio rs.  To get
    #qs from rs, use the formula qs = (1 - qt)*rs, where qt is the total
    #water mass fraction


    omeg = omega(T,T1,T2)
    if iceflag==0:
        term1=(cpv-cpl)/Rv
        term2=(xlv-ttrip*(cpv-cpl))/Rv
        esl=np.exp((T-ttrip)*term2/(T*ttrip))*eref*(T/ttrip)**(term1)
        qsat=epsilon*esl/(p-esl)
    elif iceflag==1: #give linear combination of mixing ratio with respect to liquid and ice (eq. 20 in Peters et al. 2022)
        term1=(cpv-cpl)/Rv
        term2=(xlv-ttrip*(cpv-cpl))/Rv
        esl_l=np.exp((T-ttrip)*term2/(T*ttrip))*eref*(T/ttrip)**(term1)
        qsat_l=epsilon*esl_l/(p-esl_l);
        term1=(cpv-cpi)/Rv
        term2=( xls-ttrip*(cpv-cpi))/Rv
        esl_i=np.exp((T-ttrip)*term2/(T*ttrip))*eref*(T/ttrip)**(term1);
        qsat_i=epsilon*esl_i/(p-esl_i)
        qsat=(1-omeg)*qsat_l + (omeg)*qsat_i
    elif iceflag==2: #only give mixing ratio with respect to ice
        term1=(cpv-cpi)/Rv
        term2=( xls-ttrip*(cpv-cpi))/Rv
        esl=np.exp((T-ttrip)*term2/(T*ttrip))*eref*(T/ttrip)**(term1)
        esl = min( esl , p*0.5 )
        qsat=epsilon*esl/(p-esl);
    return qsat
#==============================================================================



#==============================================================================
#FUNCTION THAT CALCULATES THE SATURATION MIXING RATIO
def compute_esat(T,iceflag):
    
    #THIS FUNCTION COMPUTES THE SATURATION MIXING RATIO, USING THE INTEGRATED
    #CLAUSIUS CLAPEYRON EQUATION (eq. 7-12 in Peters et al. 2022).
    #https://doi-org.ezaccess.libraries.psu.edu/10.1175/JAS-D-21-0118.1 

    #input arguments
    #T temperature (in K)
    #p pressure (in Pa)
    #iceflag (give mixing ratio with respect to liquid (0), combo liquid and
    #ice (2), or ice (3)
    #T1 warmest mixed-phase temperature
    #T2 coldest mixed-phase temperature
    
    #NOTE: most of my scripts and functions that use this function need
    #saturation mass fraction qs, not saturation mixing ratio rs.  To get
    #qs from rs, use the formula qs = (1 - qt)*rs, where qt is the total
    #water mass fraction


    if iceflag==0:
        term1=(cpv-cpl)/Rv
        term2=(xlv-ttrip*(cpv-cpl))/Rv
        esl=np.exp((T-ttrip)*term2/(T*ttrip))*eref*(T/ttrip)**(term1)
    elif iceflag==1: #only give mixing ratio with respect to ice
        term1=(cpv-cpi)/Rv
        term2=( xls-ttrip*(cpv-cpi))/Rv
        esl=np.exp((T-ttrip)*term2/(T*ttrip))*eref*(T/ttrip)**(term1)
    return esl
#==============================================================================


#==============================================================================
#LAPSE RATE FOR AN UNSATURATED PARCEL
def drylift(T,qv,T0,qv0,fracent):
    #CONSTANTS

    
    cpmv = (1 - qv)*cp + qv*cpv
    B = g*( (T-T0)/T0 + (Rv/Rd - 1)*(qv - qv0) )
    eps = -fracent*(T - T0)
    gamma_d = - (g + B)/cpmv + eps
    return gamma_d
#==============================================================================


#==============================================================================
#LIFTED CONDENSATION LEVEL USING THE ROMPS 2017 FORMULA
def compute_LCL(T,qv,p):
    #CONSTANTS


    cpm = (1 - qv)*cp + qv*cpv
    Rm = (1 - qv)*Rd + qv*Rv
    
    a = cpm/Rm + ( cpl - cpv )/Rv
    b = -(xlv - (cpv - cpl)*ttrip)/(Rv*T)
    c = b/a
    
    r_sat = compute_rsat(T,p,0,273.15,253.15)
    q_sat = r_sat/(1 + r_sat)
    RH = qv/q_sat
    arg1 = RH**(1/a)
    arg2 = c*np.exp(1)**c
    arg3 = lambertw(arg1*arg2,k=-1)
    T_LCL = c*T/arg3
    P_LCL = p*(T_LCL/T)**(cpm/Rm)
    Z_LCL = (cpm/g)*(T - T_LCL)
    
    return Z_LCL
#==============================================================================





#==============================================================================
#LIFTED CONDENSATION LEVEL USING NUMERICAL INTEGRATION
def compute_LCL_NUMERICAL(T,qv,p,dz):
    #CONSTANTS
    #NOTE, WE ARE ASSUMING ZERO BUOYANCY (I.E., WELL MIXED PBL)


    nfound_LCL = True
    
    zon = 0
    ind_hgt = 0
    Ton = T
    Qon = qv
    Pon = p
    while nfound_LCL:
        ind_hgt = ind_hgt+1
        Ton = Ton + dz*drylift(Ton,Qon,Ton,Qon,0)
        Pon = Pon - dz*(Pon*g)/(Rd*(1 + (Rv/Rd - 1)*Qon )*Ton )
        rsat = compute_rsat(Ton,Pon,0,273.15,253.15)
        qsat = rsat/(1 + rsat)
        if Qon >= qsat:
            nfound_LCL = False
    Z_LCL = ind_hgt*dz
    
    return Z_LCL
#==============================================================================








#==============================================================================
#LAPSE RATE FOR A SATURATED PARCEL
def moislif(T,qv,qvv,qvi,p0,T0,q0,qt,fracent,prate,T1,T2):
    
    #CONSTANTS

 
    qt=max(qt,0.0)
    qv=max(qv,0.0)
    
    OMEGA = omega(T,T1,T2)
    dOMEGA = domega(T,T1,T2)
    
    
    cpm = (1 - qt)*cp + qv*cpv + (1 - OMEGA)*(qt-qv)*cpl + OMEGA*(qt-qv)*cpi
    Lv = xlv + (T - ttrip)*(cpv - cpl)
    Li = (xls-xlv) + (T - ttrip)*(cpl - cpi);
    Rm0 = (1 - q0)*Rd + q0*Rv
    

    T_rho=T*(1 - qt + qv/epsilon)
    T_rho0=T0*( 1 - q0 + q0/epsilon )
    B = g*(T_rho - T_rho0)/(T_rho0)
    
    Qvsl = qvv/( epsilon - epsilon*qt + qv)
    Qvsi = qvi/( epsilon - epsilon*qt + qv)
    Q_M = (1 - OMEGA)*qvv/(1 - Qvsl) + OMEGA*qvi/(1 - Qvsi)
    L_M = Lv*(1 - OMEGA)*qvv/(1 - Qvsl) + (Lv + Li)*OMEGA*qvi/(1 - Qvsi)

    
    eps_T = -fracent*(T - T0)
    eps_qv = -fracent*(qv - q0)
    eps_qt = -fracent*(qt - q0)-prate*(qt-qv)
    term1 = -B
    
    term2 = - Q_M*(Lv + Li*OMEGA)*g/(Rm0*T0)
    
    term3 = -g
    term4 = (cpm - Li*(qt - qv)*dOMEGA)*eps_T
    term5 = (Lv + Li*OMEGA)*(eps_qv + (qv/(1-qt))*eps_qt)

    term6 = cpm
    term7 = -Li*(qt - qv)*dOMEGA
    term8 = (Lv + Li*OMEGA)*(-dOMEGA*(qvv - qvi) + (1/(Rv*(T**2)))*(L_M))
    gamma_m =( term1 + term2 + term3 + term4 + term5)/(term6 + term7 + term8)
    return gamma_m
#==============================================================================


#==============================================================================
#FUNCTION THAT LIFTS A PARCEL
def lift_parcel_adiabatic(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2):
    #[T_lif,Qv_lif,Qt_lif,B_lif]

    #this function computes lifted parcel properties using the unsaturated
    #and saturated lapse rate formulas from (Peters et al. 2022)
    #https://doi-org.ezaccess.libraries.psu.edu/10.1175/JAS-D-21-0118.1 
    
    #input arguments
    #T0: sounding profile of temperature (in K)
    #p0: sounding profile of pressure (in Pa)
    #q0: sounding profile of water vapor mass fraction (in kg/kg)
    #start_loc: index of the parcel starting location (set to 1 for the
    #lowest: level in the sounding)
    #fracent: fractional entrainment rate (in m^-1)
    
    #output arguments
    #T_lif: lifted parcel temperature
    #Qv_lif: lifted parcel water vapor mass fraction
    #Qt_lif: lifted parcel total water mass fraction
    #B_lif: Lifted parcel buoyancy, computed using Eq. B6 in (Peters et al.
    #2022) (accounts for virtual temperature and loading effects)
    
    #prate: precipitation rate (in m^-1) large values make parcel more
    #pseudoadiabatic, small values make parcel more adiabatic.  I usually
    #just set it to 0 to get an adiabatic parce
    
    #z0: sounding profile of height above ground level (first height should
    #be 0 m)
    #T1 warmest mixed-phase temperature
    #T2 coldest mixed-phase temperature

    #CONSTANTS

    
    #ESTIMATE THE MOIST STATIC ENERGY (MSE)
    MSE = cp*T0 + xlv*q0 + g*z0
    mn_hgt = np.min(np.where(MSE==np.nanmin(MSE))) #FIND THE INDEX OF THE HEIGHT OF MINIMUM MSE
    
    #descriminator function between liquid and ice (i.e., omega defined in the
    #beginning of section 2e in Peters et al. 2022)

    
    T_lif=np.zeros(T0.shape)*np.nan #temperature of the lifted parcel
    Qv_lif=np.zeros(T0.shape)*np.nan #water vapor mass fraction of the lifted parcel
    Qt_lif=np.ones(T0.shape)*np.nan #total water mass fraction of the lifted parcel

    if start_loc>0:
        T_lif[0:start_loc+1]=T0[0:start_loc+1] #set initial values to that of the environment
        Qv_lif[0:start_loc+1]=q0[0:start_loc+1] #set initial values to that of the environment
        Qt_lif[0:start_loc+1]=Qv_lif[0:start_loc+1] #set initial values to that of the environment
    else:
        T_lif[0]=T0[0] #set initial values to that of the environment
        Qv_lif[0]=q0[0] #set initial values to that of the environment
        Qt_lif[0]=Qv_lif[0] #set initial values to that of the environment


    q_sat_prev=0
    B_run = 0
    iz=start_loc
    #
    #for iz in np.arange(start_loc+1,z0.shape[0]):
    #
    #
    #I REVISED THIS A BIT.  TO MAKE THE CODE FASTER, I HAVE THE CALCULATION CUT OUT WHEN THE INTEGRATED NEGATIVE BUOYANCY ("BRUN") 
    #BECOMES MORE NEGATIVE THAN THAN THE TOTAL INTEGRATED POSITIVE BUOYANCY.  I RESTRICT THIS TO ONLY HAPPEN AFTER WE HAVE PASSED 
    #THE HEIGHT OF MINIMUM MSE.  UNCOMMENT THE FOR LOOP ABOVE AND COMMENT OUT THE WHILE LOOP IF YOU JUST WANT TO INTEGRATE TO THE TOP OF THE SOUNDING.
    #THE +25 PART IN THE WHILE STATEMENT IS A PAD ON B_RUN (THE NEGATIVE CAPE HAS TO BE 25 J/KG LESS THAN THE POSITIVE CAPE TO KILL THE LOOP)
    #while iz<(z0.shape[0])-1 and (z0[iz]<z0[mn_hgt] or (B_run+25)>0):
    while iz<(z0.shape[0])-1 and (z0[iz]<z0[mn_hgt] or (B_run+25000)>0):
        iz = iz + 1
        q_sat=(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],1,T1,T2)
        if Qv_lif[iz-1]<q_sat: #if we are unsaturated, go up at the unsaturated adiabatic lapse rate (eq. 19 in Peters et al. 2022)
            
        
        
            T_lif[iz] = T_lif[iz-1] + (z0[iz] - z0[iz-1])*drylift(T_lif[iz-1],Qv_lif[iz-1],T0[iz-1],q0[iz-1],fracent)
            Qv_lif[iz] = Qv_lif[iz-1] - (z0[iz] - z0[iz-1])*fracent*( Qv_lif[iz-1] - q0[iz-1] )
            Qt_lif[iz] = Qv_lif[iz]
            q_sat=(1-Qt_lif[iz])*compute_rsat(T_lif[iz],p0[iz],1,T1,T2)
            
            if Qv_lif[iz]>=q_sat: #if we hit saturation, split the vertical step into two stages.  The first stage advances at the saturated lapse rate to the saturation point, and the second stage completes the grid step at the moist lapse rate
                OMEGA = omega(T_lif[iz-1],T1,T2)
                dOMEGA = domega(T_lif[iz-1],T1,T2)
                satrat=(Qv_lif[iz]-q_sat_prev)/(q_sat-q_sat_prev)
                dz_dry=satrat*(z0[iz]-z0[iz-1])
                dz_wet=(1-satrat)*(z0[iz]-z0[iz-1])


                
                T_halfstep = T_lif[iz-1] + dz_dry*drylift(T_lif[iz-1],Qv_lif[iz-1],T0[iz-1],q0[iz-1],fracent)
                Qv_halfstep = Qv_lif[iz-1] - dz_dry*fracent*( Qv_lif[iz-1] - q0[iz-1] )
                Qt_halfstep = Qv_lif[iz]
                p_halfstep=p0[iz-1]*satrat + p0[iz]*(1-satrat)
                T0_halfstep=T0[iz-1]*satrat + T0[iz]*(1-satrat)
                Q0_halfstep=q0[iz-1]*satrat + q0[iz]*(1-satrat)

                T_lif[iz] = T_halfstep + dz_wet*moislif(T_halfstep,Qv_halfstep,(1-Qt_halfstep)*compute_rsat(T_halfstep,p_halfstep,0,T1,T2),(1-Qt_halfstep)*compute_rsat(T_halfstep,p_halfstep,2,T1,T2),p_halfstep,T0_halfstep,Q0_halfstep,Qt_halfstep,fracent,prate,T1,T2)
                
                
                Qt_lif[iz] = Qt_lif[iz-1] - (z0[iz] - z0[iz-1])*fracent*( Qt_halfstep - Q0_halfstep )
                Qv_lif[iz] = (1-Qt_lif[iz])*compute_rsat(T_lif[iz],p0[iz],1,T1,T2)

                if Qt_lif[iz]<Qv_lif[iz]:
                    Qv_lif[iz]=Qt_lif[iz]

            q_sat_prev=q_sat;
            
        else: #if we are already at saturation, just advance upward using the saturated lapse rate (eq. 24 in Peters et al. 2022)
            OMEGA = omega(T_lif[iz-1],T1,T2)
            dOMEGA = domega(T_lif[iz-1],T1,T2)

            T_lif[iz] = T_lif[iz-1] + (z0[iz] - z0[iz-1])*moislif(T_lif[iz-1],Qv_lif[iz-1],(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],0,T1,T2),(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],2,T1,T2),p0[iz-1],T0[iz-1],q0[iz-1],Qt_lif[iz-1],fracent,prate,T1,T2);
                     
             
            Qt_lif[iz] = Qt_lif[iz-1] - (z0[iz] - z0[iz-1])*(fracent*( Qt_lif[iz-1] - q0[iz-1] )  + prate*( Qt_lif[iz-1]-Qv_lif[iz-1]) )
            Qv_lif[iz] = (1-Qt_lif[iz])*compute_rsat(T_lif[iz],p0[iz],1,T1,T2)
            
            if Qt_lif[iz]<Qv_lif[iz]:
                Qv_lif[iz]=Qt_lif[iz]

        B_run = B_run + (g*T_lif[iz]*(1 + (Rv/Rd)*Qv_lif[iz] - Qt_lif[iz])/(T0[iz]*(1 + (Rv/Rd)*q0[iz] - q0[iz])) - g)*(z0[iz]-z0[iz-1])

    T_rho_lif = T_lif*(1 + (Rv/Rd)*Qv_lif - Qt_lif)
    T_0_lif = T0*(1 + (Rv/Rd - 1)*q0)
    #T_rho_lif=T_lif*(1 - Qt_lif + Qv_lif)/( 1 + (epsilon - 1)/( ( epsilon*(1 - Qt_lif)/Qv_lif - 1) ) )
    #T_0_lif=T0/( 1 + (epsilon - 1)/( ( epsilon*(1 - q0)/q0 - 1) ) )
    
    B_lif=g*(T_rho_lif - T_0_lif)/T_0_lif
   
    
    return T_lif,Qv_lif,Qt_lif,B_lif

#==============================================================================
#FUNCTION THAT COMPUTES CAPE, CIN, EL, LFC
def compute_CAPE_AND_CIN(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2,useenergy=True,useleapfrog=False):
#[CAPE,CIN,LFC,EL]

    #this function computes CAPE and CIN
    
    #input arguments
    #T0: sounding profile of temperature (in K)
    #p0: sounding profile of pressure (in Pa)
    #q0: sounding profile of water vapor mass fraction (in kg/kg)
    #start_loc: index of the parcel starting location (set to 1 for the
    #lowest: level in the sounding)
    #fracent: fractional entrainment rate (in m^-1)
    
    #CONSTANTS

    
    #compute lifted parcel buoyancy
    if useenergy:
        
        if useleapfrog:
            T_lif,Qv_lif,Qt_lif,B_lif=lift_parcel_adiabatic_leapfrog(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2)
        else:
            T_lif,Qv_lif,Qt_lif,B_lif=lift_parcel_adiabatic(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2)       
        #T_lif,Qv_lif,Qt_lif,B_lif=lift_parcel_adiabatic_leapfrog(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2)

    else:
        T_lif,Qv_lif,Qt_lif,B_lif=lift_parcel_ME(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2)
    
    if np.nanmax(B_lif)>0:
        #CAPE will be the total integrated positive buoyancy
        B_pos = np.zeros(B_lif.shape)
        B_pos[:] = B_lif[:]
        B_pos[np.where(B_pos<0)]=0
        dz = z0[1:z0.shape[0]] - z0[0:z0.shape[0]-1]
        CAPE = np.nansum( 0.5*B_pos[0:z0.shape[0]-1]*dz + 0.5*B_pos[1:z0.shape[0]]*dz )
        
        #CIN will be the total negative buoyancy below the height of maximum
        #buoyancy
        B_neg = np.zeros(B_lif.shape)
        B_neg[:] = B_lif[:]
        mx = np.nanmax(B_lif)
        imx = np.where(B_lif==mx)
        imx=imx[0][0]
        B_neg[0:imx]=np.minimum( B_neg[0:imx], 0 )
        B_neg[imx:z0.shape[0]]= 0
        CIN = np.nansum( 0.5*B_neg[0:z0.shape[0]-1]*dz + 0.5*B_neg[1:z0.shape[0]]*dz )
        
        #LFC will be the last instance of negative buoyancy before the
        #continuous interval that contains the maximum in buoyancy
        fneg = np.where(B_lif<0)
        fneg=fneg[0]
        inn = np.where(fneg<imx)
        inn = inn[0]
        fneg = fneg[inn]
        if len(fneg)>0:
            LFC = 0.5*z0[np.max(fneg)] + 0.5*z0[min(np.max(fneg)+1,z0.shape[0]-1)]
        else:
            LFC = z0[start_loc]
        
        #EL will be last instance of positive buoyancy
        fpos = np.where(B_lif>0)
        fpos=fpos[0]
        EL = 0.5*z0[np.max(fpos)] + 0.5*z0[min(np.max(fpos)+1,z0.shape[0])]
    else:
        CAPE = 0
        CIN = 0
        LFC = np.nan
        EL = np.nan

    return CAPE,CIN,LFC,EL


#==============================================================================
#FUNCTION THAT COMPUTES NCAPEf
def compute_NCAPE(T0,p0,q0,z0,T1,T2,LFC,EL):


    
    #COMPUTE THE MOIST STATIC ENERGY
    MSE0 = cp*T0 + xlv*q0 + g*z0
    
    #COMPUTE THE SATURATED MOIST STATIC ENERGY
    rsat = compute_rsat(T0,p0,0,T1,T2)
    qsat = (1 - rsat)*rsat
    MSE0_star = cp*T0 + xlv*qsat + g*z0
    
    #COMPUTE MSE0_BAR
    MSE0bar=np.zeros(MSE0.shape)
    #for iz in np.arange(0,MSE0bar.shape[0],1):
     #   MSE0bar[iz]=np.mean(MSE0[1:iz])
        
    MSE0bar[0]=MSE0[0]
    for iz in np.arange(1,MSE0bar.shape[0],1):
        MSE0bar[iz] = 0.5*np.sum( (MSE0[0:iz] + MSE0[1:iz+1])*(z0[1:iz+1]-z0[0:iz]) )/(z0[iz]-z0[0])
    
    int_arg = - ( g/(cp*T0) )*( MSE0bar - MSE0_star)
    ddiff = abs(z0-LFC)
    mn = np.min(ddiff)
    ind_LFC = np.where(ddiff==mn)[0][0]
    ddiff = abs(z0-EL)
    mn = np.min(ddiff)
    ind_EL = np.where(ddiff==mn)[0][0]
    #ind_LFC=max(ind_LFC);
    #ind_EL=max(ind_EL);
    
    NCAPE = np.maximum(np.nansum( (0.5*int_arg[ind_LFC:ind_EL-1] + 0.5*int_arg[ind_LFC+1:ind_EL] )*(z0[ind_LFC+1:ind_EL] - z0[ind_LFC:ind_EL-1] ) ),0)
    return NCAPE,MSE0_star,MSE0bar


#==============================================================================
#FUNCTION THAT COMPUTES BUNKERS SR MOTION
def compute_VSR(z0,u0,v0,propfac_coef=1):
    #compute 0-1 km storm-relative flow (V_SR) using the storm motion
    #estimate of Bunkers et al. (2000)
    #https://doi.org/10.1175/1520-0434(2000)015<0061:PSMUAN>2.0.CO;2
    
    
    f6000 = np.where(z0<=6000)[0]
    f6001 = f6000+1
    meanx=np.sum(u0[f6000]*(z0[f6001]-z0[f6000]))/np.sum(z0[f6001]-z0[f6000])
    meany=np.sum(v0[f6000]*(z0[f6001]-z0[f6000]))/np.sum(z0[f6001]-z0[f6000])
    
    f0500 = np.where(z0<=500)[0]
    lowx=np.mean(u0[f0500])
    lowy=np.mean(v0[f0500])
    
    f560 = np.where(np.logical_and(z0<=6000,z0>=5500))[0]
    highx=np.mean(u0[f560])
    highy=np.mean(v0[f560])
    BK_SHRx=highx-lowx
    BK_SHRy=highy-lowy
    BK_mag=np.sqrt(BK_SHRx**2 + BK_SHRy**2)
    BK_dirx=BK_SHRx/BK_mag
    BK_diry=BK_SHRy/BK_mag
    BK_orthx=BK_diry*7.5
    BK_orthy=-BK_dirx*7.5


    SR_mean_u= u0 - meanx
    SR_mean_v= v0 - meany
    dudz=np.zeros(u0.shape)
    dvdz=np.zeros(v0.shape)
    dudz[1:dudz.shape[0]-1]= ( u0[2:dudz.shape[0]]-u0[0:dudz.shape[0]-2] )/( z0[2:dudz.shape[0]]-z0[0:dudz.shape[0]-2] )
    dudz[0]=2*dudz[1]-dudz[2]
    dvdz[1:dudz.shape[0]-1]= ( v0[2:dudz.shape[0]]-v0[0:dudz.shape[0]-2] )/( z0[2:dudz.shape[0]]-z0[0:dudz.shape[0]-2] )
    dvdz[0]=2*dvdz[1]-dvdz[2]
    f1000 = np.where(z0<=1000)[0]
    SRH_mean = abs(np.mean(-SR_mean_u[f1000]*dvdz[f1000] + SR_mean_v[f1000]*dudz[f1000])*1000.0)
    
    
    if SRH_mean>0:
        sign_SRH = SRH_mean/abs(SRH_mean)
        propfac= sign_SRH*propfac_coef*min(abs(SRH_mean)/150,1)
        #propfac = 1
        if propfac_coef ==2:
            propfac = 1
    else:
        propfac = 0


    C_x=meanx+propfac*BK_orthx
    C_y=meany+propfac*BK_orthy
    
    u_sr = u0 - C_x
    v_sr = v0 - C_y
    
    f1000 = np.where(z0<=1000)[0]
    V_SR = np.nanmean(np.sqrt(  u_sr[f1000]**2 + v_sr[f1000]**2  ))
    return V_SR,C_x,C_y


#==============================================================================
def compute_ETILDE(CAPE,NCAPE,V_SR,EL,L,squall_line=False,dynamfac=1):
    #THESE ARE A BUNCH OF CONSTANT PARAMTERS SET FOR THE ECAPE CALCULATION
    H=EL
    l=L/H
    sigma = 1.1
    alpha=0.8

    if squall_line:
        pitchfork=ksq*(alpha**2)*L/(Pr*(sigma**2)*H)
    else:
        pitchfork=ksq*(alpha**2)*(np.pi**2)*L/(4*Pr*(sigma**2)*H)
    vsr_tilde = V_SR/np.sqrt(2*CAPE)
    N_tilde = NCAPE/CAPE
    
    #EQUATION SOLVES FOR THE NONDIMENSIONAL ECAPE (E_TILDE_A IN THE PAPER)
    E_tilde = dynamfac*vsr_tilde**2 + ( -1 - pitchfork - (pitchfork/(vsr_tilde**2 ))*N_tilde + \
                              np.sqrt((1 + pitchfork + (pitchfork/(vsr_tilde**2 ))*N_tilde)**2 + \
                                      (4*(pitchfork/(vsr_tilde**2 ))*(1 - pitchfork*N_tilde) ) ) )/( 2*pitchfork/(vsr_tilde**2) )
        
    E_tilde_ = E_tilde - dynamfac*vsr_tilde**2
        
   # varepsilon = 0.65*2*((1 - E_tilde_)/(E_tilde_ + N_tilde))/(EL)  
    varepsilon = 2*((1 - E_tilde_)/(E_tilde_ + N_tilde))/(EL)
    

    #eps = 2*ksq*L/(EL*Pr)
    
    #Rm2 = ( (alpha*np.pi/(sigma) )**2 )*( E_tilde/vsr_tilde + 1)
    #Radius =  EL*Rm2**(-1/2)
    #varepsilon = 2*ksq*L/(Pr*Radius**2 )
    
    #Radius=Radius/2
    
    #varepsilon = 0.65*eps*(alpha**2)*(np.pi**2)*E_tilde/(4*(sigma**2)*EL*(vsr_tilde**2 ) ) #THIS IS THE FRACTIONAL ENTRAINMENT RATE
    Radius = np.sqrt(2*ksq*L/(Pr*varepsilon))

    return E_tilde,varepsilon,Radius

#==============================================================================
def CI_model(T0,p0,q0,z0,u0,v0,T1,T2,radrng,itmax,L,prate_global):
    
    #THIS FUNCTION EXECUTES THE "PROGRESSIVE ROOTING" TOY MODEL DESCRIBED BY PETERS ET AL. 2022A
    #https://journals.ametsoc.org/view/journals/atsc/79/6/JAS-D-21-0145.1.xml
    
    #NOTE, A VAREITY OF THINGS HAVE CHANGED SINCE THAT PUBLICATION.  I WILL TRY TO 
    #POINT SPECIFIC EQUATIONS HERE TO EQUATION NUMBERS IN THE PUBLCIATION.  I WILL PROBABLY
    #CREATE A TECHNICAL DOCUMENT TO DESCRIBE THESE CHANGES SOMETIME SOON.  STAY TUNED...
    
    #THE FUNCTION TAKES AS INPUT:
        #T0, profile of temperature (K)
        #p0, profile of pressure (Pa)
        #q0, profile of specific humidity (kg/kg)
        #z0, profile of height above ground level (m)
        #u0, profile of u wind (m/s)
        #v0, profile of v wind (m/s)
        #T1, temperature at which freezing begins in parcel calculations (I usually set to 273.15 K)
        #T2, temperature at which freezing ends in the parcel calculation (K).  This will control the temperature
            #range over which mixed-phase occurs.  I usually set to 253.15 k
        #radrng, a vector containing the initial radii we are going to test.  A reasonable
            #choice here is a range from 100 m to 6000 m at intervals of 100 m (np.arange(100,6000,100))
        #itmax, the number of iterations (I usually set to 20)
        #L, the mixing length (I usually set to 250 m)
        #prate_global, the precipitation loss inverse length scale (km^(-1)).  Larger values make the
            #parcel more pseudoadiabatic, smaller values make it more adiabatic.
    
    #STANDARD THERMODYNAMIC CONSTANTS

    

    #PARAMTERS UNIQUE TO THE CI MODEL
    alpha=0.8 #ASSUMED RATIO OF HORIZONTALLY AVERAGED W TO HORIZONTAL MAX OF W AT A GIVEN LEVEL
    Pr=1/3 #PRANDTL NUMBER
    ksq=0.18 #VON KARMAN CONSTANT
    start_loc = 0 #STARTING HEIGHT OF THE AIR PARCEL WE ARE LIFTING
    sig = 0.5 #RATIO OF THE HEIGHT OF WMAX TO EQUILBIRIUM LEVEL HEIGHT (SHOULD PROBABLY SET THIS TO 1)
    rfac = 1/4 #RELAXATION FACTOR FOR MODEL INTEGRATION.  SMALLER VALUE GIVES A SMOOTHER SOLUTION
    
    #WE WILL NEED THE DENSITY PROFILE TO COMPUTE THE STORM-RELATIVE WIND LATER
    rho0 = p0/(Rd*T0*(1 + (Rv/Rd - 1)*q0))

    #TIME SERIES OF QUANTITIES OUTPUTTED FROM THE CI MODEL
    R_TS = np.zeros((radrng.shape[0],itmax)) #RADIUS OF THE UPDRAFT
    H_TS = np.zeros((radrng.shape[0],itmax)) #EL HEIGHT
    W_TS = np.zeros((radrng.shape[0],itmax)) #MAX VERTICAL VELOCITY
    VSR_TS = np.zeros((radrng.shape[0],itmax)) #STORM-RELATIVE FLOW
    
    #INITIAL CONDITION ON RADIUS: SET TO R0
    R_TS[:,0]=radrng 
    
    
    #dudz=np.zeros(u0.shape)
    #dvdz=np.zeros(v0.shape)
    #dudz[1:dudz.shape[0]-1]= ( u0[2:dudz.shape[0]]-u0[0:dudz.shape[0]-2] )/( z0[2:dudz.shape[0]]-z0[0:dudz.shape[0]-2] )
    #dudz[0]=2*dudz[1]-dudz[2]
    #dvdz[1:dudz.shape[0]-1]= ( v0[2:dudz.shape[0]]-v0[0:dudz.shape[0]-2] )/( z0[2:dudz.shape[0]]-z0[0:dudz.shape[0]-2] )
    #dvdz[0]=2*dvdz[1]-dvdz[2]
    #SHR_mag = np.sqrt(dudz**2 + dvdz**2)
    
    #IN THE FUTRE, WE'LL PROBABLY WANT TO COMPUTE THE DENSITY WEIGHTED STORM-RELATIVE FLOW, LIKE IN THE ECAPE THEORY
    #UDCAPE,UDCIN,UDLFC,UDEL=compute_CAPE_AND_CIN(T0,p0,q0,start_loc,0,prate_global,z0,T1,T2)

    #PARAMETERS FOR CI MODEL
    for it in np.arange(0,itmax-1,1): #LOOP THROUGH THE SPECIFIED NUMBER OF ITERATIONS
        for ir in np.arange(0,radrng.shape[0],1): #LOOP THROUGH EACH OF THE STARTING RADII
            R_on = R_TS[ir,it] #STORE THE RADIUS (IN M)
            #
            fracent = 2*ksq*L/(Pr*(R_on**2)) #USE RADIUS TO COMPUTE FRACTION ENTRAINMENT RATE WITH EQ. XX IN XX
            
            #WHEN COMPUTING THE VERTICAL PROFILE OF KINETIC ENERGY, THE LOWER BOUNDARY CONDITION IS THAT A PARCEL
            #BEGINS WITH THE KINETIC ENERGY OF THE INFLOW.  THIS MEANS WE HAVE TO GIVE THE VERTICAl VELOCITY
            #FUNCTION THE STORM RELATIVE WIND.
            if it == 0:
                #AT THE FIRST TIME STEP, WE WONT HAVE THE STORM RELATIVE WIND YET SO WE'LL MAKE AN AD-HOC ESTIMATE
                #V_SR = 15*R_on/5000 #5.0  
                V_SR = 20*R_on/5000 #5.0  
            else:
                #AT LATER TIMES, WE JUST USE THE STORM RELATIVE FLOW FROM THE PREVIOUS TIME STEP
                V_SR = VSR_TS[ir,it-1]
                    
            #GET THE MAXIMUM VERTICAL VELOCITY PROFILE FOR A RISING CLOUD THERMAL
            CAPE,LFC,EL,B_pos=compute_w(T0,p0,q0,start_loc,fracent,prate_global,z0,T1,T2,R_on,u0,v0,V_SR)
            #NOW THE VERTICAL VELOCITY AT THE BASE OF THE THERMAL, WHICH WILL EXPERIENCE A HIGHER ENTRAINMENT RATE
            CAPE2,LFC2,EL2,null=compute_w(T0,p0,q0,start_loc,fracent*9/4,prate_global,z0,T1,T2,R_on,u0,v0,V_SR)
            
            #WE WILL NEED THE PROFILE OF POSITIVE BUOYANCY TO ESTIMATE STORM MOTION LATER.  
            #ZERO OUT THE NEGATIVE BUOYANCY
            B_pos = np.maximum(B_pos,0)

            #IF WE ACTUALLY HAVE ANY POSITIVE BUOYANCY, WE'LL ADVANCE THE MODEL
            if ~np.isnan(EL):
                
                #GET THE 0-1 KM STORM-RELATIVE FLOW
                V_SR = compute_VSR_DIFF(z0,u0,v0,rho0,EL,B_pos)
                #V_SR = compute_VSR(z0,u0,v0)
                     
                #ADVANCE TO THE NEXT RADIUS USING EQ XX IN XX
                R_next = 1.7*( (EL2/EL)**2 )*2*V_SR*(EL-LFC)*sig/(np.pi*alpha*np.sqrt(2*CAPE))
                R_next = (rfac)*R_next + (1-rfac)*R_on #RELAXATION PROCEEDURE
                
            else: #OTHERWISE SET THE RADIUS AT THE NEXT TIME TO ZERO
                R_next = 0
                
            if EL<LFC: #THIS HAPPENS SOMETIMES.  SET TO ZERO IF EL IS LESS THAN LFC
                R_next = 0
                
            #STORE TIME SERIES'
            R_TS[ir,it+1] = R_next
            H_TS[ir,it]=EL
            W_TS[ir,it]=np.sqrt(2*CAPE)
            VSR_TS[ir,it]=V_SR
        R_TS[np.where(np.isnan(R_TS))]=0
    
    return R_TS,H_TS,W_TS,VSR_TS
            
    
#==============================================================================
#FUNCTION THAT COMPUTES CAPE, CIN, EL, LFC
def compute_w(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2,Radius,u0,v0,V_SR,include_arg=False):
#[CAPE,CIN,LFC,EL]

    #this function computes CAPE and CIN
    
    #input arguments
    #T0: sounding profile of temperature (in K)
    #p0: sounding profile of pressure (in Pa)
    #q0: sounding profile of water vapor mass fraction (in kg/kg)
    #start_loc: index of the parcel starting location (set to 1 for the
    #lowest: level in the sounding)
    #fracent: fractional entrainment rate (in m^-1)
    
    #CONSTANTS
    c_d = 0.2 #DRAG COEFICIENT ON A SPHERE
    Lambda=0.6 #RATIO OF ASCENT RATE OF THERMAL TO ITS MAX W
    alpha=0.8 #ASSUMED RATIO OF HORIZONTALLY AVERAGED W TO HORIZONTAL MAX OF W AT A GIVEN LEVEL
    
    #COMPUTE A VERTICAL PROFILE OF THE MAGNITUDE OF VERTICAL WIND SHEAR
    dz = np.zeros(u0.shape)
    dz[0:u0.shape[0]-1]=z0[1:u0.shape[0]]-z0[0:u0.shape[0]-1]
    dudz = np.zeros(u0.shape)
    dvdz = np.zeros(u0.shape)
    dudz[0:dudz.shape[0]-1]=(u0[1:dudz.shape[0]]-u0[0:dudz.shape[0]-1])/dz[0:dudz.shape[0]-1]
    dvdz[0:dudz.shape[0]-1]=(v0[1:dudz.shape[0]]-v0[0:dudz.shape[0]-1])/dz[0:dudz.shape[0]-1]  
    S = np.sqrt( dudz**2 + dvdz**2)                                            
    
    #COMPUTE THE LIFTED PARCEL BUOYANCY
    T_lif,Qv_lif,Qt_lif,B_lif=lift_parcel_adiabatic(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2)
    
    #CALCULATE THE LIFTED CONDENSATION LEVEL
    qdiff = abs(Qt_lif - Qv_lif) #FIGURE OUT THE FIRST HEIGHT WHERE QV STARTS DEVIATING FROM QT, IMPLYING CONDENSATION
    if np.logical_and(~np.isnan(qdiff[1]),np.nanmax(qdiff)>0):
        lcl_ind = np.where(qdiff>0)[0][0]
        LCL = z0[lcl_ind]
    else:
        LCL = 1000
        lcl_ind = np.where(abs(LCL-z0)==np.amin(abs(LCL-z0)))[0][0]
    
    #IF WE HAVE SOME POSITIVE BUOYANCY, PROCEED
    if np.nanmax(B_lif)>0:
        #MAKE A NEW MATRIX THAT WILL ONLY CONTAIN THE POSITIVE PART OF BUOYANCY
        B_pos = np.zeros(B_lif.shape)
        B_pos[:] = B_lif[:]

        #GET RID OF ALL NEGATIVE BUOYANCY BELOW THE LCL
        B_pos[0:lcl_ind]=0
        wpos = np.where(B_pos>0)[0]
        if len(wpos)>0:
            wpos=wpos[0] #WPOS CONTAINS INDEX OF LCL.  SET TO 0 IF THERE IS NO POSTIVE BUOYANCY
        else:
            wpos=lcl_ind
        B_pos[0:wpos]=0
        dz = z0[1:z0.shape[0]] - z0[0:z0.shape[0]-1]      
  
        #LFC WILL BE THE LAST INSTANCE OF NEGATIVE BUOYANCY BEFORE THE PARCEL REACHES ITS CONTINUOUS INTERVAL OF POSITIVE BUOY
        mx = np.nanmax(B_lif)
        imx = np.where(B_lif==mx)
        imx=imx[0][0]
        
        fneg = np.where(B_lif<0)
        fneg=fneg[0]
        inn = np.where(fneg<imx)
        
        inn = inn[0]
        fneg = fneg[inn]
        if len(inn)>0:
            LFC = 0.5*z0[np.max(fneg)] + 0.5*z0[np.max(fneg)+1]
        else:
            LFC = z0[start_loc]
        
        #EL WILL BE THE LAST INSTANCE OF POSITIVE BUOYANCY
        fpos = np.where(B_lif>0)
        fpos=fpos[0]
        EL = 0.5*z0[np.max(fpos)] + 0.5*z0[np.max(fpos)+1]
        
        #INTIALIZE PROFILE OF SQUARED VERTICAL VELOCITY (I.E., VERTICAL KINETIC ENERGY)
        WSQ_prof = np.zeros(B_pos.shape[0])
        WSQ_prof[start_loc]=(V_SR**2)/2 #LOWER BOUNADRY CONDITION ON VERTICAL KE IS THE KE OF INFLOW
        uprime_prof = np.zeros(B_pos.shape[0]) #INITIALIZE UPRIME PROFILE
        for iz in np.arange(0,WSQ_prof.shape[0]-1,1): #VERTICALLY INTEGRATE
            B_on = B_pos[iz] #STORE THE CURRENT BUOYANCY
            ebuoy_fac = 1/(1 + 2*(alpha**2)*(Radius**2)/((EL-LFC)**2 ) ) #SCALE FACTOR THAT ACCOUNTS FOR EFFECITVE BUOYANCY
            #ns_drag = -2.5*c_d*(3/8)/Radius #COEFICIENT ON THE NON-SHEARED PART OF DRAG
            ns_drag = -c_d*(3/8)/Radius #COEFICIENT ON THE NON-SHEARED PART OF DRAG
            s_drag = -( c_d/Radius )*(1 - Lambda)/(Lambda**2) #COEFICIENT ON THE SHEARED PART OF DRAG
            sh_drag = (  1/(0.5*np.sqrt(2*WSQ_prof[iz-1])) )*(3*c_d/(8*Radius)) #SHEARED DRAG TERM
            
            
            if np.sqrt(2*WSQ_prof[iz-1])<1: #IF WE HAVE VERY SMALL VERTICAL VELOICTY (LESS THAN 1 M/S, WE NEED TO ZERO OUT THE SHEAR DRAG TERM OR THINGS BLOW UP)
                sh_drag = 0
                
            #NOW VERTICALLY INTEGRATE THE UPRIME AND WSQ EQUATIONS TOGETHER, FOLLOWING EQ. XX AND XX IN XX RESPECTIVELY
            uprime_prof[iz+1] = uprime_prof[iz-1] + ( z0[iz+1]-z0[iz] )*(-sh_drag*uprime_prof[iz]**2 + S[iz] )
            WSQ_prof[iz+1] = WSQ_prof[iz] + ( z0[iz+1]-z0[iz] )*(ebuoy_fac*B_on + ns_drag*WSQ_prof[iz] + s_drag*uprime_prof[iz]*np.sqrt(2*WSQ_prof[iz]))
          
        #WE WILL OUTPUT THE MAXIMUM KE AS THE "CAPE" ARGUMENT                                                                                                          
        CAPE = np.nanmax(WSQ_prof)
        
        #SET THE EL TO THE HEIGHT OF MAXIMUM VERTICAL VELOCITY
        mxval = np.nanmax(WSQ_prof)
        fnval=np.where(WSQ_prof==mxval)
        LFC = LCL
        if fnval[0].shape[0]>0:
            EL = z0[fnval[0][0]]
            ELtemp = T0[fnval[0][0]]
            imaxheight = np.max(np.where(WSQ_prof>0)[0])
            zmax = z0[imaxheight]
            TOPtemp = T_lif[imaxheight]
        else:
            EL = np.nan
            ELtemp = np.nan
            zmax = np.nan
            TOPtemp = np.nan
        
    else:
        #IF WE HAVE NO POSITIVE BUOYANCY, SET EVERYTHING TO 0S AND NANS
        CAPE = 0      
        LFC = np.nan
        EL = np.nan
        B_pos = np.zeros(T0.shape)
        ELtemp = np.nan
        zmax = np.nan
        TOPtemp = np.nan
            
    if include_arg:
        return CAPE,LFC,EL,B_pos,ELtemp,zmax,TOPtemp
    else:
        return CAPE,LFC,EL,B_pos


#==============================================================================
#FUNCTION THAT COMPUTES BUNKERS SR MOTION
def compute_VSR_DIFF(z0,u0,v0,rho0,EL,B_pos):
    #compute 0-1 km storm-relative flow (V_SR) using the storm motion
    #estimate of Bunkers et al. (2000)
    #https://doi.org/10.1175/1520-0434(2000)015<0061:PSMUAN>2.0.CO;2
    
    zdiff = ( z0 - EL )**2
    ind_top = np.where(zdiff==np.min(zdiff))[0][0]
    inds_avg=np.arange(0,ind_top,1)
    
    meanx = np.nanmean(B_pos[inds_avg]*rho0[inds_avg]*u0[inds_avg])/np.nanmean(B_pos[inds_avg]*rho0[inds_avg])
    meany = np.nanmean(B_pos[inds_avg]*rho0[inds_avg]*v0[inds_avg])/np.nanmean(B_pos[inds_avg]*rho0[inds_avg])
    
    #meanx = np.nanmean(rho0[inds_avg]*u0[inds_avg])/np.nanmean(rho0[inds_avg])
    #meany = np.nanmean(rho0[inds_avg]*v0[inds_avg])/np.nanmean(rho0[inds_avg])
    
    f6000 = np.where(z0<=6000)[0]
    #meanx=np.mean(u0[f6000])
    #meany=np.mean(v0[f6000])
    
    f0500 = np.where(z0<=500)[0]
    lowx=np.mean(u0[f0500])
    lowy=np.mean(v0[f0500])
    
    f560 = np.where(np.logical_and(z0<=6000,z0>=5500))[0]
    highx=np.mean(u0[f560])
    highy=np.mean(v0[f560])
    BK_SHRx=highx-lowx
    BK_SHRy=highy-lowy
    BK_mag=np.sqrt(BK_SHRx**2 + BK_SHRy**2)
    BK_dirx=BK_SHRx/BK_mag
    BK_diry=BK_SHRy/BK_mag
    BK_orthx=BK_diry*7.5
    BK_orthy=-BK_dirx*7.5


    SR_mean_u= u0 - meanx
    SR_mean_v= v0 - meany
    dudz=np.zeros(u0.shape)
    dvdz=np.zeros(v0.shape)
    dudz[1:dudz.shape[0]-1]= ( u0[2:dudz.shape[0]]-u0[0:dudz.shape[0]-2] )/( z0[2:dudz.shape[0]]-z0[0:dudz.shape[0]-2] )
    dudz[0]=2*dudz[1]-dudz[2]
    dvdz[1:dudz.shape[0]-1]= ( v0[2:dudz.shape[0]]-v0[0:dudz.shape[0]-2] )/( z0[2:dudz.shape[0]]-z0[0:dudz.shape[0]-2] )
    dvdz[0]=2*dvdz[1]-dvdz[2]
    f1000 = np.where(z0<=1000)[0]
    SRH_mean = abs(np.mean(-SR_mean_u[f1000]*dvdz[f1000] + SR_mean_v[f1000]*dudz[f1000])*1000.0)
    
    
    propfac=min(SRH_mean/150,1)


    C_x=meanx+propfac*BK_orthx
    C_y=meany+propfac*BK_orthy
    
    u_sr = u0 - C_x
    v_sr = v0 - C_y
    
    f1000 = np.where(z0<=1000)[0]
    V_SR = np.nanmean(np.sqrt(  u_sr[f1000]**2 + v_sr[f1000]**2  ))
    return V_SR



#==============================================================================
#FUNCTION THAT COMPUTES BUNKERS SR MOTION
def compute_OMEGA_AND_SRH(z0,u0,v0,C_x,C_y,rho0,EL):
    #compute 0-1 km storm-relative flow (V_SR) using the storm motion
    #estimate of Bunkers et al. (2000)
    #https://doi.org/10.1175/1520-0434(2000)015<0061:PSMUAN>2.0.CO;2
    
    zdiff = ( z0 - EL )**2
    ind_top = np.where(zdiff==np.min(zdiff))[0][0]
    inds_avg=np.arange(0,ind_top,1)
    
    meanx = np.nanmean(rho0[inds_avg]*u0[inds_avg])/np.nanmean(rho0[inds_avg])
    meany = np.nanmean(rho0[inds_avg]*v0[inds_avg])/np.nanmean(rho0[inds_avg])
    
    f6000 = np.where(z0<=6000)[0]
    #meanx=np.mean(u0[f6000])
    #meany=np.mean(v0[f6000])
    
    f0500 = np.where(z0<=500)[0]
    lowx=np.mean(u0[f0500])
    lowy=np.mean(v0[f0500])
    
    f560 = np.where(np.logical_and(z0<=6000,z0>=5500))[0]
    highx=np.mean(u0[f560])
    highy=np.mean(v0[f560])
    BK_SHRx=highx-lowx
    BK_SHRy=highy-lowy
    BK_mag=np.sqrt(BK_SHRx**2 + BK_SHRy**2)
    BK_dirx=BK_SHRx/BK_mag
    BK_diry=BK_SHRy/BK_mag
    BK_orthx=BK_diry*7.5
    BK_orthy=-BK_dirx*7.5


    SR_mean_u= u0 - meanx
    SR_mean_v= v0 - meany
    dudz=np.zeros(u0.shape)
    dvdz=np.zeros(v0.shape)
    dudz[1:dudz.shape[0]-1]= ( u0[2:dudz.shape[0]]-u0[0:dudz.shape[0]-2] )/( z0[2:dudz.shape[0]]-z0[0:dudz.shape[0]-2] )
    dudz[0]=2*dudz[1]-dudz[2]
    dvdz[1:dudz.shape[0]-1]= ( v0[2:dudz.shape[0]]-v0[0:dudz.shape[0]-2] )/( z0[2:dudz.shape[0]]-z0[0:dudz.shape[0]-2] )
    dvdz[0]=2*dvdz[1]-dvdz[2]
    f1000 = np.where(z0<=1000)[0]
    SRH_mean = abs(np.mean(-SR_mean_u[f1000]*dvdz[f1000] + SR_mean_v[f1000]*dudz[f1000])*1000.0)
    
    
    #propfac=min(SRH_mean/150,2)
    propfac=min(SRH_mean/250,2)
    #propfac=1

    if C_x<0:
        C_x=meanx+propfac*BK_orthx
        C_y=meany+propfac*BK_orthy
        
    u_sr = u0 - C_x
    v_sr = v0 - C_y
    sr_mag = np.sqrt(  u_sr**2 + v_sr**2 )
    dudz=np.zeros(u0.shape)
    dvdz=np.zeros(v0.shape)
    dudz[1:dudz.shape[0]-1]= ( u0[2:dudz.shape[0]]-u0[0:dudz.shape[0]-2] )/( z0[2:dudz.shape[0]]-z0[0:dudz.shape[0]-2] )
    dudz[0]=2*dudz[1]-dudz[2]
    dvdz[1:dudz.shape[0]-1]= ( v0[2:dudz.shape[0]]-v0[0:dudz.shape[0]-2] )/( z0[2:dudz.shape[0]]-z0[0:dudz.shape[0]-2] )
    dvdz[0]=2*dvdz[1]-dvdz[2]
    f1000 = np.where(z0<=1000)[0]
    SRH = abs(np.mean(-u_sr[f1000]*dvdz[f1000] + v_sr[f1000]*dudz[f1000])*1000.0)
    OMEGA = np.mean( (-u_sr[f1000]*dvdz[f1000] + v_sr[f1000]*dudz[f1000])/sr_mag[f1000] )
        
    
    
    f1000 = np.where(z0<=1000)[0]
    V_SR = np.nanmean(np.sqrt(  u_sr[f1000]**2 + v_sr[f1000]**2  ))
    return V_SR,C_x,C_y,SRH,OMEGA




#==============================================================================
#FUNCTION THAT COMPUTES CAPE, CIN, EL, LFC
def compute_CAPE_CONTS(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2):
#[CAPE,CIN,LFC,EL]

    #this function computes CAPE and CIN
    
    #input arguments
    #T0: sounding profile of temperature (in K)
    #p0: sounding profile of pressure (in Pa)
    #q0: sounding profile of water vapor mass fraction (in kg/kg)
    #start_loc: index of the parcel starting location (set to 1 for the
    #lowest: level in the sounding)
    #fracent: fractional entrainment rate (in m^-1)
    
    #CONSTANTS

    
    #compute lifted parcel buoyancy
    T_lif,Qv_lif,Qt_lif,B_lif=lift_parcel_adiabatic(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2)
    
    #compute lifted parcel buoyancy
    T_lif_p,Qv_lif_p,Qt_lif_p,B_lif_p=lift_parcel_adiabatic(T0,p0,q0,start_loc,fracent,0.01,z0,T1,T2)
    
    #compute thermal buoyancy 
    B_therm = g*(T_lif - T0)/T0
    
    B_therm_p = g*(T_lif_p - T0)/T0
    B_cond = B_therm - B_therm_p
    
    #water vapor buoyancy contribution
    B_vap = g*(Rv/Rd - 1)*(Qv_lif - q0)
    
    #condensate loading contribution
    B_load = -g*(Qt_lif - Qv_lif)
    
    
    
    if np.nanmax(B_lif)>0:
        #CAPE will be the total integrated positive buoyancy
        B_pos = np.zeros(B_lif.shape)
        B_pos[:] = B_lif[:]
        
        B_therm[np.where(B_pos<0)]=0
        B_cond[np.where(B_pos<0)]=0
        B_vap[np.where(B_pos<0)]=0
        B_load[np.where(B_pos<0)]=0
        
        B_pos[np.where(B_pos<0)]=0
        dz = z0[1:z0.shape[0]] - z0[0:z0.shape[0]-1]
        CAPE = np.nansum( 0.5*B_pos[0:z0.shape[0]-1]*dz + 0.5*B_pos[1:z0.shape[0]]*dz )
        CAPE_therm = np.nansum( 0.5*B_therm[0:z0.shape[0]-1]*dz + 0.5*B_therm[1:z0.shape[0]]*dz )
        CAPE_cond = np.nansum( 0.5*B_cond[0:z0.shape[0]-1]*dz + 0.5*B_cond[1:z0.shape[0]]*dz )
        CAPE_vap = np.nansum( 0.5*B_vap[0:z0.shape[0]-1]*dz + 0.5*B_vap[1:z0.shape[0]]*dz )
        CAPE_load = np.nansum( 0.5*B_load[0:z0.shape[0]-1]*dz + 0.5*B_load[1:z0.shape[0]]*dz )
        
        #CIN will be the total negative buoyancy below the height of maximum
        #buoyancy
        B_neg = np.zeros(B_lif.shape)
        B_neg[:] = B_lif[:]
        mx = np.nanmax(B_lif)
        imx = np.where(B_lif==mx)
        imx=imx[0][0]
        B_neg[0:imx]=np.minimum( B_neg[0:imx], 0 )
        B_neg[imx:z0.shape[0]]= 0
        CIN = np.nansum( 0.5*B_neg[0:z0.shape[0]-1]*dz + 0.5*B_neg[1:z0.shape[0]]*dz )
        
        #LFC will be the last instance of negative buoyancy before the
        #continuous interval that contains the maximum in buoyancy
        fneg = np.where(B_lif<0)
        fneg=fneg[0]
        inn = np.where(fneg<imx)
        inn = inn[0]
        fneg = fneg[inn]
        if len(fneg)>0:
            LFC = 0.5*z0[np.max(fneg)] + 0.5*z0[np.max(fneg)+1]
        else:
            LFC = z0[start_loc]
        
        #EL will be last instance of positive buoyancy
        fpos = np.where(B_lif>0)
        fpos=fpos[0]
        EL = 0.5*z0[np.max(fpos)] + 0.5*z0[np.max(fpos)+1]
    else:
        CAPE = 0
        CIN = 0
        LFC = np.nan
        EL = np.nan
        CAPE_therm = 0
        CAPE_cond = 0
        CAPE_vap = 0
        CAPE_load = 0

    return CAPE,CAPE_therm,CAPE_cond,CAPE_vap,CAPE_load



#==============================================================================
#FUNCTION THAT COMPUTES CAPE, CIN, EL, LFC
def compute_CAPES_DRAG(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2,Radius,V_SR):
#[CAPE,CIN,LFC,EL]

    #this function computes CAPE and CIN
    
    #input arguments
    #T0: sounding profile of temperature (in K)
    #p0: sounding profile of pressure (in Pa)
    #q0: sounding profile of water vapor mass fraction (in kg/kg)
    #start_loc: index of the parcel starting location (set to 1 for the
    #lowest: level in the sounding)
    #fracent: fractional entrainment rate (in m^-1)
    
    #CONSTANTS

    alpha = 0.8
    #c_d = 0.2
    c_d = 0.2
    
    th0 = T0*(1000*100/p0)**(Rd/cp)
    
    #compute lifted parcel buoyancy
    T_lif,Qv_lif,Qt_lif,B_lif=lift_parcel_adiabatic(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2)
    
    #CALCULATE THE LIFTED CONDENSATION LEVEL
    qdiff = abs(Qt_lif - Qv_lif) #FIGURE OUT THE FIRST HEIGHT WHERE QV STARTS DEVIATING FROM QT, IMPLYING CONDENSATION
    if np.logical_and(~np.isnan(qdiff[1]),np.nanmax(qdiff)>0):
        lcl_ind = np.where(qdiff>0)[0][0]
        LCL = z0[lcl_ind]
    else:
        LCL = 1000
        lcl_ind = np.where(abs(LCL-z0)==np.amin(abs(LCL-z0)))[0][0]
    
    #IF WE HAVE SOME POSITIVE BUOYANCY, PROCEED
    if np.nanmax(B_lif)>0:
        #MAKE A NEW MATRIX THAT WILL ONLY CONTAIN THE POSITIVE PART OF BUOYANCY
        B_pos = np.zeros(B_lif.shape)
        B_pos[:] = B_lif[:]

        #GET RID OF ALL NEGATIVE BUOYANCY BELOW THE LCL
        B_pos[0:lcl_ind]=0
        wpos = np.where(B_pos>0)[0]
        if len(wpos)>0:
            wpos=wpos[0] #WPOS CONTAINS INDEX OF LCL.  SET TO 0 IF THERE IS NO POSTIVE BUOYANCY
        else:
            wpos=lcl_ind
        B_pos[0:wpos]=0
        dz = z0[1:z0.shape[0]] - z0[0:z0.shape[0]-1]      
  
        #LFC WILL BE THE LAST INSTANCE OF NEGATIVE BUOYANCY BEFORE THE PARCEL REACHES ITS CONTINUOUS INTERVAL OF POSITIVE BUOY
        mx = np.nanmax(B_lif)
        imx = np.where(B_lif==mx)
        imx=imx[0][0]
        
        fneg = np.where(B_lif<0)
        fneg=fneg[0]
        inn = np.where(fneg<imx)
        
        inn = inn[0]
        fneg = fneg[inn]
        if len(inn)>0:
            LFC = 0.5*z0[np.max(fneg)] + 0.5*z0[np.max(fneg)+1]
        else:
            LFC = z0[start_loc]
        
        #EL WILL BE THE LAST INSTANCE OF POSITIVE BUOYANCY
        fpos = np.where(B_lif>0)
        fpos=fpos[0]
        EL = 0.5*z0[np.max(fpos)] + 0.5*z0[np.max(fpos)+1]
        
        #INTIALIZE PROFILE OF SQUARED VERTICAL VELOCITY (I.E., VERTICAL KINETIC ENERGY)
        WSQ_prof = np.zeros(B_pos.shape[0])
        WSQ_prof[start_loc]=(V_SR**2)/2 #LOWER BOUNADRY CONDITION ON VERTICAL KE IS THE KE OF INFLOW
        for iz in np.arange(0,WSQ_prof.shape[0]-1,1): #VERTICALLY INTEGRATE
            B_on = B_pos[iz] #STORE THE CURRENT BUOYANCY
            ebuoy_fac = 1/(1 + 2*(alpha**2)*(Radius**2)/((EL-LFC)**2 ) ) #SCALE FACTOR THAT ACCOUNTS FOR EFFECITVE BUOYANCY
            #ns_drag = -2.5*c_d*(3/8)/Radius #COEFICIENT ON THE NON-SHEARED PART OF DRAG
            
            N = np.max( ( g/th0[iz] )*( th0[iz+1]-th0[iz] )/( z0[iz+1]-z0[iz] ), 0 )
            N = ( np.minimum(z0[iz]/Radius,1)*1/2 + np.minimum( np.maximum(( z0[iz] - Radius)/Radius,0), 1)*(2/3-1/2))*N
            
            F = np.sqrt( np.max( WSQ_prof/2 ,0) )/( np.sqrt(N)*Radius )
            
            if N>0:
                ns_drag = -(c_d + comp_cdwave(F))*(3/8)/Radius #COEFICIENT ON THE NON-SHEARED PART OF DRAG
            else:
                ns_drag = -(c_d)*(3/8)/Radius #COEFICIENT ON THE NON-SHEARED PART OF DRAG
            
            if np.sqrt(2*WSQ_prof[iz-1])<1: #IF WE HAVE VERY SMALL VERTICAL VELOICTY (LESS THAN 1 M/S, WE NEED TO ZERO OUT THE SHEAR DRAG TERM OR THINGS BLOW UP)
                sh_drag = 0
                
            #NOW VERTICALLY INTEGRATE THE UPRIME AND WSQ EQUATIONS TOGETHER, FOLLOWING EQ. XX AND XX IN XX RESPECTIVELY
            WSQ_prof[iz+1] = WSQ_prof[iz] + ( z0[iz+1]-z0[iz] )*(ebuoy_fac*B_on + ns_drag*WSQ_prof[iz])
          
        #WE WILL OUTPUT THE MAXIMUM KE AS THE "CAPE" ARGUMENT                                                                                                          
        CAPE = np.nanmax(WSQ_prof)
        
        #SET THE EL TO THE HEIGHT OF MAXIMUM VERTICAL VELOCITY
        mxval = np.nanmax(WSQ_prof)
        fnval=np.where(WSQ_prof==mxval)
        LFC = LCL
        if fnval[0].shape[0]>0:
            EL = z0[fnval[0][0]]
        else:
            EL = np.nan
    else:
        #IF WE HAVE NO POSITIVE BUOYANCY, SET EVERYTHING TO 0S AND NANS
        CAPE = 0      
        LFC = np.nan
        EL = np.nan
        B_pos = np.zeros(T0.shape)

    return CAPE



#==============================================================================
#FUNCTION THAT COMPUTES BUNKERS SR MOTION
def compute_VSR_squall(z0,u0,v0,T0,q0):
    #compute 0-1 km storm-relative flow (V_SR) using the storm motion
    #estimate of Bunkers et al. (2000)
    #https://doi.org/10.1175/1520-0434(2000)015<0061:PSMUAN>2.0.CO;2
    

    
    dz = np.zeros(z0.shape)
    dz[0:dz.shape[0]-1] = z0[1:dz.shape[0]] - z0[0:dz.shape[0]-1]
    
    #COMPUTE THE MOIST STATIC ENERGY
    MSE0 = cp*T0 + xlv*q0 + g*z0
    
    f3000 = np.where(z0<3000)[0]
    minind = np.argmin(MSE0[f3000])
    if minind!=0:
        
        intarg = ( g/(cp*T0[0:minind]) )*(MSE0[minind]-MSE0[0:minind])
        csq = z0[minind]*np.sum(intarg[0:minind]*dz[0:minind])/np.sum(dz[0:minind])
        c = np.sqrt(-2*csq)
    
        
        du = np.sqrt( (u0[minind]-u0[0])**2 + (v0[minind]-v0[0])**2)

        V_SR = (du**2)/c
    else:
        V_SR = 0
    return V_SR,c,du

#==============================================================================
#FUNCTION THAT COMPUTES THE WATER VAPOR PRESSURE FROM THE MIXING RATIO
def vaporpressure_from_mixrat(qv,p):
    return qv*p/epsilon


#==============================================================================
#FUNCTION THAT COMPUTES THE MOIST ENTROPY FROM STATE VARIABLES
def moist_entropy(T,p,qv,ql,qi):
    #NOTE THE FORMULA IN HERE FOLLOWS ROMPS AND KUANG (2010) https://journals.ametsoc.org/view/journals/atsc/67/2/2009jas3184.1.xml
    #THIS FORMULA COMPUTES MOIST ENTROPY RELATIVE TO THE ENTROPY OF WATER VAPOR AT THE TRIPLE POINT TEMPERATURE 
    #INPUTS FOR WATER SPECIES ARE MASS FRACTIONS (denoted by q), BECAUSE THIS IS HOW THEY ARE TYPICALLY STORED IN MY CODE
    e = vaporpressure_from_mixrat(qv,p) #WATER VAPOR PRESSURE
    pd = p - e #DRY AIR PRESSURE
    qt = qv+ql+qi #TOTAL WATER MASS FRACTION
    
    s_d = cp*np.log(T/ttrip) - Rd*np.log(pd/eref) #DRY AIR ENTROPY
    s_v = cpv*np.log(T/ttrip) - Rv*np.log(e/eref) + xlv/ttrip #WATER VAPOR ENTROPY       
    s_l = cpl*np.log(T/ttrip) #LIQUID WATER ENTROPY
    s_i = cpi*np.log(T/ttrip) - (xls-xlv)/ttrip #SOLID WATER ENTROPY
    
    #COMPUTE MIXING RATIOS FROM MASS FRACTIONS
    r_v = qv/(1 - qt)
    r_l = ql/(1 - qt)
    r_i = qi/(1 - qt) 

    #THE MOIST ENTROPY IS A MASS WEIGHTED SUM OF THE IDIVIDUAL COMPONENTS
    s_m = s_d + r_v*s_v + r_l*s_l + r_i*s_i
    #COMPUTE THETA_E FROM ENTROPY FOLLOWING ROMPS AND KUANG (2010)
    theta_e = np.exp(s_m/cp)*ttrip*(1000*100/eref)**(Rd/cp)
    
    return s_m,theta_e

#==============================================================================
#WRAPPER FUNCTION FOR USE IN NONLINEAR SOLVER, FOR SATURATED PARCEL
def me_inner_wrapper(s_m,T,p,qt,T1,T2):
    rv = compute_rsat(T,p,1,T1,T2)
    qvs = (1 - qt)*rv
    
    #if qvs<qv:
    ql = (1 - omega(T,T1,T2))*(qt - qvs)
    qi = omega(T,T1,T2)*(qt - qvs)
    smguess,null = moist_entropy(T,p,qvs,ql,qi)
    #else:
    #    smguess,null = moist_entropy(T,p,qv,0,0)
    return s_m - smguess

#==============================================================================
#WRAPPER FUNCTION FOR USE IN NONLINEAR SOLVER, UNSATURATED PARCEL
def me_inner_wrapper_unsat(s_m,T,p,qv):
    
    smguess,null = moist_entropy(T,p,qv,0,0)
    return s_m - smguess

 

#==============================================================================
#FUNCTION THAT COMPUTES THE MOIST ENTROPY FROM STATE VARIABLES
def get_hydrostatic_pressure(T0,q0,z0,psfc):
  
    p0 = np.zeros(T0.shape)
    p0[0] = psfc
    theta_rho = T0*(1 + (Rv/Rd - 1)*q0)
    intarg = g/(Rd*theta_rho)
    dz = z0[1: ] - z0[0: -1]
    
    
    for iz in np.arange(1,p0.shape[0],1):
        intsum = - 0.5*np.sum( ( intarg[0:iz] + intarg[1:iz+1] )*dz[0:iz] )
        p0[iz] = psfc*np.exp(intsum)
    
    return p0




#==============================================================================
#FUNCTION THAT LIFTS A PARCEL
def lift_parcel_adiabatic_leapfrog(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2):
    #[T_lif,Qv_lif,Qt_lif,B_lif]

    #this function computes lifted parcel properties using the unsaturated
    #and saturated lapse rate formulas from (Peters et al. 2022)
    #https://doi-org.ezaccess.libraries.psu.edu/10.1175/JAS-D-21-0118.1 
    
    #input arguments
    #T0: sounding profile of temperature (in K)
    #p0: sounding profile of pressure (in Pa)
    #q0: sounding profile of water vapor mass fraction (in kg/kg)
    #start_loc: index of the parcel starting location (set to 1 for the
    #lowest: level in the sounding)
    #fracent: fractional entrainment rate (in m^-1)
    
    #output arguments
    #T_lif: lifted parcel temperature
    #Qv_lif: lifted parcel water vapor mass fraction
    #Qt_lif: lifted parcel total water mass fraction
    #B_lif: Lifted parcel buoyancy, computed using Eq. B6 in (Peters et al.
    #2022) (accounts for virtual temperature and loading effects)
    
    #prate: precipitation rate (in m^-1) large values make parcel more
    #pseudoadiabatic, small values make parcel more adiabatic.  I usually
    #just set it to 0 to get an adiabatic parce
    
    #z0: sounding profile of height above ground level (first height should
    #be 0 m)
    #T1 warmest mixed-phase temperature
    #T2 coldest mixed-phase temperature

    #CONSTANTS

    
    #ESTIMATE THE MOIST STATIC ENERGY (MSE)
    MSE = cp*T0 + xlv*q0 + g*z0
    mn_hgt = np.min(np.where(MSE==np.nanmin(MSE))) #FIND THE INDEX OF THE HEIGHT OF MINIMUM MSE
    
    #descriminator function between liquid and ice (i.e., omega defined in the
    #beginning of section 2e in Peters et al. 2022)

    
    T_lif=np.zeros(T0.shape)*np.nan #temperature of the lifted parcel
    Qv_lif=np.zeros(T0.shape)*np.nan #water vapor mass fraction of the lifted parcel
    Qt_lif=np.ones(T0.shape)*np.nan #total water mass fraction of the lifted parcel

    if start_loc>0:
        T_lif[0:start_loc+1]=T0[0:start_loc+1] #set initial values to that of the environment
        Qv_lif[0:start_loc+1]=q0[0:start_loc+1] #set initial values to that of the environment
        Qt_lif[0:start_loc+1]=Qv_lif[0:start_loc+1] #set initial values to that of the environment
    else:
        T_lif[0]=T0[0] #set initial values to that of the environment
        Qv_lif[0]=q0[0] #set initial values to that of the environment
        Qt_lif[0]=Qv_lif[0] #set initial values to that of the environment


    q_sat_prev=0
    B_run = 0
    iz=start_loc
    #
    #for iz in np.arange(start_loc+1,z0.shape[0]):
    #
    #
    #I REVISED THIS A BIT.  TO MAKE THE CODE FASTER, I HAVE THE CALCULATION CUT OUT WHEN THE INTEGRATED NEGATIVE BUOYANCY ("BRUN") 
    #BECOMES MORE NEGATIVE THAN THAN THE TOTAL INTEGRATED POSITIVE BUOYANCY.  I RESTRICT THIS TO ONLY HAPPEN AFTER WE HAVE PASSED 
    #THE HEIGHT OF MINIMUM MSE.  UNCOMMENT THE FOR LOOP ABOVE AND COMMENT OUT THE WHILE LOOP IF YOU JUST WANT TO INTEGRATE TO THE TOP OF THE SOUNDING.
    #THE +25 PART IN THE WHILE STATEMENT IS A PAD ON B_RUN (THE NEGATIVE CAPE HAS TO BE 25 J/KG LESS THAN THE POSITIVE CAPE TO KILL THE LOOP)
    satflag = True
    #while iz<(z0.shape[0])-1 and (z0[iz]<z0[mn_hgt] or (B_run+25)>0):
    while iz<(z0.shape[0])-1 and (z0[iz]<z0[mn_hgt] or (B_run+250)>0):
        iz = iz + 1
        q_sat=(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],1,T1,T2)
        if Qv_lif[iz-1]<q_sat: #if we are unsaturated, go up at the unsaturated adiabatic lapse rate (eq. 19 in Peters et al. 2022)
            
        
            if iz==start_loc+1:
                T_lif[iz] = T_lif[iz-1] + (z0[iz] - z0[iz-1])*drylift(T_lif[iz-1],Qv_lif[iz-1],T0[iz-1],q0[iz-1],fracent)
                Qv_lif[iz] = Qv_lif[iz-1] - (z0[iz] - z0[iz-1])*fracent*( Qv_lif[iz-1] - q0[iz-1] )
            else:
                T_lif[iz] = T_lif[iz-2] + (z0[iz] - z0[iz-2])*drylift(T_lif[iz-1],Qv_lif[iz-1],T0[iz-1],q0[iz-1],fracent)
                Qv_lif[iz] = Qv_lif[iz-2] - (z0[iz] - z0[iz-2])*fracent*( Qv_lif[iz-1] - q0[iz-1] )
                
            Qt_lif[iz] = Qv_lif[iz]
            q_sat=(1-Qt_lif[iz])*compute_rsat(T_lif[iz],p0[iz],1,T1,T2)
            
            if Qv_lif[iz]>=q_sat: #if we hit saturation, split the vertical step into two stages.  The first stage advances at the saturated lapse rate to the saturation point, and the second stage completes the grid step at the moist lapse rate
                OMEGA = omega(T_lif[iz-1],T1,T2)
                dOMEGA = domega(T_lif[iz-1],T1,T2)
                satrat=(Qv_lif[iz]-q_sat_prev)/(q_sat-q_sat_prev)
                dz_dry=satrat*(z0[iz]-z0[iz-1])
                dz_wet=(1-satrat)*(z0[iz]-z0[iz-1])


                
                T_halfstep = T_lif[iz-1] + dz_dry*drylift(T_lif[iz-1],Qv_lif[iz-1],T0[iz-1],q0[iz-1],fracent)
                Qv_halfstep = Qv_lif[iz-1] - dz_dry*fracent*( Qv_lif[iz-1] - q0[iz-1] )
                Qt_halfstep = Qv_lif[iz]
                p_halfstep=p0[iz-1]*satrat + p0[iz]*(1-satrat)
                T0_halfstep=T0[iz-1]*satrat + T0[iz]*(1-satrat)
                Q0_halfstep=q0[iz-1]*satrat + q0[iz]*(1-satrat)

                T_lif[iz] = T_halfstep + dz_wet*moislif(T_halfstep,Qv_halfstep,(1-Qt_halfstep)*compute_rsat(T_halfstep,p_halfstep,0,T1,T2),(1-Qt_halfstep)*compute_rsat(T_halfstep,p_halfstep,2,T1,T2),p_halfstep,T0_halfstep,Q0_halfstep,Qt_halfstep,fracent,prate,T1,T2)
                
                
                Qt_lif[iz] = Qt_lif[iz-1] - (z0[iz] - z0[iz-1])*fracent*( Qt_halfstep - Q0_halfstep )
                Qv_lif[iz] = (1-Qt_lif[iz])*compute_rsat(T_lif[iz],p0[iz],1,T1,T2)

                if Qt_lif[iz]<Qv_lif[iz]:
                    Qv_lif[iz]=Qt_lif[iz]

            q_sat_prev=q_sat;
            
        else: #if we are already at saturation, just advance upward using the saturated lapse rate (eq. 24 in Peters et al. 2022)
            OMEGA = omega(T_lif[iz-1],T1,T2)
            dOMEGA = domega(T_lif[iz-1],T1,T2)

            if iz==start_loc+1 or satflag:
                T_lif[iz] = T_lif[iz-1] + (z0[iz] - z0[iz-1])*moislif(T_lif[iz-1],Qv_lif[iz-1],(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],0,T1,T2),(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],2,T1,T2),p0[iz-1],T0[iz-1],q0[iz-1],Qt_lif[iz-1],fracent,prate,T1,T2);
                Qt_lif[iz] = Qt_lif[iz-1] - (z0[iz] - z0[iz-1])*(fracent*( Qt_lif[iz-1] - q0[iz-1] )  + prate*( Qt_lif[iz-1]-Qv_lif[iz-1]) )
                satflag = False
            else:
                T_lif[iz] = T_lif[iz-2] + (z0[iz] - z0[iz-2])*moislif(T_lif[iz-1],Qv_lif[iz-1],(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],0,T1,T2),(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],2,T1,T2),p0[iz-1],T0[iz-1],q0[iz-1],Qt_lif[iz-1],fracent,prate,T1,T2);
                Qt_lif[iz] = Qt_lif[iz-2] - (z0[iz] - z0[iz-2])*(fracent*( Qt_lif[iz-1] - q0[iz-1] )  + prate*( Qt_lif[iz-1]-Qv_lif[iz-1]) )
            Qv_lif[iz] = (1-Qt_lif[iz])*compute_rsat(T_lif[iz],p0[iz],1,T1,T2)
            
            if Qt_lif[iz]<Qv_lif[iz]:
                Qv_lif[iz]=Qt_lif[iz]

        B_run = B_run + (g*T_lif[iz]*(1 + (Rv/Rd)*Qv_lif[iz] - Qt_lif[iz])/(T0[iz]*(1 + (Rv/Rd)*q0[iz] - q0[iz])) - g)*(z0[iz]-z0[iz-1])

    T_rho_lif = T_lif*(1 + (Rv/Rd)*Qv_lif - Qt_lif)
    T_0_lif = T0*(1 + (Rv/Rd - 1)*q0)
    #T_rho_lif=T_lif*(1 - Qt_lif + Qv_lif)/( 1 + (epsilon - 1)/( ( epsilon*(1 - Qt_lif)/Qv_lif - 1) ) )
    #T_0_lif=T0/( 1 + (epsilon - 1)/( ( epsilon*(1 - q0)/q0 - 1) ) )
    
    B_lif=g*(T_rho_lif - T_0_lif)/T_0_lif
    
    
    return T_lif,Qv_lif,Qt_lif,B_lif

#==============================================================================

#==============================================================================
#LIFT A PARCEL BY CONSERVING THE MOIST ENTROPY
def lift_parcel_ME_old(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2):
    
    #INITIALIZE ARRAYS
    ME_0,null = moist_entropy(T0[start_loc],p0[start_loc],q0[start_loc],0,0)
    ME = np.zeros(z0.shape)
    qt = np.zeros(z0.shape)
    qv = np.copy(qt)
    ql = np.copy(qt)
    qi = np.copy(qt)
    T_lif = np.zeros(z0.shape)
    Qv_lif = np.zeros(z0.shape)
    Rv_lif = np.zeros(z0.shape)
    
    ME0,null = moist_entropy(T0,p0,q0,0,0)
    dluga = np.where((ME0-100)<ME_0)
    if len(dluga[0])>0:
        bluga = np.minimum(np.max(dluga[0])+10,ME0.shape[0])
        wherepos = np.arange(0,bluga,1)
        whereneg = np.arange(bluga,z0.shape,1)
        
        #IF WE AREN'T ENTRAINING, JUST HOLD ENTROPY AND TOTAL WATER CONSTANT ALONG PARCEL PATH
        if fracent==0:
            ME[:] = ME_0
            qt[:] = q0[start_loc]
        else:        
            intarg_ME = np.exp(fracent*(z0-z0[start_loc]))*fracent*ME0
            intarg_qt = np.exp(fracent*(z0-z0[start_loc]))*fracent*q0
            meint = np.copy(qt)
            q0int = np.copy(qt)
            for iz in np.arange(start_loc+1,z0.shape[0],1):
                meint[iz] = 0.5*np.nansum( ( intarg_ME[start_loc:iz] + intarg_ME[start_loc+1:iz+1] )*(z0[start_loc+1:iz+1]-z0[start_loc:iz]) )
                q0int[iz] = 0.5*np.nansum( ( intarg_qt[start_loc:iz] + intarg_qt[start_loc+1:iz+1] )*(z0[start_loc+1:iz+1]-z0[start_loc:iz]) )
            ME = np.exp(-fracent*(z0-z0[start_loc]))*(ME_0 + meint)
            qt = np.exp(-fracent*(z0-z0[start_loc]))*(q0[start_loc] + q0int)
        
        #def me_outer_wrapper(T):
        #    return me_inner_wrapper(ME[wherepos],T[wherepos],p0[wherepos],qt[wherepos],T1,T2)
        
        #T_lif = np.zeros(z0.shape)
        #tmp = fsolve(me_outer_wrapper,T0[wherepos],col_deriv=True)
        #======================================================================
        #======================================================================
        #======================================================================

        T_lif[0] = T0[0]
        Qv_lif[0] = q0[0]
        unsat = True
        satinds = []
        for iz in np.arange(1,len(wherepos),1):
        #for iz in np.arange(1,5,1):
            if unsat:
                def me_outer_wrapper(T):
                    return me_inner_wrapper_unsat(ME[wherepos[iz]],T,p0[wherepos[iz]],qt[wherepos[iz]])
                T_lif[wherepos[iz]] = fsolve(me_outer_wrapper,T_lif[wherepos[iz-1]],col_deriv=True)
                Qv_lif[wherepos[iz]] = qt[wherepos[iz]]

                rstar = compute_rsat(T_lif[wherepos[iz]],p0[wherepos[iz]],1,T1,T2)
                qstar = (1 - qt[wherepos[iz]])*rstar
                if Qv_lif[wherepos[iz]]>qstar:
                    unsat = False
            if not unsat:
                satinds.append(wherepos[iz])
                def me_outer_wrapper(T):
                    return me_inner_wrapper(ME[wherepos[iz]],T,p0[wherepos[iz]],qt[wherepos[iz]],T1,T2)
                T_lif[wherepos[iz]] = fsolve(me_outer_wrapper,T_lif[wherepos[iz-1]],col_deriv=True)
        #======================================================================
        #======================================================================
        #======================================================================
        
        T_lif[whereneg] = np.nan
        Rv_lif[satinds] = compute_rsat(T_lif[satinds],p0[satinds],1,T1,T2)
        Qv_lif[satinds] = (1 - qt[satinds])*Rv_lif[satinds]
        Qt_lif = qt
        
        
        T_rho_lif = T_lif*(1 + (Rv/Rd)*Qv_lif - Qt_lif)
        T_0_lif = T0*(1 + (Rv/Rd - 1)*q0)
    
        B_lif=g*(T_rho_lif - T_0_lif)/T_0_lif
        
        B_lif[start_loc] = 0
        B_lif[whereneg] = np.nan
    else:
        T_lif = np.zeros(z0.shape)
        Qv_lif = np.zeros(z0.shape)
        Qt_lif = np.zeros(z0.shape)
        B_lif = np.zeros(z0.shape)
    
    return T_lif,Qv_lif,Qt_lif,B_lif




#==============================================================================
#LIFT A PARCEL BY CONSERVING THE MOIST ENTROPY
def lift_parcel_ME(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2):
    
    #INITIALIZE ARRAYS
    ME_0,null = moist_entropy(T0[start_loc],p0[start_loc],q0[start_loc],0,0)
    ME = np.zeros(z0.shape)
    qt = np.zeros(z0.shape)
    qv = np.copy(qt)
    ql = np.copy(qt)
    qi = np.copy(qt)
    
    ME0,null = moist_entropy(T0,p0,q0,0,0)
    dluga = np.where((ME0-100)<ME_0)
    if len(dluga[0])>0:
        bluga = np.minimum(np.max(dluga[0])+10,ME0.shape[0])
        wherepos = np.arange(0,bluga,1)
        whereneg = np.arange(bluga,z0.shape,1)
        
        #IF WE AREN'T ENTRAINING, JUST HOLD ENTROPY AND TOTAL WATER CONSTANT ALONG PARCEL PATH
        if fracent==0:
            ME[:] = ME_0
            qt[:] = q0[start_loc]
        else:        
            intarg_ME = np.exp(fracent*(z0-z0[start_loc]))*fracent*ME0
            intarg_qt = np.exp(fracent*(z0-z0[start_loc]))*fracent*q0
            meint = np.copy(qt)
            q0int = np.copy(qt)
            for iz in np.arange(start_loc+1,z0.shape[0],1):
                meint[iz] = 0.5*np.nansum( ( intarg_ME[start_loc:iz] + intarg_ME[start_loc+1:iz+1] )*(z0[start_loc+1:iz+1]-z0[start_loc:iz]) )
                q0int[iz] = 0.5*np.nansum( ( intarg_qt[start_loc:iz] + intarg_qt[start_loc+1:iz+1] )*(z0[start_loc+1:iz+1]-z0[start_loc:iz]) )
            ME = np.exp(-fracent*(z0-z0[start_loc]))*(ME_0 + meint)
            qt = np.exp(-fracent*(z0-z0[start_loc]))*(q0[start_loc] + q0int)
        
        def me_outer_wrapper(T):
            return me_inner_wrapper(ME[wherepos],T[wherepos],p0[wherepos],qt[wherepos],T1,T2)
        
        T_lif = np.zeros(z0.shape)
        #tmp = fsolve(me_outer_wrapper,T0[wherepos],col_deriv=True,xtol=1.49012e-03)
        tmp = fsolve(me_outer_wrapper,T0[wherepos],col_deriv=True)
        # #======================================================================
        # #======================================================================
        # #======================================================================
        # tmp = np.zeros(len(wherepos))
        # tmp[0] = fsolve(me_outer_wrapper,T0[wherepos[0]],col_deriv=True)
        # for iz in np.arange(1,len(wherepos),1):
        #     def me_outer_wrapper(T):
        #         return me_inner_wrapper(ME[wherepos[iz]],T,p0[wherepos[iz]],qt[wherepos[iz],T1,T2)
        #     tmp[iz] = fsolve(me_outer_wrapper,tmp[iz-1],col_deriv=True)
        # #======================================================================
        # #======================================================================
        # #======================================================================
        
        T_lif[wherepos] = tmp
        T_lif[whereneg] = np.nan
        Rv_lif = compute_rsat(T_lif,p0,1,T1,T2)
        Qv_lif = (1 - qt)*Rv_lif
        Qt_lif = qt
        
        usatinds = np.where(Qv_lif>q0[start_loc])[0]
        if len(usatinds)>0:
            def me_outer_wrapper(T):
                return me_inner_wrapper_unsat(ME[usatinds],T,p0[usatinds],qt[usatinds])
                #return me_inner_wrapper_unsat(ME[usatinds],T[usatinds],p0[usatinds],qt[usatinds])
                #for some reason the code was breaking so i had to make this fix.  Not sure what it means
            
            T_lif[usatinds] = fsolve(me_outer_wrapper,T0[usatinds],col_deriv=True)
            #T_lif[usatinds] = fsolve(me_outer_wrapper,T0[usatinds],col_deriv=True)
            Qv_lif[usatinds] = qt[usatinds]
            Qt_lif[usatinds] = qt[usatinds]
        
        T_rho_lif = T_lif*(1 + (Rv/Rd)*Qv_lif - Qt_lif)
        T_0_lif = T0*(1 + (Rv/Rd - 1)*q0)
    
        B_lif=g*(T_rho_lif - T_0_lif)/T_0_lif
        
        B_lif[start_loc] = 0
        B_lif[whereneg] = np.nan
    else:
        T_lif = np.zeros(z0.shape)
        Qv_lif = np.zeros(z0.shape)
        Qt_lif = np.zeros(z0.shape)
        B_lif = np.zeros(z0.shape)
    
    return T_lif,Qv_lif,Qt_lif,B_lif




#==============================================================================
#FUNCTION THAT COMPUTES CAPE, CIN, EL, LFC
def compute_CAPE(B_lif,z0,start_loc=0):
#[CAPE,CIN,LFC,EL]

    #this function computes CAPE and CIN
    
    #input arguments
    #T0: sounding profile of temperature (in K)
    #p0: sounding profile of pressure (in Pa)
    #q0: sounding profile of water vapor mass fraction (in kg/kg)
    #start_loc: index of the parcel starting location (set to 1 for the
    #lowest: level in the sounding)
    #fracent: fractional entrainment rate (in m^-1)
    
    #CONSTANTS

    
    if np.nanmax(B_lif)>0:
        #CAPE will be the total integrated positive buoyancy
        B_pos = np.zeros(B_lif.shape)
        B_pos[:] = B_lif[:]
        B_pos[np.where(B_pos<0)]=0
        dz = z0[1:z0.shape[0]] - z0[0:z0.shape[0]-1]
        CAPE = np.nansum( 0.5*B_pos[0:z0.shape[0]-1]*dz + 0.5*B_pos[1:z0.shape[0]]*dz )
        
        #CIN will be the total negative buoyancy below the height of maximum
        #buoyancy
        B_neg = np.zeros(B_lif.shape)
        B_neg[:] = B_lif[:]
        mx = np.nanmax(B_lif)
        imx = np.where(B_lif==mx)
        imx=imx[0][0]
        B_neg[0:imx]=np.minimum( B_neg[0:imx], 0 )
        B_neg[imx:z0.shape[0]]= 0
        CIN = np.nansum( 0.5*B_neg[0:z0.shape[0]-1]*dz + 0.5*B_neg[1:z0.shape[0]]*dz )
        
        #LFC will be the last instance of negative buoyancy before the
        #continuous interval that contains the maximum in buoyancy
        fneg = np.where(B_lif<0)
        fneg=fneg[0]
        inn = np.where(fneg<imx)
        inn = inn[0]
        fneg = fneg[inn]
        if len(fneg)>0:
            LFC = 0.5*z0[np.max(fneg)] + 0.5*z0[min(np.max(fneg)+1,z0.shape[0]-1)]
        else:
            LFC = z0[start_loc]
        
        #EL will be last instance of positive buoyancy
        fpos = np.where(B_lif>0)
        fpos=fpos[0]
        EL = 0.5*z0[np.max(fpos)] + 0.5*z0[min(np.max(fpos)+1,z0.shape[0])]
    else:
        CAPE = 0
        CIN = 0
        LFC = np.nan
        EL = np.nan

    return CAPE,CIN,LFC,EL










#==============================================================================
#FUNCTION THAT COMPUTES ENVIRONMENTAL PROFILE USING ZERO BUOYANCY MODEL
def zero_buoyancy_model_leapfrog(T_sfc,p_sfc,RH_sfc,fracent,prate,z0,T1,T2,RH_ft,z_trop):
    start_loc = 0
    #[T_lif,Qv_lif,Qt_lif,B_lif]

    #this function computes lifted parcel properties using the unsaturated
    #and saturated lapse rate formulas from (Peters et al. 2022)
    #https://doi-org.ezaccess.libraries.psu.edu/10.1175/JAS-D-21-0118.1 
    
    #input arguments
    #T0: sounding profile of temperature (in K)
    #p0: sounding profile of pressure (in Pa)
    #q0: sounding profile of water vapor mass fraction (in kg/kg)
    #start_loc: index of the parcel starting location (set to 1 for the
    #lowest: level in the sounding)
    #fracent: fractional entrainment rate (in m^-1)
    
    #output arguments
    #T_lif: lifted parcel temperature
    #Qv_lif: lifted parcel water vapor mass fraction
    #Qt_lif: lifted parcel total water mass fraction
    #B_lif: Lifted parcel buoyancy, computed using Eq. B6 in (Peters et al.
    #2022) (accounts for virtual temperature and loading effects)
    
    #prate: precipitation rate (in m^-1) large values make parcel more
    #pseudoadiabatic, small values make parcel more adiabatic.  I usually
    #just set it to 0 to get an adiabatic parce
    
    #z0: sounding profile of height above ground level (first height should
    #be 0 m)
    #T1 warmest mixed-phase temperature
    #T2 coldest mixed-phase temperature

    #CONSTANTS

    
    
    #descriminator function between liquid and ice (i.e., omega defined in the
    #beginning of section 2e in Peters et al. 2022)

    
    T_lif=np.zeros(z0.shape)*np.nan #temperature of the lifted parcel
    Qv_lif=np.zeros(z0.shape)*np.nan #water vapor mass fraction of the lifted parcel
    Qt_lif=np.ones(z0.shape)*np.nan #total water mass fraction of the lifted parcel
    T0 = np.copy(T_lif)
    q0 = np.copy(Qv_lif)
    p0 = np.copy(Qv_lif)


    T_lif[0]=T_sfc #set initial values to that of the environment
    rsat_sfc = compute_rsat(T_sfc,p_sfc,0,T1,T2)
    qsat_sfc = rsat_sfc/(1 + rsat_sfc)
    Qv_lif[0]=qsat_sfc*RH_sfc #set initial values to that of the environment
    Qt_lif[0]=qsat_sfc*RH_sfc #set initial values to that of the environment
    
    T0[0] = T_sfc
    q0[0] = qsat_sfc*RH_sfc
    p0[0] = p_sfc


    q_sat_prev=0
    B_run = 0
    iz=0
    #
    #for iz in np.arange(start_loc+1,z0.shape[0]):
    #
    #
    #I REVISED THIS A BIT.  TO MAKE THE CODE FASTER, I HAVE THE CALCULATION CUT OUT WHEN THE INTEGRATED NEGATIVE BUOYANCY ("BRUN") 
    #BECOMES MORE NEGATIVE THAN THAN THE TOTAL INTEGRATED POSITIVE BUOYANCY.  I RESTRICT THIS TO ONLY HAPPEN AFTER WE HAVE PASSED 
    #THE HEIGHT OF MINIMUM MSE.  UNCOMMENT THE FOR LOOP ABOVE AND COMMENT OUT THE WHILE LOOP IF YOU JUST WANT TO INTEGRATE TO THE TOP OF THE SOUNDING.
    #THE +25 PART IN THE WHILE STATEMENT IS A PAD ON B_RUN (THE NEGATIVE CAPE HAS TO BE 25 J/KG LESS THAN THE POSITIVE CAPE TO KILL THE LOOP)
    satflag = True
    #while iz<(z0.shape[0])-1 and (z0[iz]<z0[mn_hgt] or (B_run+25)>0):
    while iz<(z0.shape[0])-1:
        iz = iz + 1
        q_sat=(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],1,T1,T2)
        
        if z0[iz]<z_trop:
            if Qv_lif[iz-1]<q_sat: #if we are unsaturated, go up at the unsaturated adiabatic lapse rate (eq. 19 in Peters et al. 2022)
                
            
                if iz==start_loc+1:
                    T_lif[iz] = T_lif[iz-1] + (z0[iz] - z0[iz-1])*drylift(T_lif[iz-1],Qv_lif[iz-1],T0[iz-1],q0[iz-1],fracent)
                    Qv_lif[iz] = Qv_lif[iz-1] - (z0[iz] - z0[iz-1])*fracent*( Qv_lif[iz-1] - q0[iz-1] )
                    T0[iz] = T_lif[iz]
                    q0[iz] = Qv_lif[iz]
                    #INTEGRATE UPWARD TO GET PRESSURE HYDROSTATICALLY
                    p0[iz] = p0[iz-1] + (z0[iz]-z0[iz-1])*( -g*p0[iz-1]/(Rd*T0[iz-1]*(1 + (Rd/Rv-1)*q0[iz-1])  ) )
                else:
                    T_lif[iz] = T_lif[iz-2] + (z0[iz] - z0[iz-2])*drylift(T_lif[iz-1],Qv_lif[iz-1],T0[iz-1],q0[iz-1],fracent)
                    Qv_lif[iz] = Qv_lif[iz-2] - (z0[iz] - z0[iz-2])*fracent*( Qv_lif[iz-1] - q0[iz-1] )
                    T0[iz] = T_lif[iz]
                    q0[iz] = Qv_lif[iz]
                    #INTEGRATE UPWARD TO GET PRESSURE HYDROSTATICALLY
                    p0[iz] = p0[iz-2] + (z0[iz]-z0[iz-2])*( -g*p0[iz-1]/(Rd*T0[iz-1]*(1 + (Rd/Rv-1)*q0[iz-1])  ) )
                    
                Qt_lif[iz] = Qv_lif[iz]
                q_sat=(1-Qt_lif[iz])*compute_rsat(T_lif[iz],p0[iz],1,T1,T2)
                
                if Qv_lif[iz]>=q_sat: #if we hit saturation, split the vertical step into two stages.  The first stage advances at the saturated lapse rate to the saturation point, and the second stage completes the grid step at the moist lapse rate
    
                    OMEGA = omega(T_lif[iz-1],T1,T2)
                    dOMEGA = domega(T_lif[iz-1],T1,T2)
                    satrat=(Qv_lif[iz]-q_sat_prev)/(q_sat-q_sat_prev)
                    dz_dry=satrat*(z0[iz]-z0[iz-1])
                    dz_wet=(1-satrat)*(z0[iz]-z0[iz-1])
    
                    
                    
                    T_halfstep = T_lif[iz-1] + dz_dry*drylift(T_lif[iz-1],Qv_lif[iz-1],T0[iz-1],q0[iz-1],fracent)
                    Qv_halfstep = Qv_lif[iz-1] - dz_dry*fracent*( Qv_lif[iz-1] - q0[iz-1] )
                    Qt_halfstep = Qv_lif[iz]
                    p_halfstep=p0[iz-1]*satrat + p0[iz]*(1-satrat)
                    T0_halfstep=T0[iz-1]*satrat + T0[iz]*(1-satrat)
                    Q0_halfstep=q0[iz-1]*satrat + q0[iz]*(1-satrat)
    
                    T_lif[iz] = T_halfstep + dz_wet*moislif(T_halfstep,Qv_halfstep,(1-Qt_halfstep)*compute_rsat(T_halfstep,p_halfstep,0,T1,T2),(1-Qt_halfstep)*compute_rsat(T_halfstep,p_halfstep,2,T1,T2),p_halfstep,T0_halfstep,Q0_halfstep,Qt_halfstep,fracent,prate,T1,T2)
                    
                    #INTEGRATE UPWARD TO GET PRESSURE HYDROSTATICALLY
                    p0[iz] = p0[iz-1] + (z0[iz]-z0[iz-1])*( -g*p0[iz-1]/(Rd*T0[iz-1]*(1 + (Rd/Rv-1)*q0[iz-1])  ) )
                    
                    Qt_lif[iz] = Qt_lif[iz-1] - (z0[iz] - z0[iz-1])*fracent*( Qt_halfstep - Q0_halfstep )
                    Qv_lif[iz] = (1-Qt_lif[iz])*compute_rsat(T_lif[iz],p0[iz],1,T1,T2)
                    
                    T0[iz] = T_lif[iz]*( 1 + (Rv/Rd)*Qv_lif[iz] - Qt_lif[iz])/(1 + RH_ft*(Rv/Rd-1)*Qv_lif[iz] )
                    q0[iz] = RH_ft*Qv_lif[iz]
                    
    
                    if Qt_lif[iz]<Qv_lif[iz]:
                        Qv_lif[iz]=Qt_lif[iz]
    
                q_sat_prev=q_sat;
                
            else: #if we are already at saturation, just advance upward using the saturated lapse rate (eq. 24 in Peters et al. 2022)
                OMEGA = omega(T_lif[iz-1],T1,T2)
                dOMEGA = domega(T_lif[iz-1],T1,T2)
    
                #if iz==start_loc+1 or satflag:
                #    T_lif[iz] = T_lif[iz-1] + (z0[iz] - z0[iz-1])*moislif(T_lif[iz-1],Qv_lif[iz-1],(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],0,T1,T2),(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],2,T1,T2),p0[iz-1],T0[iz-1],q0[iz-1],Qt_lif[iz-1],fracent,prate,T1,T2);
                #    Qt_lif[iz] = Qt_lif[iz-1] - (z0[iz] - z0[iz-1])*(fracent*( Qt_lif[iz-1] - q0[iz-1] )  + prate*( Qt_lif[iz-1]-Qv_lif[iz-1]) )
                #    satflag = False
                #else:
                T_lif[iz] = T_lif[iz-2] + (z0[iz] - z0[iz-2])*moislif(T_lif[iz-1],Qv_lif[iz-1],(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],0,T1,T2),(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],2,T1,T2),p0[iz-1],T0[iz-1],q0[iz-1],Qt_lif[iz-1],fracent,prate,T1,T2);
                Qt_lif[iz] = Qt_lif[iz-2] - (z0[iz] - z0[iz-2])*(fracent*( Qt_lif[iz-1] - q0[iz-1] )  + prate*( Qt_lif[iz-1]-Qv_lif[iz-1]) )
                
                #INTEGRATE UPWARD TO GET PRESSURE HYDROSTATICALLY
                p0[iz] = p0[iz-2] + (z0[iz]-z0[iz-2])*( -g*p0[iz-1]/(Rd*T0[iz-1]*(1 + (Rd/Rv-1)*q0[iz-1])  ) )        
                Qv_lif[iz] = (1-Qt_lif[iz])*compute_rsat(T_lif[iz],p0[iz],1,T1,T2)
                T0[iz] = T_lif[iz]*( 1 + (Rv/Rd)*Qv_lif[iz] - Qt_lif[iz])/(1 + RH_ft*(Rv/Rd-1)*Qv_lif[iz] )
                q0[iz] = RH_ft*Qv_lif[iz]
                
                if Qt_lif[iz]<Qv_lif[iz]:
                    Qv_lif[iz]=Qt_lif[iz]
        else:
            T0[iz] = T0[iz-1]
            p0[iz] = p0[iz-2] + (z0[iz]-z0[iz-2])*( -g*p0[iz-1]/(Rd*T0[iz-1]*(1 + (Rd/Rv-1)*q0[iz-1])  ) )  
            rsat = compute_rsat(T0[iz],p0[iz],1,T1,T2)
            q0[iz] = RH_ft*(rsat/(1 + rsat))


    T_rho_lif = T_lif*(1 + (Rv/Rd)*Qv_lif - Qt_lif)
    T_0_lif = T0*(1 + (Rv/Rd - 1)*q0)

    
    B_lif=g*(T_rho_lif - T_0_lif)/T_0_lif
    
    
    return T0,q0,p0,B_lif






#==============================================================================
#FUNCTION THAT COMPUTES ENVIRONMENTAL PROFILE USING ZERO BUOYANCY MODEL
def zero_buoyancy_model(T_sfc,p_sfc,RH_sfc,fracent,prate,z0,T1,T2,RH_ft,T_fat):
    start_loc = 0
    #[T_lif,Qv_lif,Qt_lif,B_lif]
    

    #this function computes lifted parcel properties using the unsaturated
    #and saturated lapse rate formulas from (Peters et al. 2022)
    #https://doi-org.ezaccess.libraries.psu.edu/10.1175/JAS-D-21-0118.1 
    
    #input arguments
    #T0: sounding profile of temperature (in K)
    #p0: sounding profile of pressure (in Pa)
    #q0: sounding profile of water vapor mass fraction (in kg/kg)
    #start_loc: index of the parcel starting location (set to 1 for the
    #lowest: level in the sounding)
    #fracent: fractional entrainment rate (in m^-1)
    
    #output arguments
    #T_lif: lifted parcel temperature
    #Qv_lif: lifted parcel water vapor mass fraction
    #Qt_lif: lifted parcel total water mass fraction
    #B_lif: Lifted parcel buoyancy, computed using Eq. B6 in (Peters et al.
    #2022) (accounts for virtual temperature and loading effects)
    
    #prate: precipitation rate (in m^-1) large values make parcel more
    #pseudoadiabatic, small values make parcel more adiabatic.  I usually
    #just set it to 0 to get an adiabatic parce
    
    #z0: sounding profile of height above ground level (first height should
    #be 0 m)
    #T1 warmest mixed-phase temperature
    #T2 coldest mixed-phase temperature

    #CONSTANTS

    
    
    #descriminator function between liquid and ice (i.e., omega defined in the
    #beginning of section 2e in Peters et al. 2022)

    
    T_lif=np.zeros(z0.shape)*np.nan #temperature of the lifted parcel
    Qv_lif=np.zeros(z0.shape)*np.nan #water vapor mass fraction of the lifted parcel
    Qt_lif=np.ones(z0.shape)*np.nan #total water mass fraction of the lifted parcel
    T0 = np.copy(T_lif)
    q0 = np.copy(Qv_lif)
    p0 = np.copy(Qv_lif)


    T_lif[0]=T_sfc #set initial values to that of the environment
    rsat_sfc = compute_rsat(T_sfc,p_sfc,0,T1,T2)
    qsat_sfc = rsat_sfc/(1 + rsat_sfc)
    Qv_lif[0]=qsat_sfc*RH_sfc #set initial values to that of the environment
    Qt_lif[0]=qsat_sfc*RH_sfc #set initial values to that of the environment
    
    T0[0] = T_sfc
    q0[0] = qsat_sfc*RH_sfc
    p0[0] = p_sfc


    q_sat_prev=0
    B_run = 0
    iz=0
    #
    #for iz in np.arange(start_loc+1,z0.shape[0]):
    #
    #
    #I REVISED THIS A BIT.  TO MAKE THE CODE FASTER, I HAVE THE CALCULATION CUT OUT WHEN THE INTEGRATED NEGATIVE BUOYANCY ("BRUN") 
    #BECOMES MORE NEGATIVE THAN THAN THE TOTAL INTEGRATED POSITIVE BUOYANCY.  I RESTRICT THIS TO ONLY HAPPEN AFTER WE HAVE PASSED 
    #THE HEIGHT OF MINIMUM MSE.  UNCOMMENT THE FOR LOOP ABOVE AND COMMENT OUT THE WHILE LOOP IF YOU JUST WANT TO INTEGRATE TO THE TOP OF THE SOUNDING.
    #THE +25 PART IN THE WHILE STATEMENT IS A PAD ON B_RUN (THE NEGATIVE CAPE HAS TO BE 25 J/KG LESS THAN THE POSITIVE CAPE TO KILL THE LOOP)
    satflag = True
    #while iz<(z0.shape[0])-1 and (z0[iz]<z0[mn_hgt] or (B_run+25)>0):
    while iz<(z0.shape[0])-1:
        iz = iz + 1
        q_sat=(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],1,T1,T2)
        
        if T0[iz-1]>T_fat:
            if Qv_lif[iz-1]<q_sat: #if we are unsaturated, go up at the unsaturated adiabatic lapse rate (eq. 19 in Peters et al. 2022)
                
            

                T_lif[iz] = T_lif[iz-1] + (z0[iz] - z0[iz-1])*drylift(T_lif[iz-1],Qv_lif[iz-1],T0[iz-1],q0[iz-1],fracent)
                Qv_lif[iz] = Qv_lif[iz-1] - (z0[iz] - z0[iz-1])*fracent*( Qv_lif[iz-1] - q0[iz-1] )
                T0[iz] = T_lif[iz]
                q0[iz] = Qv_lif[iz]
                #INTEGRATE UPWARD TO GET PRESSURE HYDROSTATICALLY
                p0[iz] = p0[iz-1] + (z0[iz]-z0[iz-1])*( -g*p0[iz-1]/(Rd*T0[iz-1]*(1 + (Rd/Rv-1)*q0[iz-1])  ) )

                    
                Qt_lif[iz] = Qv_lif[iz]
                q_sat=(1-Qt_lif[iz])*compute_rsat(T_lif[iz],p0[iz],1,T1,T2)
                
                if Qv_lif[iz]>=q_sat: #if we hit saturation, split the vertical step into two stages.  The first stage advances at the saturated lapse rate to the saturation point, and the second stage completes the grid step at the moist lapse rate
    
                    OMEGA = omega(T_lif[iz-1],T1,T2)
                    dOMEGA = domega(T_lif[iz-1],T1,T2)
                    satrat=(Qv_lif[iz]-q_sat_prev)/(q_sat-q_sat_prev)
                    dz_dry=satrat*(z0[iz]-z0[iz-1])
                    dz_wet=(1-satrat)*(z0[iz]-z0[iz-1])
    
                    
                    
                    T_halfstep = T_lif[iz-1] + dz_dry*drylift(T_lif[iz-1],Qv_lif[iz-1],T0[iz-1],q0[iz-1],fracent)
                    Qv_halfstep = Qv_lif[iz-1] - dz_dry*fracent*( Qv_lif[iz-1] - q0[iz-1] )
                    Qt_halfstep = Qv_lif[iz]
                    p_halfstep=p0[iz-1]*satrat + p0[iz]*(1-satrat)
                    T0_halfstep=T0[iz-1]*satrat + T0[iz]*(1-satrat)
                    Q0_halfstep=q0[iz-1]*satrat + q0[iz]*(1-satrat)
    
                    T_lif[iz] = T_halfstep + dz_wet*moislif(T_halfstep,Qv_halfstep,(1-Qt_halfstep)*compute_rsat(T_halfstep,p_halfstep,0,T1,T2),(1-Qt_halfstep)*compute_rsat(T_halfstep,p_halfstep,2,T1,T2),p_halfstep,T0_halfstep,Q0_halfstep,Qt_halfstep,fracent,prate,T1,T2)
                    
                    #INTEGRATE UPWARD TO GET PRESSURE HYDROSTATICALLY
                    p0[iz] = p0[iz-1] + (z0[iz]-z0[iz-1])*( -g*p0[iz-1]/(Rd*T0[iz-1]*(1 + (Rd/Rv-1)*q0[iz-1])  ) )
                    
                    Qt_lif[iz] = Qt_lif[iz-1] - (z0[iz] - z0[iz-1])*fracent*( Qt_halfstep - Q0_halfstep )
                    Qv_lif[iz] = (1-Qt_lif[iz])*compute_rsat(T_lif[iz],p0[iz],1,T1,T2)
                    
                    T0[iz] = T_lif[iz]*( 1 + (Rv/Rd)*Qv_lif[iz] - Qt_lif[iz])/(1 + RH_ft[iz]*(Rv/Rd-1)*Qv_lif[iz] )
                    q0[iz] = RH_ft[iz]*Qv_lif[iz]
                    
    
                    if Qt_lif[iz]<Qv_lif[iz]:
                        Qv_lif[iz]=Qt_lif[iz]
    
                q_sat_prev=q_sat;
                
            else: #if we are already at saturation, just advance upward using the saturated lapse rate (eq. 24 in Peters et al. 2022)
                OMEGA = omega(T_lif[iz-1],T1,T2)
                dOMEGA = domega(T_lif[iz-1],T1,T2)
    
                #if iz==start_loc+1 or satflag:
                #    T_lif[iz] = T_lif[iz-1] + (z0[iz] - z0[iz-1])*moislif(T_lif[iz-1],Qv_lif[iz-1],(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],0,T1,T2),(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],2,T1,T2),p0[iz-1],T0[iz-1],q0[iz-1],Qt_lif[iz-1],fracent,prate,T1,T2);
                #    Qt_lif[iz] = Qt_lif[iz-1] - (z0[iz] - z0[iz-1])*(fracent*( Qt_lif[iz-1] - q0[iz-1] )  + prate*( Qt_lif[iz-1]-Qv_lif[iz-1]) )
                #    satflag = False
                #else:
                T_lif[iz] = T_lif[iz-1] + (z0[iz] - z0[iz-1])*moislif(T_lif[iz-1],Qv_lif[iz-1],(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],0,T1,T2),(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],2,T1,T2),p0[iz-1],T0[iz-1],q0[iz-1],Qt_lif[iz-1],fracent,prate,T1,T2);
                Qt_lif[iz] = Qt_lif[iz-1] - (z0[iz] - z0[iz-1])*(fracent*( Qt_lif[iz-1] - q0[iz-1] )  + prate*( Qt_lif[iz-1]-Qv_lif[iz-1]) )
                
                #INTEGRATE UPWARD TO GET PRESSURE HYDROSTATICALLY
                p0[iz] = p0[iz-1] + (z0[iz]-z0[iz-1])*( -g*p0[iz-1]/(Rd*T0[iz-1]*(1 + (Rd/Rv-1)*q0[iz-1])  ) )        
                Qv_lif[iz] = (1-Qt_lif[iz])*compute_rsat(T_lif[iz],p0[iz],1,T1,T2)
                T0[iz] = T_lif[iz]*( 1 + (Rv/Rd)*Qv_lif[iz] - Qt_lif[iz])/(1 + RH_ft[iz]*(Rv/Rd-1)*Qv_lif[iz] )
                q0[iz] = RH_ft[iz]*Qv_lif[iz]
                
                if Qt_lif[iz]<Qv_lif[iz]:
                    Qv_lif[iz]=Qt_lif[iz]
        else:
            T0[iz] = T0[iz-1]
            p0[iz] = p0[iz-1] + (z0[iz]-z0[iz-1])*( -g*p0[iz-1]/(Rd*T0[iz-1]*(1 + (Rd/Rv-1)*q0[iz-1])  ) )  
            rsat = compute_rsat(T0[iz],p0[iz],1,T1,T2)
            q0[iz] = RH_ft[iz]*(rsat/(1 + rsat))


    T_rho_lif = T_lif*(1 + (Rv/Rd)*Qv_lif - Qt_lif)
    T_0_lif = T0*(1 + (Rv/Rd - 1)*q0)

    
    B_lif=g*(T_rho_lif - T_0_lif)/T_0_lif
    
    
    return T0,q0,p0,B_lif

#==============================================================================
#==============================================================================
#FUNCTION THAT COMPUTES ENVIRONMENTAL PROFILE USING ZERO BUOYANCY MODEL
def calculate_dewpoint(T,p,iceflag,T1,T2,RH):
    esat = compute_esat_zb(T,p,iceflag,T1,T2)
    
    pa = (RH/1)*esat
    c = 243
    a=eref
    b=17
    Td = c*np.log(pa/a)/(b - np.log(pa/a))
    return Td
    
#==============================================================================
#FUNCTION THAT CALCULATES THE SATURATION VAPOR PRESSURE
def compute_esat_zb(T,p,iceflag,T1,T2):
    
    #THIS FUNCTION COMPUTES THE SATURATION MIXING RATIO, USING THE INTEGRATED
    #CLAUSIUS CLAPEYRON EQUATION (eq. 7-12 in Peters et al. 2022).
    #https://doi-org.ezaccess.libraries.psu.edu/10.1175/JAS-D-21-0118.1 

    #input arguments
    #T temperature (in K)
    #p pressure (in Pa)
    #iceflag (give mixing ratio with respect to liquid (0), combo liquid and
    #ice (2), or ice (3)
    #T1 warmest mixed-phase temperature
    #T2 coldest mixed-phase temperature
    
    #NOTE: most of my scripts and functions that use this function need
    #saturation mass fraction qs, not saturation mixing ratio rs.  To get
    #qs from rs, use the formula qs = (1 - qt)*rs, where qt is the total
    #water mass fraction

    rsat = compute_rsat(T,p,iceflag,T1,T2)

    esat = rsat*p/(Rd/Rv + rsat )
    return esat
#==============================================================================


#==============================================================================
#FUNCTION THAT COMPUTES ENVIRONMENTAL PROFILE USING ZERO BUOYANCY MODEL
def ROMPS_analytic_profile(T_sfc,p_sfc,RH_sfc,fracent,fracdet,z0,T1,T2,T_fat,alpha_ratio = 0.5):
    start_loc = 0
    #[T_lif,Qv_lif,Qt_lif,B_lif]
    varepsilon = fracent
    prate = 0
    RH_ft = 0.1
    delta = fracdet
    

    #this function computes lifted parcel properties using the unsaturated
    #and saturated lapse rate formulas from (Peters et al. 2022)
    #https://doi-org.ezaccess.libraries.psu.edu/10.1175/JAS-D-21-0118.1 
    
    #input arguments
    #T0: sounding profile of temperature (in K)
    #p0: sounding profile of pressure (in Pa)
    #q0: sounding profile of water vapor mass fraction (in kg/kg)
    #start_loc: index of the parcel starting location (set to 1 for the
    #lowest: level in the sounding)
    #fracent: fractional entrainment rate (in m^-1)
    
    #output arguments
    #T_lif: lifted parcel temperature
    #Qv_lif: lifted parcel water vapor mass fraction
    #Qt_lif: lifted parcel total water mass fraction
    #B_lif: Lifted parcel buoyancy, computed using Eq. B6 in (Peters et al.
    #2022) (accounts for virtual temperature and loading effects)
    
    #prate: precipitation rate (in m^-1) large values make parcel more
    #pseudoadiabatic, small values make parcel more adiabatic.  I usually
    #just set it to 0 to get an adiabatic parce
    
    #z0: sounding profile of height above ground level (first height should
    #be 0 m)
    #T1 warmest mixed-phase temperature
    #T2 coldest mixed-phase temperature

    #CONSTANTS

    
    
    #descriminator function between liquid and ice (i.e., omega defined in the
    #beginning of section 2e in Peters et al. 2022)

    
    T_lif=np.zeros(z0.shape)*np.nan #temperature of the lifted parcel
    Qv_lif=np.zeros(z0.shape)*np.nan #water vapor mass fraction of the lifted parcel
    Qt_lif=np.ones(z0.shape)*np.nan #total water mass fraction of the lifted parcel
    T0 = np.copy(T_lif)
    q0 = np.copy(Qv_lif)
    p0 = np.copy(Qv_lif)


    T_lif[0]=T_sfc #set initial values to that of the environment
    rsat_sfc = compute_rsat(T_sfc,p_sfc,0,T1,T2)
    qsat_sfc = rsat_sfc/(1 + rsat_sfc)
    Qv_lif[0]=qsat_sfc*RH_sfc #set initial values to that of the environment
    Qt_lif[0]=qsat_sfc*RH_sfc #set initial values to that of the environment
    
    T0[0] = T_sfc
    q0[0] = qsat_sfc*RH_sfc
    p0[0] = p_sfc


    q_sat_prev=0
    B_run = 0
    iz=0
    #
    #for iz in np.arange(start_loc+1,z0.shape[0]):
    #
    #
    #I REVISED THIS A BIT.  TO MAKE THE CODE FASTER, I HAVE THE CALCULATION CUT OUT WHEN THE INTEGRATED NEGATIVE BUOYANCY ("BRUN") 
    #BECOMES MORE NEGATIVE THAN THAN THE TOTAL INTEGRATED POSITIVE BUOYANCY.  I RESTRICT THIS TO ONLY HAPPEN AFTER WE HAVE PASSED 
    #THE HEIGHT OF MINIMUM MSE.  UNCOMMENT THE FOR LOOP ABOVE AND COMMENT OUT THE WHILE LOOP IF YOU JUST WANT TO INTEGRATE TO THE TOP OF THE SOUNDING.
    #THE +25 PART IN THE WHILE STATEMENT IS A PAD ON B_RUN (THE NEGATIVE CAPE HAS TO BE 25 J/KG LESS THAN THE POSITIVE CAPE TO KILL THE LOOP)
    satflag = True
    #while iz<(z0.shape[0])-1 and (z0[iz]<z0[mn_hgt] or (B_run+25)>0):
    while iz<(z0.shape[0])-1:
        iz = iz + 1
        q_sat=(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],1,T1,T2)
        
        if T0[iz-1]>T_fat:
            
            if Qv_lif[iz-1]>=q_sat:
                satflag = False
            if satflag: #if we are unsaturated, go up at the unsaturated adiabatic lapse rate (eq. 19 in Peters et al. 2022)
                
            

                T_lif[iz] = T_lif[iz-1] + (z0[iz] - z0[iz-1])*drylift(T_lif[iz-1],Qv_lif[iz-1],T0[iz-1],q0[iz-1],fracent)
                Qv_lif[iz] = Qv_lif[iz-1] - (z0[iz] - z0[iz-1])*fracent*( Qv_lif[iz-1] - q0[iz-1] )
                T0[iz] = T_lif[iz]
                q0[iz] = Qv_lif[iz]
                #INTEGRATE UPWARD TO GET PRESSURE HYDROSTATICALLY
                p0[iz] = p0[iz-1] + (z0[iz]-z0[iz-1])*( -g*p0[iz-1]/(Rd*T0[iz-1]*(1 + (Rd/Rv-1)*q0[iz-1])  ) )

                    
                Qt_lif[iz] = Qv_lif[iz]
                q_sat=(1-Qt_lif[iz])*compute_rsat(T_lif[iz],p0[iz],1,T1,T2)
                
                if Qv_lif[iz]>=q_sat: #if we hit saturation, split the vertical step into two stages.  The first stage advances at the saturated lapse rate to the saturation point, and the second stage completes the grid step at the moist lapse rate
    
                    #compute_qsat(1)
                    OMEGA = omega(T_lif[iz-1],T1,T2)
                    dOMEGA = domega(T_lif[iz-1],T1,T2)
                    satrat=(Qv_lif[iz]-q_sat_prev)/(q_sat-q_sat_prev)
                    dz_dry=satrat*(z0[iz]-z0[iz-1])
                    dz_wet=(1-satrat)*(z0[iz]-z0[iz-1])
    
                    
                    
                    T_halfstep = T_lif[iz-1] + dz_dry*drylift(T_lif[iz-1],Qv_lif[iz-1],T0[iz-1],q0[iz-1],fracent)
                    Qv_halfstep = Qv_lif[iz-1] - dz_dry*fracent*( Qv_lif[iz-1] - q0[iz-1] )
                    Qt_halfstep = Qv_lif[iz]
                    p_halfstep=p0[iz-1]*satrat + p0[iz]*(1-satrat)
                    T0_halfstep=T0[iz-1]*satrat + T0[iz]*(1-satrat)
                    Q0_halfstep=q0[iz-1]*satrat + q0[iz]*(1-satrat)
    
                    Gamma,mgamma = mgamma_romps(T_halfstep,p_halfstep,T1,T2,varepsilon,varepsilon)
                    T_lif[iz] = T_halfstep + dz_wet*moislif(T_halfstep,Qv_halfstep,(1-Qt_halfstep)*compute_rsat(T_halfstep,p_halfstep,0,T1,T2),(1-Qt_halfstep)*compute_rsat(T_halfstep,p_halfstep,2,T1,T2),p_halfstep,T0_halfstep,Q0_halfstep,Qt_halfstep,fracent,prate,T1,T2)
                    
                    #INTEGRATE UPWARD TO GET PRESSURE HYDROSTATICALLY
                    p0[iz] = p0[iz-1] + (z0[iz]-z0[iz-1])*( -g*p0[iz-1]/(Rd*T0[iz-1]*(1 + (Rd/Rv-1)*q0[iz-1])  ) )
                    
                    Qt_lif[iz] = Qt_lif[iz-1] - (z0[iz] - z0[iz-1])*fracent*( Qt_halfstep - Q0_halfstep )
                    Qv_lif[iz] = (1-Qt_lif[iz])*compute_rsat(T_lif[iz],p0[iz],1,T1,T2)
                    
                    T0[iz] = T0_halfstep + dz_wet*Gamma
                    #comp_qsat(1)
                    RH_on = compute_RH_romps(varepsilon,delta,-mgamma,alpha_ratio)
                    #RH_on = ( delta + alpha_ratio*delta  - alpha_ratio*varepsilon)/(delta - mgamma - alpha_ratio*varepsilon)
                    rsat = compute_rsat( T0[iz] , p0[iz],1,T1,T2)
                    q0[iz] = (rsat/(1+rsat))*RH_on

    
                    if Qt_lif[iz]<Qv_lif[iz]:
                        Qv_lif[iz]=Qt_lif[iz]
    
                q_sat_prev=q_sat;
                
            else: #if we are already at saturation, just advance upward using the saturated lapse rate (eq. 24 in Peters et al. 2022)
                OMEGA = omega(T_lif[iz-1],T1,T2)
                dOMEGA = domega(T_lif[iz-1],T1,T2)
    
                #if iz==start_loc+1 or satflag:
                #    T_lif[iz] = T_lif[iz-1] + (z0[iz] - z0[iz-1])*moislif(T_lif[iz-1],Qv_lif[iz-1],(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],0,T1,T2),(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],2,T1,T2),p0[iz-1],T0[iz-1],q0[iz-1],Qt_lif[iz-1],fracent,prate,T1,T2);
                #    Qt_lif[iz] = Qt_lif[iz-1] - (z0[iz] - z0[iz-1])*(fracent*( Qt_lif[iz-1] - q0[iz-1] )  + prate*( Qt_lif[iz-1]-Qv_lif[iz-1]) )
                #    satflag = False
                #else:
                T_lif[iz] = T_lif[iz-1] + (z0[iz] - z0[iz-1])*moislif(T_lif[iz-1],Qv_lif[iz-1],(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],0,T1,T2),(1-Qt_lif[iz-1])*compute_rsat(T_lif[iz-1],p0[iz-1],2,T1,T2),p0[iz-1],T0[iz-1],q0[iz-1],Qt_lif[iz-1],fracent,prate,T1,T2);
                Qt_lif[iz] = Qt_lif[iz-1] - (z0[iz] - z0[iz-1])*(fracent*( Qt_lif[iz-1] - q0[iz-1] )  + prate*( Qt_lif[iz-1]-Qv_lif[iz-1]) )
                
                #INTEGRATE UPWARD TO GET PRESSURE HYDROSTATICALLY
                p0[iz] = p0[iz-1] + (z0[iz]-z0[iz-1])*( -g*p0[iz-1]/(Rd*T0[iz-1]*(1 + (Rd/Rv-1)*q0[iz-1])  ) )        
                Qv_lif[iz] = (1-Qt_lif[iz])*compute_rsat(T_lif[iz],p0[iz],1,T1,T2)
                
                Gamma,mgamma = mgamma_romps(T0[iz-1],p0[iz-1],T1,T2,varepsilon,varepsilon,alphaRH = RH_on*alpha_ratio)
                
                T0[iz] = T0[iz-1] + (z0[iz] - z0[iz-1])*Gamma
                RH_on = compute_RH_romps(varepsilon,delta,-mgamma,alpha_ratio)
                #RH_on = ( delta + alpha_ratio*delta  - alpha_ratio*varepsilon)/(delta - mgamma - alpha_ratio*varepsilon)
                rsat = compute_rsat( T0[iz] , p0[iz],1,T1,T2)
                q0[iz] = (rsat/(1+rsat))*RH_on
                
                if Qt_lif[iz]<Qv_lif[iz]:
                    Qv_lif[iz]=Qt_lif[iz]
        else:
            T0[iz] = T0[iz-1]
            p0[iz] = p0[iz-1] + (z0[iz]-z0[iz-1])*( -g*p0[iz-1]/(Rd*T0[iz-1]*(1 + (Rd/Rv-1)*q0[iz-1])  ) )  
            rsat = compute_rsat(T0[iz],p0[iz],1,T1,T2)
            q0[iz] = RH_ft*(rsat/(1 + rsat))


    T_rho_lif = T_lif*(1 + (Rv/Rd)*Qv_lif - Qt_lif)
    T_0_lif = T0*(1 + (Rv/Rd - 1)*q0)

    
    B_lif=g*(T_rho_lif - T_0_lif)/T_0_lif
    
    rsat = compute_rsat(T0,p0,1,273.15,253.15)

    qsat = rsat/(1 + rsat)
    
    RH0 = q0/qsat
    
    return T0,q0,p0,RH0

#==============================================================================

#==============================================================================

def Gamma_romps(T,p,T1,T2,varepsilon,delta,alphaRH = 0):
    
    rsat = compute_rsat(T,p,1,T1,T2)
    qsat = rsat/(1 + rsat)
    a_1 =  ( Rv*cp*T**2 )/( xlv ) + qsat*xlv
    a_2 = ( ( Rv*cp*T**2 )/xlv )*(delta -alphaRH*varepsilon + g/(Rd*T)) + qsat*xlv*(delta - varepsilon) - g
    a_3 = (Rv*cp*T/(Rd*xlv) - 1 )*g*(delta - alphaRH*varepsilon)
    
    Gamma = -( ( Rv*T**2 )/xlv )*( ( -a_2 + np.sqrt(a_2**2 - 4*a_1*a_3) )/(2*a_1) + g/(Rd*T) )
    
    return Gamma

#==============================================================================

#==============================================================================

def mgamma_romps(T,p,T1,T2,varepsilon,delta,alphaRH = 0):
    Gamma = Gamma_romps(T,p,T1,T2,varepsilon,delta,alphaRH = alphaRH)
    mgamma = -( -xlv*Gamma/(Rv*T**2) - g/(Rd*T) )
    
    return Gamma,mgamma



#==============================================================================

#==============================================================================

def get_all_romps(T_sfc,p_sfc,RH_sfc,fracent,fracdet,z0,T1,T2,T_fat,second_ent,alpha_ratio = 0.5):
    T0,q0,p0,RH0 = ROMPS_analytic_profile(T_sfc,p_sfc,RH_sfc,fracent,fracdet,z0,T1,T2,T_fat,alpha_ratio = alpha_ratio)
    
    
    T_UD,Qv_UD,Qt_UD,B_UD = lift_parcel_adiabatic_leapfrog(T0,p0,q0,0,0,0,z0,T1,T2)
    T_D,Qv_D,Qt_D,B_D = lift_parcel_adiabatic_leapfrog(T0,p0,q0,0,second_ent,0,z0,T1,T2)
    
    Trho_UD = T_UD*(1 + (Rd/Rv)*Qv_UD - Qt_UD)
    
    
    CAPE,CIN,LFC,EL = compute_CAPE_AND_CIN(T0,p0,q0,0,0,0,z0,T1,T2,useleapfrog=True)
    ECAPE,CIN,LFC,EL = compute_CAPE_AND_CIN(T0,p0,q0,0,second_ent,0,z0,T1,T2,useleapfrog=True)
    
    return T0,p0,q0,RH0,Trho_UD,CAPE,ECAPE


#==============================================================================

#==============================================================================

def compute_RH_romps(epsilon,delta,mgamma,beta):
    if beta>0:    
        a = 1
        b = -(delta + mgamma - beta*mgamma + beta*epsilon)/(beta*epsilon)
        c = delta/(beta*epsilon)
        RH = ( -b - np.sqrt(b**2 - 4*c) )/2
    else:
        RH = delta/(delta + mgamma)
    return RH


#==============================================================================

def ALL_ECAPE_PARAMS(T0,p0,z0,q0,u0,v0,useenergy=True,useleapfrog=True):
    L_mix = 120
    T1 = 273.15
    T2 = 253.15
    CAPE_ud,CIN,LFC,EL=compute_CAPE_AND_CIN(T0,p0,q0,0,0,0,z0,T1,T2,useenergy=useenergy,useleapfrog=useleapfrog)
    #GET NCAPE, WHICH IS NEEDED FOR ECAPE CALULATION
    NCAPE,MSE0_star,MSE0bar=compute_NCAPE(T0,p0,q0,z0,T1,T2,LFC,EL)

    #GET THE 0-1 KM MEAN STORM-RELATIVE WIND, ESTIMATED USING BUNKERS METHOD FOR RIGHT-MOVER STORM MOTION
    V_SR,C_x,C_y = compute_VSR(z0,u0,v0)
    

    EL_on = []
    ECAPE_on = []
    EL_ud = EL

    #THIS IS A PROCEEDURE THAT ITERATIVELY SETTLES ONTO THE ENTRAINING EL
    for iit in range(0,10):
        E_tilde,varepsilon,Radius  = compute_ETILDE(CAPE_ud,NCAPE,V_SR,EL,L_mix,dynamfac=1)

        ECAPE,ECIN,ELFC,EL=compute_CAPE_AND_CIN(T0,p0,q0,0,varepsilon,0,z0,T1,T2,useleapfrog=True)
        #CAPE,CIN,LFC,EL=compute_CAPE_AND_CIN(T0,p0,q0,0,varepsilon,0,z0,T1,T2)
        
        EL_on.append(EL)
        ECAPE_on.append(ECAPE)

    E_tilde = ECAPE/CAPE_ud
    
    EL_tilde = EL/EL_ud
    
    return ECAPE,E_tilde,ECIN,EL,EL_tilde

#==============================================================================
#==============================================================================
#=============================================================================
#END FUNCTION DEFINITIONS======================================================
#==============================================================================
#==============================================================================
#==============================================================================
