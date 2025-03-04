Module cal_cape
IMPLICIT NONE

REAL*8, PARAMETER    :: alpha=0.5,& 
                        sigma=1.,&  
                        k2=0.18,&  
                        pr=0.33333333,&  
                        lmix=120.,&
                        pi=3.14159265359,&  
                        fraclen=1.,&  
                        t1=273.15,&  
                        t2=253.15,&
                        prate=0.,&  
                        rd=287.04,&  
                        rv=461.5,&  
                        cp=1005.,&  
                        hv=2.26e6,&
                        g=9.81,&  
                        xlv=2.501e6,&  
                        th0r=273.15,&  
                        xls=2.834e6,&
                        cpv=1870.,& 
                        cpl=4190.,& 
                        cpi=2106.,& 
                        pref=611.65,&
                        ttrip=273.15,& 
                        eref=611.2,& 
                        cvv=cpv-rv,& 
                        epsilon=rd/rv
INTEGER,  DIMENSION(12)  :: dlist
DATA dlist/31,28,31,30,31,30,31,31,30,31,30,31/

contains

real*8 function qsat(T,p,iceflag)
implicit none
integer   :: iceflag
real*8    :: T, p
!local
real*8    :: term1, term2, esl, esl_l, qsat_l, esl_i, qsat_i, omeg, &
             h1, h2, h3

  h1=(T - t1)/(t2-t1)
  h2=(1 - (T - t1)/(t2-t1)) 
  h3=-(1 - (T - t1)/(t2-t1))
  if (h1 .lt. 0.) then
   h1=0.
  elseif(h1 .eq. 0.) then
   h1=0.5
  elseif(h1 .gt. 0.) then
   h1=1.
  endif
  if (h2 .lt. 0.) then
   h2=0.
  elseif(h2 .eq. 0.) then
   h2=0.5
  elseif(h2 .gt. 0.) then
   h2=1.
  endif
  if (h3 .lt. 0.) then
   h3=0.
  elseif(h3 .eq. 0.) then
   h3=0.5
  elseif(h3 .gt. 0.) then
   h3=1.
  endif

  omeg = ((T - t1)/(t2-t1))*h1*h2 + h3
  if (iceflag .eq. 0) then
   term1=(cpv-cpl)/rv
   term2=(xlv-ttrip*(cpv-cpl))/rv
   esl = dexp((T-ttrip)*term2/(T*ttrip))*eref*(T/ttrip)**term1
   esl = min( esl , p*0.5 )
   qsat=epsilon*esl/(p-esl)
  elseif (iceflag .eq. 1) then
   term1=(cpv-cpl)/rv
   term2=(xlv-ttrip*(cpv-cpl))/rv
   esl_l=dexp((T-ttrip)*term2/(T*ttrip))*eref*(T/ttrip)**term1
   qsat_l=epsilon*esl_l/(p-esl_l)

   term1=(cpv-cpi)/rv
   term2=(xls-ttrip*(cpv-cpi))/rv
   esl_i=dexp((T-ttrip)*term2/(T*ttrip))*eref*(T/ttrip)**term1
   qsat_i=epsilon*esl_i/(p-esl_i)
   qsat=(1-omeg)*qsat_l +omeg*qsat_i
  elseif (iceflag .eq. 2) then
   term1=(cpv-cpi)/rv
   term2=(xls-ttrip*(cpv-cpi))/rv
   esl=dexp((T-ttrip)*term2/(T*ttrip))*eref*(T/ttrip)**term1
   esl=min( esl , p*0.5 )
   qsat=epsilon*esl/(p-esl)
  endif

end function qsat

subroutine lift_parcel_adiabatic(lev,t0,p0,q0,z0,buoc)
implicit none
integer,intent(in)                :: lev
real*8,intent(in),dimension(lev)  :: t0, p0, q0, z0
real*8,intent(out),dimension(lev) :: buoc
!local
integer                           :: iz
real*8,dimension(lev)             :: T_lif, Qv_lif, Qt_lif, T_rho_lif, T_0_lif
real*8                            :: fracent, q_sat_prev, q_sat, satrat, dz_dry
real*8                            :: dz_wet, T_halfstep, Qv_halfstep, Qt_halfstep
real*8                            :: p_halfstep, T0_halfstep, Q0_halfstep
real*8                            :: qent_ref


  T_lif=0.
  Qv_lif=0.
  Qt_lif=0.

  T_lif(1)=t0(1)
  Qv_lif(1)=q0(1)
  Qt_lif=Qv_lif

  q_sat_prev=0.
  do iz = 2, lev
   q_sat=(1-Qt_lif(iz-1))*qsat(T_lif(iz-1),p0(iz-1),1)
   if (Qv_lif(iz-1) .lt. q_sat) then
    T_lif(iz) = T_lif(iz-1) + (z0(iz) - z0(iz-1))*gamma_d(T_lif(iz-1),Qv_lif(iz-1),T0(iz-1),q0(iz-1),fracent)
    Qv_lif(iz) = Qv_lif(iz-1) - (z0(iz) - z0(iz-1))*fracent*( Qv_lif(iz-1) - q0(iz-1) )
    Qt_lif(iz) = Qv_lif(iz)
    q_sat=(1-Qt_lif(iz))*qsat(T_lif(iz),p0(iz),1)

          if (Qv_lif(iz) .gt. q_sat) then
              satrat=(Qv_lif(iz)-q_sat_prev)/(q_sat-q_sat_prev)
              dz_dry=satrat*(z0(iz)-z0(iz-1))
              dz_wet=(1-satrat)*(z0(iz)-z0(iz-1))
              T_halfstep = T_lif(iz-1) + dz_dry*gamma_d(T_lif(iz-1),Qv_lif(iz-1),T0(iz-1),q0(iz-1),fracent)
              Qv_halfstep = Qv_lif(iz-1) - dz_dry*fracent*( Qv_lif(iz-1) - q0(iz-1) )
              Qt_halfstep = Qv_lif(iz)
              p_halfstep=p0(iz-1)*satrat + p0(iz)*(1-satrat)
              T0_halfstep=T0(iz-1)*satrat + T0(iz)*(1-satrat)
              Q0_halfstep=q0(iz-1)*satrat + q0(iz)*(1-satrat)

              T_lif(iz) = T_halfstep + dz_wet*gamma_m(T_halfstep,Qv_halfstep,&
                  (1-Qt_halfstep)*qsat(T_halfstep,p_halfstep,0),(1-Qt_halfstep)*qsat(T_halfstep,p_halfstep,2),p_halfstep,T0_halfstep,Q0_halfstep&
                  ,omega(T_halfstep),domega(T_halfstep),Qt_halfstep,fracent)
              qent_ref=Qt_halfstep
              
              Qt_lif(iz) = Qt_lif(iz-1) - (z0(iz) - z0(iz-1))*fracent*( qent_ref - Q0_halfstep )
              Qv_lif(iz) = (1-Qt_lif(iz))*qsat(T_lif(iz),p0(iz),1)

              if (Qt_lif(iz) .lt. Qv_lif(iz)) then
                  Qv_lif(iz)=Qt_lif(iz)
              endif
          endif
          q_sat_prev=q_sat
   else
          T_lif(iz) = T_lif(iz-1) + (z0(iz) - z0(iz-1))*gamma_m(T_lif(iz-1),Qv_lif(iz-1),&
              (1-Qt_lif(iz-1))*qsat(T_lif(iz-1),p0(iz-1),0),(1-Qt_lif(iz-1))*qsat(T_lif(iz-1),p0(iz-1),2),p0(iz-1),T0(iz-1),q0(iz-1)&
              ,omega(T_lif(iz-1)),domega(T_lif(iz-1)),Qt_lif(iz-1),fracent)
          
          qent_ref=Qt_lif(iz-1)
          Qt_lif(iz) = Qt_lif(iz-1) - (z0(iz) - z0(iz-1))*(fracent*( qent_ref - q0(iz-1) )  + max(prate*( Qt_lif(iz-1)-Qv_lif(iz-1)),0.) )
          Qv_lif(iz) = (1-Qt_lif(iz))*qsat(T_lif(iz),p0(iz),1)
          if (Qt_lif(iz) .lt. Qv_lif(iz)) then
              Qv_lif(iz)=Qt_lif(iz)
          endif
   endif
  enddo
  
  T_rho_lif=T_lif*(1 - Qt_lif + Qv_lif)/( 1 + (epsilon - 1)/( ( epsilon*(1 - Qt_lif)/Qv_lif - 1) ) )
  T_0_lif=T0/( 1 + (epsilon - 1)/( ( epsilon*(1 - q0)/q0 - 1) ) )
  
  buoc=g*(T_rho_lif - T_0_lif)/T_0_lif

end subroutine lift_parcel_adiabatic

real*8 function omega(T)
implicit none
real*8   :: T
  omega=min(max((T - t1)/(t2-t1),0.),1.)
end function omega

real*8 function domega(T)
implicit none
real*8   :: T
real*8   :: h1, h2
  h1=(t1-T)
  h2=(t2-T)
  if (h1 .lt. 0.) then
   h1=0.
  elseif(h1 .eq. 0.) then
   h1=0.5
  elseif(h1 .gt. 0.) then
   h1=1.
  endif
  if (h2 .lt. 0.) then
   h2=0.
  elseif(h2 .eq. 0.) then
   h2=0.5
  elseif(h2 .gt. 0.) then
   h2=1.
  endif
  domega=(h1-h2)/(t2-t1)
end function domega


real*8 function gamma_m(T,qv,qvv,qvi,p0,T0,q0,omega,domega,qt,fracent)
implicit none
real*8  :: T, qv, qvv, qvi, p0, T0, q0, omega, domega, qt, fracent
real*8  :: zeroB, Lv, Li, Rm0, pseudofac, T_rho, T_rho0, B
real*8  :: Qvsl, Qvsi, Q_M, L_M, eps_T, eps_qv, eps_qt, cpm, eps
real*8  :: term1, term2, term3, term4, term5, term6, term7, term8
logical :: noT
  if (qt .lt. 0.) then
   zeroB=0.
   qt=-qt
  else
   zeroB=1.
  endif
  
  if (qv .lt. 0.) then
   noT=.true.
   qv=-qv
  else
   noT=.false.
  endif
    
  cpm = (1 - qt)*cp + qv*cpv + (1 - omega)*(qt-qv)*cpl + omega*(qt-qv)*cpi
  Lv = xlv + (T - ttrip)*(cpv - cpl)
  Li = (xls-xlv) + (T - ttrip)*(cpl - cpi);
  Rm0 = (1 - q0)*Rd + q0*Rv
    
  if (qt .eq. qv) then
   pseudofac=(1-qv)
  else
   pseudofac=1
  endif
  
  T_rho=T*(1 - qt + qv/epsilon)
  T_rho0=T0*( 1 - q0 + q0/epsilon )
  B=zeroB*g*(T_rho - T_rho0)/(T_rho0)
    
  Qvsl = qvv/( epsilon - epsilon*qt + qv)
  Qvsi = qvi/( epsilon - epsilon*qt + qv)
  Q_M = (1 - omega)*qvv/(1 - Qvsl) + omega*qvi/(1 - Qvsi)
  L_M = Lv*(1 - omega)*qvv/(1 - Qvsl) + (Lv + Li)*omega*qvi/(1 - Qvsi)
    
  eps_T = -fracent*(T - T0)
  eps_qv = -fracent*(qv - q0)
  eps_qt = -fracent*(qt - q0)
  term1 = -B
  if (noT) then
   term2 = - pseudofac*Q_M*(Lv + Li*omega)*g/(Rd*T_rho)
  else
   term2 = - pseudofac*Q_M*(Lv + Li*omega)*g/(Rm0*T0)
  endif
  term3 = -g
  term4 = (cpm - Li*(qt - qv)*domega)*eps_T
  term5 = (Lv + Li*omega)*(eps_qv + (qv/(1-qt))*eps_qt)
  term6 = cpm
  term7 = -Li*(qt - qv)*domega
  term8 = pseudofac*(Lv + Li*omega)*(-domega*(qvv - qvi) + (1/(Rv*T**2))*(L_M))
  eps = -fracent*(T - T0)
  gamma_m =( term1 + term2 + term3 + term4 + term5)/(term6 + term7 + term8)
end function gamma_m

real*8 function gamma_d(T,qv,T0,qv0,fracent)
implicit none
real*8  :: T, qv, T0, qv0, fracent
real*8  :: cpmv, B, eps
  cpmv = (1 - qv)*cp + qv*cpv
  B = g*( (T-T0)/T0 + (Rv/Rd - 1)*(qv - qv0) )
  eps = -fracent*(T - T0)
  gamma_d = - (g + B)/cpmv + eps
end function gamma_d

end module cal_cape
