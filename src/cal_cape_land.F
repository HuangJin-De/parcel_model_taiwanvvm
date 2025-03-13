program cal_mean
use netcdf
use mpi
use cal_cape
implicit none

integer, parameter :: nx=1024,ny=1024,nz=70,nt=72
real, parameter :: Lv=2.5e6
integer :: i,j,k,m,n,o,t
integer :: ierr, ncid1, varid1
integer :: myid,nproc,ts,te,tt
integer :: idum1,idum2
integer :: mlfc, lfc, el
integer :: is,ie,js,je,ks,ke,k1
integer, dimension(:), allocatable :: mpi_s,mpi_n
real :: vsr,ncape,fork, term1, term2, term3, term4, term5,&
        term6, term7,duma,dumb,dumc,dumd,et,nncape,v2, ent, dume
real :: qv_sat
real :: dum1,dum2,dum3,dum4,dum5
real, dimension(nx,ny,nz) :: th,qv,the,u,v
real, dimension(nz) :: thm,qvm,um,vm,np
real*4, dimension(:,:,:), allocatable :: cape,ecape1,ecape2
real, dimension(nz+1) :: zz,zc,rho,rhoz,thbar,thvbar,pbar,pibar,q1,q2
real, dimension(nz) :: dz
real, dimension(:), allocatable :: t0,q0,p0,z0,h0,hs0,h0h,ht,buoc
real :: start, finish
integer, dimension(nx,ny) :: topo,lu
character(200) :: path,run,filename,workpath

call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world,nproc,ierr)
call mpi_comm_rank(mpi_comm_world,myid,ierr)

allocate(mpi_s(nproc),mpi_n(nproc))

path="/data2/VVM/taiwanvvm_tpe/"
workpath="/data/der0318/parcel_model_taiwanvvm/"
run="runname"
!run="tpe20110816nor"

idum2=nt/nproc
ts=idum2*myid+1
if (ts==1) then
  ts=ts-1
  idum2=idum2+1
endif
if (myid<mod(nt,nproc)) then
  ts=ts+myid
  idum2=idum2+1
else
  ts=ts+mod(nt,nproc)
endif
te=ts+idum2-1
idum1=ts

write(*,*) myid, ts, te

call mpi_allgather(idum2,1,mpi_int,mpi_n,1,mpi_int,mpi_comm_world,ierr)
call mpi_allgather(ts,1,mpi_int,mpi_s,1,mpi_int,mpi_comm_world,ierr)

if (myid==0) then
  allocate(cape(nx,ny,nt+1),ecape1(nx,ny,nt+1),ecape2(nx,ny,nt+1))
else
  allocate(cape(nx,ny,mpi_n(myid+1)),ecape1(nx,ny,mpi_n(myid+1)),ecape2(nx,ny,mpi_n(myid+1)))
endif

!write(*,*) myid, ts, te, idum2
!call mpi_barrier(mpi_comm_world,ierr)

WRITE(filename,111) trim(path),"/",trim(run),"/fort.98"
111 FORMAT(A,A,A,A)

! read mean profiles
OPEN(10,FILE=filename)
DO i=1,188
  READ(10,*)
ENDDO
DO i=1,nz+1
  READ(10,*) j, zz(i), zc(i)
ENDDO
DO i=1,3
  READ(10,*)
ENDDO
DO i=1,nz+1
  READ(10,*) j, rho(i), thbar(i), pbar(i), pibar(i)
ENDDO
DO i=1,3
  READ(10,*)
ENDDO
DO i=1,nz+1
  READ(10,*) j, rhoz(i), rhoz(i), q1(i), q2(i)
ENDDO
DO i=1,3
  READ(10,*)
ENDDO
DO i=2,nz+1
  READ(10,*) j, rhoz(i)
ENDDO
DO k=2,nz
  dz(k)=zz(k)-zz(k-1)
ENDDO
CLOSE(10)

call mpi_barrier(mpi_comm_world,ierr)

cape=-999.
ecape1=-999.
ecape2=-999.

do t=ts,te
  
tt=t-ts+1

! read 3-D data field
WRITE(filename,112) trim(path), "/", trim(run),"/archive/",&
trim(run),".L.Thermodynamic-",t,".nc"
112 FORMAT(6A,I6.6,A)

ierr = nf90_open(trim(filename),NF90_NOWRITE,ncid1)
IF (ierr/=nf90_noerr) WRITE(*,*) "open fail"
ierr = nf90_inq_varid(ncid1,"th",varid1)
IF (ierr/=nf90_noerr) WRITE(*,*) "var_inq fail"
ierr = nf90_get_var(ncid1,varid1,th,start=(/ 1,1,1,1 /),count=(/ nx,ny,nz,1 /))
IF (ierr/=nf90_noerr) WRITE(*,*) "read fail"
ierr = nf90_inq_varid(ncid1,"qv",varid1)
IF (ierr/=nf90_noerr) WRITE(*,*) "var_inq fail"
ierr = nf90_get_var(ncid1,varid1,qv,start=(/ 1,1,1,1 /),count=(/ nx,ny,nz,1 /))
IF (ierr/=nf90_noerr) WRITE(*,*) "read fail"
ierr = nf90_close(ncid1)

WRITE(filename,112) trim(path), "/", trim(run),"/archive/",&
trim(run),".L.Dynamic-",t,".nc"

ierr = nf90_open(trim(filename),NF90_NOWRITE,ncid1)
IF (ierr/=nf90_noerr) WRITE(*,*) "open fail"
ierr = nf90_inq_varid(ncid1,"u",varid1)
IF (ierr/=nf90_noerr) WRITE(*,*) "var_inq fail"
ierr = nf90_get_var(ncid1,varid1,u,start=(/ 1,1,1,1 /),count=(/ nx,ny,nz,1 /))
IF (ierr/=nf90_noerr) WRITE(*,*) "read fail"
ierr = nf90_inq_varid(ncid1,"v",varid1)
IF (ierr/=nf90_noerr) WRITE(*,*) "var_inq fail"
ierr = nf90_get_var(ncid1,varid1,v,start=(/ 1,1,1,1 /),count=(/ nx,ny,nz,1 /))
IF (ierr/=nf90_noerr) WRITE(*,*) "read fail"
ierr = nf90_close(ncid1)

WRITE(filename,113) trim(path), "/", trim(run),"/TOPO.nc"
113 FORMAT(4A)

ierr = nf90_open(trim(filename),NF90_NOWRITE,ncid1)
IF (ierr/=nf90_noerr) WRITE(*,*) "open fail"
ierr = nf90_inq_varid(ncid1,"topo",varid1)
IF (ierr/=nf90_noerr) WRITE(*,*) "var_inq fail"
ierr = nf90_get_var(ncid1,varid1,topo,start=(/ 1,1,1,1 /),count=(/ nx,ny,1,1 /))
IF (ierr/=nf90_noerr) WRITE(*,*) "read fail"
ierr = nf90_inq_varid(ncid1,"lu",varid1)
IF (ierr/=nf90_noerr) WRITE(*,*) "var_inq fail"
ierr = nf90_get_var(ncid1,varid1,lu,start=(/ 1,1,1,1 /),count=(/ nx,ny,1,1 /))
IF (ierr/=nf90_noerr) WRITE(*,*) "read fail"
ierr = nf90_close(ncid1)

where(topo<1) topo=1 

thm=0.
qvm=0.
um=0.
vm=0.
np=0.
do j=1,ny
do i=1,nx
  if (lu(i,j)==16) goto 755
  do k=topo(i,j)+1,nz
    thm(k)=thm(k)+th(i,j,k)
    qvm(k)=qvm(k)+qv(i,j,k)
    um(k)=um(k)+u(i,j,k)
    vm(k)=vm(k)+v(i,j,k)
    np(k)=np(k)+1.
  enddo
  755 continue
enddo
enddo
thm=thm/np
qvm=qvm/np
um=um/np
vm=vm/np

do k=1,nz
do j=1,ny
do i=1,nx
  dum1=th(i,j,k)*pibar(k)
  the(i,j,k)=th(i,j,k)*exp(Lv*qv(i,j,k)/(1004.5*dum1))
enddo
enddo
enddo

n=nz-1
allocate(t0(n),q0(n),p0(n),z0(n),buoc(n))
t0=thm(2:nz)*pibar(2:nz)
q0=qvm(2:nz)
p0=pbar(2:nz)
z0=zc(2:nz)
call lift_parcel_adiabatic(nz-1,t0,p0,q0,z0,buoc)
mlfc=-1

do k=1,nz
  if (buoc(k)>0.) then
    mlfc=k
    exit
  endif
enddo
deallocate(t0,q0,p0,z0,buoc)

!write(*,*) "LFC=",lfc,zc(lfc)

!! CAPE
do j=1,ny
do i=1,nx
  if (lu(i,j)==16 .or. topo(i,j)>=mlfc ) goto 756
  
  !m=minloc(abs(zc-zz(topo(i,j))-1000.),1)
  !o=maxloc(the(i,j,topo(i,j)+1:m),1)  
  
  o=topo(i,j)+1

  n=nz-topo(i,j)
  allocate(t0(n),q0(n),p0(n),z0(n),h0(n),hs0(n),h0h(n),ht(n),buoc(n))
  t0=thm(o:nz)*pibar(o:nz)
  q0=qvm(o:nz)
  p0=pbar(o:nz)
  z0=zc(o:nz)
  h0=cp*t0+g*z0+xlv*q0
  do k=1,n-1
    hs0(k)=cp*t0(k)+g*z0(k)+xlv*qsat(t0(k),p0(k),0)
  enddo

  ! replace with local parcel
  t0(1)=th(i,j,o)*pibar(o)
  q0(1)=qv(i,j,o)
  call lift_parcel_adiabatic(nz-topo(i,j),t0,p0,q0,z0,buoc)
 
  lfc=0
  el=n-1
  lfc=findloc(buoc>0.,.True.,dim=1)
  el=findloc(buoc>0.,.True.,dim=1,back=.True.)
  if (lfc==0 .or. el==n-1) goto 744

  dum1=0.
  do k=lfc,el-1
    if (buoc(k)>0. .and. buoc(k+1)>0.) then
      dum1=dum1+0.5*(buoc(k+1)+buoc(k))*(z0(k+1)-z0(k))
    endif
  enddo  
  cape(i,j,tt)=dum1


  if (cape(i,j,tt)>0.) then
    dum1=0.5*( u(i-1,j,max(topo(i-1,j),topo(i,j))+1) &
              +u(i,j,max(topo(i,j),topo(i+1,j))+1) )
    dum2=0.5*( v(i,j-1,max(topo(i,j-1),topo(i,j))+1) &
              +v(i,j,max(topo(i,j),topo(i,j+1))+1) )
    
    vsr=0.5*sqrt((um(46)-dum1)**2.+(vm(46)-dum2)**2)
  
    h0h=0.
    do k=lfc,el
      do k1=1,k-1
        h0h(k)=h0h(k)+0.5*(h0(k1)+h0(k1+1))*(z0(k1+1)-z0(k1))
      enddo
      h0h(k)=h0h(k)/z0(k)
      ht(k)=g/cp*(h0h(k)-hs0(k))/t0(k)
    enddo 
     
    ncape=0.
    do k=lfc,el-1
      ncape=ncape-0.5*(ht(k)+ht(k+1))*(z0(k+1)-z0(k))
    enddo
    
    fork=k2*alpha**2*pi**2*lmix/pr/sigma**2/(z0(el)-z0(lfc))/4.
    term1=vsr**2/2.
    term2=-1.-fork-2.*fork*ncape/vsr**2
    term3=(1.+fork+2.*fork*ncape/vsr**2)**2
    term4=8.*fork*(cape(i,j,tt)-fork*ncape)/vsr**2
    term5=4.*fork/vsr**2
    ecape1(i,j,tt)=term1+(term2+sqrt(term3+term4))/term5

    v2=vsr**2/2./cape(i,j,tt)
    nncape=ncape/cape(i,j,tt)
    duma=-1.-fork-fork*nncape/v2
    dumb=(1.+fork+fork*nncape/v2)**2
    dumc=4.*fork*(1.-fork*nncape)/v2
    dumd=2*fork/v2
    et=v2+(duma+sqrt(dumb+dumc))/dumd
    if (et .le. 0.) et=0.
    ecape2(i,j,tt)=et*cape(i,j,tt)
  endif

  744 continue
  deallocate(t0,q0,p0,z0,h0,hs0,h0h,ht,buoc)
  756 continue
enddo
enddo

write(*,*) t
enddo

call mpi_barrier(mpi_comm_world,ierr)

if (myid==0) then
  call mpi_gatherv(mpi_in_place,mpi_n(myid+1)*nx*ny,mpi_real,cape,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(mpi_in_place,mpi_n(myid+1)*nx*ny,mpi_real,ecape1,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(mpi_in_place,mpi_n(myid+1)*nx*ny,mpi_real,ecape2,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
else
  call mpi_gatherv(cape,mpi_n(myid+1)*nx*ny,mpi_real,cape,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(ecape1,mpi_n(myid+1)*nx*ny,mpi_real,ecape1,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(ecape1,mpi_n(myid+1)*nx*ny,mpi_real,ecape2,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
endif

call mpi_barrier(mpi_comm_world,ierr)


if (myid==0) then
  open(10,file=trim(workpath)//"/data/old_cape_land_"//trim(run)//".dat",access="direct",recl=nx*ny*3)
  do t=1,nt+1
    write(10,rec=t) cape(:,:,t),ecape1(:,:,t),ecape2(:,:,t)
  enddo
  close(10)
endif

deallocate(mpi_s,mpi_n)

call mpi_finalize(ierr)

end program cal_mean


! qv_sat from P3
 real function qv_sat(t_atm,p_atm,i_wrt)

!------------------------------------------------------------------------------------
! Calls polysvp1 to obtain the saturation vapor pressure, and then computes
! and returns the saturation mixing ratio, with respect to either liquid or ice,
! depending on value of 'i_wrt'
!------------------------------------------------------------------------------------

 implicit none

 !Calling parameters:
 real    :: t_atm  !temperature [K]
 real    :: p_atm  !pressure    [Pa]
 integer :: i_wrt  !index, 0 = w.r.t. liquid, 1 = w.r.t. ice

 !Local variables:
 real    :: e_pres         !saturation vapor pressure [Pa]
 real    :: polysvp1

  e_pres = polysvp1(t_atm,i_wrt)
  qv_sat = 0.622*e_pres/max(1.e-3,(p_atm-e_pres))

 return
 end function qv_sat

 real function polysvp1(T,i_type)

!-------------------------------------------
!  COMPUTE SATURATION VAPOR PRESSURE
!  POLYSVP1 RETURNED IN UNITS OF PA.
!  T IS INPUT IN UNITS OF K.
!  i_type REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)
!-------------------------------------------

      implicit none

      real    :: T
      integer :: i_type

! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

! ice
      real a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i
      data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
        6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/

! liquid
      real a0,a1,a2,a3,a4,a5,a6,a7,a8

! V1.7
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
        6.11239921, 0.443987641, 0.142986287e-1, &
        0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
        0.640689451e-10,-0.952447341e-13,-0.976195544e-15/
      real dt

!-------------------------------------------

      if (i_type.EQ.1 .and. T.lt.273.15) then
! ICE

! use Goff-Gratch for T < 195.8 K and Flatau et al. equal or above 195.8 K
         if (t.ge.195.8) then
            dt=t-273.15
            polysvp1 = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+a8i*dt)))))))
            polysvp1 = polysvp1*100.
         else
            polysvp1 = 10.**(-9.09718*(273.16/t-1.)-3.56654* &
                alog10(273.16/t)+0.876793*(1.-t/273.16)+ &
                alog10(6.1071))*100.
         end if

      elseif (i_type.EQ.0 .or. T.ge.273.15) then
! LIQUID

! use Goff-Gratch for T < 202.0 K and Flatau et al. equal or above 202.0 K
         if (t.ge.202.0) then
            dt = t-273.15
            polysvp1 = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
            polysvp1 = polysvp1*100.
         else
! note: uncertain below -70 C, but produces physical values (non-negative) unlike flatau
            polysvp1 = 10.**(-7.90298*(373.16/t-1.)+ &
                5.02808*alog10(373.16/t)- &
                1.3816e-7*(10**(11.344*(1.-t/373.16))-1.)+ &
                8.1328e-3*(10**(-3.49149*(373.16/t-1.))-1.)+ &
                alog10(1013.246))*100.
         end if

         endif

 end function polysvp1
! qv_sat from P3



