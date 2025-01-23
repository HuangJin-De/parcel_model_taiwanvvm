program cal_mean
use netcdf
use mpi
implicit none

integer, parameter :: nx=1024,ny=1024,nz=70,nt=144
real, parameter :: Lv=2.5e6
integer :: i,j,k,m,n,o,t
integer :: ierr, ncid1, varid1
integer :: myid,nproc,ts,te,tt
integer :: idum1,idum2
integer :: is,ie,js,je,ks,ke
integer, dimension(:), allocatable :: mpi_s,mpi_n
real :: qv_sat
real :: dum1,dum2,dum3,dum4,dum5
real, dimension(nx,ny,nz) :: th,qv,u,v
real, dimension(nz) :: thm,qvm,um,vm,np
real*4, dimension(:,:), allocatable :: tm,tdm,wsm,wdm,pm
real, dimension(nz+1) :: zz,zc,rho,rhoz,thbar,thvbar,pbar,pibar,q1,q2
real, dimension(nz) :: dz
real, dimension(:), allocatable :: t0,q0,p0,z0, buoc
real :: start, finish
integer, dimension(nx,ny) :: topo
character(200) :: path,run,filename,workpath

call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world,nproc,ierr)
call mpi_comm_rank(mpi_comm_world,myid,ierr)

allocate(mpi_s(nproc),mpi_n(nproc))

path="/data/der0318/taiwanVVM/"
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
  allocate(tm(nz,nt+1),tdm(nz,nt+1),wsm(nz,nt+1),wdm(nz,nt+1),pm(nz,nt+1))
else
  allocate(tm(nz,mpi_n(myid+1)),tdm(nz,mpi_n(myid+1)),wsm(nz,mpi_n(myid+1)),wdm(nz,mpi_n(myid+1)),pm(nz,mpi_n(myid+1)))
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
ierr = nf90_close(ncid1)

where(topo<1) topo=1 

thm=0.
qvm=0.
um=0.
vm=0.
np=0.
do j=1,ny
do i=1,nx
  do k=topo(i,j)+1,nz
    thm(k)=thm(k)+th(i,j,k)
    qvm(k)=qvm(k)+qv(i,j,k)
    um(k)=um(k)+u(i,j,k)
    vm(k)=vm(k)+v(i,j,k)
    np(k)=np(k)+1.
  enddo
enddo
enddo
thm=thm/np
qvm=qvm/np
um=um/np
vm=vm/np


do k=1,nz
  tm(k,tt)=thm(k)*pibar(k)-273.15
  
  dum1 = 100.
  dum2 = thm(k)*pibar(k)
  dum3 = qvm(k)
  do while (dum1>=5e-2)
    i = 0
    if (dum2<273.15) i = 1
    dum4 = qv_sat(dum2,pbar(k),i) - dum3
    !dum1 = abs(dum4) / dum3
    dum1 = 0.7 * dum4 * 2.5e6 / 1004.5 
    dum2 = dum2 - dum1
    !write(*,*) k, dum1,dum2
  enddo
  tdm(k,tt) = dum2 - 273.15
  pm(k,tt) = pbar(k) / 100. 

  wsm(k,tt) = sqrt(um(k)**2. + vm(k)**2.)
  dum1 = atan2(um(k),vm(k))/3.14159265*180.
  wdm(k,tt) = 270. - dum1 
  if (wdm(k,tt)<0.) wdm(k,tt) = wdm(k,tt) + 360.
  if (wdm(k,tt)>=360.) wdm(k,tt) = wdm(k,tt) - 360.
enddo

write(*,*) t
enddo

call mpi_barrier(mpi_comm_world,ierr)

if (myid==0) then
  call mpi_gatherv(mpi_in_place,mpi_n(myid+1)*nz,mpi_real,tm,mpi_n*nz,mpi_s*nz,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(mpi_in_place,mpi_n(myid+1)*nz,mpi_real,tdm,mpi_n*nz,mpi_s*nz,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(mpi_in_place,mpi_n(myid+1)*nz,mpi_real,wsm,mpi_n*nz,mpi_s*nz,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(mpi_in_place,mpi_n(myid+1)*nz,mpi_real,wdm,mpi_n*nz,mpi_s*nz,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(mpi_in_place,mpi_n(myid+1)*nz,mpi_real,pm,mpi_n*nz,mpi_s*nz,mpi_real,0,mpi_comm_world,ierr)
else
  call mpi_gatherv(tm,mpi_n(myid+1)*nz,mpi_real,tm,mpi_n*nz,mpi_s*nz,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(tdm,mpi_n(myid+1)*nz,mpi_real,tdm,mpi_n*nz,mpi_s*nz,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(wsm,mpi_n(myid+1)*nz,mpi_real,wsm,mpi_n*nz,mpi_s*nz,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(wdm,mpi_n(myid+1)*nz,mpi_real,wdm,mpi_n*nz,mpi_s*nz,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(pm,mpi_n(myid+1)*nz,mpi_real,pm,mpi_n*nz,mpi_s*nz,mpi_real,0,mpi_comm_world,ierr)
endif

call mpi_barrier(mpi_comm_world,ierr)


if (myid==0) then
  open(10,file=trim(workpath)//"/data/mean_"//trim(run)//".dat",access="direct",recl=nz*5)
  do t=1,nt+1
    write(10,rec=t) tm(:,t),tdm(:,t),wsm(:,t),wdm(:,t),pm(:,t)
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



