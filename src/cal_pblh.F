program cal_pblh
use netcdf
use mpi
use cal_cape
implicit none

integer, parameter :: nx=1024,ny=1024,nz=70,nt=144
real, parameter :: Lv=2.5e6
integer :: i,j,k,m,n,o,t
integer :: ierr, ncid1, varid1
integer :: myid,nproc,ts,te,tt
integer :: idum1,idum2
integer :: is,ie,js,je,ks,ke
integer, dimension(:), allocatable :: mpi_s,mpi_n
real :: dum1,dum2,dum3,dum4,dum5
real :: qv_sat
real, dimension(nx,ny,nz) :: th,thv,qv,the,w,dth
real, dimension(nz) :: thep,qvp
real, dimension(nz+1) :: zz,zc,rho,rhoz,thbar,thvbar,pbar,pibar,q1,q2
real, dimension(nz) :: dz
real, dimension(:), allocatable :: t0,q0,p0,z0, buoc
real*4, dimension(:,:,:), allocatable :: pblh1, pblh2, the_pblh1, the_pblh2, the_sur, cape, cin, wmax, rh5k
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
!run="tpe20110616nor"


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
  allocate(pblh1(nx,ny,nt+1),pblh2(nx,ny,nt+1))
  allocate(the_pblh1(nx,ny,nt+1),the_pblh2(nx,ny,nt+1),the_sur(nx,ny,nt+1) &
          ,cape(nx,ny,nt+1),cin(nx,ny,nt+1),wmax(nx,ny,nt+1),rh5k(nx,ny,nt+1))
else
  allocate(pblh1(nx,ny,mpi_n(myid+1)),pblh2(nx,ny,mpi_n(myid+1)))
  allocate(the_pblh1(nx,ny,mpi_n(myid+1)),the_pblh2(nx,ny,mpi_n(myid+1)),the_sur(nx,ny,mpi_n(myid+1)) &
          ,cape(nx,ny,mpi_n(myid+1)),cin(nx,ny,mpi_n(myid+1)),wmax(nx,ny,mpi_n(myid+1)),rh5k(nx,ny,mpi_n(myid+1)))
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


do k=1,nz
  do j=1,ny
    do i=1,nx 
      dum1=th(i,j,k)*pibar(k)
      the(i,j,k)=th(i,j,k)*exp(Lv*qv(i,j,k)/(1004.5*dum1))
      thv(i,j,k)=th(i,j,k)*(1.+0.608*qv(i,j,k))
     
      if (k==46) then 
        dum2=th(i,j,k)*pibar(k)
        dum1=qv_sat(dum2,pbar(k),0)
        rh5k(i,j,tt)=qv(i,j,k)/dum1 
      endif
    enddo
  enddo
enddo

call cpu_time(start)
do j=1,ny
  do i=1,nx
    n=1
    do k=topo(i,j)+1,nz-1
      dth(i,j,k)=(th(i,j,k+1)-th(i,j,k-1))/(zc(k+1)-zc(k-1))
      if (th(i,j,topo(i,j)+1)+0.5>=th(i,j,k)) n=k
      if (zc(k)-zz(topo(i,j))<=1000.) m=k 
    enddo
    !m=maxloc(dth(i,j,topo(i,j)+1:45),1)+topo(i,j)
    pblh1(i,j,tt)=zc(m)!-zz(topo(i,j))
    pblh2(i,j,tt)=zc(n)!-zz(topo(i,j))

    the_pblh1(i,j,tt)=maxval(the(i,j,topo(i,j)+1:m),1)
    the_pblh2(i,j,tt)=maxval(the(i,j,topo(i,j)+1:n),1)
    the_sur(i,j,tt)=the(i,j,topo(i,j)+1)

    !if (i==512 .and. j==512) write(*,*) i,j, pblh1(i,j,tt),pblh2(i,j,tt),zz(topo(i,j))

    !! CAPE
    n=nz-topo(i,j)
    allocate(t0(n),q0(n),p0(n),z0(n),buoc(n))
    t0=th(i,j,topo(i,j)+1:nz)*pibar(topo(i,j)+1:nz)
    q0=qv(i,j,topo(i,j)+1:nz)
    p0=pbar(topo(i,j)+1:nz)
    z0=zc(topo(i,j)+1:nz)
    call lift_parcel_adiabatic(nz-topo(i,j),t0,p0,q0,z0,buoc)
    do k=1,n-1
      if (buoc(k)>0.) then
        cape(i,j,tt)=cape(i,j,tt)+0.5*(buoc(k+1)+buoc(k))*(z0(k+1)-z0(k))
      elseif (buoc(k)<=0. .and. z0(k)-zz(topo(i,j))<3000.) then
        cin(i,j,tt)=cin(i,j,tt)+0.5*(buoc(k+1)+buoc(k))*(z0(k+1)-z0(k))
      endif
    enddo
    wmax(i,j,tt)=sqrt(2.*cape(i,j,tt))
    deallocate(t0,q0,p0,z0,buoc)
  enddo
enddo
call cpu_time(finish)

write(*,*) t,"time: ",(finish-start)
enddo

call mpi_barrier(mpi_comm_world,ierr)

if (myid==0) then
  call mpi_gatherv(mpi_in_place,mpi_n(myid+1)*nx*ny,mpi_real,pblh1,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(mpi_in_place,mpi_n(myid+1)*nx*ny,mpi_real,pblh2,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(mpi_in_place,mpi_n(myid+1)*nx*ny,mpi_real,the_pblh1,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(mpi_in_place,mpi_n(myid+1)*nx*ny,mpi_real,the_pblh2,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(mpi_in_place,mpi_n(myid+1)*nx*ny,mpi_real,the_sur,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(mpi_in_place,mpi_n(myid+1)*nx*ny,mpi_real,cape,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(mpi_in_place,mpi_n(myid+1)*nx*ny,mpi_real,wmax,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(mpi_in_place,mpi_n(myid+1)*nx*ny,mpi_real,rh5k,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
else
  call mpi_gatherv(pblh1,mpi_n(myid+1)*nx*ny,mpi_real,pblh1,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(pblh2,mpi_n(myid+1)*nx*ny,mpi_real,pblh2,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(the_pblh1,mpi_n(myid+1)*nx*ny,mpi_real,the_pblh1,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(the_pblh2,mpi_n(myid+1)*nx*ny,mpi_real,the_pblh2,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(the_sur,mpi_n(myid+1)*nx*ny,mpi_real,the_sur,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(cape,mpi_n(myid+1)*nx*ny,mpi_real,cape,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(wmax,mpi_n(myid+1)*nx*ny,mpi_real,wmax,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(rh5k,mpi_n(myid+1)*nx*ny,mpi_real,rh5k,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
endif

call mpi_barrier(mpi_comm_world,ierr)


if (myid==0) then
  open(10,file=trim(workpath)//"/data/pblh_"//trim(run)//".dat",access="direct",recl=nx*ny*9)
  do t=1,nt+1
    write(10,rec=t) pblh1(:,:,t),pblh2(:,:,t),the_pblh1(:,:,t),the_pblh2(:,:,t),the_sur(:,:,t)&
                   ,cape(:,:,t),cin(:,:,t),wmax(:,:,t),rh5k(:,:,t)
  enddo
  close(10)
endif

deallocate(mpi_s,mpi_n)

call mpi_finalize(ierr)

end program cal_pblh


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
 real    :: polysvp1

 !Local variables:
 real    :: e_pres         !saturation vapor pressure [Pa]

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


