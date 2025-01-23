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
real, dimension(nx,ny,nz) :: th,thv,qv,the,w,dth
real, dimension(nz) :: thep,qvp
real, dimension(nz+1) :: zz,zc,rho,rhoz,thbar,thvbar,pbar,pibar,q1,q2
real, dimension(nz) :: dz
real, dimension(:), allocatable :: t0,q0,p0,z0, buoc
real*4, dimension(:,:,:), allocatable :: pblh1, pblh2, the_pblh1, the_pblh2, the_sur, cape, cin, wmax
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
          ,cape(nx,ny,nt+1),cin(nx,ny,nt+1),wmax(nx,ny,nt+1))
else
  allocate(pblh1(nx,ny,mpi_n(myid+1)),pblh2(nx,ny,mpi_n(myid+1)))
  allocate(the_pblh1(nx,ny,mpi_n(myid+1)),the_pblh2(nx,ny,mpi_n(myid+1)),the_sur(nx,ny,mpi_n(myid+1)) &
          ,cape(nx,ny,mpi_n(myid+1)),cin(nx,ny,mpi_n(myid+1)),wmax(nx,ny,mpi_n(myid+1)))
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
    enddo
    m=maxloc(dth(i,j,topo(i,j)+1:45),1)+topo(i,j)
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
      elseif (buoc(k)<=0. .and. z0(k)<3000.) then
        cin(i,j,tt)=cin(i,j,tt)+0.5*(buoc(k+1)+buoc(k))*(z0(k+1)-z0(k))
      endif
    enddo
    wmax(i,j,tt)=sqrt(2.*cape(i,j,tt))
    deallocate(t0,q0,p0,z0,buoc)
  enddo
enddo
call cpu_time(finish)

write(*,*) "time: ",(finish-start)

write(*,*) t
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
else
  call mpi_gatherv(pblh1,mpi_n(myid+1)*nx*ny,mpi_real,pblh1,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(pblh2,mpi_n(myid+1)*nx*ny,mpi_real,pblh2,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(the_pblh1,mpi_n(myid+1)*nx*ny,mpi_real,the_pblh1,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(the_pblh2,mpi_n(myid+1)*nx*ny,mpi_real,the_pblh2,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(the_sur,mpi_n(myid+1)*nx*ny,mpi_real,the_sur,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(cape,mpi_n(myid+1)*nx*ny,mpi_real,cape,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
  call mpi_gatherv(wmax,mpi_n(myid+1)*nx*ny,mpi_real,wmax,mpi_n*nx*ny,mpi_s*nx*ny,mpi_real,0,mpi_comm_world,ierr)
endif

call mpi_barrier(mpi_comm_world,ierr)


if (myid==0) then
  open(10,file=trim(workpath)//"/data/pblh_"//trim(run)//".dat",access="direct",recl=nx*ny*8)
  do t=1,nt+1
    write(10,rec=t) pblh1(:,:,t),pblh2(:,:,t),the_pblh1(:,:,t),the_pblh2(:,:,t),the_sur(:,:,t)&
                   ,cape(:,:,t),cin(:,:,t),wmax(:,:,t)
  enddo
  close(10)
endif

deallocate(mpi_s,mpi_n)

call mpi_finalize(ierr)

end program cal_pblh

