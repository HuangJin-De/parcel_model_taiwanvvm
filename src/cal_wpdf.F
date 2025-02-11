program cal_wpdf
use netcdf
implicit none

integer, parameter :: nx=1024,ny=1024,nz=70,nt=144
integer :: i,j,k,m,n,o,t
integer :: idum1,idum2
integer :: is,ie,js,je,ks,ke
integer :: ierr,ncid1,varid1
real :: dum1,dum2,dum3,dum4,dum5
real :: start, finish
real, dimension(nz+1) :: zz,zc,rho,rhoz,thbar,thvbar,pbar,pibar,q1,q2
real, dimension(nx,ny) :: pbl14, pbl24, the14, the24, cape4, cin4, wmax4
real*4, dimension(nx,ny,nt+1) :: pblh1,pblh2,the_pblh1,the_pblh2,the_sur,cape,cin,wmax,rh5k
integer, dimension(nx,ny) :: lu
character(200) :: path,run,filename,workpath

integer, parameter :: ncri=100
real, dimension(ncri+1) :: cri
real*4, dimension(ncri,nt+1) :: wcon


path="/data/der0318/taiwanVVM/"
workpath="/data/der0318/parcel_model_taiwanvvm/"
run="runname"
!run="tpe20050702nor"


open(10,file=trim(workpath)//"/data/pblh_"//trim(run)//".dat",access="direct",recl=nx*ny*9)
do t=1,nt+1
  read(10,rec=t) pblh1(:,:,t),pblh2(:,:,t),the_pblh1(:,:,t),the_pblh2(:,:,t),the_sur(:,:,t)&
                ,cape(:,:,t),cin(:,:,t),wmax(:,:,t),rh5k(:,:,t)
enddo
close(10)

WRITE(filename,113) trim(path), "/", trim(run),"/TOPO.nc"
113 FORMAT(4A)

ierr = nf90_open(trim(filename),NF90_NOWRITE,ncid1)
IF (ierr/=nf90_noerr) WRITE(*,*) "open fail"
ierr = nf90_inq_varid(ncid1,"lu",varid1)
IF (ierr/=nf90_noerr) WRITE(*,*) "var_inq fail"
ierr = nf90_get_var(ncid1,varid1,lu,start=(/ 1,1,1,1 /),count=(/ nx,ny,1,1 /))
IF (ierr/=nf90_noerr) WRITE(*,*) "read fail"
ierr = nf90_close(ncid1)

wcon=0.
cri=0.
do i=1,ncri 
  cri(i)=real(i-1)
enddo
cri(ncri+1)=200.


do t=1,nt+1
  where(lu==16) wmax(:,:,t)=-20.
  do j=1,ny
  do i=1,nx
    dum1=wmax(i,j,t)
    do n=1,ncri
      if (dum1<cri(n)) goto 999
      if (cri(n)<=dum1 .and. dum1<cri(n+1)) then
        wcon(n,t)=wcon(n,t)+1.
        goto 999
      endif
    enddo
    999 continue
  enddo
  enddo
  !write(*,*) t,"here",sum(wcon(:,t),1)
enddo


open(10,file=trim(workpath)//"/data/wpdf_"//trim(run)//".dat",access="direct",recl=ncri)
do t=1,nt+1
  write(10,rec=t) wcon(:,t)
enddo
close(10)

end program cal_wpdf
