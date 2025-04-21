program main
use share
use netcdf
IMPLICIT NONE
INTEGER                           :: time_num, time_bin
INTEGER                           :: ierr, cpun, cpuid, err, ncid,&
                                     varid, i ,j, m, t, inv, k, k1,&
                                     x_start,x_end,y_start,y_end, loc0, loc1, ik
INTEGER                           :: yyyy, mm, dd, hh, dd_end, lfc, el
REAL                              :: lon(NX), lat_yinv(NY), lat(NY)
REAL*8                            :: sf, ao, var_miss, dum, &
                                     fork, term1, term2, term3, term4, term5,&
                                     term6, term7,&
                                     duma,dumb,dumc,dumd,et,nncape,v2, ent, dume
INTEGER*2,DIMENSION(NX,NY)        :: var_yinv
INTEGER*2,DIMENSION(NX,NY,NZ)     :: var3d_yinv
REAL*8,DIMENSION(NX,NY)           :: topo,ps,u10,v10,t2m,d2m
REAL*8,DIMENSION(NX,NY,NZ)        :: tl,ql,phi,u,v
REAL*8,DIMENSION(NZ+1)            :: t0,qv0,z0,p0,h0,buoc,h0s,h0h,ht
REAL,  DIMENSION(NX,NY)           :: ecape_a,cape,vsr,ncape,ecape_a1
CHARACTER                         :: time_ini*10, time_now*10, out_dtg*8, &
                                     inpath*200, outpath*200

  print *,'========================================'
  print *,' Program starts'
  print *,'========================================'
  namelist/cntrl/time_ini,time_num,time_bin,inpath,outpath

! initialize MPI
! call mpi_init(ierr)
! call mpi_comm_size(mpi_comm_world, cpun, ierr)
! call mpi_comm_rank(mpi_comm_world, cpuid, ierr)

  open(11,file='namelist',status='old',iostat=err)
   if (err .ne. 0) then
    print *,'namelist not found or problem occurs, program stops'
    stop
   endif

  read(11,nml=cntrl,iostat=err)
   if (err .ne. 0) then
    print *,'domain in namelist is wrong, program stops'
    stop
   endif
  close(11)

! read lat lon from any .nc
  err = nf90_open(trim(inpath)//'2001/era5_sfc_025d_2001010100.nc',nf90_nowrite,ncid)
  if (err/=nf90_noerr) print *,'open input data fail'
  err = nf90_inq_varid(ncid,'longitude',varid)
  if (err/=nf90_noerr) print *,'q variable problem:','longitude'
  err = nf90_get_var(ncid,varid,lon,start=(/ 1 /),count=(/NX/))
  if (err/=nf90_noerr) print *,'read variable problem:','longitude'
  err = nf90_inq_varid(ncid,'latitude',varid)
  if (err/=nf90_noerr) print *,'q variable problem:','latitude'
  err = nf90_get_var(ncid,varid,lat_yinv,start=(/ 1 /),count=(/NY/))
  if (err/=nf90_noerr) print *,'read variable problem:','latitude'
  do j=1,NY
   inv=NY-j+1
   lat(j)=lat_yinv(inv)
  enddo
  err = nf90_close(ncid)

! set analysis domain x-y indices
  x_start=minloc(abs(lon-230.),dim=1)
  x_end=minloc(abs(lon-295.),dim=1)
  y_start=minloc(abs(lat-25.),dim=1)
  y_end=minloc(abs(lat-49.),dim=1)
  print *,' analysis domain:'
  print *,'  longitude ', lon(x_start), lon(x_end),' x ',x_start,x_end
  print *,'   latitude ', lat(y_start), lat(y_end),' y ',y_start,y_end

!==============================================================
  read(time_ini(1:4),*)  yyyy
  read(time_ini(5:6),*)  mm
  read(time_ini(7:8),*)  dd
  read(time_ini(9:10),*) hh
!==============================================================
!  time loop
   do t = 1, time_num
    if (t .gt. 1) hh=hh+time_bin
    dd_end=dlist(mm)
    if ((mod(yyyy,4) .eq. 0) .and. (mm .eq. 2)) dd_end=dd_end+1
    if (hh .gt. 23) then
     dd=dd+1
     hh=0
     if (dd .gt. dd_end) then
      dd=1
      mm=mm+1
      if (mm .gt. 12) then
       mm=1
       yyyy=yyyy+1
      endif
     endif
    endif
   write (time_now(1:4),'(I4.4)')  yyyy
   write (time_now(5:6),'(I2.2)')  mm
   write (time_now(7:8),'(I2.2)')  dd
   write (time_now(9:10),'(I2.2)') hh
   out_dtg=time_now(3:10)
!   print *,'========================================'
    print *,' time  : ',time_now
!   print *,'========================================'


!  print *,' reading input...'
!  read surface data

   err = nf90_open(trim(inpath)//time_now(1:4)//"/era5_sfc_025d_"//time_now//".nc",nf90_nowrite,ncid)
   if (err/=nf90_noerr) print *,'open input data (surface) fail'

   do m = 1, sfcvarn
   err = nf90_inq_varid(ncid,trim(sfcvar(m)),varid)
   if (err/=nf90_noerr) print *,'q variable problem:',trim(sfcvar(m))
   err = nf90_get_att(ncid, varid, "scale_factor", sf)
   if (err /= nf90_noerr) print *,'scale factor wrong:',trim(sfcvar(m))
   err = nf90_get_att(ncid, varid, "add_offset", ao)
   if (err /= nf90_noerr) print *,'offset wrong:',trim(sfcvar(m))
   err = nf90_get_var(ncid,varid,var_yinv,start=(/ 1,1,1/),count=(/NX,NY,1/))
   if (err/=nf90_noerr) print *,'read variable problem:',trim(sfcvar(m))

   if (sfcvar(m) .eq. 'z') then
   do j=1,NY
    inv=NY-j+1
    topo(:,j)=var_yinv(:,inv)*sf+ao!surface geopotential (m2/s2)
   enddo
   endif
   if (sfcvar(m) .eq. 'sp') then
   do j=1,NY
    inv=NY-j+1
    ps(:,j)=var_yinv(:,inv)*sf+ao!surface pressure (Pa)
   enddo
   endif
   if (sfcvar(m) .eq. 'u10') then
   do j=1,NY
    inv=NY-j+1
    u10(:,j)=var_yinv(:,inv)*sf+ao!10 m u component wind speed
   enddo
   endif
   if (sfcvar(m) .eq. 'v10') then
   do j=1,NY
    inv=NY-j+1
    v10(:,j)=var_yinv(:,inv)*sf+ao!10 m v component wind speed
   enddo
   endif
   if (sfcvar(m) .eq. 't2m') then
   do j=1,NY
    inv=NY-j+1
    t2m(:,j)=var_yinv(:,inv)*sf+ao!2 m temperature K
   enddo
   endif
   if (sfcvar(m) .eq. 'd2m') then
   do j=1,NY
    inv=NY-j+1
    d2m(:,j)=var_yinv(:,inv)*sf+ao!2 m dew point temperature K
   enddo
   endif
   
   enddo!m

   err = nf90_close(ncid)

!  read profile at pressure level
!  vertical index starts from top
   err = nf90_open(trim(inpath)//time_now(1:4)//"/era5_tqzw_025d_"//time_now//".nc",nf90_nowrite,ncid)
   if (err/=nf90_noerr) print *,'open input data (pressure essentials) fail'

   do m = 1, pevarn
!  print *,'reading field - ',pevar(m)
   err = nf90_inq_varid(ncid,trim(pevar(m)),varid)
   if (err/=nf90_noerr) print *,'q variable problem:',trim(pevar(m))
   err = nf90_get_att(ncid, varid, "scale_factor", sf)
   if (err /= nf90_noerr) print *,'scale factor wrong:',trim(pevar(m))
   err = nf90_get_att(ncid, varid, "add_offset", ao)
   if (err /= nf90_noerr) print *,'offset wrong:',trim(pevar(m))
   err = nf90_get_var(ncid,varid,var3d_yinv,start=(/ 1,1,1,1/),count=(/NX,NY,NZ,1/))
   if (err/=nf90_noerr) print *,'read variable problem:',trim(pevar(m))

   if (pevar(m) .eq. 't') then
    do k=1,NZ
     do j=1,NY
      inv=NY-j+1
      tl(:,j,k)=var3d_yinv(:,inv,k)*sf+ao!K
     enddo
    enddo
   endif
   if (pevar(m) .eq. 'q') then
    do k=1,NZ
     do j=1,NY
      inv=NY-j+1
      ql(:,j,k)=var3d_yinv(:,inv,k)*sf+ao!kg/kg
     enddo
    enddo
   endif
   if (pevar(m) .eq. 'z') then
    do k=1,NZ
     do j=1,NY
      inv=NY-j+1
      phi(:,j,k)=var3d_yinv(:,inv,k)*sf+ao!m2/s2
     enddo
    enddo
   endif
   if (pevar(m) .eq. 'u') then
    do k=1,NZ
     do j=1,NY
      inv=NY-j+1
      u(:,j,k)=var3d_yinv(:,inv,k)*sf+ao!m2/s2
     enddo
    enddo
   endif
   if (pevar(m) .eq. 'v') then
    do k=1,NZ
     do j=1,NY
      inv=NY-j+1
      v(:,j,k)=var3d_yinv(:,inv,k)*sf+ao!m2/s2
     enddo
    enddo
   endif

   enddo!m

   err = nf90_close(ncid)

!  print *,' start x-y loop...'
  
   cape=0.
   ecape_a=0.
   ecape_a1=0.
   vsr=0.
   ncape=0.
   do j = y_start, y_end
    do i = x_start, x_end
     vsr(i,j)=sqrt((u(i,j,17)-u10(i,j))**2+(v(i,j,17)-v10(i,j))**2)/2.!500hPa-10m
!
     t0=0.
     z0=0.
     qv0=0.
     p0=0.
     h0=0.
     h0s=0.
     buoc=0.
!    loc0=z-index(top->bottom) closiest to topo
     loc0=NZ
     do k = 1, NZ
      if (topo(i,j) .ge. phi(i,j,k)) then
       loc0=k-1
       exit
      endif
     enddo 
     
!    var(1) is surface
     t0(1)=t2m(i,j)
     z0(1)=0.
     p0(1)=ps(i,j)
     qv0(1)=qsat(d2m(i,j),p0(1),1)
     h0(1)=cp*t0(1)+topo(i,j)+xlv*qv0(1)
     h0s(1)=cp*t0(1)+topo(i,j)+xlv*qsat(t2m(i,j),p0(1),0)
  
!    change vertical indexing to bottom->top
!    atmosphere starts from var(2) (z-index closiest to the topo)
     ik=2
     do k = loc0, 1, -1
     t0(ik)=tl(i,j,k)
     z0(ik)=(phi(i,j,k)-topo(i,j))/g
     qv0(ik)=ql(i,j,k)
     p0(ik)=100.*pl(k)
     h0(ik)=cp*t0(ik)+phi(i,j,k)+xlv*qv0(ik)
     h0s(ik)=cp*t0(ik)+phi(i,j,k)+xlv*qsat(t0(ik),p0(ik),0)
     ik=ik+1
     enddo
!    z dimension is 1:ik-1, ik-1 depends on topo index

!    loc1=level lower than 5km, including surface, with highest mse
     loc1=1
     dum=h0(loc1)
     do k = 2, ik-1
      if (z0(k) .gt. 5000.) exit
      if (h0(k) .gt. dum) then
       loc1=k
       dum=h0(k)
      endif
     enddo

!    z dimension is loc1:ik-1
     call lift_parcel_adiabatic(ik-1-loc1+1,t0(loc1:ik-1),p0(loc1:ik-1),qv0(loc1:ik-1),&
                                z0(loc1:ik-1),buoc(loc1:ik-1))

     lfc=loc1+1
     do k = loc1+1, ik-1
     if (buoc(k) .gt. 0.) then
      lfc=k
      exit
     endif
     enddo

     el=ik-1
     do k = ik-1, lfc, -1
     if (buoc(k) .gt. 0.) then
      el=k
      exit
     endif
     enddo
!    print *,lfc,el,z0(lfc),z0(el)

     cape(i,j)=0.
     do k = lfc, el-1
      if ((buoc(k) .gt. 0.) .and. (buoc(k+1) .gt. 0.)) then
      cape(i,j)=cape(i,j)+(buoc(k)+buoc(k+1))*(z0(k+1)-z0(k))/2.
      endif
     enddo

 
     if (cape(i,j) .le. 0.) cycle

     h0h=0.
     do k = lfc, el
       do k1 = 1, k-1
        h0h(k)=h0h(k)+(h0(k1)+h0(k1+1))*(z0(k1+1)-z0(k1))/2.
       enddo
       h0h(k)=h0h(k)/z0(k)
     ht(k)=g/cp*(h0h(k)-h0s(k))/t0(k)
     enddo

     ncape(i,j)=0.
     do k = lfc, el-1
     ncape(i,j)=ncape(i,j)+(ht(k)+ht(k+1))*(z0(k+1)-z0(k))/2.
     enddo
     ncape(i,j)=-1.*ncape(i,j)
!
     fork=k2*alpha**2*pi**2*lmix/pr/sigma**2/(z0(el)-z0(lfc))/4.
     term1=vsr(i,j)**2/2.
     term2=-1.-fork-2.*fork*ncape(i,j)/vsr(i,j)**2
     term3=(1.+fork+2.*fork*ncape(i,j)/vsr(i,j)**2)**2
     term4=8.*fork*(cape(i,j)-fork*ncape(i,j))/vsr(i,j)**2
     term5=4.*fork/vsr(i,j)**2
     ecape_a1(i,j)=term1+(term2+sqrt(term3+term4))/term5

     v2=vsr(i,j)**2/2./cape(i,j)
     nncape=ncape(i,j)/cape(i,j)
     duma=-1.-fork-fork*nncape/v2
     dumb=(1.+fork+fork*nncape/v2)**2
     dumc=4.*fork*(1.-fork*nncape)/v2
     dumd=2*fork/v2
     et=v2+(duma+sqrt(dumb+dumc))/dumd
     if (et .le. 0.) et=0.
     ecape_a(i,j)=et*cape(i,j)


!    print *,cape(i,j),ncape,ecape_a(i,j)
!    print *,ncape-vsr**2/(2.*fork)+sqrt(vsr**4/(4.*fork)*(1.-fork*ncape/vsr**2)**2+cape(i,j)*vsr**2/fork)
!    print *,z0(el)-z0(lfc)
!    print *,cape(i,j),vsr,ncape,ecape_a(i,j)
!    print *,et,ecape_a(i,j),nncape
!

!     ent=2.*k2*lmix/pr/2.25e6
!     dume=cape(i,j)*(1.-dexp(-(z0(el)-z0(lfc))*ent))/(z0(el)-z0(lfc))/ent+&
!          ncape(i,j)*(1.-(1.-dexp(-(z0(el)-z0(lfc))*ent))/(z0(el)-z0(lfc))/ent)
!     dume=max(dume,0.)
!     ecape_a(i,j)=max(ecape_a(i,j),dume)
!    print *,cape(i,j),ecape_a(i,j),ecape_a1(i,j)

!    endif!positive cape

    enddo
   enddo

 
   open(unit=700,file=trim(outpath)//"ecapef_"//time_now//".dat",form='unformatted',access='direct',&
        recl=(x_end-x_start+1)*(y_end-y_start+1)*5)
    write(700,rec=1) cape(x_start:x_end,y_start:y_end),&
                     ecape_a(x_start:x_end,y_start:y_end),&
                     ecape_a1(x_start:x_end,y_start:y_end),&
                     ncape(x_start:x_end,y_start:y_end),&
                     vsr(x_start:x_end,y_start:y_end)
   close(700)
   
   
   enddo!time loop

! call mpi_finalize(ierr)
END program main
