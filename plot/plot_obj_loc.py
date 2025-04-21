import os
import numpy as np
import netCDF4 as nc
import metpy.calc as mc
import matplotlib.pyplot as plt

nx=1024
ny=1024
nz=70

ts=48  # 8 LST
te=54  # 9 LST

path='/data/der0318/parcel_model_taiwanvvm/'
caselist=os.listdir(path+'/taiwanVVM/')
#caselist=['tpe20050702nor','tpe20050712nor','tpe20050723nor','tpe20050826nor','tpe20060508nor']
ncase=len(caselist)
print(ncase)

fname1=path+'/taiwanVVM/'+caselist[0]+'/TOPO.nc'
fid=nc.Dataset(fname1,'r')
height=fid.variables['height'][:]

vvm_w=np.zeros([10000,ncase])
vvm_w_x=np.zeros([10000,ncase])
vvm_w_y=np.zeros([10000,ncase])
vvm_w_z=np.zeros([10000,ncase])
n=0
for casename in caselist:
  fname1=path+'/data/old_cape_land_'+casename+'.dat'
  fname2=path+'/obj_wmax_data/'+casename+'.dat'
  
  w_vvm=np.fromfile(fname2,dtype=np.float32)[1:].reshape(-1,7)
  nw=w_vvm.shape[0]
  vvm_w[0:nw,n]=w_vvm[:,0]
  vvm_w_x[0:nw,n]=w_vvm[:,-2] 
  vvm_w_y[0:nw,n]=w_vvm[:,-1] 
  vvm_w_z[0:nw,n]=w_vvm[:,-3] 
 
  w_ecape=np.fromfile(fname1,dtype=np.float32,count=(te-ts)*3*ny*nx,offset=ts*3*ny*nx*4).reshape(-1,3,ny,nx)[:,2,:,:].flatten()
  w_cape=np.compress(w_ecape>0.,w_ecape)
 
  w_cape90=np.percentile(w_cape,99)
  w_capemax=w_cape.max()
  #print(w_cape90,w_capemax)
  
  w_vvm90=np.percentile(vvm_w[0:nw,n],99)
  vvm_wx=vvm_w_x[vvm_w[:,n]>w_vvm90,n]
  vvm_wy=vvm_w_y[vvm_w[:,n]>w_vvm90,n]

  w_ecape=np.fromfile(fname1,dtype=np.float32,count=(te-ts)*3*ny*nx,offset=ts*3*ny*nx*4).reshape(-1,3,ny,nx)[:,2,:,:]
 
  wecape_loc=np.argmax(w_ecape)
  dnt=np.int(wecape_loc/nx/ny)
  dny=np.int((wecape_loc-dnt*nx*ny)/nx)
  dnx=wecape_loc-dnt*nx*ny-dny*nx

  w_ecape=np.where(w_ecape>0.,w_ecape,np.nan)

  loc=np.sum(np.where(w_ecape>w_cape90,1.,0.),0)
  loc=np.where(loc>0.5,1.,0.)
 
  fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(6,6),dpi=300)
  ax.pcolormesh(loc,vmax=2.0,vmin=0.5,cmap='Reds')
  #im=ax.pcolormesh(w_ecape[0,:,:],vmax=1000,vmin=0,cmap='Blues')
  ax.contour(height,levels=[0.05,0.5],colors=['k','k'],linewidths=0.5)
  ax.plot(vvm_wx,vvm_wy,'bo',markersize=1,alpha=0.7)
  ax.plot(dnx,dny,'c*',markersize=4,alpha=0.7)

  #cbar=fig.colorbar(im)
  #cbar.set_ticks(np.arange(0,1000+.1,200.))

  ax.set_title(casename)
  #plt.show()
  plt.savefig("./figure/obj_"+casename+".png")
  plt.close()

  n=n+1
