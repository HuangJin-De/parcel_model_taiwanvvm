import os
import numpy as np
import metpy.calc as mc
import matplotlib.pyplot as plt


nx=1024
ny=1024
nz=70

ts=48  # 8 LST
te=54  # 9 LST

path='/data/der0318/parcel_model_taiwanvvm/'
caselist=os.listdir(path+'/taiwanVVM/')
#caselist=['tpe20050702nor']
#caselist=['tpe20050702nor','tpe20050712nor','tpe20050723nor','tpe20050826nor','tpe20060508nor']
ncase=len(caselist)
print(ncase)

w_vvm_m=np.zeros(ncase)
w_cape_m=np.zeros(ncase)
w_vvm_r=np.zeros([2,ncase])
w_cape_r=np.zeros([2,ncase])
vvm_w=np.zeros([10000,ncase])
vvm_w_time=np.zeros([10000,ncase])
vvm_w_size=np.zeros([10000,ncase])
vvm_w_prec=np.zeros([10000,ncase])
cape=np.zeros(ncase)
dcape_max=np.zeros(ncase)
ecape=np.zeros(ncase)
rh_mid=np.zeros(ncase)
qv_sur=np.zeros(ncase)
hm_sur=np.zeros(ncase)
ws_mid=np.zeros(ncase)
wd_mid=np.zeros(ncase)
n=0
for casename in caselist:
  fname1=path+'/obj_wmax_data/'+casename+'.dat'
  fname2=path+'/data/old_cape_land_'+casename+'.dat'
  fname3=path+'/data/mean_'+casename+'.dat'

  w_vvm=np.fromfile(fname1,dtype=np.float32)[1:].reshape(-1,4)
  #w_vvm=w_vvm[w_vvm[:,3]>20.,:]
  nw=w_vvm.shape[0]
  vvm_w[0:nw,n]=w_vvm[:,0]
  vvm_w_time[0:nw,n]=w_vvm[:,1]/6.
  vvm_w_size[0:nw,n]=w_vvm[:,2]
  vvm_w_prec[0:nw,n]=w_vvm[:,3]
  w_vvm=w_vvm[:,0]
  w_ecape=np.fromfile(fname2,dtype=np.float32,count=(te-ts)*3*ny*nx,offset=ts*3*ny*nx*4).reshape(-1,3,ny,nx)[:,2,:,:].flatten()
  w_ocape=np.fromfile(fname2,dtype=np.float32,count=(te-ts)*3*ny*nx,offset=ts*3*ny*nx*4).reshape(-1,3,ny,nx)[:,0,:,:].flatten()

  dcape=w_ecape/w_ocape

  dcape=np.compress(w_ecape>0.,dcape)
  w_cape=np.compress(w_ecape>0.,w_ecape)
  w_cape=np.sqrt(2.*w_cape)

  #print(w_vvm.shape,w_cape.shape)
 
  w_vvm_m[n]=w_vvm.max() #w_vvm.mean()
  w_vvm_r[:,n]=np.abs([np.percentile(w_vvm,25),np.percentile(w_vvm,75)]-w_vvm_m[n])
  #w_vvm_r[:,n]=[np.percentile(w_vvm,25),np.percentile(w_vvm,75)]

  dcape_max[n]=dcape[w_cape.argmax()]
  w_cape_m[n]=w_cape.mean()
  #w_cape_r[:,n]=np.abs([np.percentile(w_cape,25),np.percentile(w_cape,75)]-w_cape_m[n])
  w_cape_r[:,n]=0. #np.abs([w_cape.max()-w_cape.min(),0.])
  #w_cape_r[:,n]=[np.percentile(w_cape,25),np.percentile(w_cape,75)]

  ## env vars
  #env=np.fromfile(fname3,dtype=np.float32,count=(nz*5),offset=ts*(nz*5)*4)
  #pro=env.reshape(5,nz)
  #rh_mid[n]=mc.temperature_from_potential_temperature(pro[4,46]*units('K'),,pro[0,46]*units(''),)
  #rh_mid[n]=mc.relative_humidity_from_mixing_ratio(pro[4,46],rh_mid[n],pro[1,46])
  #qv_sur[n]=pro[1,2]*1000.
  #hm_sur[n]=mc.temperature_from_potential_temperature(pro[4,2],pro[0,2])
  #hm_sur[n]=mc.moist_static_energy(100.,hm_sur[n],qv_sur[n])
  #ws_mid[n]=mc.wind_speed(pro[2,46]*units('m/s'),pro[2,46]*units('m/s'))
  #wd_mid[n]=mc.wind_direction(pro[2,46]*units('m/s'),pro[2,46]*units('m/s'),)
 
  if np.abs(w_vvm_m[n]-w_cape_m[n])<2:
    print(casename,w_vvm_m[n],w_cape_m[n],w_vvm_m[n]-w_cape_m[n],ecape[n]) 
  
  n=n+1 

print(wd_mid)
#print(w_vvm_r,w_cape_r)

fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(4,4),dpi=300)
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,yerr=w_cape_r,fmt='bo',alpha=0.1,markersize=1,ecolor='c',elinewidth=0.5)
ax.plot(np.arange(-0,300,50),np.arange(-0,300,50),'k-',lw=1)
#im=ax.scatter(vvm_w,w_cape_m*np.ones(vvm_w.shape),s=0.1,alpha=0.6,c=vvm_w_prec,vmax=50,vmin=0,cmap='jet',edgecolors=None,linewidths=0.)
im=ax.scatter(w_vvm_m,w_cape_m,s=8,c=dcape_max,vmax=0.6,vmin=0,cmap='turbo',edgecolors='gray',linewidths=0.3)

cbar=fig.colorbar(im)
cbar.set_ticks(np.arange(0,1+.1,.1))

[bar.set_alpha(0.2) for bar in bars]
[cap.set_alpha(0.2) for cap in caps]

ax.set_xlim([-0,80])
ax.set_ylim([-0,80])
ax.set_xticks(np.arange(-0,80.1,20))
ax.set_yticks(np.arange(-0,80.1,20))
ax.set_xlabel('W in VVM [m/s]',fontsize=8)
ax.set_ylabel('W from CAPE [m/s]',fontsize=8)
ax.tick_params(labelsize=5)

plt.show()
#plt.savefig('./figure/w_cape_max_obj_prec.png')

exit()

# plot
fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(3,3),dpi=300)
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,yerr=w_cape_r,fmt='bo',alpha=0.01,markersize=0.5,ecolor='c',elinewidth=0.5)
ax.plot(np.arange(-0,300,50),np.arange(-0,300,50),'k-',lw=1)
im=ax.scatter(w_vvm_m,w_cape_m,s=8,c=rh_mid,vmax=80,vmin=10,cmap='turbo',edgecolors='gray',linewidths=0.3)

fig.colorbar(im)

[bar.set_alpha(0.3) for bar in bars]
[cap.set_alpha(0.3) for cap in caps]

ax.set_xlim([-0,80])
ax.set_ylim([-0,80])
ax.set_xticks(np.arange(-0,80.1,20))
ax.set_yticks(np.arange(-0,80.1,20))
ax.set_xlabel('W in VVM [m/s]',fontsize=6)
ax.set_ylabel('W from CAPE [m/s]',fontsize=6)
ax.set_title('5-km RH',fontsize=10)
ax.tick_params(labelsize=5)

plt.savefig('./figure/rh_mid.png')


fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(3,3),dpi=300)
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,yerr=w_cape_r,fmt='bo',alpha=0.01,markersize=0.5,ecolor='c',elinewidth=0.5)
ax.plot(np.arange(-0,300,50),np.arange(-0,300,50),'k-',lw=1)
im=ax.scatter(w_vvm_m,w_cape_m,s=8,c=hm_sur,vmax=330,vmin=350,cmap='turbo',edgecolors='gray',linewidths=0.3)

fig.colorbar(im)

[bar.set_alpha(0.3) for bar in bars]
[cap.set_alpha(0.3) for cap in caps]

ax.set_xlim([-0,80])
ax.set_ylim([-0,80])
ax.set_xticks(np.arange(-0,80.1,20))
ax.set_yticks(np.arange(-0,80.1,20))
ax.set_xlabel('W in VVM [m/s]',fontsize=6)
ax.set_ylabel('W from CAPE [m/s]',fontsize=6)
ax.set_title('Surface Hm',fontsize=10)
ax.tick_params(labelsize=5)
ax.tick_params(labelsize=5)

plt.savefig('./figure/hm_sur.png')


fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(3,3),dpi=300)
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,yerr=w_cape_r,fmt='bo',alpha=0.01,markersize=0.5,ecolor='c',elinewidth=0.5)
ax.plot(np.arange(-0,300,50),np.arange(-0,300,50),'k-',lw=1)
im=ax.scatter(w_vvm_m,w_cape_m,s=8,c=qv_sur,vmax=16,vmin=20,cmap='turbo',edgecolors='gray',linewidths=0.3)

fig.colorbar(im)

[bar.set_alpha(0.3) for bar in bars]
[cap.set_alpha(0.3) for cap in caps]

ax.set_xlim([-0,80])
ax.set_ylim([-0,80])
ax.set_xticks(np.arange(-0,80.1,20))
ax.set_yticks(np.arange(-0,80.1,20))
ax.set_xlabel('W in VVM [m/s]',fontsize=6)
ax.set_ylabel('W from CAPE [m/s]',fontsize=6)
ax.set_title('Surface Qv',fontsize=10)
ax.tick_params(labelsize=5)
ax.tick_params(labelsize=5)

plt.savefig('./figure/qv_sur.png')


fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(3,3),dpi=300)
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,yerr=w_cape_r,fmt='bo',alpha=0.01,markersize=0.5,ecolor='c',elinewidth=0.5)
ax.plot(np.arange(-0,300,50),np.arange(-0,300,50),'k-',lw=1)
im=ax.scatter(w_vvm_m,w_cape_m,s=8,c=cape,vmax=500,vmin=2500,cmap='turbo',edgecolors='gray',linewidths=0.3)

fig.colorbar(im)

[bar.set_alpha(0.3) for bar in bars]
[cap.set_alpha(0.3) for cap in caps]

ax.set_xlim([-0,80])
ax.set_ylim([-0,80])
ax.set_xticks(np.arange(-0,80.1,20))
ax.set_yticks(np.arange(-0,80.1,20))
ax.set_xlabel('W in VVM [m/s]',fontsize=6)
ax.set_ylabel('W from CAPE [m/s]',fontsize=6)
ax.set_title('CAPE',fontsize=10)
ax.tick_params(labelsize=5)
ax.tick_params(labelsize=5)

plt.savefig('./figure/cape.png')


fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(3,3),dpi=300)
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,yerr=w_cape_r,fmt='bo',alpha=0.01,markersize=0.5,ecolor='c',elinewidth=0.5)
ax.plot(np.arange(-0,300,50),np.arange(-0,300,50),'k-',lw=1)
im=ax.scatter(w_vvm_m,w_cape_m,s=8,c=ecape,vmax=0,vmin=500,cmap='turbo',edgecolors='gray',linewidths=0.3)

fig.colorbar(im)

[bar.set_alpha(0.3) for bar in bars]
[cap.set_alpha(0.3) for cap in caps]

ax.set_xlim([-0,80])
ax.set_ylim([-0,80])
ax.set_xticks(np.arange(-0,80.1,20))
ax.set_yticks(np.arange(-0,80.1,20))
ax.set_xlabel('W in VVM [m/s]',fontsize=6)
ax.set_ylabel('W from CAPE [m/s]',fontsize=6)
ax.set_title('ECAPE',fontsize=10)
ax.tick_params(labelsize=5)
ax.tick_params(labelsize=5)

plt.savefig('./figure/ecape.png')


fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(3,3),dpi=300)                                                          
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,yerr=w_cape_r,fmt='bo',alpha=0.01,markersize=0.5,ecolor='c',elinewidth=0.5)
ax.plot(np.arange(-0,300,50),np.arange(-0,300,50),'k-',lw=1)
im=ax.scatter(w_vvm_m,w_cape_m,s=8,c=ws_mid,vmax=0,vmin=8,cmap='turbo',edgecolors='gray',linewidths=0.3)

fig.colorbar(im)

[bar.set_alpha(0.3) for bar in bars]
[cap.set_alpha(0.3) for cap in caps]

ax.set_xlim([-0,80])
ax.set_ylim([-0,80])
ax.set_xticks(np.arange(-0,80.1,20))
ax.set_yticks(np.arange(-0,80.1,20))
ax.set_xlabel('W in VVM [m/s]',fontsize=6)
ax.set_ylabel('W from CAPE [m/s]',fontsize=6)
ax.set_title('5-km WS',fontsize=10)
ax.tick_params(labelsize=5)
ax.tick_params(labelsize=5)

plt.savefig('./figure/ws_mid.png')

fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(3,3),dpi=300)                                                          
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,yerr=w_cape_r,fmt='bo',alpha=0.01,markersize=0.5,ecolor='c',elinewidth=0.5)
ax.plot(np.arange(-0,300,50),np.arange(-0,300,50),'k-',lw=1)
im=ax.scatter(w_vvm_m,w_cape_m,s=8,c=wd_mid,vmax=0,vmin=360,cmap='hsv',edgecolors='gray',linewidths=0.3)

fig.colorbar(im,ticks=[0,45,90,135,180,225,270,315,360])

[bar.set_alpha(0.3) for bar in bars]
[cap.set_alpha(0.3) for cap in caps]

ax.set_xlim([-0,80])
ax.set_ylim([-0,80])
ax.set_xticks(np.arange(-0,80.1,20))
ax.set_yticks(np.arange(-0,80.1,20))
ax.set_xlabel('W in VVM [m/s]',fontsize=6)
ax.set_ylabel('W from CAPE [m/s]',fontsize=6)
ax.set_title('5-km WD',fontsize=10)
ax.tick_params(labelsize=5)
ax.tick_params(labelsize=5)

plt.savefig('./figure/wd_mid.png')


plt.show()  
