import os
import numpy as np
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
cape=np.zeros(ncase)
ecape=np.zeros(ncase)
rh_mid=np.zeros(ncase)
qv_sur=np.zeros(ncase)
hm_sur=np.zeros(ncase)
ws_mid=np.zeros(ncase)
wd_mid=np.zeros(ncase)
n=0
for casename in caselist:
  fname1=path+'/obj_wmax_data/'+casename+'.dat'
  fname2=path+'/data/cape_land_'+casename+'.dat'
  fname3=path+'/data/mean_'+casename+'.dat'

  w_vvm=np.fromfile(fname1,dtype=np.float32)[1:].reshape(-1,2)
  nw=w_vvm.shape[0]
  vvm_w[0:nw,n]=w_vvm[:,0]
  vvm_w_time[0:nw,n]=w_vvm[:,1]/6.
  w_vvm=w_vvm[:,0]
  w_cape=np.fromfile(fname2,dtype=np.float32,count=(te-ts)*3*ny*nx,offset=ts*3*ny*nx*4).reshape(-1,3,ny,nx)[:,2,:,:].flatten()
  
  w_cape=np.compress(w_cape>=0.,w_cape)
  w_cape=np.sqrt(2.*w_cape)

  #print(w_vvm.shape,w_cape.shape)
 
  w_vvm_m[n]=w_vvm.mean()
  w_vvm_r[:,n]=np.abs([np.percentile(w_vvm,25),np.percentile(w_vvm,75)]-w_vvm_m[n])
  #w_vvm_r[:,n]=[np.percentile(w_vvm,25),np.percentile(w_vvm,75)]

  w_cape_m[n]=w_cape.mean()
  w_cape_r[:,n]=np.abs([np.percentile(w_cape,25),np.percentile(w_cape,75)]-w_cape_m[n])
  #w_cape_r[:,n]=[np.percentile(w_cape,25),np.percentile(w_cape,75)]

  ## env vars
  #env=np.fromfile(fname3,dtype=np.float32,count=(nz*6+3),offset=ts*(nz*6+3)*4)
  #pro=env[0:420].reshape(6,nz)
  #cape[n]=env[420]
  #ecape[n]=env[-1]
  #rh_mid[n]=pro[2,46]
  #qv_sur[n]=pro[5,2]*1000.
  #hm_sur[n]=pro[1,2]/1004.5
  #ws_mid[n]=pro[3,46]
  #wd_mid[n]=pro[4,46]
 
  if np.abs(w_vvm_m[n]-w_cape_m[n])<2:
    print(casename,w_vvm_m[n],w_cape_m[n],w_vvm_m[n]-w_cape_m[n],ecape[n]) 
  
  n=n+1 


#print(w_vvm_r,w_cape_r)

fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(4,4),dpi=300)
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,yerr=w_cape_r,fmt='bo',alpha=0.9,markersize=1,ecolor='c',elinewidth=0.5)
ax.plot(np.arange(-0,60,10),np.arange(-0,60,10),'k-',lw=1)
im=ax.scatter(vvm_w,w_cape_m*np.ones(vvm_w.shape),s=0.1,alpha=0.6,c=vvm_w_time,vmax=24,vmin=8,cmap='jet',edgecolors=None,linewidths=0.)

cbar=fig.colorbar(im)
cbar.set_ticks(np.arange(0,24.1,3))

[bar.set_alpha(0.2) for bar in bars]
[cap.set_alpha(0.2) for cap in caps]

ax.set_xlim([-0,40])
ax.set_ylim([-0,40])
ax.set_xticks(np.arange(-0,40.1,5))
ax.set_yticks(np.arange(-0,40.1,5))
ax.set_xlabel('W in VVM [m/s]',fontsize=8)
ax.set_ylabel('W from CAPE [m/s]',fontsize=8)
ax.tick_params(labelsize=5)

#plt.show()
plt.savefig('./figure/w_cape_occ_time_shear_var.png')

exit()

# plot
fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(3,3),dpi=300)
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,xerr=w_vvm_r,yerr=w_cape_r,fmt='bo',alpha=0.01,markersize=0.5,ecolor='c',elinewidth=0.5)
ax.plot(np.arange(-0,40,10),np.arange(-0,40,10),'k-',lw=1)
im=ax.scatter(w_vvm_m,w_cape_m,s=8,c=rh_mid,vmax=80,vmin=10,cmap='turbo',edgecolors='gray',linewidths=0.3)

fig.colorbar(im)

[bar.set_alpha(0.3) for bar in bars]
[cap.set_alpha(0.3) for cap in caps]

ax.set_xlim([-0,25])
ax.set_ylim([-0,40])
ax.set_xticks(np.arange(-0,25.1,5))
ax.set_yticks(np.arange(-0,40.1,5))
ax.set_xlabel('W in VVM [m/s]',fontsize=6)
ax.set_ylabel('W from CAPE [m/s]',fontsize=6)
ax.set_title('5-km RH',fontsize=10)
ax.tick_params(labelsize=5)

plt.savefig('./figure/rh_mid.png')


fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(3,3),dpi=300)
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,xerr=w_vvm_r,yerr=w_cape_r,fmt='bo',alpha=0.01,markersize=0.5,ecolor='c',elinewidth=0.5)
ax.plot(np.arange(-0,40,10),np.arange(-0,40,10),'k-',lw=1)
im=ax.scatter(w_vvm_m,w_cape_m,s=8,c=hm_sur,vmax=330,vmin=350,cmap='turbo',edgecolors='gray',linewidths=0.3)

fig.colorbar(im)

[bar.set_alpha(0.3) for bar in bars]
[cap.set_alpha(0.3) for cap in caps]

ax.set_xlim([-0,25])
ax.set_ylim([-0,40])
ax.set_xticks(np.arange(-0,25.1,5))
ax.set_yticks(np.arange(-0,40.1,5))
ax.set_xlabel('W in VVM [m/s]',fontsize=6)
ax.set_ylabel('W from CAPE [m/s]',fontsize=6)
ax.set_title('Surface Hm',fontsize=10)
ax.tick_params(labelsize=5)
ax.tick_params(labelsize=5)

plt.savefig('./figure/hm_sur.png')


fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(3,3),dpi=300)
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,xerr=w_vvm_r,yerr=w_cape_r,fmt='bo',alpha=0.01,markersize=0.5,ecolor='c',elinewidth=0.5)
ax.plot(np.arange(-0,40,10),np.arange(-0,40,10),'k-',lw=1)
im=ax.scatter(w_vvm_m,w_cape_m,s=8,c=qv_sur,vmax=16,vmin=20,cmap='turbo',edgecolors='gray',linewidths=0.3)

fig.colorbar(im)

[bar.set_alpha(0.3) for bar in bars]
[cap.set_alpha(0.3) for cap in caps]

ax.set_xlim([-0,25])
ax.set_ylim([-0,40])
ax.set_xticks(np.arange(-0,25.1,5))
ax.set_yticks(np.arange(-0,40.1,5))
ax.set_xlabel('W in VVM [m/s]',fontsize=6)
ax.set_ylabel('W from CAPE [m/s]',fontsize=6)
ax.set_title('Surface Qv',fontsize=10)
ax.tick_params(labelsize=5)
ax.tick_params(labelsize=5)

plt.savefig('./figure/qv_sur.png')


fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(3,3),dpi=300)
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,xerr=w_vvm_r,yerr=w_cape_r,fmt='bo',alpha=0.01,markersize=0.5,ecolor='c',elinewidth=0.5)
ax.plot(np.arange(-0,40,10),np.arange(-0,40,10),'k-',lw=1)
im=ax.scatter(w_vvm_m,w_cape_m,s=8,c=cape,vmax=500,vmin=2500,cmap='turbo',edgecolors='gray',linewidths=0.3)

fig.colorbar(im)

[bar.set_alpha(0.3) for bar in bars]
[cap.set_alpha(0.3) for cap in caps]

ax.set_xlim([-0,25])
ax.set_ylim([-0,40])
ax.set_xticks(np.arange(-0,25.1,5))
ax.set_yticks(np.arange(-0,40.1,5))
ax.set_xlabel('W in VVM [m/s]',fontsize=6)
ax.set_ylabel('W from CAPE [m/s]',fontsize=6)
ax.set_title('CAPE',fontsize=10)
ax.tick_params(labelsize=5)
ax.tick_params(labelsize=5)

plt.savefig('./figure/cape.png')


fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(3,3),dpi=300)
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,xerr=w_vvm_r,yerr=w_cape_r,fmt='bo',alpha=0.01,markersize=0.5,ecolor='c',elinewidth=0.5)
ax.plot(np.arange(-0,40,10),np.arange(-0,40,10),'k-',lw=1)
im=ax.scatter(w_vvm_m,w_cape_m,s=8,c=ecape,vmax=0,vmin=500,cmap='turbo',edgecolors='gray',linewidths=0.3)

fig.colorbar(im)

[bar.set_alpha(0.3) for bar in bars]
[cap.set_alpha(0.3) for cap in caps]

ax.set_xlim([-0,25])
ax.set_ylim([-0,40])
ax.set_xticks(np.arange(-0,25.1,5))
ax.set_yticks(np.arange(-0,40.1,5))
ax.set_xlabel('W in VVM [m/s]',fontsize=6)
ax.set_ylabel('W from CAPE [m/s]',fontsize=6)
ax.set_title('ECAPE',fontsize=10)
ax.tick_params(labelsize=5)
ax.tick_params(labelsize=5)

plt.savefig('./figure/ecape.png')


fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(3,3),dpi=300)                                                          
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,xerr=w_vvm_r,yerr=w_cape_r,fmt='bo',alpha=0.01,markersize=0.5,ecolor='c',elinewidth=0.5)
ax.plot(np.arange(-0,40,10),np.arange(-0,40,10),'k-',lw=1)
im=ax.scatter(w_vvm_m,w_cape_m,s=8,c=ws_mid,vmax=0,vmin=8,cmap='turbo',edgecolors='gray',linewidths=0.3)

fig.colorbar(im)

[bar.set_alpha(0.3) for bar in bars]
[cap.set_alpha(0.3) for cap in caps]

ax.set_xlim([-0,25])
ax.set_ylim([-0,40])
ax.set_xticks(np.arange(-0,25.1,5))
ax.set_yticks(np.arange(-0,40.1,5))
ax.set_xlabel('W in VVM [m/s]',fontsize=6)
ax.set_ylabel('W from CAPE [m/s]',fontsize=6)
ax.set_title('5-km WS',fontsize=10)
ax.tick_params(labelsize=5)
ax.tick_params(labelsize=5)

plt.savefig('./figure/ws_mid.png')

fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(3,3),dpi=300)                                                          
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,xerr=w_vvm_r,yerr=w_cape_r,fmt='bo',alpha=0.01,markersize=0.5,ecolor='c',elinewidth=0.5)
ax.plot(np.arange(-0,40,10),np.arange(-0,40,10),'k-',lw=1)
im=ax.scatter(w_vvm_m,w_cape_m,s=8,c=wd_mid,vmax=0,vmin=360,cmap='hsv',edgecolors='gray',linewidths=0.3)

fig.colorbar(im,ticks=[0,45,90,135,180,225,270,315,360])

[bar.set_alpha(0.3) for bar in bars]
[cap.set_alpha(0.3) for cap in caps]

ax.set_xlim([-0,25])
ax.set_ylim([-0,40])
ax.set_xticks(np.arange(-0,25.1,5))
ax.set_yticks(np.arange(-0,40.1,5))
ax.set_xlabel('W in VVM [m/s]',fontsize=6)
ax.set_ylabel('W from CAPE [m/s]',fontsize=6)
ax.set_title('5-km WD',fontsize=10)
ax.tick_params(labelsize=5)
ax.tick_params(labelsize=5)

plt.savefig('./figure/wd_mid.png')


plt.show()  
