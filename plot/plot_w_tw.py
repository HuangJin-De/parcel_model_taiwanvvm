import os
import numpy as np
import matplotlib.pyplot as plt


nx=1024
ny=1024

ts=48  # 8 LST
te=54  # 9 LST

path='/data/der0318/parcel_model_taiwanvvm/'
caselist=os.listdir(path+'/taiwanVVM/')
#caselist=['tpe20050702nor']
ncase=len(caselist)
print(ncase)

w_vvm_m=np.zeros(ncase)
w_cape_m=np.zeros(ncase)
w_vvm_r=np.zeros([2,ncase])
w_cape_r=np.zeros([2,ncase])
n=0
for casename in caselist:
  fname1=path+'/obj_wmax_data/'+casename+'.dat'
  fname2=path+'/data/cape_land_'+casename+'.dat'

  w_vvm=np.fromfile(fname1,dtype=np.float32)[1:]
  w_cape=np.fromfile(fname2,dtype=np.float32,count=(te-ts)*3*ny*nx,offset=ts*3*ny*nx*4).reshape(-1,3,ny,nx)[:,2,:,:].flatten()
  
  w_cape=np.compress(w_cape>=0.,w_cape)
  w_cape=np.sqrt(2.*w_cape)

  #print(w_vvm.shape,w_cape.shape)
 
  w_vvm_m[n]=w_vvm.mean()
  w_vvm_r[:,n]=[np.percentile(w_vvm,25),np.percentile(w_vvm,75)]

  w_cape_m[n]=w_cape.mean()
  w_cape_r[:,n]=[np.percentile(w_cape,25),np.percentile(w_cape,75)]

  if np.abs(w_vvm_m[n]-w_cape_m[n])<2:
    print(casename,w_vvm_m[n],w_cape_m[n],w_vvm_m[n]-w_cape_m[n])
   
  n=n+1 



fif,ax=plt.subplots(nrows=1,ncols=1,figsize=(3,3),dpi=300)
markers,caps,bars=ax.errorbar(w_vvm_m,w_cape_m,xerr=w_vvm_r,yerr=w_cape_r,fmt='bo',markersize=0.5,ecolor='c',elinewidth=0.5)

[bar.set_alpha(0.3) for bar in bars]
[cap.set_alpha(0.3) for cap in caps]


ax.plot(np.arange(-0,40,10),np.arange(-0,40,10),'k-',lw=1)

#ax.set_xlim([-0,30])
#ax.set_ylim([-0,30])
#ax.set_xticks(np.arange(-0,30.1,10))
#ax.set_yticks(np.arange(-0,30.1,10))


plt.show()


  
