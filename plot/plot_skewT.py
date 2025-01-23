import os
import numpy as np
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.plots import add_metpy_logo, SkewT
from metpy.units import units

nz=70 
nt=145

path='/data/der0318/parcel_model_taiwanvvm/'
#caselist=os.listdir(path+'/taiwanVVM/')
caselist=['tpe20110816nor']

for casename in caselist:
  filename=path+'/data/mean_'+casename+'.dat'
  var=np.fromfile(filename,dtype=np.float32).reshape(nt,-1,nz)

  fig=plt.figure(figsize=(9,9))
  for t in np.arange(0,49+1,49):
    p=var[t,4,1:] 
    ptd=p
    ta=var[t,0,1:] *units.degC
    td=var[t,1,1:] *units.degC
    ws=var[t,2,1:] *units('m/s')
    wd=var[t,3,1:] *units.deg
    p=p *units.hPa
    u,v=mpcalc.wind_components(ws,wd)

    skew=SkewT(fig,rotation=45,rect=[0.05,0.1,0.8,0.8])
    
    if t==0:
      linestyle=':'
      shift=0.15
    else:
      linestyle='-'
      shift=0.3
    
    skew.plot(p,ta,'r'+linestyle,lw=1.5)
    skew.plot(ptd[ptd>=200.],td[ptd>=200.],'g'+linestyle,lw=1.5)
    skew.plot_barbs(ptd[ptd>=200.][0::2],u[ptd>=200.][0::2],v[ptd>=200.][0::2],xloc=1+shift)
    skew.ax.set_ylim(1000, 100)
    skew.ax.set_xlim(-15, 35)

    skew.ax.set_xlabel(f'Temperature ({ta.units:~P})')
    skew.ax.set_ylabel(f'Pressure ({p.units:~P})')

    lcl_pressure,lcl_temperature=mpcalc.lcl(p[0],ta[0],td[0])
    skew.plot(lcl_pressure,lcl_temperature,'ko',markerfacecolor='black')

    prof=mpcalc.parcel_profile(p,ta[0],td[0]).to('degC')
    skew.plot(p,prof,'k'+linestyle,linewidth=2)

    skew.shade_cin(p,ta,prof,td)
    skew.shade_cape(p,ta,prof)

  skew.ax.axvline(0,color='k',linewidth=1)

  skew.plot_dry_adiabats(lw=0.8)
  skew.plot_moist_adiabats(lw=0.8)
  skew.plot_mixing_lines(lw=0.8) 

  plt.show()

