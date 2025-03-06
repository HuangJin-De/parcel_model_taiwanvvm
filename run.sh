#!/usr/bin/bash

path=/data/der0318/parcel_model_taiwanvvm/taiwanVVM/

export list=$(ls $path)

for fname in $list
do
  echo $fname
  #sed -e 's/runname/'$fname'/g' src/cal_cape.F > cal_cape.f
  #sed -e 's/runname/'$fname'/g' src/cal_pblh.F > cal_pblh.f
  
  sed -e 's/runname/'$fname'/g' src/cal_cape.F > cal_cape.f
  sed -e 's/runname/'$fname'/g' src/cal_mean.F > cal_mean.f

  #sed -e 's/runname/'$fname'/g' src/cal_wpdf.F > cal_wpdf.f

  #sed -e 's/runname/'$fname'/g' src/cal_cape.F > cal_cape.f
  #sed -e 's/runname/'$fname'/g' src/cal_cape_land.F > cal_cape_land.f

  make

  mpirun -np 24 ./a.out
done
