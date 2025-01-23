#!/usr/bin/bash

path=/data/der0318/parcel_model/taiwanVVM/

export list=$(ls $path)

for fname in $list
do
  echo $fname
  #sed -e 's/runname/'$fname'/g' src/cal_cape.F > cal_cape.f
  #sed -e 's/runname/'$fname'/g' src/cal_pblh.F > cal_pblh.f
  
  sed -e 's/runname/'$fname'/g' src/cal_mean.F > cal_mean.f

  make

  mpirun -np 24 ./a.out
done
