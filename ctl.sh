#!/usr/bin/bash

path=/data/der0318/parcel_model_taiwanvvm/taiwanVVM/

export list=$(ls $path)

for fname in $list
do
  #echo $fname
  cp template/pblh.ctl gs_ctl/pblh_$fname.ctl
  sed -i 's/runname/'$fname'/g' gs_ctl/pblh_$fname.ctl

  cp template/wpdf.ctl gs_ctl/wpdf_$fname.ctl
  sed -i 's/runname/'$fname'/g' gs_ctl/wpdf_$fname.ctl
done

