#!/usr/bin/bash

path=/data/der0318/parcel_model_taiwanvvm/taiwanVVM/

export list=$(ls $path)

for fname in $list
do
  #echo $fname
  cd $path/$fname/
  ./make_ctl.sh

done

