#!/usr/bin/bash

# remove the cases without rainfall
caselist=$(ls -I tpe20100802nor -I tpe20110615nor -I tpe20110723nor -I tpe20130703nor -I tpe20130825nor /data2/VVM/taiwanvvm_tpe/)


cd ./taiwanVVM

for casename in $caselist
do
  echo $casename
  ln -s /data2/VVM/taiwanvvm_tpe/$casename .
done


