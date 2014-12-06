#!/bin/sh
#WRAP=gw
wufolder=files/smallexpx/ser1
for i in `seq 0 40000`
do
  wuname=smallexp_x_seria_1_120_1_5_ss7_n$i
  inpf=inp$wuname
  stdinpf=stdinp$wuname
  echo $(($wuname + 1))>./files/wucnt
  #echo Copiing files
  cp $wufolder/in.txt `dir_hier_path $inpf`
  cp stdin `dir_hier_path $stdinpf` 
  echo Submitting $wuname $inpf $stdinpf
  create_work -appname smallexpx -wu_name $wuname -wu_template templates/itmpl.smallexpx -result_template templates/rtmpl.smallexpx $inpf  $stdinpf
done
