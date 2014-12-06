#!/bin/sh
wufolder=input_wu_files
for i in `seq 0 40000`
do
  wuname=bivium8$i
  inpf=inp$wuname
  echo $(($wuname + 1))>./files/wucnt
  #echo Copiing files
  cp $wufolder/in `dir_hier_path $inpf`
  echo Submitting $wuname $inpf
  create_work -appname smallexpx -wu_name $wuname -wu_template templates/itmpl.smallexpx -result_template templates/rtmpl.smallexpx $inpf
done
