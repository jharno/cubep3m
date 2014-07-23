export OUTPUT_PATH='/vol-th/home/bnu_ztj_haoran/haoran/tides/LOS12/'

let irank=0
let nrank=7

while [ $irank -le $nrank ]; do

   ##echo $irank

   dirname=$OUTPUT_PATH'node'$irank
   mkdir -vp $dirname
   let irank=$irank+1
done
