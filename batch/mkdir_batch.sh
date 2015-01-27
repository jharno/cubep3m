export OUTPUT_PATH='../output/coarseproj/'

let irank=0
let nrank=511

while [ $irank -le $nrank ]; do

   ##echo $irank

   dirname=$OUTPUT_PATH'node'$irank
   mkdir -vp $dirname
   let irank=$irank+1
done
