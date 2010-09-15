#!/bin/csh

set init_firstSim='1'
set init_lastSim='1'
set numSims_per_Jobs='4'
set MAXSIM='200'

set old_firstSim=$init_firstSim
set old_lastSim=$init_lastSim
set firstSim='141' #or any number
@ lastSim = ( $firstSim + $numSims_per_Jobs - 1 )


while ( $lastSim <= $MAXSIM ) 

   echo 'sending [' $firstSim 'to' $lastSim ']'

   sed '14 s/'"$old_firstSim"'/'"$firstSim"'/' cubep3m_loop.pbs -i 
   sed '15 s/'"$old_lastSim"'/'"$lastSim"'/' cubep3m_loop.pbs -i

   #grep MIN cubep3m_loop.pbs 
   #grep MAX cubep3m_loop.pbs 



   qsub cubep3m_loop.pbs
   #echo qsubed

   set old_firstSim=$firstSim
   set old_lastSim=$lastSim

   @ firstSim = ( $lastSim + 1 )
   @ lastSim = ( $lastSim + $numSims_per_Jobs )

end

   sed '14 s/'"$old_firstSim"'/'"$init_firstSim"'/' cubep3m_loop.pbs -i
   sed '15 s/'"$old_lastSim"'/'"$init_lastSim"'/' cubep3m_loop.pbs -i


