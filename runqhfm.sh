#!/bin/bash

# running version 5.8


# run the program using the relative path


cpus=1
jobid=script_test

# delete mcr cache
rm -fr ~/.mcrCache9.12/
nohup ./skyrmionseffective.sh.app /usr/remote/apps/matlab/current $size $anneal $static $dyn $initcond $sk $skdist $impx $impy $read $bc $Eprotocol $Bprotocol $cpus $U $w $B $J $damping $ext $Bamp $Bfreq $tmax $jobid $texture &
