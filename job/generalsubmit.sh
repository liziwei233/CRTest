#!/bin/bash
echo "starting configure condor environment"

#source /opt/geant4-10.2.2_condor/bin/geant4.sh 
#source /opt/root/root-6.06.06/bin/thisroot.sh 
source /data2/R710liziwei/.bashrc
#mkdir -p /data2/R710liziwei/g4/output/log
#mkg4
#cd ..
condor_submit geant4.condo


