#!/bin/bash
echo "starting configure condor environment"

#source /opt/geant4-10.2.2_condor/bin/geant4.sh 
#source /opt/root/root-6.06.06/bin/thisroot.sh 
source /data2/R710liziwei/.bashrc
#mkdir -p /data2/R710liziwei/g4/output/log
#mkg4
#cd ..
for pos in -12 -8 -4 0 
do
macname="$gwk/mac/cdor_pos$pos.mac"
filename="pos$pos"
condoname="g4batch.condo"
echo $macname

cp $gwk/mac/cdor_pc.mac $macname
sed -i "/centre/s/z/$pos/" $macname

cp g4batch_template.condo $condoname
sed -i "/Arguments/s/filename/$filename/" $condoname
sed -i "/Arguments/s#macname#$macname#" $condoname

condor_submit $condoname
done

