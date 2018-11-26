#!/bin/bash
echo "starting configure condor environment"

#source /opt/geant4-10.2.2_condor/bin/geant4.sh 
#source /opt/root/root-6.06.06/bin/thisroot.sh 
source /data2/R710liziwei/.bashrc
#mkdir -p /data2/R710liziwei/g4/output/log
#mkg4
#cd ..
for pos in -20 -10 0 10 20
do
macname="$gwk/mac/cdor_pos$pos.mac"
filename="pos$pos"
condoname="g4batch.condo"
echo $macname

cp $gwk/mac/cdor_poschange.mac $macname
sed -i "/position/s/z/$pos/" $macname

cp g4batch_template.condo $condoname
sed -i "/Arguments/s/filename/$filename/" $condoname
sed -i "/Arguments/s#macname#$macname#" $condoname

condor_submit $condoname
done

