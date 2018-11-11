#!/bin/bash
echo "starting configure condor environment"

#source /opt/geant4-10.2.2_condor/bin/geant4.sh 
#source /opt/root/root-6.06.06/bin/thisroot.sh 
source /data2/R710liziwei/.bashrc
#mkdir -p /data2/R710liziwei/g4/output/log
#mkg4
#cd ..
#for theta in 0 0.27 0.58 1 1.73 3.73
for theta in -75 -50 -25 0 30 75
do
macname="$gwk/mac/cdor_angle$theta.mac"
filename="theta$theta"
condoname="g4batch.condo"
echo $macname
cp $gwk/mac/cdor_anglechange.mac $macname

rad=`echo | awk -v theta=$theta '{print theta/180*3.1415926}'`
echo $rad
z=`echo | awk -v rad=$rad '{print sin(rad)/cos(rad)}'`
echo $z
sed -i "/direction/s/z/$z/" $macname

cp geant4.condo $condoname
sed -i "/Arguments/s/filename/$filename/" $condoname
sed -i "/Arguments/s#macname#$macname#" $condoname

condor_submit $condoname
done

