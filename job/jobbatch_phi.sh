#!/bin/bash
echo "starting configure condor environment"

#source /opt/geant4-10.2.2_condor/bin/geant4.sh 
#source /opt/root/root-6.06.06/bin/thisroot.sh 
source /data2/R710liziwei/.bashrc
#mkdir -p /data2/R710liziwei/g4/output/log
#mkg4
#cd ..
#for phi in 0 0.27 0.58 1 1.73 3.73
#for phi in -75 -50 -25 0 30 75
for phi in 0 30 60 90 120 150 180

do
macname="$gwk/mac/cdor_phi$phi.mac"
filename="phi$phi"
condoname="g4batch.condo"
echo $macname
cp $gwk/mac/cdor_phi.mac $macname

rad=`echo | awk -v phi=$phi '{print phi/180*3.1415926}'`
echo $rad
x=`echo | awk -v rad=$rad '{print -1/sqrt(3)*sin(rad)}'`
z=`echo | awk -v rad=$rad '{print 1/sqrt(3)*cos(rad)}'`
echo $z
sed -i "/direction/s/z/$z/" $macname
sed -i "/direction/s/x/$x/" $macname

cp g4batch_template.condo $condoname
sed -i "/Arguments/s/filename/$filename/" $condoname
sed -i "/Arguments/s#macname#$macname#" $condoname

condor_submit $condoname
done

