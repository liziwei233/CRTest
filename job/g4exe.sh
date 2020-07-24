#!/bin/bash
echo "starting configure condor environment"
echo " CWD: " $PWD
echo " ENV: "
source /home/lizw/.bashrc
path=$1
NAME=${path}/$2
gdmlname=$3
macname=$4
process=$5
cluster=$6
rootname=${NAME}/${process}
mkdir ${NAME}
echo ${path}/CRTest ${path}/../mac/${gdmlname} ${path}/../mac/${macname} ${rootname} ${process} 
#${path}/CRTest ${path}/../mac/${gdmlname} ${path}/../mac/${macname} ${rootname} ${process} 

#echo $NAME > $PWD/name.log
#root -b -q "MultiTiersOutputfun_SiPM.C(\"$NAME\")"

#echo $NAME > ../job/name.log
echo `date` >> ${path}/cluster.log
echo -e "$cluster\n\n" >> ${path}/cluster.log
#root -b -q "../job/ME
root -b -q "$path/../job/Output.C(\"$NAME\",\"$process\")"
