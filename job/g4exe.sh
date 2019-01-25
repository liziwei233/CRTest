#!/bin/bash
echo "starting configure condor environment"
echo " CWD: " $PWD
echo " ENV: "
which root
which geant4-config
mkdir -p `dirname $3`
./CRTest $1 $2 $3 $4
NAME=$3
cluster=$5
#echo $NAME > $PWD/name.log
#root -b -q "MultiTiersOutputfun_SiPM.C(\"$NAME\")"

echo $NAME > ../job/name.log
echo `date` >> $(dirname $3)/cluster.log
echo -e "$cluster\n\n" >> $(dirname $3)/cluster.log
#root -b -q "../job/MultiTiersOutputfun_SiPM.C(\"$NAME\",$thrd)"
root -b -q "../job/Outputfun_MCP.C(\"$NAME\")"
