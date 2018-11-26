#!/bin/bash

FULLNAME=$(cat name.log)
echo $FULLNAME
POS=$(dirname $FULLNAME)
echo $POS
NAME=$(basename $POS)
echo $NAME
thrd=1
#NAME='/data2/R710liziwei/g4/output/data/quartz/Test1/0'
#root -b -q "MultiTiersOutputfun_SiPM.C(\"$FULLNAME\",$thrd)"
#root -b -q "Outputfun_MCP.C(\"$FULLNAME\")"
root -b -q "Outputfun_MCP.C(\"/data2/R710liziwei/g4/FTOF_server/CRTest/build/CRTest\")"
