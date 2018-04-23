#!/bin/bash
echo "starting configure condor environment"
echo " CWD: " $PWD
echo " ENV: "
which root
which geant4-config
mkdir -p `dirname $3`
./CRTest $@
NAME=$3
echo $NAME > $PWD/name.log
root -b -q "Outputfun_MCP.C(\"$NAME\")"
