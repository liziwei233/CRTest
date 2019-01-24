#!/bin/bash
#NAME=/baby/one/more/time
#echo $NAME>name.log

echo "  ./tr_batch.sh yourfilepath"


path=$1
if [ ${path:0-1} == "/" ]
then
path=${path%/*}
fi
echo "Your Path is: $path"

subdirN=$(ls -l $path | grep "^d" | wc -l)
flag=0
if [ $subdirN -gt $flag ]
then
    for dirname in `ls $path`
    do
    echo "Open the direct : $dirname"
    if [ -d $path/$dirname ]
    then
        rootname=$path/$dirname/${dirname}
        hadd -f ${rootname}data.root $path/$dirname/{0..1599}data.root
        #hadd -f ${rootname}dataFIX.root $path/$dirname/{0..1599}dataFIX.root
        #hadd -f ${rootname}dataCFD.root $path/$dirname/{0..1599}dataCFD.root
        hadd -f ${rootname}.root $path/$dirname/{0..1599}.root
        #root -b -q "MultiTiers_TACor_sep.C(\"${rootname}data\")"
    fi
    done
else
    rootname=$path/$(basename $path)
    hadd -f ${rootname}data.root $path/{0..1599}data.root
    #hadd -f ${rootname}dataFIX.root $path/{0..1599}dataFIX.root
    #hadd -f ${rootname}dataCFD.root $path/{0..1599}dataCFD.root
    hadd -f ${rootname}.root $path/{0..1599}.root
    #root -b -q "MultiTiers_TACor_sep.C(\"${rootname}data\")"
fi



