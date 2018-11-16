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
        hadd ${rootname}data.root $path/$dirname/{0..39}data.root
        hadd ${rootname}.root $path/$dirname/{0..39}.root
        root -b -q "MultiTiers_TACor_sep.C(\"${rootname}data\")"
    fi
    done
else
    rootname=$path/$(basename $path)
    hadd ${rootname}data.root $path/{0..39}data.root
    hadd ${rootname}.root $path/{0..39}.root
    root -b -q "MultiTiers_TACor_sep.C(\"${rootname}data\")"
fi



