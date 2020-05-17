
#!/bin/bash



#theta: 0-34, phi: 45-135
#position 0 0 0
macdir=/Users/liziwei/learning/CRTest/build/mac/angle
outdir=/Users/liziwei/learning/CRTest/build/angleop
exedir=/Users/liziwei/learning/CRTest/build
macname=AngleTemplate.mac

#for i in 5 
for i in 5 10 15 20 25 30 35
do
#for j in 45 
for j in 45 60 75 90 105 120 135 
do
    
        
        theta=`awk -v i=$i  'BEGIN {printf("%.3f",i/180*3.1415926)}'`
        phi=`awk -v j=$j 'BEGIN {printf("%.3f",j/180*3.1415926)}'`
        #x=r*cos(theta)
        #z=r*sin(theta)*cos(phi)
        #y=r*sin(theta)*sin(phi)
        x=1
        y=`echo | awk -v x=$x -v phi=$phi -v theta=$theta '{printf("%.3f",x/cos(theta)*sin(theta)*sin(phi));}'`
        z=`echo | awk -v x=$x -v phi=$phi -v theta=$theta '{printf("%.3f",x/cos(theta)*sin(theta)*cos(phi));}'`
        echo " theta & phi is: $theta, $phi"
        echo " x & y & z is: $x, $y, $z"

        themac="runtheta${i}phi${j}.mac"
        echo $themac
        cp $macdir/$macname $macdir/$themac
        sed -i "" "/direction/s/y/$y/g" $macdir/$themac
        sed -i "" "/direction/s/z/$z/g" $macdir/$themac
        $exedir/CRTest $exedir/mac/quartz.gdml $macdir/$themac $outdir/model1_theta${i}phi${j} >$outdir/model1_theta${i}phi${j}.log
        echo "$exedir/CRTest $exedir/mac/quartz.gdml $macdir/$themac $outdir/model1_theta${i}phi${j} >$outdir/model1_theta${i}phi${j}.log"

        #string="/gps/direction $x $y $z"
        #constantstr="/run/beamOn 1000"
        #echo $string
        #echo $constantstr
        #echo $string >> $macname
        #echo $constantstr >> $macname

    done
done

