This text is meaning to explain the function and using method of every file among this achive.


===================EXCUATE FILES==============
==============================================
>env.sh         ----set up environment, and submit your condo file.

>job.sh         ----set up environment, have a loop to produce batch condo files according to your requrement,for example, different thresholds. And submit these condo files.

>jobbatch_pos.sh ----set up environment, have a loop to produce batch condo files according to different positon. and submit these condo files.

>jobbatch_angle.sh  ----set up environment, have a loop to produce batch condo files according to different angle. and submit these condo files.

>g4exe.sh       ----write for your condo file. Its function is to excuate your geant4 exe.


>rtexe.sh       ----write for your condo file. Its function is to excuate your root code. MultiTiersOutputfun_SiPM.C.

>test.sh        ----write for your condo file. Its function is to test if the job your submit is work normaly.

>tr.sh          ----analysis data genarate from gent4 exe. Its have a loop to analysis data that have regular rule, for example, have different thershold.

>debug.sh       ----test the MultiTiersOutputfun_SiPM.C.If you have some modification in this file, you can excuate debug.sh to debug the code.



================CONDOR FILES===============
===========================================
>geant4.condo   ----write for excuating your jobs. It set up a exe (g4exe.sh) you want to excuate, the numbers of job (N_JOBS), and the directory to output results( RUN = your dir)

>rt.condo       ----write for excuating your jobs. It set up a exe (rtexe.sh) you want to excuate, the numbers of job(N_JOBS), and the directory to output results (RUN).

>g4batch.condo  ----write for excuate your jobs. It's a template, your can use .sh to modify the parameter of it and submit batch jobs.

===============ROOT CODE=================
=========================================
MultiTiersOutputfun_SiPM.C  ----convert the time information to waveform of SiPM or MCP, and discriminate it to get the time stamp. (For different detector construction, you should modify this file)

MultiTiers_TACor_sep.C      ----Correct the time stamp by TACor method.(For different detector construction,you should modify this file)


how to submit your job?

First of all, change the parameter of g4exe under the achive mac/,
then, change the dirname in the condor file geant.condo,
finally, excuate env.sh to submit your job.

how to submit your batch jobs?

First of all, prepare a .mac template under the achive mac/,
then,change the jobbatch.sh depends on your requirement,
finally, excuate it.



