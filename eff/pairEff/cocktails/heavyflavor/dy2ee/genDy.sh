#!/bin/bash
date

 dir="/gpfs/mnt/gpfs01/star/pwg/syang/run12/uu/minibias/eff/pairEff/cocktails_new/heavyflavor/dy2ee"
 
 if [ ! -d $dir/output ]; then
      mkdir -p $dir/output
 fi
 
 if [ ! -d $dir/script ]; then
      mkdir -p $dir/script
 fi 
 
 if [ ! -d  $dir/log ]; then
      mkdir -p $dir/log
 fi
 
 if [ ! -d output ]; then
      ln -s $dir/output ./
 fi
 
 if [ ! -d script ]; then
      ln -s $dir/script ./
 fi
 
 if [ ! -d log ]; then
      ln -s $dir/log ./
 fi
 
 echo $PWD

 echo "Universe     = vanilla">runDy.job
 echo "Notification = never">>runDy.job
 echo "Requirements = (CPU_Type != \"crs\") && (CPU_Experiment == \"star\")">>runDy.job
 echo "Initialdir   = $PWD">>runDy.job
 echo "GetEnv       = True">>runDy.job
 echo "+Experiment  = \"star\"">>runDy.job
 echo "+Job_Type    = \"cas\"">>runDy.job
 echo "">>runDy.job
 echo "">>runDy.job
   
 cenArray=(080 010 1040 4080 4060 6080 6070 7080)

#	num=${#cenArray[@]}
#	for ((i=0;i<num;i++))
#    {
#	    echo ${cenArray[i]}
#    }

 for cen in ${cenArray[*]}
 do
     echo $cen;
 
     if [ $cen != "080" -a $cen != "010" -a $cen != "1040" -a $cen != "4080" -a $cen != "4060" -a $cen != "6080" -a $cen != "6070" -a $cen != "7080" ]; then
          echo "Please check the centality string, which should be \"080\", \"010\", \"1040\", \"4080\", \"4060\", \"6080\", \"6070\", or \"7080\" !"
          exit 0
     fi
 
     echo "#!/bin/bash">script/runDy_Cen${cen}.sh
     echo "./sampledataeff $cen Scale dydata.list">>script/runDy_Cen${cen}.sh
#     echo "./sampledataeff $cen Individual dydata.list">>script/runDy_Cen${cen}.sh
     chmod 755 script/runDy_Cen${cen}.sh
     
     echo "Executable       = script/runDy_Cen${cen}.sh">>runDy.job
     echo "Output           = log/runDy_Cen${cen}.out">>runDy.job
     echo "Error            = log/runDy_Cen${cen}.err">>runDy.job
     echo "Log              = log/runDy_Cen${cen}.olog">>runDy.job
     echo  "Queue" >>runDy.job
     echo  "     " >>runDy.job

 done
      
 condor_submit runDy.job
