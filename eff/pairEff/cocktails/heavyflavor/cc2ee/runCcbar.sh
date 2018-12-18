#!/bin/bash
date

 if [ $# -ne 1 ]; then
    echo "Please input the centrality parameter, which should be \"080\", \"010\", \"1040\", \"4080\", \"4060\", \"6080\", \"6070\", or \"7080\""
    exit 0
 fi
 
 cen=$1
 
 if [ $cen != "080" -a $cen != "010" -a $cen != "1040" -a $cen != "4080" -a $cen != "4060" -a $cen != "6080" -a $cen != "6070" -a $cen != "7080" ]; then
    echo "Please check the centality string, which should be \"080\", \"010\", \"1040\", \"4080\", \"4060\", \"6080\", \"6070\", or \"7080\" !"
    exit 0
 fi

 dir="/gpfs/mnt/gpfs01/star/pwg/syang/run12/uu/minibias/eff/pairEff/cocktails_new/heavyflavor/cc2ee"
 
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
 
 echo "Universe     = vanilla">runCcbar_Cen${cen}.job
 echo "Notification = never">>runCcbar_Cen${cen}.job
 echo "Requirements = (CPU_Type != \"crs\") && (CPU_Experiment == \"star\")">>runCcbar_Cen${cen}.job
 echo "Initialdir   = $PWD">>runCcbar_Cen${cen}.job
 echo "GetEnv       = True">>runCcbar_Cen${cen}.job
 echo "+Experiment  = \"star\"">>runCcbar_Cen${cen}.job
 echo "+Job_Type    = \"cas\"">>runCcbar_Cen${cen}.job
 echo "">>runCcbar_Cen${cen}.job
 echo "">>runCcbar_Cen${cen}.job
   
 ifile=0
 for FILE in `cat datalist_all`
 do
      echo $FILE

	  echo "#!/bin/bash">script/Cen${cen}_${ifile}.sh
      echo "./sampledataeff $cen Scale $FILE $ifile">>script/Cen${cen}_${ifile}.sh
#      echo "./sampledataeff $cen Individual $FILE $ifile">>script/Cen${cen}_${ifile}.sh
	  chmod 755 script/Cen${cen}_${ifile}.sh

      echo "Executable       = script/Cen${cen}_${ifile}.sh">>runCcbar_Cen${cen}.job
      echo "Output           = log/Cen${cen}_${ifile}.out">>runCcbar_Cen${cen}.job
      echo "Error            = log/Cen${cen}_${ifile}.err">>runCcbar_Cen${cen}.job
      echo "Log              = log/Cen${cen}_${ifile}.olog">>runCcbar_Cen${cen}.job
      echo  "Queue" >>runCcbar_Cen${cen}.job
      echo  "     " >>runCcbar_Cen${cen}.job
       
      let "ifile+=1";
 done
condor_submit runCcbar_Cen${cen}.job
