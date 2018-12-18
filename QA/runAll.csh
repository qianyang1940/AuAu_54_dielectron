#!/bin/bash
date

if [ $# -ne 1 ]; then
     echo "Please input one arguement!"
	 exit 1
fi

dir="/star/u/tc88qy/AuAu/run17/54GeV/QA"
echo $dir

if [ ! -d $dir/output_$1 ]; then
     mkdir -p $dir/output_$1
fi

if [ ! -d $dir/script_$1 ]; then
     mkdir -p $dir/script_$1
fi 

if [ ! -d  $dir/log_$1 ]; then
     mkdir -p $dir/log_$1
fi

if [ ! -d output_$1 ]; then
     ln -s $dir/output_$1 ./
fi

if [ ! -d script_$1 ]; then
     ln -s $dir/script_$1 ./
fi

if [ ! -d log_$1 ]; then
     ln -s $dir/log_$1 ./
fi

cp run.con runAll_$1.job
  
ifile=0
for FILE in `cat /star/u/tc88qy/AuAu/run17/54GeV/QA/datalist_$1`
do
     echo $FILE
     cp ./run.csh script_$1/$1_$ifile.csh
 
     echo "./produceQA $FILE output_$1/$ifile">>script_$1/$1_$ifile.csh

     echo "Executable       = script_$1/$1_$ifile.csh">>runAll_$1.job
     echo "Output           = log_$1/$1_$ifile.out">>runAll_$1.job
     echo "Error            = log_$1/$1_$ifile.err">>runAll_$1.job
     echo "Log              = log_$1/$1_$ifile.olog">>runAll_$1.job
     echo  "Queue" >>runAll_$1.job
     echo  "     " >>runAll_$1.job
      
     let "ifile+=1";
done
condor_submit runAll_$1.job
