#!/bin/bash
date

if [ $# -ne 1 ]; then
     echo "Please input one arguement --- '#job'(should be <1000) !"
	 exit 1
fi

dir="/star/u/tc88qy/AuAu/run17/54GeV/eff/pairEff/cocktails/heavyflavor/dy2ee/genTree"

if [ ! -d $dir/output ]; then
     mkdir -p $dir/output
fi

if [ ! -d $dir/script ]; then
     mkdir -p $dir/script
fi

if [ ! -d $dir/log ]; then
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

rm -rf script/*
rm -rf log/*

cp runMc.con runMcAll.con

njobs=$1
count=0
while [ $count -le $njobs ] 
do

  echo $count
  cp ./run.csh ./run_$count.csh
  
  echo "root -l -b <<EOF ">>run_$count.csh
  echo ".O2" >> run_$count.csh 
  echo ".L runhfevent_C.so">> run_$count.csh
  echo -n "runhfevent(">>run_$count.csh
  echo -n $count>>run_$count.csh
  echo -n ",1000000,">>run_$count.csh
  echo -n $RANDOM>>run_$count.csh
  echo -n ",0">>run_$count.csh
  echo ")">>run_$count.csh 
  echo ".q">>run_$count.csh      
  echo "EOF">>run_$count.csh 

  #qsub -o log/run_$count.log -e log/run_$count.err run_$count.csh

  mv run_$count.csh script/ 
	
  echo "Executable     = script/run_$count.csh">>runMcAll.con
  echo "Output         = log/run_${count}.out" >> runMcAll.con
  echo "Error          = log/run_${count}.err" >> runMcAll.con
  echo "Log            = log/run_${count}.log" >> runMcAll.con
  echo "Queue"         >> runMcAll.con
  echo  "     " >>runMcAll.con
	
  let "count=count+1" 

done   

condor_submit runMcAll.con
