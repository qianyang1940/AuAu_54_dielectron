#!/bin/bash
date

echo submitting pre-crashed jobs for "$1" 
cp $1.session.xml $1.session.xml.old$3
tag=0
ifile=0
while [ "$ifile" -le $2 ];
do  
    grep Goodbye /star/u/tc88qy/AuAu/run17/54GeV/picoDst/submitjob/submitStdout/$1_$ifile.out
    if [ $? -eq 0 ]; then
      echo "/star/u/tc88qy/AuAu/run17/54GeV/picoDst/miniTree/rootfiles_PicoDst/$1_$ifile.root exit, skip ..."
    else
      echo "re-submit crashed job for miniTree/rootfiles_PicoDst/$1_$ifile.root"
      cp -f ./submitjob/submitStdout/$1_$ifile.out ./submitjob/submitStdout/$1_$ifile.out.old$3
      cp -f ./submitjob/submitStdout/$1_$ifile.err ./submitjob/submitStdout/$1_$ifile.err.old$3
      star-submit -r $ifile $1.session.xml
      echo -n
     fi

    let "ifile+=1";
done
