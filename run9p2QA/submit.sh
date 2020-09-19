#!/bin/bash

dir="/star/data01/pwg/yjzhou19/Local/run9p2QA"

echo $dir

if [ ! -d $dir/submitjob ]; then
     echo "NO \"submitjob\" directory !!!"
     exit
fi

if [ ! -d $dir/submitjob/submitErrInfo ]; then
     mkdir $dir/submitjob/submitErrInfo
fi

if [ ! -d $dir/submitjob/rootfiles ]; then
     mkdir $dir/submitjob/rootfiles
fi

mkdir -p submitReport submitScript submitList

rm -rf $dir/submitjob/submitErrInfo/*
rm -rf submitReport/*
rm -rf submitScript/*
rm -rf submitList/*

ln -fs $dir/Run20_9p2_AllTrigQA.xml ./

star-submit-template -template Run20_9p2_AllTrigQA.xml -entities dir=$dir

#for list in `ls $dir/picoDataList/sublist`
#do
#   echo $list
#   star-submit-template -template Run14_AuAu200_miniDst.xml -entities dir=$dir,list=$list
#done
