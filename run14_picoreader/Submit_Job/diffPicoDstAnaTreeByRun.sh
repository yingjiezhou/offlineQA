#!/bin/bash

#: << 'MYCOMMENT'
#####################
# obtain the difference between the current PicoDst List and produced Pico List
#####################
picoDstList=$1
anaTreeList=$2
copiedanaTreeList=$3
submittedList=$4
logList=$5
outputList=$6

if [ -e $outputList ]; then
  rm -f $outputList
fi
touch $outputList

for picoDst in `cat $picoDstList`
do
  anaTree=`basename ${picoDst/picoDst/anaTree}`
  del=`grep -l $anaTree $anaTreeList` 
  copiedl=`grep -l $anaTree $copiedanaTreeList` 
  picoName=`basename ${picoDst}`
  logl=`grep -l $picoName $logList` 
  submitl=`grep -l $picoName $submittedList`

  #echo "nAnaTree = ${#del}, nCopied = ${#copiedl}"
  #Oct 11, 2016 remove this line to resubmit lost jobs.
  #if [ ${#submitl} -gt 0 ] && [ ${#logl} -le 0 ]; then
  #  continue
  #fi
  if [ ${#del} -lt 1 ] && [ ${#copiedl} -le 0 ]; then
    #echo "add to list"
    echo $picoDst >> $outputList
  else
    anaTreeWithPath=`grep $anaTree $anaTreeList`
    #echo "------------------"
    #echo "anaTree = $anaTree anaTreeList = $anaTreeList"
    #echo "anaTreeWithPath = ${anaTreeWithPath}"
    if [ ${#anaTreeWithPath} -gt 1 ]; then
      size=`ls -l ${anaTreeWithPath} | awk '{print $5}'`
      #echo "size = $size"
      if [ ${size} -lt 10000 ]; then
        #echo "add to list"
        echo $picoDst >> $outputList
      fi
    fi
  fi
done

