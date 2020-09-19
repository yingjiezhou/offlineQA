#!/bin/bash
date

if [ $# -lt 3 ]; then
  echo "not enough arguements"
  #echo "./submitAnaTreePDSF.sh <outName> <mode: 0-mb, 1-ht, 2-mtd> <prod type:0-low and mid, 1-high> <start day> <end day>"
  echo "./submitAnaTreePDSF.sh <outName> <mode: 0-mb, 1-ht, 2-mtd> <prod type:0-low and mid, 1-high>"
  exit
fi

if [ "$STAR_LEVEL" != "SL16d" ]; then
  echo "wrong library: $STAR_LEVEL, please use SL16d"
  exit
fi

#echo "./submitAnaTreePDSF.sh $1 $2 $3 $4 $5"
echo "./submitAnaTreePDSF.sh $1 $2 $3"

outName=$1
prodType=$3
#startDay=$4
#endDay=$5
max=1300

# cp src
scriptDir=/global/homes/x/xiao00/pwg_disk/AuAu14/Run14_HT_Anatree_production


#cp -r $scriptDir/StRoot .
#cp -r $scriptDir/StRoot/script/AuAu200.xml .
#cp -r $scriptDir/runNumberList_run14AuAu200mb .
#cp -r $scriptDir/recenter_correction.root .
#cp -r $scriptDir/runNumberList_run14AuAu200ht_high .
#cp -r $scriptDir/recenter_correction_ht_high.root .

#pathTmp=$(pwd)
#sed -i "s#${scriptDir}#${pathTmp}#g" AuAu200.xml

#echo $NDHOSTNAME
#starver SL16d
#cons

if [ ! -d .sl64_gcc447 ]; then
  ln -s .sl64_gcc482 .sl64_gcc447
fi

if [ ! -d out ]; then
  mkdir out
fi
if [ ! -d log ]; then
  mkdir log
fi

dir=$(pwd)
echo $dir
if [ ! -d ${dir}/out/out_$outName ]; then
  mkdir -p ${dir}/out/out_$outName
fi
if [ ! -d ${dir}/log/log_$outName ]; then
  mkdir -p ${dir}/log/log_$outName
fi

if [ ! -d submit ]; then
  mkdir submit
fi

if [ ! -d submit/$outName ]; then
  mkdir submit/$outName
fi

if [ ! -d fileList/$outName ]; then
  mkdir fileList/$outName
fi

#while true; do
  ./generateDayList.sh
#  for day in `seq 163 1 163`
  for day in `seq 121 1 145`
#  for day in `seq 163 1 163`
  do
    if [ $day -gt 0 ] && [ $day -lt 10 ]; then
      dayxxx="00${day}"
    fi
    if [ $day -gt 10 ] && [ $day -lt 100 ]; then
      dayxxx="0${day}"
    fi
    if [ $day -gt 100 ]; then
      dayxxx=${day}
    fi

   
    echo $dayxxx
    #loop runs 
    for runId in `cat $dir/List/day${dayxxx}.lis`
    do
      if [ ! -d /project/projectdirs/starprod/picodsts/Run14/AuAu/200GeV/physics2/P16id/${dayxxx}/${runId} ]; then
        echo "can't find picoDst for run $runId"
        continue
      fi
      cd $dir

      while :
      do
        jobs="$(qstat -u xiao00 | wc -l)" # check how many files are running under my name currently
        echo "jobs = $jobs"
        if [ $jobs -gt 0 ]; then
          jobs=$(($jobs-2))
        fi

        if [ $jobs -gt $max ]; then
          echo "current number of jobs = $jobs"
          echo "sleep for 20 m"
          date
          sleep 60m # sleep for x minutes before checking the job list
        else
          break;
        fi
      done

      echo "day = $dayxxx run = $runId"
      if [ ! -d $dir/submit/$1/$runId ]; then
        mkdir -p $dir/submit/$1/$runId
      fi
      if [ -e $dir/fileList/${outName}/run_PicoDst_${runId}.lis ]; then
        rm -rf $dir/fileList/${outName}/run_PicoDst_${runId}.lis
      fi
      cp $dir/blank $dir/fileList/${outName}/run_PicoDst_${runId}.lis

      ls -1 /project/projectdirs/starprod/picodsts/Run14/AuAu/200GeV/physics2/P16id/${dayxxx}/$runId/*.picoDst.root >> $dir/fileList/${outName}/run_PicoDst_$runId.lis 

      if [ -e $dir/fileList/${outName}/run_anaTree_${runId}.lis ]; then
        rm -rf $dir/fileList/${outName}/run_anaTree_${runId}.lis
      fi
      cp $dir/blank $dir/fileList/${outName}/run_anaTree_${runId}.lis
      ls -1 $dir/out/out_$1/$runId/*.anaTree.root >> $dir/fileList/${outName}/run_anaTree_${runId}.lis

      if [ -e $dir/fileList/${outName}/run_log_${runId}.lis ]; then
        rm -rf $dir/fileList/${outName}/run_log_${runId}.lis
      fi
      cp $dir/blank $dir/fileList/${outName}/run_log_${runId}.lis
      ls -1 $dir/log/log_$1/$runId/*.root.log >> $dir/fileList/${outName}/run_log_${runId}.lis

      if [ -e $dir/fileList/${outName}/run_tobeSubmitted_${runId}.lis ]; then
        rm -rf $dir/fileList/${outName}/run_tobeSubmitted_${runId}.lis
      fi
      cp $dir/blank $dir/fileList/${outName}/run_tobeSubmitted_${runId}.lis
      
      if [ ! -e $dir/fileList/${outName}/run_submitted_${runId}.lis ]; then
        cp $dir/blank $dir/fileList/${outName}/run_submitted_${runId}.lis
      fi
      
      if [ ! -e $dir/fileList/${outName}/run_transferred_${runId}.lis ]; then
        cp $dir/blank $dir/fileList/${outName}/run_transferred_${runId}.lis
      fi

      
      #generate filelist to be submitted
      $dir/diffPicoDstAnaTreeByRun.sh $dir/fileList/${outName}/run_PicoDst_${runId}.lis $dir/fileList/${outName}/run_anaTree_${runId}.lis $dir/fileList/${outName}/run_transferred_${runId}.lis $dir/fileList/${outName}/run_submitted_${runId}.lis $dir/fileList/${outName}/run_log_${runId}.lis $dir/fileList/${outName}/run_tobeSubmitted_${runId}.lis

      nLines=`cat $dir/fileList/$outName/run_tobeSubmitted_${runId}.lis | wc -l`
      echo ${nLines}
      
      if [ ${nLines} -le 0 ]; then
        echo "no input filelist for run $runId, skip..."
        continue
      fi

      cp $dir/AuAu200.xml $dir/submit/$1/temp.xml
      cd $dir/submit/$1/
      if [ ! -d ${dir}/out/out_$1/$runId ]; then
        mkdir -p ${dir}/out/out_$1/$runId
      fi
      if [ ! -d ${dir}/log/log_$1/$runId ]; then
        mkdir -p ${dir}/log/log_$1/$runId
      fi

      if [ ${nLines} -le 50 ]
      then
	  minfile_perJob=${nLines}
	  maxfile_perJob=${nLines}
      else
	  minfile_perJob=50
	  maxfile_perJob=100
      fi
      echo ${minfile_perJob}  ${maxfile_perJob}       
      sed 's/minFilesPerProcess="1" maxFilesPerProcess="1"/minFilesPerProcess="'${minfile_perJob}'" maxFilesPerProcess="'${maxfile_perJob}'"/g' temp.xml > AuAu200.xml
      star-submit-template -template AuAu200.xml -entities path=${dir},ver=$1,mode=$2,run=$runId,ptype=$prodType,listOfFiles=$dir/fileList/$outName/run_tobeSubmitted_${runId}.lis

      #     Minfile_perJob=${minfile_perJob} Maxfile_perJob=${maxfile_perJob}
      echo "submit!!!"
      if [ -f temp ]; then
        rm -f temp
      fi
      touch temp
      cat $dir/fileList/$outName/run_submitted_${runId}.lis $dir/fileList/$outName/run_tobeSubmitted_${runId}.lis | sort | uniq >> temp
      cp temp $dir/fileList/$outName/run_submitted_${runId}.lis

    done
    echo "Job submission for day $day finished!"
  done

#done
