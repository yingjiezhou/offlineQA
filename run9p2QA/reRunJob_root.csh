#!/bin/csh

if ( $#argv == 0 ) then 
    echo 
	echo "Please input 2 arguments!" 
	echo "--------- first: The number of last root file"
	echo "--------- second: JobId"
    echo 
    exit
else if ( $#argv == 1 ) then
    echo 
    echo "Please input one more arguments!"
    echo 
	exit
else 
    echo 
    echo "Checking the lost root files, then resubmit job to produce them ..."
    echo 
endif	

set rootDir = $PWD/rootfiles

set lastIndex = $1
set name = st_physics_QA_$2

if ( -e reSubmit_$2.csh ) then
    rm -f reSubmit_$2.csh
	touch reSubmit_$2.csh
else 
    touch reSubmit_$2.csh
endif

#echo "#\!/bin/csh" >> reSubmit_$2.csh
echo -n "star-submit -r  " >> reSubmit_$2.csh

set index = 0
set tag = 0

while ( $index <= ${lastIndex} )
	if ( -e ${rootDir}/${name}_${index}.root ) then
	     echo "${name}_${index}.root is exist, skip..."
		 @ index++
	else 
	     echo "Resubmit job to produce ${name}_${index}.root !"
		 @ tag++
		 if( $tag == 1 ) then
		     echo -n "${index}" >> reSubmit_$2.csh
		 else
		     echo -n ",${index}" >> reSubmit_$2.csh
		 endif
		 @ index++
    endif
end

echo -n "  $2.session.xml" >> reSubmit_$2.csh 
chmod 755 reSubmit_$2.csh

if ( ${tag} == 0 ) then
    echo 
    echo "All root files are exist. Great!"
    echo 
	rm -f reSubmit_$2.csh
else 
    echo
    echo "There are ${tag} jobs need to resubmit ..."
    echo 
endif 

#./reSubmit_$2.csh
