#!/bin/tcsh


date

set BaseDir=/global/u1/x/xiao00/WWW

set FileName=$1

scp -r ${FileName} xiao00@pdsf.nersc.gov:${BaseDir}
