#!/bin/tcsh


date

set BaseDir=/star/starlib/doc_protected/pwg/heavy

set FileName=$1

scp -r ${FileName} xiao00@rftpexp.rhic.bnl.gov:~/

