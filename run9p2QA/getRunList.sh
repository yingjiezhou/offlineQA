#!/bin/csh

rm -rf run19_stPhysics_runnumber_DD.dat

get_file_list.pl -keys 'runnumber' -cond 'filetype=daq_reco_PicoDst,filename~st_physics,trgsetupname=production_9p2GeV_2020' -limit 0 >& run19_stPhysics_runnumber_DD.dat
