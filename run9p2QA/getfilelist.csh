#!/bin/csh

rm -rf micro.list 

# in range use
#get_file_list.pl -keys 'path,filename' -cond 'storage=local,filetype=daq_reco_MuDst,filename~st_ht,production=P11id,trgsetupname=AuAu200_production_2011' -limit 0 -distinct -delim '/' >& micro.list

get_file_list.pl -keys 'path,filename' -cond 'filetype=daq_reco_MuDst,filename~st_physics,trgsetupname=production_9p2GeV_2020' -limit 0 -distinct -delim '/' >& micro.list


## outside range use
#get_file_list.pl -keys 'path,filename' -cond 'storage=hpss,filetype=daq_reco_MuDst,filename~st_mtd,production=P15ie,trgsetupname=AuAu_200_production_high_2014,runnumber][15166023-15166028,tpx=1,sanity=1' -limit 0 -distinct -delim '/' >& micro.list

## or logic 
#get_file_list.pl -keys 'path,filename' -cond 'storage=hpss,filetype=daq_reco_MuDst,filename~st_mtd,production=P15ie,trgsetupname=AuAu_200_production_high_2014||AuAu_200_production_mid_2014||AuAu_200_production_low_2014,runnumber[]15166001-15166060,tpx=1,sanity=1' -limit 0 -distinct -delim '/' >& micro.list

##day select
#get_file_list.pl -keys 'path,filename' -cond 'storage=local,filetype=daq_reco_MuDst,filename~st_physics,trgsetupname=production_pp201long3_2015||production_pp200long2_2015||production_pp200long_2015||production_pp200trans_2015,daynumber~41,tpx=1,sanity=1' -limit 0 -distinct -delim '/' >& micro.list

##output daynumber
#get_file_list.pl -keys 'daynumber' -cond 'storage=local,filetype=daq_reco_MuDst,filename~st_physics,trgsetupname=production_pp201long3_2015||production_pp200long2_2015||production_pp200long_2015||production_pp200trans_2015,tpx=1,sanity=1' -limit 0 -distinct -delim '/' >& micro.list

## day select 
#get_file_list.pl -keys 'path,filename' -cond 'storage=hpss,filetype=daq_reco_MuDst,filename~st_mtd,production=P15ie,trgsetupname=AuAu_200_production_high_2014||AuAu_200_production_mid_2014||AuAu_200_production_low_2014,daynumber~151' -limit 0 -distinct -delim '/' >& micro.list
