#include "/star/u/syang/Macro/headers.h"
#include "/star/u/syang/Macro/function.C"
//#include "/Users/syang/Tools/Macro/headers.h"
//#include "/Users/syang/Tools/Macro/function.C"

vector<Int_t> DDRunIds;
vector<Int_t> JamieRunIds;
vector<Int_t> zeroHT2EvtsRunIds;

vector<Int_t> run11OffcialBadRunIds;

void checkRunIds(){
	DDRunIds.clear();
	JamieRunIds.clear();
	zeroHT2EvtsRunIds.clear();

	run11OffcialBadRunIds.clear();

	Int_t runId;

	ifstream indata;
	indata.open("run11_stHT_runnumber_DD.dat");
	if(indata.is_open()){
		cout<<"read in total run number from DD ...";
		while(indata>>runId){
			DDRunIds.push_back(runId);
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the total run number list from DD !!!"<<endl;
		return kFALSE;
	}
	indata.close();

	indata.open("run11_stHT_runnumber_Jamie.dat");
	if(indata.is_open()){
		cout<<"read in total run number from Jamie's webpage ...";
		while(indata>>runId){
			JamieRunIds.push_back(runId);
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the total run number list from Jamie's webpage !!!"<<endl;
		return kFALSE;
	}
	indata.close();

	indata.open("run11_zeroHT2EvtsInDD_runnumber.dat");
	if(indata.is_open()){
		cout<<"read in total run number from zeroHT2EvtsInDD ...";
		while(indata>>runId){
			zeroHT2EvtsRunIds.push_back(runId);
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the total run number list from zeroHT2EvtsInDD !!!"<<endl;
		return kFALSE;
	}
	indata.close();

	indata.open("bad_runs_refmult_year2011.txt");
	if(indata.is_open()){
		cout<<"read in total run number from bad_runs_refmult_year2011 ...";
		while(indata>>runId){
			run11OffcialBadRunIds.push_back(runId);
		}
		cout<<" [OK]"<<endl;
	}else{
		cout<<"Failed to load the total run number list bad_runs_refmult_year2011 !!!"<<endl;
		return kFALSE;
	}
	indata.close();

	// these runs do not have either EMC detector nor HT related triggers
	ofstream outdata;
	outdata.open("runIds_InDD_NotInJaimeWebpage.dat");
	for(Int_t i=0; i<DDRunIds.size(); i++){
		vector<Int_t>::iterator it = find(JamieRunIds.begin(), JamieRunIds.end(), DDRunIds[i]);
		if(it == JamieRunIds.end()) outdata<<DDRunIds[i]<<endl;
	}
	outdata.close();

	outdata.open("runIds_InJaimeWebpage_NoInDD.dat");
	for(Int_t i=0; i<JamieRunIds.size(); i++){
		vector<Int_t>::iterator it = find(DDRunIds.begin(), DDRunIds.end(), JamieRunIds[i]);
		if(it == DDRunIds.end()) outdata<<JamieRunIds[i]<<endl;
	}
	outdata.close();

	outdata.open("run11_offical_badruns_contain_HT2trigger.dat");
	for(Int_t i=0; i<run11OffcialBadRunIds.size(); i++){
		vector<Int_t>::iterator it = find(JamieRunIds.begin(), JamieRunIds.end(), run11OffcialBadRunIds[i]);
		if(it != JamieRunIds.end()) outdata<<run11OffcialBadRunIds[i]<<endl;
	}
	outdata.close();

	cout<<"Can find the run number in DD, however this run has no HT2 triggered events !!!"<<endl;
	for(Int_t i=0; i<zeroHT2EvtsRunIds.size(); i++){
		vector<Int_t>::iterator it = find(DDRunIds.begin(), DDRunIds.end(), zeroHT2EvtsRunIds[i]);
		if(it == DDRunIds.end()) cout<<zeroHT2EvtsRunIds[i]<<endl;
	}
	outdata.close();
}
