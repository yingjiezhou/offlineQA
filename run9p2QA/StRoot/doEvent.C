#include <TSystem>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void load();
void doEvent(Int_t nEvents=-1, const Char_t *inputFile="test.list", const TString outputFile="test/test.root", const Bool_t debug = kFALSE)
{
	load();

  string SL_version = "SL19b";
  string env_SL = getenv("STAR");
  if (env_SL.find(SL_version) == string::npos)
  {
    cout << "Environment Star Library does not match the requested library in runPicoD0EventMaker.C. Exiting..." << endl;
    exit(1);
  }

	StChain *chain = new StChain();
	chain->SetDebug(0);

	StMuTimer timer;
	timer.start();
  
  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(2, inputFile, "picoDstMaker");
	StQAMaker *qaMaker = new StQAMaker("StQAMaker",picoDstMaker);
	qaMaker->setOutFileName(outputFile);
	//qaMaker->setUseDefaultVtx(1);
	qaMaker->setUseDefaultVtx(0);
	qaMaker->setSelectVtxRanking(0);
	qaMaker->setMaxVtxR(2.);
//	qaMaker->setMaxVtxZ(100.);// default
  qaMaker->setMaxVtxZ(200.);
	qaMaker->setMaxVzDiff(3.);
//	qaMaker->setMinTrackPt(0.2); // default
  qaMaker->setMinTrackPt(0.15);

	qaMaker->setMaxTrackEta(1.0); // not used
	qaMaker->setMinNHitsFit(15);
	qaMaker->setMinNHitsFitRatio(0.52);
//	qaMaker->setMinNHitsDedx(10); // default
  qaMaker->setMinNHitsDedx(15);

	qaMaker->setMaxDca(3.); // default
//  qaMaker->setMaxDca(1.5);

	qaMaker->setMaxnSigmaE(2.);
	qaMaker->setMaxBeta2TOF(0.03);
	qaMaker->setMinBemcPt(3.5);
	qaMaker->setMinAdc0(290);
	qaMaker->setMinMaxPoverE(0.3, 1.5);
	qaMaker->setMaxZDist(3);
	qaMaker->setMaxPhiDist(0.02);
	qaMaker->setMinNEta(1);
	qaMaker->setMinNPhi(1);
	//qaMaker->setPrintMemory(1);
	//qaMaker->setPrintCpu(1);
	//qaMaker->setPrintConfig(1);
//	if(debug)
//		qaMaker->SetDebug(1);

	if(chain->Init()==kStERR) return;
	cout<<"chain->Init();"<<endl;

	if(nEvents<0){
		nEvents = picoDstMaker->chain()->GetEntries();
	}

	cout << "****************************************** " << endl;
	cout << "total number of events  " << nEvents << endl;
	cout << "****************************************** " << endl;

	for(Int_t i=0; i<nEvents; i++) {
		if(debug) {
			cout<<endl;
			cout<<"Working on eventNumber "<< i <<endl;
		} else {
			if(i%1000==0)
				cout << "Working on eventNumber " << i << endl;
		}

		chain->Clear();
		int iret = chain->Make(i);

		if(iret) { cout << "Bad return code!" << iret << endl; break;}
	}

	chain->Finish();
	delete chain;

	timer.stop();
	cout << "Total time = " << timer.elapsedTime() << " s" << endl;

	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
}

void load(){
  //Load all the System libraries
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StTpcDb");
  gSystem->Load("StEvent");
  gSystem->Load("StMcEvent");
  gSystem->Load("StMcEventMaker");
  gSystem->Load("StDaqLib");
  gSystem->Load("libgen_Tables");
  gSystem->Load("libsim_Tables");
  gSystem->Load("libglobal_Tables");
  gSystem->Load("StMagF");

  gSystem->Load("St_g2t.so");
  gSystem->Load("St_geant_Maker.so");
  gSystem->Load("StAssociationMaker");
  gSystem->Load("StMcAnalysisMaker");
  gSystem->Load("libgeometry_Tables");
  gSystem->Load("StTriggerUtilities");

  gSystem->Load("StEmcUtil");
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StPreEclMaker");
  gSystem->Load("StEpcMaker");
  gSystem->Load("StEmcSimulatorMaker");

  gSystem->Load("StDbLib");
  gSystem->Load("StDbUtilities");
  gSystem->Load("StDbBroker");
  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("St_db_Maker");

  gSystem->Load("StMtdHitMaker");
  gSystem->Load("StMtdUtil");
  gSystem->Load("StMtdMatchMaker");
  gSystem->Load("StMtdCalibMaker");
  gSystem->Load("StBTofUtil");
  gSystem->Load("StVpdCalibMaker");

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StQAMaker");
}
