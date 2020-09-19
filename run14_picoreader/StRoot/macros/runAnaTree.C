#include <TSystem>


class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoAnaTreeMaker;


StChain *chain;
void runAnaTree(Int_t nEvents = 10, const Char_t *inputFile="test.list", const Char_t *outputFile="test.anaTree.root")
{
	
	//  Int_t nEvents = 1000;	
	//Load all the System libraries

	Int_t trigType = 0; // 0 = minbias; 1 = ht; 2 = st_mtd;
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();

	gSystem->Load("StPicoDstMaker");
	gSystem->Load("StPicoAnaTreeMaker");

	chain = new StChain();

	StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst");

	StPicoMtdCalibMaker *calibMaker;
	if(trigType==0 || trigType==1){
		calibMaker = new StPicoMtdCalibMaker();
		calibMaker->setInitFromFile(kTRUE);
		calibMaker->setApplyT0(kTRUE);
		calibMaker->setCalibFileT0("StRoot/StPicoDstMaker/Run14_AuAu200_CalibDtof.offline.dat");
	}

	StPicoAnaTreeMaker *treeMaker = new StPicoAnaTreeMaker(1,outputFile,picoMaker);
	treeMaker->setTriggerSelection(trigType);
	if(trigType==0){
		treeMaker->setVzCut(-8,8);
		treeMaker->setVzDiffCut(-4,4);
	}	
	if(trigType==1){
		treeMaker->setVzCut(-100,100);
	}	
	if(trigType==2){
		treeMaker->setVzCut(-100,100);
	}	

	/*
	treeMaker->setVzCut(-6,6);
	treeMaker->setVzDiffCut(-3,3);
	treeMaker->setPtCut(-0.15,30);
	treeMaker->setEtaCut(-1,1);
	treeMaker->setDcaCut(0,1);
	treeMaker->setnHitsFitCut(15,50);
	treeMaker->setnHitsDedxCut(10,50);
	treeMaker->setRatioCut(0.52,1);
	treeMaker->setEnSigECut(-2.5,2.5);
	treeMaker->setEInvBetaCut(0.97,1.03);
	*/

	chain->Init();
	cout<<"chain->Init();"<<endl;
	int total = picoMaker->chain()->GetEntries();
	cout << " Total entries = " << total << endl;
	if(nEvents>total) nEvents = total;



	for (Int_t i=0; i<nEvents; i++){
		if(i%1000==0)
			cout << "Working on eventNumber " << i << endl;

		chain->Clear();
		int iret = chain->Make(i);

		if (iret) { cout << "Bad return code!" << iret << endl; break;}

		total++;

	}

	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
	chain->Finish();
	cout << "****************************************** " << endl;
	cout << "total number of events  " << nEvents << endl;
	cout << "****************************************** " << endl;

	delete chain;


}
