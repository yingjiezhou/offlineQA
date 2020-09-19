#include <TSystem>


class StMaker;
class StChain;
class StPicoDstMaker;


StChain *chain;
void readPicoDst(const Char_t *inputFile="test.list", const Char_t *outputFile="test.root")
{
	Int_t nEvents = 200000000;
	//  Int_t nEvents = 1000;	
	//Load all the System libraries

	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();

	gSystem->Load("StPicoDstMaker");
	gSystem->Load("StMyAnalysisMaker");
	gSystem->Load("StMyTreeMaker");

	chain = new StChain();

	StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst");
	StPicoAnaTreeMaker *treeMaker = new StPicoAnaTreeMaker(1,inputFile,picoMaker);
	chain->Init();
	cout<<"chain->Init();"<<endl;
	int total = picoMaker->chain()->GetEntries();
	cout << " Total entries = " << total << endl;
	if(nEvents>total) nEvents = total;


	for (Int_t i=0; i<nEvents; i++){
		// if(i%1000==0)
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
