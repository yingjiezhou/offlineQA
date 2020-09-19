
#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;
class StMuDstMaker;


StChain *chain;
void makeAnaTreeFromPicoDst(Int_t nEvents = 10000000, const Int_t runnumber=15107008,
    const Char_t *inputFile="/project/projectdirs/starprod/picodsts/Run14/AuAu/200GeV/physics2/P16id/107/15107008/st_physics_15107008_raw_3000056.picoDst.root",
    const bool creatingPhiWgt = kFALSE, const int prodMod = 0, const int emcMode=1, const int prodType = 0
    ){
  //Int_t nEvents = 10000000;
  //Int_t nEvents = 1000;	
  //Load all the System libraries

  gSystem->Load("libTable");
  gSystem->Load("libPhysics");
  gSystem->Load("St_base");
  gSystem->Load("StChain");
  gSystem->Load("St_Tables");
  gSystem->Load("StUtilities");        // new addition 22jul99
  gSystem->Load("StTreeMaker");
  gSystem->Load("StIOMaker");
  gSystem->Load("StarClassLibrary");
  gSystem->Load("StTriggerDataMaker"); // new starting from April 2003
  gSystem->Load("StBichsel");
  gSystem->Load("StEvent");
  gSystem->Load("StEventUtilities");
  gSystem->Load("StDbLib");
  gSystem->Load("StEmcUtil");
  gSystem->Load("StTofUtil");
  gSystem->Load("StPmdUtil");
  gSystem->Load("StPreEclMaker");
  gSystem->Load("StStrangeMuDstMaker");
  gSystem->Load("StMuDSTMaker");

  if(!creatingPhiWgt&&emcMode) {
    gSystem->Load("StTpcDb");
    gSystem->Load("StMcEvent");
    gSystem->Load("StMcEventMaker");
    gSystem->Load("StDaqLib");
    gSystem->Load("libgen_Tables");
    gSystem->Load("libsim_Tables");
    gSystem->Load("libglobal_Tables");
    gSystem->Load("StEmcTriggerMaker");
    gSystem->Load("StEmcUtil");//mine
    gSystem->Load("StEmcRawMaker");
    gSystem->Load("StEmcADCtoEMaker");
    gSystem->Load("StPreEclMaker");
    gSystem->Load("StEpcMaker");
    gSystem->Load("StEmcSimulatorMaker");
    gSystem->Load("StEmcUtil");
    gSystem->Load("StDbBroker");
    gSystem->Load("StDetectorDbMaker");
    gSystem->Load("StDbUtilities");
    gSystem->Load("StEEmcUtil");
    gSystem->Load("StEEmcDbMaker");
    gSystem->Load("St_db_Maker");
    gSystem->Load("StTriggerUtilities");
  }

  // New additions for MTD fixes in Run14
  gSystem->Load("StMagF");
  gSystem->Load("StMtdUtil");
  gSystem->Load("StMtdMatchMaker");
  gSystem->Load("StMtdCalibMaker");

  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StPicoAnaTreeMaker");
  gSystem->Load("StPicoQAMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StPicoElecPurityMaker");

  chain = new StChain();

  //StMuDstMaker *MuDstMaker = new StMuDstMaker(0,0,"",inputFile,"MuDst",100);
  //MuDstMaker->SetStatus("*",0);
  //MuDstMaker->SetStatus("MuEvent",1);
  //MuDstMaker->SetStatus("PrimaryVertices",1);
  //MuDstMaker->SetStatus("PrimaryTracks",1);
  //MuDstMaker->SetStatus("GlobalTracks",1);
  //MuDstMaker->SetStatus("CovGlobTrack",1);
  //MuDstMaker->SetStatus("BTof*",1);
  //MuDstMaker->SetStatus("Emc*",1);
  //MuDstMaker->SetStatus("MTD*",1);


  //StMagFMaker *magfMk = new StMagFMaker; 
  //StMtdMatchMaker *mtdMatchMaker = new StMtdMatchMaker();
  //StMtdCalibMaker *mtdCalibMaker = new StMtdCalibMaker("mtdcalib"); 

  /*
     if(!creatingPhiWgt&&emcMode) {
     St_db_Maker *dbMk = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","StarDb");

  // Endcap database
  StEEmcDbMaker* eemcDb = new StEEmcDbMaker;

  StEmcADCtoEMaker *adc2e = new StEmcADCtoEMaker();
  adc2e->setPrint(false);
  //adc2e->setFillHisto(false);
  //adc2e->setDebug(false); //more histograms
  //adc2e->setSMDRmsCut(0,0);
  adc2e->saveAllStEvent(true);
  //adc2e->setRemoveGhostEvent(false);
  //adc2e->setTallyHist(mTally);
  //adc2e->setDbName("Calibrations/emc/y3");

  StPreEclMaker *pre_ecl=new StPreEclMaker();
  pre_ecl->setPrint(kFALSE);
  StEpcMaker *epc=new StEpcMaker();
  epc->setPrint(kFALSE);

  if(prodMod==1){
  // Trigger simulator
  StTriggerSimuMaker* trigSimu = new StTriggerSimuMaker;
  trigSimu->setMC(false);
  trigSimu->useBemc();
  trigSimu->useEemc();
  trigSimu->useOnlineDB();
  trigSimu->bemc->setConfig(StBemcTriggerSimu::kOffline);
  }
  }
  */

  StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst"); //0-read, 1-write
  picoMaker->setRunNumber(runnumber);
  picoMaker->setProdMode(prodMod); // 0-mb, 1-ht, 4-only e or mu  
  picoMaker->setEmcMode(emcMode); // 0-No EMC, 1-EMC ON
  //picoMaker->SetDebug(1);

  TString inputFileName = inputFile;
  Int_t index = inputFileName.Index("st_");
  TString mInputFileName="";
  for(int i=index;i<(int)inputFileName.Length();i++) {
    mInputFileName.Append(inputFileName(i));
  }

  TString outputFile,outQAFile;
  outQAFile=mInputFileName;
  outQAFile.ReplaceAll("picoDst.root","qa.root");
  StPicoQAMaker *qaMaker = new StPicoQAMaker("ana",picoMaker,outQAFile);

  TString outPurityFile=mInputFileName;
  outPurityFile.ReplaceAll("picoDst.root","purity.root");
  StPicoElecPurityMaker *ePurMaker = new StPicoElecPurityMaker("purity",picoMaker,outPurityFile);
  ePurMaker->setRunMode(1); // 0 - pp (no centrality), 1- AuAu (centrality on)
  if(prodMod==0){
   ePurMaker->addTrigger(450050,4); //0 -BHT0, 1-BHT1, 2-BHT2, 3-BHT3, 4-MB
   ePurMaker->addTrigger(450060,4);
   ePurMaker->addTrigger(450005,4);
   ePurMaker->addTrigger(450009,4);
   ePurMaker->addTrigger(450015,4);
   ePurMaker->addTrigger(450025,4);
  }
  if(prodMod==1){
   ePurMaker->addTrigger(450201,1); //0 -BHT0, 1-BHT1, 2-BHT2, 3-BHT3, 4-MB
   ePurMaker->addTrigger(450211,1);
   ePurMaker->addTrigger(450202,2);
   ePurMaker->addTrigger(450212,2);
   ePurMaker->addTrigger(450203,3);
   ePurMaker->addTrigger(450213,3);
  }
	
  outputFile=mInputFileName;
  outputFile.ReplaceAll("picoDst.root","anaTree.root");

  StPicoAnaTreeMaker *treeMaker = new StPicoAnaTreeMaker(1,outputFile,picoMaker);
  treeMaker->setTriggerSelection(prodMod); //0-mb, 1-ht, 2-mtd
  if(prodMod==0){
    treeMaker->setVzCut(-8,8);
    treeMaker->setVzDiffCut(-4,4);
    treeMaker->setInputRunList("./runNumberList_run14AuAu200mb");
    treeMaker->setInputRecenterFile("./recenter_correction.root");
    treeMaker->setPhoEPairMassCut(0.2);
    treeMaker->addTrigger(450050); // vpdmb-5-p-nobsmd-hlt (production_mid_2014, production_low_2014)
    treeMaker->addTrigger(450060); // vpdmb-5-p-nobsmd-hlt (production_mid_2014, production_low_2014)
    treeMaker->addTrigger(450009); // vPDMB-5-p-nobsmd-ssd-hlt (production_mid_2014, production_low_2014)
    treeMaker->addTrigger(450005); // vpdmb-5-p-nobsmd (production_2014)
    treeMaker->addTrigger(450015); // vpdmb-5-p-nobsmd (production_2014, production_mid_2014, production_low_2014)
    treeMaker->addTrigger(450025); // vpdmb-5-p-nobsmd (production_mid_2014, production_low_2014)
    treeMaker->addTrigger(450010); // VPDMB-30 (production_2014, production_mid_2014, production_low_2014)
    treeMaker->addTrigger(450020); // VPDMB-30 (production_2014, production_mid_2014, production_low_2014)
    treeMaker->addTrigger(450011); // MB-mon (production_2014, production_mid_2014, production_low_2014)
    treeMaker->addTrigger(450021); // MB-mon (production_2014, production_mid_2014, production_low_2014)
  }
  if(prodMod==1){
    treeMaker->setVzCut(-100,100);
    //treeMaker->setDoEvtPlane(false); //default is true
    if(prodType==0){ // prod low and mid
      treeMaker->setInputRunList("./runNumberList_run14AuAu200mb");
      treeMaker->setInputRecenterFile("./recenter_correction.root");
      treeMaker->setMaxRunId(1700);
    }
    if(prodType==1){ // prod high
      treeMaker->setInputRunList("./runNumberList_run14AuAu200ht_high");
      treeMaker->setInputRecenterFile("./recenter_correction_ht_high.root");
      treeMaker->setMaxRunId(1000); //need to check
    }
    treeMaker->setPhoEPairMassCut(0.24);
    treeMaker->setSaveHadron(true);
    treeMaker->addTrigger(450201);    // BHT1*VPDMB-30 (production_2014, production_mid_2014, production_low_2014)
    treeMaker->addTrigger(450211);    // BHT1*VPDMB-30 (production_2014, production_mid_2014, production_low_2014)
    treeMaker->addTrigger(450202);    // BHT2*VPDMB-30 (production_2014, production_mid_2014, production_low_2014)
    treeMaker->addTrigger(450212);    // BHT2*VPDMB-30 (production_2014, production_mid_2014, production_low_2014)
    treeMaker->addTrigger(450203);    // BHT3 (production_2014, production_mid_2014, production_low_2014)
    treeMaker->addTrigger(450213);   // BHT3 (production_2014, production_mid_2014, production_low_2014)
  }
  if(prodMod==2){
    treeMaker->setVzCut(-100,100);
    //treeMaker->setDoEvtPlane(false); //default is true
    //treeMaker->setInputRunList("./runNumberList_run14AuAu200");
    //treeMaker->setInputRecenterFile("./recenter_correction.root");
  }
  //treeMaker->setDoCalcRecenter(false);  //default is false
  treeMaker->setPartEnSigECut(-3.5,3);
  treeMaker->setPhoEPairDcaCut(1);

  chain->Init();
  cout<<"chain->Init();"<<endl;

  int total = picoMaker->chain()->GetEntries();
	cout << " Total entries = " << total << endl;
	if(nEvents>total) nEvents = total;

  int counts = 0;
  for (Int_t i=0; i<nEvents; i++){
    if(i%100==0)
      cout << "Working on eventNumber " << i << endl;

    chain->Clear();
    int iret = chain->Make(i);
    if (iret) { cout << "Bad return code!" << iret << endl; break;}
    counts++;
  }

  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!"<< endl;
  cout << "****************************************** " << endl;
  chain->Finish();
  cout << "****************************************** " << endl;
  cout << "total number of events  " << counts << endl;
  cout << "****************************************** " << endl;

  //	delete chain;


}
