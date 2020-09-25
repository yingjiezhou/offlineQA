/* **************************************************
 *
 *  Authors: Yuanjing Ji
 Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StEvent/StDcaGeometry.h"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoMtdPidTraits.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoDstarMixedMaker.h"
#include "StAnaCuts.h"
#include "StMemStat.h"
#include "calmean.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#ifndef DEBUG
#define DEBUG 1
#endif
ClassImp(StPicoDstarMixedMaker)
  StPicoDstarMixedMaker::StPicoDstarMixedMaker(char const * name, TString const inputFilesList, TString const outFileBaseName, StPicoDstMaker* picoDstMaker):
    StMaker(name), mPicoDstMaker(picoDstMaker),
    mInputFilesList(inputFilesList), mOutFileBaseName(outFileBaseName),mRunbyRunQA(true)
{}

Int_t StPicoDstarMixedMaker::Init()
{
  mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
  // -------------- USER VARIABLES -------------------------
  mFile = new TFile(mOutFileBaseName+".QA.root", "RECREATE");
  //mFile_RunID = new TFile(mOutFileBaseName+".RunID.root","RECREATE");
  //initialize trees
  initHists();

  return kStOK;
}
//-----------------------------------------------------------------------------
StPicoDstarMixedMaker::~StPicoDstarMixedMaker()
{}
//-----------------------------------------------------------------------------
void StPicoDstarMixedMaker::initHists(){
  //int totalNum = 170;//9.2GeV 
  //int totalNum = 1950;//11.5GeV 
  //int totalNum = 29;//19.5GeV fixed target
  int totalNum = 15;//26.5GeV fixed target
  
  //char name_RunID[100];

  if (mRunbyRunQA){
    ifstream readnum;
    readnum.open(mRunNumList);
    int tmp=0;
    if (DEBUG) cout<<"start initial run number..."<<endl;
    for (int i=0;i<totalNum;i++){
      readnum>>tmp;
      runnum.insert(pair<int,int>(tmp,i));
      
      // sprintf(name_RunID,"hinvBetavsP_%d",tmp);
      // hinvBetavsP_RunID[i] = new TH2F(name_RunID,"#frac{1}{#beta} vs P;P(GeV/c);#frac{1}{#beta}",300,0,3,200,0.5,1.5);
       //memset(name_RunID,0,sizeof(name_RunID));

      if (DEBUG) cout <<"run number : " << tmp <<" id :" <<runnum[tmp] <<endl;
    }
    readnum.close();
  }
  // event level QA
  hevt = new TH1D("hevt","hevt",totalNum,0,totalNum);
  hevtcut = new TH1D("hevtcut","hevtcut",totalNum,0,totalNum);
  hevtbadcut = new TH1D("hevtbadcut","Events after remove bad run;Run;Counts",totalNum,0,totalNum);
  hpassevtcut = new TH1D("hpassevtcut","pass event cut",6  , -0.5 , 5.5 );
  //run by run QA
  if (mRunbyRunQA){ 
    pVpdVz = new TProfile("VpdVz","VpdVz vs runId;runId;VpdVz(cm)",totalNum,0,totalNum);
    pVzVpdVz = new TProfile("VzVpdVz","VzVpdVz vs runId;runId;VpdVz-Vz(cm)",totalNum,0,totalNum);
    pRefmult = new TProfile("Refmult","Refmult vs runId;runId;Refmult",totalNum,0,totalNum);
    pVpdHitEast = new TProfile("pVpdHitEast","pVpdHitEast;runId;pVpdHitEast",totalNum,0,totalNum);
    pVpdHitWest = new TProfile("pVpdHitWest","pVpdHitWest;runId;pVpdHitWest",totalNum,0,totalNum);
    pVz  = new TProfile("pVz","pVz;runId;pVz(cm)",totalNum,0,totalNum);
    pVx = new TProfile("pVx","pVx;runId;pVx(cm)",totalNum,0,totalNum);
    pVy = new TProfile("pVy","pVy;runId;pVy(cm)",totalNum,0,totalNum);
    pVr = new TProfile("pVr","pVr;runId;pVr(cm)",totalNum,0,totalNum);
    //track level QA
    pTof = new TProfile("Tof","1/#beta vs runId;runId;1/#beta",totalNum,0,totalNum);
    pDedx = new TProfile("Dedx","dEdx vs runId;runId;dEdx",totalNum,0,totalNum);
    //pRMSDedx = new TProfile("RMSDedx","RMSdEdx vs runId;runId;RMSdEdx",totalNum,0,totalNum);
    pgDCA = new TProfile("gDCA","gDCA vs runId;runId;global DCA(cm)",totalNum,0,totalNum);
    pgPt = new TProfile("gPt","global Pt vs runId;runId;global p_{T}(GeV/c)",totalNum,0,totalNum);
    pgPhi = new TProfile("gPhi","global Phi vs runId;runId;gPhi",totalNum,0,totalNum);
    pgEta = new TProfile("gEta","global Eta vs runId;runId;Eta",totalNum,0,totalNum);
    pNFits = new TProfile("NFits","NHitsFit vs runId;runId;nHitsFit",totalNum,0,totalNum);
    ppPt = new TProfile("pPt","primary Pt vs runId;runId;primary p_{T}(GeV/c)",totalNum,0,totalNum);
    ppEta = new TProfile("pEta","primary Eta vs runId;runId;pEta",totalNum,0,totalNum);
    ppPhi = new TProfile("pPhi","primary Phi vs runId;runId;pPhi",totalNum,0,totalNum);
  }
   //event QA
    hVxVyVz = new TH3F("hVxVyVz","VxVyVz;Vx(cm);Vy(cm);Vz(cm)",100,-0.5,0.5,100,-0.5,0.5,240,-60,60);
    hVz = new TH1F("hVz","Vz;Vz(cm);Counts",1600,-800,800);
    hVpdVz = new TH1F("hVpdVz","VpdVz;VpdVz(cm);Counts",3000,-1500,1500);
    hVr = new TH1F("hVr","Vr;Vr(cm);Counts",1000,0,10);
    hVzVpdVz = new TH1F("hVzVpdVz","Vz-VpdVz(cm)",32000,-1600,1600);
    h_Vx_Vy = new TH2F("h_Vx_Vy","Vertex XY",1200,-6,6,1200,-6,6);
    hnEvsEtavsVz = new TH2F("hnEvsEtavsVz","nElectron;#eta;Vz(cm)",40,-1.5,1.5,1600,-800,800);
    hnEvsPhivsVz = new TH2F("hnEvsPhivsVz","nElectron;#phi;Vz(cm)",100,-1*TMath::Pi(),TMath::Pi(),1600,-800,800);
    //hnTofMulvsRef = new TH2F("hnTofMulvsRef","hnTofMul vs Ref;btofTrayMultiplicity;refMult",2000,0,2000,900,0,900); 
    hnTofMulvsRef = new TH2F("hnTofMulvsRef","RefMul vs nTofMul;RefMul;btofTrayMultiplicity",900,0,900,2000,0,2000); 
    //hnTofMatvsRef= new TH2F("hnTofMatvsRef","nTofmatch VS Refmult;nTofMatch;refMult",900,0,900,900,0,900);
    hnTofMatvsRef= new TH2F("hnTofMatvsRef","RefMul VS nTofmatch;RefMul;nTofMatch",900,0,900,900,0,900);
    hnTofHitvsRef= new TH2F("hnTofHitvsRef","hnTofHit vs Ref;nTofHits;refMult",900,0,900,900,0,900);
    hrefmult = new TH1F("hrefmult","hrefmult",700,0,700);

   //tracl level QA
    hnHitsFit = new TH1F("hnHitsFit","nHitsFit;nHitsFit",90,0,90);
    h_p_E0 = new TH1F("h_p_E0","#frac{p}{E0};#frac{p}{E0}",40,0,4);
    ph_p_E0 = new TProfile("ph_p_E0","#frac{p}{E0};RunId;#frac{p}{E0}",totalNum,0,totalNum);
    h_bemcdz = new TH1F("h_bemcdz","bemcdz;bemcdz",400,-200,20);
    ph_bemcdz = new TProfile("ph_bemcdz","bemcdz;RunId;bemcdz",totalNum,0,totalNum);
    h_bemcDphi = new TH1F("h_bemcDphi","bemcDphi;bemcDphi",100,-0.5,0.5);
    ph_bemcDphi = new TProfile("ph_bemcDphi","bemcDphi;RunId;bemcDphi",totalNum,0,totalNum);
    h_nSMDphi = new TH1F("h_nSMDphi","nSMDphi;nSMDphi",10,0,10);
    ph_nSMDphi = new TProfile("ph_nSMDphi","nSMDphi;RunId;nSMDphi",totalNum,0,totalNum);
    h_nSMDeta = new TH1F("h_nSMDeta","nSMDeta;nSMDeta",10,0,10);
    ph_nSMDeta = new TProfile("ph_nSMDeta","nSMDeta;RunId;nSMDeta",totalNum,0,totalNum);
    h_btowPhiDz = new TH1F("h_btowPhiDz","btowPhiDz;btowPhiDz",500,-50,50);
    ph_btowPhiDz = new TProfile("ph_btowPhiDz","btowPhiDz;RunId;btowPhiDz",totalNum,0,totalNum);
    h_btowEtaDz = new TH1F("h_btowEtaDz","btowEtaDz;btowEtaDz",500,-50,50);
    ph_btowEtaDz = new TProfile("ph_btowEtaDz","btowEtaDz;RunId;btowEtaDz",totalNum,0,totalNum);
    
    h_mtdDeltaY = new TH1F("h_mtdDeltaY","#Delta Y(mtd);#Delta Y",1000,-50,50);
    ph_mtdDeltaY = new TProfile("ph_mtdDeltaY","#Delta Y(mtd);RunId;#Delta Y",totalNum,0,totalNum);
    h_mtdDeltaZ = new TH1F("h_mtdDeltaZ","#delta Z(mtd);#delta Z",1000,-50,50);
    ph_mtdDeltaZ = new TProfile("ph_mtdDeltaZ","#Delta Z(mtd);RunId;#Delta Z",totalNum,0,totalNum);
    h_mtdDeltaTOF = new TH1F("h_mtdDeltaTOF","#delta TOF(mtd);#delta TOF",20000,-2000,2000);
    ph_mtdDeltaTOF = new TProfile("ph_mtdDeltaTOF","#Delta TOF(mtd);RunId;#Delta TOF",totalNum,0,totalNum);
    // hpDca = new TH1F("hpDca","pDca",50,0,5);
    hgDca = new TH1F("hgDca","gDca",100,0,10);
    hinvBetavsP = new TH2F("hinvBetavsP","#frac{1}{#beta} vs p;p(GeV/c);#frac{1}{#beta}",300,0,3,200,0.5,1.5);
    // hinvBetavsY;
    hdEdx = new TH2F("hdEdx","dEdx vs p*charge;p*charge(GeV/c);dEdx",200,-2,2,400,0,25);
    h_mTpc = new TH2F("h_mTpc","mass vs p*charge;p*charge(GeV/c);mass(GeV/c^{2})",200,-2,2,100,0,2);
    hNsigEvsinvBeta = new TH3F("hNsigEvsinvBeta","nSigmaE vs 1/#beta;nSigmaE;1/#beta;p",200,-10,10,100,0.8,1.2,100,0.15,2.5);
    hpt = new TH1F("hpt","hpt;p_{T}(GeV/c)",240,0,12);
    hGpt = new TH1F("hGpt","hGpt;global p_{T}(GeV/c)",240,0,12);
    hEta = new TH1F("hEta","Eta;#eta",160,-15,1);
    hPhi = new TH1F("hPhi","Phi;#phi",80,-4,4);
    hBadTofId = new TH3F("hBadTofId","hBadTofId;tray;module;cell",125,-0.5,124.5,33,-0.5,32.5,7,-0.5,6.5);

    //invariant mass electron
    //hMeeCount = new TH1F("hMee","hMee;Count;Mee(GeV/c^{2})",400,0,4);
    hMeeCount = new TH1F("hMee","hMee;Mee(GeV/c^{2})",40000,0,4);
    hMeeCount_like1 = new TH1F("hMee_like1","hMee like sign electron;Mee(GeV/c^{2})",40000,0,4);
    hMeeCount_like2 = new TH1F("hMee_like2","hMee like sign positron;Mee(GeV/c^{2})",40000,0,4);
    hMeeCountPt = new TH2F("hMeePt","hMee vs p_{T};Mee(GeV/c^{2});p_{T}",40000,0,4,200,0,10);
    hMeeCountPt_like1 = new TH2F("hMeePt_like1","hMee vs p_{T} like sign electron;Mee(GeV/c^{2});p_{T}",40000,0,4,200,0,10);
    hMeeCountPt_like2 = new TH2F("hMeePt_like2","hMee vs p_{T} like sign positron;Mee(GeV/c^{2});p_{T}",40000,0,4,200,0,10);
    h_nSigmaElectron_P = new TH2F("h_nSigmaElectron_P","nSigmaElectron_P;P(GeV/c);nSigmaElectron",300,0,3,100,-5,5);
    h_nSigmaElectron_P_tpc = new TH2F("h_nSigmaElectron_P_tpc","nSigmaElectron_P;P(GeV/c);nSigmaElectron",300,0,3,100,-5,5);

    //tof module id
    /*ModuleId_1 = new TH1F("ModuleId 1","0.8<1/#beta<0.9 0.4<P;ModuleId",40,0,40);
    TofId_1 = new TH1F("TofId 1","0.8<1/#beta<0.9 0.4<P;TofId",23100,0,23100);
    TrayId_1 = new TH1F("TrayId 1","0.8<1/#beta<0.9 0.4<P;TrayId",130,0,130);
    ModuleId_2 = new TH1F("ModuleId 2","0.82<1/#beta<0.9 0.4<P;ModuleId",40,0,40);
    TofId_2 = new TH1F("TofId 2","0.8<1/#beta<0.9 0.4<P;TofId",23100,0,23100);
    ModuleId_3 = new TH1F("ModuleId 3","0.82<1/#beta<0.88 0.4<P;ModuleId",40,0,40);
    TofId_3 = new TH1F("TofId 3","0.8<1/#beta<0.9 0.4<P;TofId",23100,0,23100);
    ModuleId_4 = new TH1F("ModuleId 4","0.84<1/#beta<0.88 0.4<P;ModuleId",40,0,40);
    TofId_4 = new TH1F("TofId 4","0.8<1/#beta<0.9 0.4<P;TofId",23100,0,23100);
    ModuleId_5 = new TH1F("ModuleId 5","0.84<1/#beta<0.9 & 0.4<P;ModuleId",40,0,40);
    TofId_5 = new TH1F("TofId 5","0.8<1/#beta<0.9 0.4<P;TofId",23100,0,23100);*/
}
//-----------------------------------------------------------------------------
Int_t StPicoDstarMixedMaker::Finish()
{
mFile->cd();
  //write the hists
 //event QA
    hVxVyVz->Write();
    hVz->Write();
    hVpdVz->Write();
    hVr->Write();
    hVzVpdVz->Write();
    hnEvsEtavsVz->Write();
    hnEvsPhivsVz->Write();
    hnTofMulvsRef->Write(); 
    hnTofMatvsRef->Write();
    hnTofHitvsRef->Write();

    hevt->Write();
    h_Vx_Vy->Write();
    hevtcut->Write();
    hevtbadcut->Write();
    hpassevtcut->Write(); 
    hrefmult->Write();
    hNsigEvsinvBeta->Write();
   //tracl level QA
    hnHitsFit->Write();
    // hpDca->Write();
    hgDca->Write();
    hinvBetavsP->Write();
    // hinvBetavsY->Write();
    hdEdx->Write();
    h_mTpc->Write();
    hpt->Write();
    hGpt->Write();
    hEta->Write();
    hPhi->Write();
    hBadTofId->Write();
  hMeeCount->Write();
  hMeeCount_like1->Write();
  hMeeCount_like2->Write();
  hMeeCountPt->Write();
  hMeeCountPt_like1->Write();
  hMeeCountPt_like2->Write();
  h_nSigmaElectron_P->Write();
  h_nSigmaElectron_P_tpc->Write();
  h_p_E0->Write();
  h_bemcdz->Write();
  h_bemcDphi->Write();
  h_nSMDphi->Write();
  h_nSMDeta->Write();
  h_btowPhiDz->Write();
  h_btowEtaDz->Write();
  h_mtdDeltaY->Write();
  h_mtdDeltaZ->Write();
  h_mtdDeltaTOF->Write();
  ph_p_E0->Write();
  ph_bemcdz->Write();
  ph_bemcDphi->Write();
  ph_nSMDphi->Write();
  ph_nSMDeta->Write();
  ph_btowPhiDz->Write();
  ph_btowEtaDz->Write();
  ph_mtdDeltaY->Write();
  ph_mtdDeltaZ->Write();
  ph_mtdDeltaTOF->Write();

    /*ModuleId_1->Write();
    TofId_1->Write();
    TrayId_1->Write();
    ModuleId_2->Write();
    TofId_2->Write();
    ModuleId_3->Write();
    TofId_3->Write();
    ModuleId_4->Write();
    TofId_4->Write();
    ModuleId_5->Write();
    TofId_5->Write();*/
  if (mRunbyRunQA) {
    pVpdVz->Write();
    pVzVpdVz->Write();
    pRefmult->Write();
    pVpdHitEast->Write();
    pVpdHitWest->Write();
    pVz->Write();
    pVx->Write();
    pVy->Write();
    pVr->Write();
    pTof->Write(); 
    pDedx->Write();
    pgDCA->Write();
    pgPt->Write();
    pgPhi->Write();
    pgEta->Write();
    pNFits->Write();
    ppPt->Write();
    ppEta->Write();
    ppPhi->Write();
  }
  mFile->Close();

  /*mFile_RunID->cd();
  int x_RunID=0;
  for(x_RunID=0;x_RunID<78;x_RunID++)
  {
    hinvBetavsP_RunID[x_RunID]->Write();
  }
  mFile_RunID->Close();*/

  return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoDstarMixedMaker::Make()
{
  ParticleInfo_Electron particleinfo;
  vector<ParticleInfo_Electron> electroninfo;
  vector<ParticleInfo_Electron> positroninfo;
  // StMemStat mem;
  if (!mPicoDstMaker)
  {
    LOG_WARN << " StPicoDstarMixedMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }
  StPicoDst const* picoDst = mPicoDstMaker->picoDst();
  if (!picoDst)
  {
    LOG_WARN << "StPicoDstarMixedMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  // -------------- USER ANALYSIS -------------------------
  StPicoEvent const * picoEvent = picoDst->event();
  //trigger
//skf  if (!isGoodTrigger(picoEvent)) return 0;    
  mRunId = picoEvent->runId();
  hevt->Fill(runnum[mRunId]);

    TVector3 pVtx = picoEvent->primaryVertex();
  if (mRunbyRunQA && isGoodQaEvent(picoEvent) ){ 
    //primary vertex
    // StThreeVectorF pVtx = picoEvent->primaryVertex();
    mVx = pVtx.x();
    mVy = pVtx.y();
    mVz = pVtx.z();
    mVr = sqrt(mVx*mVx+mVy*mVy);
    mVpdVz = picoEvent->vzVpd();
    mRefmult = picoEvent->refMult();
    mVpdHitEast = picoEvent->nVpdHitsEast();
    mVpdHitWest = picoEvent->nVpdHitsWest();
    if (DEBUG) cout<<"start filling "<<mRunId<<" "<<runnum[mRunId]<<endl; 
    //fill event level profile
    pVpdVz->Fill(runnum[mRunId],mVpdVz);
    pVzVpdVz->Fill(runnum[mRunId],mVpdVz-mVz);
    pRefmult->Fill(runnum[mRunId],mRefmult);
    pVpdHitEast->Fill(runnum[mRunId],mVpdHitEast);
    pVpdHitWest->Fill(runnum[mRunId],mVpdHitWest);
    pVx->Fill(runnum[mRunId],mVx);
    pVy->Fill(runnum[mRunId],mVy);
    pVz->Fill(runnum[mRunId],mVz);
    pVr->Fill(runnum[mRunId],mVr);
    //track level 
    

    electroninfo.clear();
    positroninfo.clear();  
  


    int nTracks = picoDst->numberOfTracks();
    //if (DEBUG)  cout << nTracks <<endl;
    //global
    for (int itrack=0;itrack<nTracks;itrack++){
      StPicoTrack* trk = picoDst->track(itrack);
      bool goodQAtrack = isGoodQaTrack(trk); 
      if (!goodQAtrack) continue;
      if (!(fabs(trk->gDCA(pVtx.x(),pVtx.y(),pVtx.z()))<anaCuts::qaDca)) continue;
      bool isprimary = trk->isPrimary();
      double ptot = trk->gMom(pVtx, picoEvent->bField()).Mag();
      float beta = getTofBeta(trk);
      bool tofmatch =( beta!=std::numeric_limits<float>::quiet_NaN() )&& beta>0;
      bool tofQaPion = false;
      if (tofmatch) {
        float beta_pi = ptot / sqrt(ptot * ptot + M_PION_PLUS * M_PION_PLUS);
        tofQaPion = fabs(1. / beta - 1. / beta_pi) < anaCuts::qaTofPion;
      }
      bool tpcQaPion = fabs(trk->nSigmaPion()) < anaCuts::qaTpcPion;
      //global
      //runbyrunQA
      //pid performance 
      if (tofQaPion && ptot<1) pTof->Fill(runnum[mRunId],(1./beta)); 
      if (tpcQaPion && ptot<0.5) {
        pDedx->Fill(runnum[mRunId],trk->dEdx());
      }
      //global tracking performance
      pgPt->Fill(runnum[mRunId],trk->gMom().Perp());
      pgDCA->Fill(runnum[mRunId],trk->gDCA(pVtx.x(),pVtx.y(),pVtx.z()));
      pgPhi->Fill(runnum[mRunId],trk->gMom().Phi());
      pgEta->Fill(runnum[mRunId],trk->gMom().Eta());
      pNFits->Fill(runnum[mRunId],trk->nHitsFit()); 

      //primary 
      if (isprimary){
        ppPt->Fill(runnum[mRunId],trk->pMom().Perp());
        ppEta->Fill(runnum[mRunId],trk->pMom().Eta());
        ppPhi->Fill(runnum[mRunId],trk->pMom().Phi());
      int pbemcId = trk->bemcPidTraitsIndex();
      if(pbemcId>=0){
         StPicoBEmcPidTraits * pbemctrait = picoDst->bemcPidTraits(pbemcId);
         if(pbemctrait){
            float pnSMDphi = pbemctrait->bemcSmdNPhi();
            ph_nSMDphi->Fill(runnum[mRunId],pnSMDphi);
            float pnSMDeta = pbemctrait->bemcSmdNEta();
            ph_nSMDeta->Fill(runnum[mRunId],pnSMDeta);
            float pbtowPhiDz = pbemctrait->btowPhiDist();
            ph_btowPhiDz->Fill(runnum[mRunId],pbtowPhiDz);
            float pbtowEtaDz = pbemctrait->btowEtaDist();
            ph_btowEtaDz->Fill(runnum[mRunId],pbtowEtaDz);
            float pbemcdz = pbemctrait->bemcZDist();
            ph_bemcdz->Fill(runnum[mRunId],pbemcdz);
            float pbemcDphi = pbemctrait->bemcPhiDist();
            ph_bemcDphi->Fill(runnum[mRunId],pbemcDphi);
            float pE0 = pbemctrait->bemcE0(); 
            ph_p_E0->Fill(runnum[mRunId],trk->gMom().Mag()*1.0/pE0); 
           }
      }
    int pmtdId = trk->mtdPidTraitsIndex();
    if(pmtdId>=0){
       StPicoMtdPidTraits * pmtdtrait = picoDst->mtdPidTraits(pmtdId);
       if(pmtdtrait){
          float pmtdDeltaY = pmtdtrait->deltaY();
          ph_mtdDeltaY->Fill(runnum[mRunId],pmtdDeltaY);
          float pmtdDeltaZ = pmtdtrait->deltaZ();
          ph_mtdDeltaZ->Fill(runnum[mRunId],pmtdDeltaZ);
          float pmtdDeltaTOF = pmtdtrait->deltaTimeOfFlight();
          ph_mtdDeltaTOF->Fill(runnum[mRunId],pmtdDeltaTOF);
       }
    }
      }
    }
  } //runbyrun QA 

  //event and track level QA
  hpassevtcut->Fill(0);
  if (!isBadrun(mRunId)){
  hpassevtcut->Fill(1); //bad run list
  // bool vzcut = fabs(pVtx.z()) < 30;
  bool vzcut = fabs(pVtx.z()) < 60;
  bool verrcut = !(fabs(pVtx.x()) < anaCuts::Verror && fabs(pVtx.y()) < anaCuts::Verror &&
        fabs(pVtx.z()) < anaCuts::Verror);
  bool vrcut =  sqrt(TMath::Power(pVtx.x(), 2) + TMath::Power(pVtx.y(), 2)) <=  anaCuts::Vr ;
  // bool vpdvzcut = fabs(pVtx.z() - picoEvent->vzVpd()) < 3;
  bool vpdvzcut = true;
  if (vzcut) hpassevtcut->Fill(2);
  if (vzcut &&  vrcut) hpassevtcut->Fill(3);
  // if (vzcut && vrcut  &&  vpdvzcut ) hpassevtcut->Fill(4);
  if (vzcut && vrcut  &&  vpdvzcut && verrcut ) hpassevtcut->Fill(4);
  bool refusepileup = picoEvent->refMult()<picoEvent->btofTrayMultiplicity()*0.36+45;
  bool refusebadtof = picoEvent->refMult()>picoEvent->btofTrayMultiplicity()*0.28-115;
  bool passCentralityCut = refusepileup && refusebadtof  && verrcut && vrcut && fabs(pVtx.z()) < 10; 
  if (passCentralityCut) hrefmult->Fill(picoEvent->refMult());
  if (isGoodEvent(picoEvent)){
    // StThreeVectorF pVtx = picoEvent->primaryVertex();
    TVector3 pVtx = picoEvent->primaryVertex();
    mVx = pVtx.x();
    mVy = pVtx.y();
    mVz = pVtx.z();
    h_Vx_Vy->Fill(mVx,mVy);
    mVpdVz = picoEvent->vzVpd();
    hevtcut->Fill(runnum[mRunId]);
    hVz->Fill(mVz);
    hVpdVz->Fill(mVpdVz);
    hVxVyVz->Fill(mVx,mVy,mVz);
    hVr->Fill(sqrt(mVy*mVy+mVx*mVx));
    hVzVpdVz->Fill(mVpdVz-mVz);  

    //hnTofMulvsRef->Fill(picoEvent->btofTrayMultiplicity(),picoEvent->refMult());  
    hnTofMulvsRef->Fill(picoEvent->refMult(),picoEvent->btofTrayMultiplicity());  
    //hnTofMatvsRef->Fill(picoEvent->nBTOFMatch(),picoEvent->refMult());  
    hnTofMatvsRef->Fill(picoEvent->refMult(),picoEvent->nBTOFMatch());  
    double ntofhits = 0;
//    int ntrack_tof_hits =0; 
    int nTracks = picoDst->numberOfTracks();
    for (int itrack=0;itrack<nTracks;itrack++){
      StPicoTrack* trk = picoDst->track(itrack);
      hgDca->Fill(trk->gDCA(mVx,mVy,mVz));
      bool isprimary = trk->isPrimary();
      bool goodtrack = isGoodTrack(trk,trk->gDCA(mVx,mVy,mVz));
      if (!goodtrack) continue;
      if (!isprimary) continue;
      // StThreeVectorF mom = trk->pMom(); 
      int bemcId = trk->bemcPidTraitsIndex();
      if(bemcId>=0){
         StPicoBEmcPidTraits * bemctrait = picoDst->bemcPidTraits(bemcId);
         if(bemctrait){
            float nSMDphi = bemctrait->bemcSmdNPhi();
            h_nSMDphi->Fill(nSMDphi);
            float nSMDeta = bemctrait->bemcSmdNEta();
            h_nSMDeta->Fill(nSMDeta);
            float btowPhiDz = bemctrait->btowPhiDist();
            h_btowPhiDz->Fill(btowPhiDz);
            float btowEtaDz = bemctrait->btowEtaDist();
            h_btowEtaDz->Fill(btowEtaDz);
            float bemcdz = bemctrait->bemcZDist();
            h_bemcdz->Fill(bemcdz);
            float bemcDphi = bemctrait->bemcPhiDist();
            h_bemcDphi->Fill(bemcDphi);
            float E0 = bemctrait->bemcE0(); 
            h_p_E0->Fill(trk->gMom().Mag()*1.0/E0); 
           }
       
      }
     
    int mtdId = trk->mtdPidTraitsIndex();
    if(mtdId>=0){
       StPicoMtdPidTraits * mtdtrait = picoDst->mtdPidTraits(mtdId);
       if(mtdtrait){
          float mtdDeltaY = mtdtrait->deltaY();
          h_mtdDeltaY->Fill(mtdDeltaY);
          float mtdDeltaZ = mtdtrait->deltaZ();
          h_mtdDeltaZ->Fill(mtdDeltaZ);
          float mtdDeltaTOF = mtdtrait->deltaTimeOfFlight();
          h_mtdDeltaTOF->Fill(mtdDeltaTOF);
          //cout<<"mtdDeltaTOF: "<<mtdDeltaTOF<<endl;
       }
    }
      TVector3 mom = trk->pMom(); 
      hpt->Fill(mom.Perp());
      hGpt->Fill(trk->gMom().Perp());
      hPhi->Fill(mom.Phi());
      hEta->Fill(mom.Eta());
      // hpDca->Fill(trk->pDca(mVx,mVy,mVz));
      hnHitsFit->Fill(trk->nHitsFit());
      double beta = getTofBeta(trk);
      bool tofmatch = (beta!=std::numeric_limits<float>::quiet_NaN()) && beta>0;
      //choose inclusive electron
      // bool isTPCElectron =  trk->nSigmaElectron()<2 && trk->nSigmaElectron()>0.75;
      bool isTPCElectron=0;
      if (mom.Mag()>0.8) isTPCElectron =  trk->nSigmaElectron()<2 && trk->nSigmaElectron()>-0.75;
      else isTPCElectron = trk->nSigmaElectron()<2 && trk->nSigmaElectron()>(3*mom.Mag()-3.15);
      bool isTOFElectron = tofmatch?fabs(1./beta-1.)<0.025:false;
  
       h_nSigmaElectron_P_tpc->Fill(mom.Mag(),trk->nSigmaElectron());      

        if (isTOFElectron && isTPCElectron) {
        hnEvsEtavsVz->Fill(mom.Eta(),mVz); 
        hnEvsPhivsVz->Fill(mom.Phi(),mVz);
        
        if(trk->charge()<0 && tofmatch)
         {
            particleinfo.charge = trk->charge();
//            cout<<"charge: "<<trk->charge()<<endl;
            particleinfo.pt = mom.Perp();
            particleinfo.eta = mom.Eta();
            particleinfo.phi = mom.Phi();
            particleinfo.p = mom.Mag();
            particleinfo.nSigmaPi = trk->nSigmaPion();
            particleinfo.beta = beta;
            particleinfo.energy = sqrt(pow(M_electron,2.0)+pow(mom.Mag(),2.0));
            particleinfo.p1 = mom.X();
            particleinfo.p2 = mom.Y();
            particleinfo.p3 = mom.Z();
            electroninfo.push_back(particleinfo);
//            cout<<"debug01"<<endl;
            /* current_eMinus[current_nEMinus].SetPx(mom.x());
             current_eMinus[current_nEMinus].SetPy(mom.y());
             current_eMinus[current_nEMinus].SetPz(mom.z());
             current_eMinus[current_nEMinus].SetE(sqrt(pow(M_electron,2.0)+pow(mom.Mag(),2.0)));
             current_nEMinus++;*/
         }     

        if(trk->charge()>0 && tofmatch)
         {
            particleinfo.charge = trk->charge();
//            cout<<"charge: "<<trk->charge()<<endl;
            particleinfo.pt = mom.Perp();
            particleinfo.eta = mom.Eta();
            particleinfo.phi = mom.Phi();
            particleinfo.p = mom.Mag();
            particleinfo.nSigmaPi = trk->nSigmaPion();
            particleinfo.beta = beta;
            particleinfo.energy = sqrt(pow(M_electron,2.0)+pow(mom.Mag(),2.0));
            particleinfo.p1 = mom.X();
            particleinfo.p2 = mom.Y();
            particleinfo.p3 = mom.Z();
            positroninfo.push_back(particleinfo);
//            cout<<"debug02"<<endl;
             /*current_ePlus[current_nEPlus].SetPx(mom.x());
             current_ePlus[current_nEPlus].SetPy(mom.y());
             current_ePlus[current_nEPlus].SetPz(mom.z());
             current_ePlus[current_nEPlus].SetE(sqrt(pow(M_electron,2.0)+pow(mom.Mag(),2.0)));
             current_nEPlus++;*/
         }

         //current_nE++;
      }
      
      if (tofmatch) {
        ntofhits++;
        hinvBetavsP->Fill(mom.Mag(),1./beta);
        
        if(fabs(1.0/beta - 1) < 0.025)
        {
         h_nSigmaElectron_P->Fill(mom.Mag(),trk->nSigmaElectron()); 
        }       

       // hinvBetavsP_RunID[runnum[mRunId]]->Fill(mom.Mag(),1./beta);

    //    bool istoftrack = trk->isTofTrack();
    //    if(istoftrack && 0.8<1.0/beta && 1.0/beta<0.9 && mom.Mag()>0.4){
    //        ntrack_tof_hits = trk->numberOfBTofHits();
    //        for(int itrack_tof_hits=0;itrack_tof_hits<ntrack_tof_hits;itrack_tof_hits++)
    //           {
    //            StPicoBTofHit* tofhit = trk->btofHit(itrack_tof_hits);
    //            ModuleId_1->Fill(tofhit->module()); 
    //           }
    //   }

        int index2tof = trk->bTofPidTraitsIndex();
        StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
        if (1./beta<0.88 && mom.Mag()>0.8 && mom.Mag()<1.5) {
          int tofid = tofPid->btofCellId();
        //  TofId_1->Fill(tofid);
          hBadTofId->Fill(tofid/192+1,(tofid%192)/6+1,tofid%6+1);
        }

       /* if(1.0/beta>0.8 && 1.0/beta<0.9 && mom.Mag()>0.4 && mom.Mag()<3)
         {
            int tofid = tofPid->btofCellId();
            TofId_1->Fill(tofid);
            TrayId_1->Fill(tofid/192+1);
            ModuleId_1->Fill((tofid%192)/6+1);
         }

        if(1.0/beta>0.82 && 1.0/beta<0.9 && mom.Mag()>0.4 && mom.Mag()<3)
         {
           int  tofid = tofPid->btofCellId();
          TofId_2->Fill(tofid);
            ModuleId_2->Fill((tofid%192)/6+1);
         }
        if(1.0/beta>0.82 && 1.0/beta<0.88 && mom.Mag()>0.4 && mom.Mag()<3)
         {
           int  tofid = tofPid->btofCellId();
          TofId_3->Fill(tofid);
            ModuleId_3->Fill((tofid%192)/6+1);
         }
        if(1.0/beta>0.84 && 1.0/beta<0.88 && mom.Mag()>0.4 && mom.Mag()<3)
         {
           int   tofid = tofPid->btofCellId();
          TofId_4->Fill(tofid);
            ModuleId_4->Fill((tofid%192)/6+1);
         }
        if(1.0/beta>0.84 && 1.0/beta<0.9 && mom.Mag()>0.4 && mom.Mag()<3)
         {
           int  tofid = tofPid->btofCellId();
          TofId_5->Fill(tofid);
            ModuleId_5->Fill((tofid%192)/6+1);
         }*/
        hNsigEvsinvBeta->Fill(trk->nSigmaElectron(),1./beta,mom.Mag());
      }
      hdEdx->Fill(mom.Mag()*trk->charge(),trk->dEdx());
      h_mTpc->Fill(mom.Mag()*trk->charge(),mom.Mag()*sqrt(1.0-beta*beta)*1.0/beta);
    }
    hnTofHitvsRef->Fill(ntofhits,picoEvent->refMult());
    
      int x=0;
      int y=0;
      int num_electron = electroninfo.size();
      int num_positron = positroninfo.size();
      float inv_mass=0;
      TVector3 momentum_particle;
      TLorentzVector eepair(0,0,0,0);
      TLorentzVector particle1_4V(0,0,0,0);
      TLorentzVector particle2_4V(0,0,0,0);
      for(x=0;x<num_electron;x++)
         {
             particle1_4V.SetPx(electroninfo[x].p1);
             particle1_4V.SetPy(electroninfo[x].p2);
             particle1_4V.SetPz(electroninfo[x].p3);
             particle1_4V.SetE(electroninfo[x].energy);
               for(y=x+1;y<num_electron;y++)
                  {
                    particle2_4V.SetPx(electroninfo[y].p1);
                    particle2_4V.SetPy(electroninfo[y].p2);
                    particle2_4V.SetPz(electroninfo[y].p3);
                    particle2_4V.SetE(electroninfo[y].energy);
                    eepair = particle1_4V + particle2_4V;
                    
                    //if(eepair.Perp()<0.2){hMeeCount_like1->Fill(eepair.M());}
                    hMeeCount_like1->Fill(eepair.M());
                    hMeeCountPt_like1->Fill(eepair.M(),eepair.Perp());
//                    cout<<"debug03"<<endl;
                    //if(eepair.Perp()<=10){hMeelike1_Pt_Cent->Fill(eepair.Perp(),mCentrality,eepair.M());}
                  }
         }
     
      for(x=0;x<num_positron;x++)
         {
             particle1_4V.SetPx(positroninfo[x].p1);
             particle1_4V.SetPy(positroninfo[x].p2);
             particle1_4V.SetPz(positroninfo[x].p3);
             particle1_4V.SetE(positroninfo[x].energy);
               for(y=x+1;y<num_positron;y++)
                  {
                    particle2_4V.SetPx(positroninfo[y].p1);
                    particle2_4V.SetPy(positroninfo[y].p2);
                    particle2_4V.SetPz(positroninfo[y].p3);
                    particle2_4V.SetE(positroninfo[y].energy);
                    eepair = particle1_4V + particle2_4V;
                    //if(eepair.Perp()<0.2){hMeeCount_like2->Fill(eepair.M());}
                    hMeeCount_like2->Fill(eepair.M());
                    hMeeCountPt_like2->Fill(eepair.M(),eepair.Perp());
                    //if(eepair.Perp()<=10){hMeelike2_Pt_Cent->Fill(eepair.Perp(),mCentrality,eepair.M());}
                  }
         }
      for(x=0;x<num_positron;x++)
         {
             particle1_4V.SetPx(positroninfo[x].p1);
             particle1_4V.SetPy(positroninfo[x].p2);
             particle1_4V.SetPz(positroninfo[x].p3);
             particle1_4V.SetE(positroninfo[x].energy);
               for(y=0;y<num_electron;y++)
                  {
                    particle2_4V.SetPx(electroninfo[y].p1);
                    particle2_4V.SetPy(electroninfo[y].p2);
                    particle2_4V.SetPz(electroninfo[y].p3);
                    particle2_4V.SetE(electroninfo[y].energy);
                    eepair = particle1_4V + particle2_4V;
                    //if(eepair.Perp()<0.2){hMeeCount->Fill(eepair.M());}
                    hMeeCount->Fill(eepair.M());
                    hMeeCountPt->Fill(eepair.M(),eepair.Perp());
                    //if(eepair.Perp()<=10){hMee_Pt_Cent->Fill(eepair.Perp(),mCentrality,eepair.M());}
                  }
         }


  } //Good Event
 }
  return kStOK;
}
bool StPicoDstarMixedMaker::isGoodTrigger(StPicoEvent const* const picoEvent) const
{
  for (auto trg : anaCuts::triggers)
  {
    if (picoEvent->isTrigger(trg)) return true;
  }

  return false;
}
bool StPicoDstarMixedMaker::isGoodTrack(StPicoTrack const* trk, float dca) const
{
  // StThreeVectorF const vtx = mPicoDstMaker->picoDst()->event()->primaryVertex();
  return trk->gPt() > anaCuts::GPt && fabs(trk->nHitsFit()) >= anaCuts::NHitsFit && 
    fabs(trk->gMom().Eta())<anaCuts::Eta &&
    fabs(trk->nHitsDedx())>=anaCuts::NHitsDedx && fabs(dca)<=anaCuts::Dca;
    // fabs(trk->nHitsDedx())>=anaCuts::NHitsDedx &&
    // fabs( trk->gDCA(vtx.x() , vtx.y(), vtx.z() )) <= anaCuts::Dca;
}
bool StPicoDstarMixedMaker::isGoodQaTrack(StPicoTrack const* const trk) const
{
  // StThreeVectorF vtx = mPicoDstMaker->picoDst()->event()->primaryVertex();
  return trk->gPt() > anaCuts::qaGPt && fabs(trk->nHitsFit()) >= anaCuts::qaNHitsFit && 
    fabs(trk->gMom().Eta())<anaCuts::qaEta &&
    // fabs(trk->nHitsDedx())>=anaCuts::qaNHitsDedx && fabs(trk->gDCA(vtx.x(),vtx.y(),vtx.z()))<=anaCuts::qaDca;
    fabs(trk->nHitsDedx())>=anaCuts::qaNHitsDedx ;
}
bool StPicoDstarMixedMaker::isGoodQaEvent(StPicoEvent const* const picoEvent) const
{
  TVector3 pVtx = picoEvent->primaryVertex();
  // StThreeVectorF pVtx = picoEvent->primaryVertex();
  return fabs(pVtx.z()) < anaCuts::qavz &&
    // fabs(pVtx.z() - picoEvent->vzVpd()) < anaCuts::qavzVpdVz &&
    !(fabs(pVtx.x()) < anaCuts::qaVerror && fabs(pVtx.y()) < anaCuts::qaVerror &&
        fabs(pVtx.z()) < anaCuts::qaVerror) &&
    sqrt(TMath::Power(pVtx.x(), 2) + TMath::Power(pVtx.y(), 2)) <=  anaCuts::qaVr;
}
bool StPicoDstarMixedMaker::isGoodEvent(StPicoEvent const* const picoEvent) const
{
  TVector3 pVtx = picoEvent->primaryVertex();
  // StThreeVectorF pVtx = picoEvent->primaryVertex();
  return fabs(pVtx.z()) < anaCuts::vz &&
    // fabs(pVtx.z() - picoEvent->vzVpd()) < anaCuts::vzVpdVz &&
    !(fabs(pVtx.x()) < anaCuts::Verror && fabs(pVtx.y()) < anaCuts::Verror &&
        fabs(pVtx.z()) < anaCuts::Verror) &&
    sqrt(TMath::Power(pVtx.x(), 2) + TMath::Power(pVtx.y(), 2)) <=  anaCuts::Vr;
}
float StPicoDstarMixedMaker::getTofBeta(StPicoTrack const* const trk) const
{
  int index2tof = trk->bTofPidTraitsIndex();
  float beta = std::numeric_limits<float>::quiet_NaN();
  if (index2tof >= 0)
  {
    StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
    if (tofPid)
    {
      beta = tofPid->btofBeta();
      if (beta < 1e-4)
      {
        // StThreeVectorF const vtx = mPicoDstMaker->picoDst()->event()->primaryVertex();
        TVector3 const vtx3 = mPicoDstMaker->picoDst()->event()->primaryVertex();
        StThreeVectorF vtx(vtx3.x(),vtx3.y(),vtx3.z());
        TVector3 const btofHitPos3 = tofPid->btofHitPos();
        // StThreeVectorF const btofHitPos = tofPid->btofHitPos();
        StThreeVectorF btofHitPos(btofHitPos3.x(),btofHitPos3.y(),btofHitPos3.z());
        StPicoPhysicalHelix helix = trk->helix(mPicoDstMaker->picoDst()->event()->bField());
        // StPhysicalHelixD helix = trk->helix(mPicoDstMaker->picoDst()->event()->bField());
        float L = tofPathLength(&vtx, &btofHitPos, helix.curvature());
        float tof = tofPid->btof();
        if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
        else beta = std::numeric_limits<float>::quiet_NaN();
      }
    }
  } 
  return beta;
}

