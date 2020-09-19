#include "StMyAnalysisMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StThreeVectorF.hh"
#include "TFile.h"
#include "Stiostream.h"
#include <TMath.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "PhysicalConstants.h"
#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include <fstream>
#include <vector>
#include <algorithm>
#include "TProfile.h"

ClassImp(StMyAnalysisMaker)
StMyAnalysisMaker::StMyAnalysisMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName)
  : StMaker(name)
{
    mPicoDstMaker = picoMaker;
    mPicoDst = 0;
    mOutName = outName;
}

StMyAnalysisMaker::~StMyAnalysisMaker(){}

void StMyAnalysisMaker::setOutputName(Char_t* dir,Char_t* file)
{
    TString fileName = (TString) file;
    mOutputName = fileName+ ".root";
    return;
}

Float_t StMyAnalysisMaker::getEventInfo(StPicoDst* mPicoDst)
{
    if(!mPicoDst) return 0;
    mPicoEvent = (StPicoEvent*) mPicoDst->event();
    if(!mPicoEvent) return 0;

    //cout<<mPicoEvent->ZDCx()<<endl;
    hNumOfEvent->Fill(1);
    hNumOfEvent->Fill(2);
    TVector3 vertexPos = mPicoEvent->primaryVertex();
    hVzVsVpdVz_noCut->Fill(vertexPos.Z(), mPicoEvent->vzVpd());
    hVxVy_noCut->Fill(vertexPos.X(), vertexPos.Y());

    int tofMult  = mPicoEvent->btofTrayMultiplicity();
    int refMult  = mPicoEvent->refMult();
    int refMult3 = mPicoEvent->refMult3();


    if(TMath::Abs(vertexPos.Z()) > 30) return 0;
    hNumOfEvent->Fill(3);
    if(sqrt((vertexPos.X()-0)*(vertexPos.X()-0)+(vertexPos.Y()-0)*(vertexPos.Y()-0)) > 2) return 0;
    hNumOfEvent->Fill(4);
    if((abs(vertexPos.X())<1e-5) && (abs(vertexPos.Y())<1e-5) && (abs(vertexPos.Z())<1e-5) ) return 0;
    hNumOfEvent->Fill(5);

    if(TMath::Abs(vertexPos.Z())<1.0) hRef3_vz1->Fill(refMult3);
    hMeanRefMult-> Fill(1, refMult);
    hMeanRefMult3->Fill(1, refMult3);
    hMeanVz->      Fill(1, vertexPos.Z());
    hMeanVr->      Fill(1, sqrt(pow(vertexPos.X(),2)+pow(vertexPos.Y(),2)) );
    hMeanZDCx->    Fill(1, mPicoEvent->ZDCx());

    hVzVsVpdVz->Fill(vertexPos.Z(), mPicoEvent->vzVpd());
    hVxVy->Fill(vertexPos.X(), vertexPos.Y());
    hRef ->Fill(refMult);
    hRef3->Fill(refMult3);
    mZDCx->Fill(mPicoEvent->ZDCx(), refMult3);
    hRef3VsTofMult->Fill(refMult3, tofMult);
    hRefVsTofMult->Fill(refMult, tofMult);
    int Tracks = mPicoDst->numberOfTracks();                   
    float mField = mPicoEvent->bField();
    int nPrimaryTrack = 0;
    int nP=0;
    int nAP=0;
    int nTofMatch=0;
    int beta_eta1=0;
    for (int i=0;i<Tracks;i++){
      mPicoTrack = mPicoDst->track(i);
      if(!mPicoTrack) continue;
      if(!mPicoTrack->isPrimary()) continue;
      nPrimaryTrack++;
      TVector3 momentum = mPicoTrack->pMom();

      float pt = momentum.Perp();
      float phi = momentum.Phi();
      float eta= momentum.PseudoRapidity();
      int Charge = mPicoTrack->charge();
      float p = momentum.Mag();
      float pz= momentum.Z();
	    double EP = sqrt(p * p + 0.938 * 0.938);
	    float YP = 0.5 * log((EP + pz)/(EP - pz));
	    if(YP!= YP) YP = -999;
      StPicoPhysicalHelix helix = mPicoTrack->helix(mField);
      float dca = helix.geometricSignedDistance(vertexPos);
      if(dca<0) dca=fabs(dca);
      int tofIndex = mPicoTrack->bTofPidTraitsIndex();
      int   btofMatchFlag =  0;
      float btofYLocal    =  -999;
      float tof = 0;
      float beta= 0;
      float L=0;
      if(tofIndex>=0) {
        StPicoBTofPidTraits *tofPid = mPicoDst->btofPidTraits(tofIndex);
        btofMatchFlag = tofPid->btofMatchFlag();
        btofYLocal    = tofPid->btofYLocal();
        if(tofPid) {
          beta = tofPid->btofBeta();
          tof = tofPid->btof();
          if(beta<1e-4) {
            TVector3 btofHitPos_ = tofPid->btofHitPos();
            const StThreeVectorF *btofHitPos = new StThreeVectorF(btofHitPos_.X(),btofHitPos_.Y(),btofHitPos_.Z());
            const StThreeVectorF *vertexPos_ = new StThreeVectorF(vertexPos.X(), vertexPos.Y(), vertexPos.Z());
            L = tofPathLength(vertexPos_, btofHitPos, helix.curvature());
            if(tof>0) beta = L/(tof*(C_C_LIGHT/1.e9));
            else  beta = std::numeric_limits<float>::quiet_NaN();
          }
        }
      }
      float mass2=0;
      bool isGoodTof = btofMatchFlag >0 && beta > 0 && fabs(btofYLocal) < 1.8;
      if(isGoodTof) mass2 = momentum.Mag()*momentum.Mag()*(1./pow(beta,2)-1); else mass2 = -999;

      if( btofMatchFlag > 0 && fabs(eta) < 0.5 && dca < 3 && mPicoTrack->nHitsFit() > 10 ) nTofMatch++;
      if( beta > 0.1 && fabs(eta) < 1 && dca < 3 && mPicoTrack->nHitsFit() > 10 ) beta_eta1++;

      //track quality cuts
      if(momentum.Perp() < 0.2) continue;
      if(momentum.Mag() > 10.) continue;
      if((Float_t)mPicoTrack->nHitsFit()/(Float_t)mPicoTrack->nHitsMax() < 0.52 ) continue;
      if(mPicoTrack->nHitsFit()<=20) continue;
      if(mPicoTrack->nHitsDedx()<=5) continue;
      if(dca>1.0) continue;


      //--------------TOF efficiency-----------
	    noTofFlag_pt->Fill(pt);
	    if(btofMatchFlag > 0) TofFlag_pt->Fill(pt);
	    
	    noTofFlag_eta->Fill(eta);
	    if(btofMatchFlag > 0) TofFlag_eta->Fill(eta);

	    noTofFlag_phi->Fill(phi);
	    if(btofMatchFlag > 0) TofFlag_phi->Fill(phi);
	    
	    if(Charge > 0) noTofFlag_eta_pt_pos->Fill(pt,eta);
	    if(Charge < 0) noTofFlag_eta_pt_neg->Fill(pt,eta);
	    if(btofMatchFlag > 0 && Charge > 0 ) TofFlag_eta_pt_pos->Fill(pt,eta);
	    if(btofMatchFlag > 0 && Charge < 0 ) TofFlag_eta_pt_neg->Fill(pt,eta);
	    
	    if(Charge > 0) noTofFlag_eta_phi_pos->Fill(phi,eta);
	    if(Charge < 0) noTofFlag_eta_phi_neg->Fill(phi,eta);
	    if(btofMatchFlag > 0 && Charge > 0) TofFlag_eta_phi_pos->Fill(phi,eta);
	    if(btofMatchFlag > 0 && Charge < 0) TofFlag_eta_phi_neg->Fill(phi,eta);
	    
	    if(TMath::Abs(mPicoTrack->nSigmaProton()) < 2 && TMath::Abs(YP) < 0.5) {

	      if(Charge == 1)  noTofFlag_pro->Fill(pt);
	      if(Charge == -1) noTofFlag_apro->Fill(pt);
	      if(btofMatchFlag > 0){
	    	  if(Charge == 1)  TofFlag_pro->Fill(pt);
	    	  if(Charge == -1) TofFlag_apro->Fill(pt);
	      }
	      if(pt > 0.4 && pt < 0.8 && p<1){
          if(Charge ==1)  nP+=1;
          if(Charge ==-1) nAP+=1;

	    	  if(Charge == 1)  noTofFlag_pro_lowpt->Fill(pt);
	    	  if(Charge == -1) noTofFlag_apro_lowpt->Fill(pt);

	    	  if(btofMatchFlag > 0){
	    	    if(Charge == 1)  TofFlag_pro_lowpt->Fill(pt);
	    	    if(Charge == -1) TofFlag_apro_lowpt->Fill(pt);
	    	  }
	      }
	      if(pt > 0.8 && pt < 2 && p<3 && mass2>0.6 && mass2<1.2 ){
          if(Charge ==1)  nP+=1;
          if(Charge ==-1) nAP+=1;

	    	  if(Charge == 1)  noTofFlag_pro_highpt->Fill(pt);
	    	  if(Charge == -1) noTofFlag_apro_highpt->Fill(pt);
	    	  if(btofMatchFlag > 0){
	    	    if(Charge == 1)  TofFlag_pro_highpt->Fill(pt);
	    	    if(Charge == -1) TofFlag_apro_highpt->Fill(pt);
	    	  }
	      }

	    }
    
      
      hMeanProton->Fill(1,  nP);
      hMeanAProton->Fill(1, nAP);
      //=============
      hMeanPt->Fill(1, momentum.Perp());
      hMeanEta->Fill(1, momentum.PseudoRapidity());
      hMeanDca->Fill(1, dca);
      hMeanPhi->Fill(1, momentum.Phi());

      hRagiVsMass2->Fill(mPicoTrack->charge()*momentum.Mag(), mass2);

      hnHits->Fill( mPicoTrack->nHitsFit() * mPicoTrack->charge());
      hnHitsPt->Fill(momentum.Perp(), mPicoTrack->nHitsFit());
      hnHitsdEdx->Fill(mPicoTrack->nHitsDedx());
      hdEdxVsRagi->Fill(momentum.Mag()*mPicoTrack->charge(),mPicoTrack->dEdx());

      hPt->Fill(momentum.Perp());
      hEta->Fill(momentum.PseudoRapidity());
      hPhi->Fill(momentum.Phi());
      hDca->Fill(dca);
      hCharge->Fill(mPicoTrack->charge());

      hPVsBeta->Fill(mPicoTrack->charge()*momentum.Mag(), beta);
      hPVsTOF->Fill(momentum.Mag(),tof);
      if(TMath::Abs(beta)>1e-5) 
        hPVsInvBeta->Fill(momentum.Mag(),1.0/beta);
      else
        hPVsInvBeta->Fill(momentum.Mag(),0.0);
    }
    hnTrksVsnTofMatch->Fill(nPrimaryTrack, mPicoEvent->nBTOFMatch());
    hTofMatchHisto->Fill(refMult3, nTofMatch);
    hBeta_eta1Histo->Fill(refMult3, beta_eta1);
    
    return 0;
}

Int_t StMyAnalysisMaker::Init()
{
    outfile = new TFile(mOutputName,"RECREATE");

    hTofMatchHisto = new TH2D("hTofMatchHisto",";refMult3;nTofMatch",1201,-0.5,1200.5,1201,-0.5,1200.5);
    hBeta_eta1Histo = new TH2D("hBeta_etaHisto",";refMult3;beta_eta1",1201,-0.5,1200.5,1201,-0.5,1200.5);

    hVzVsVpdVz_noCut = new TH2D("hVzVsVpdVz_noCut","TPC Vz (cm);VPD Vz (cm)",1000,-200,200,1000,-200,200);
    hVxVy_noCut = new TH2D("hVxVy_noCut", ";X (cm);Y (cm)", 1000, -10, 10, 1000, -10, 10);
    hVzVsVpdVz = new TH2D("hVzVsVpdVz", ";TPC Vz (cm);VPD Vz (cm)", 1000, -200, 200, 1000, -200, 200);
    hVxVy = new TH2D("hVxVy", ";X (cm);Y (cm)", 1000, -2.5, 2.5, 1000, -2.5, 2.5);

    hRef = new TH1D("hRef", ";refMult;counts", 1201, -0.5, 1200.5);
    hRef3 = new TH1D("hRef3", ";refMult3;counts", 1201, -0.5, 1200.5);
    hRef3_vz1 = new TH1D("hRef3_vz1", ";refMult3;counts", 1201, -0.5, 1200.5);
    hRefVsTofMult = new TH2D("hRefVsTofMult",";refMult;tofMult",1201,-0.5,1200.5,2001,-0.5,2000.5);
    hRef3VsTofMult = new TH2D("hRef3VsTofMult",";refMult3;tofMult",1201,-0.5,1200.5,2001,-0.5,2000.5);

    hnTrksVsnTofMatch = new TH2D("hnTrksVsnTofMatch",";nPriTrks;nTofMatch",1201,-0.5,1200.5,1201,-0.5,1200.5);
    hPt = new TH1D("hPt",";p_{T} (GeV/c);counts",500,0.0,6.0);
    hDca = new TH1D("hDca", ";dca (cm);counts", 500, 0, 5);
    hEta = new TH1D("hEta",";#eta;counts",100,-5,5);
    hPhi = new TH1D("hPhi",";#phi;counts",80,-4,4);
    hCharge = new TH1D("hCharge",";Charge;counts",5,-2.5,2.5);

    hPVsTOF = new TH2D("hPVsTimeOfFlight", ";p (GeV/c);TOF", 1000, 0, 6, 1000, 0, 100);
    hPVsBeta = new TH2D("hPVsBeta", ";p (GeV/c);#beta;", 1000, 0, 6, 1000, 0, 5);
    hPVsInvBeta = new TH2D("hPVsInvBeta", ";p (GeV/c);1/#beta;", 1000, 0, 6, 1000, 0, 5);
    hRagiVsMass2 = new TH2D("hRagiVsMass2",";p #times q (GeV/c);m^{2} (GeV^{2}/c^{4})",1200, -6., 6, 1200, -0.2, 15.8);

    hnHits = new TH1D("hnHits",";nHits * charge;counts",201,-100.5,100.5);
    hnHitsPt = new TH2D("hnHitsPt",";p_{T};nHits",100,0,5,60,0,60);
    hnHitsdEdx = new TH1D("hnHitsdEdx",";hnHitsdEdx;counts",200,-100.5,100.5);
    hdEdxVsRagi = new TH2D("hdEdxVsRagi",";Ragidity (GeV/c);dE/dx (keV/cm)", 2000, -6., 6., 1000, 0., 50.);

    mZDCx = new TProfile("mZDCx",";ZDCx;<RefMult3>",100,0,10000);

    hMeanRefMult  = new TProfile("hMeanRefMult", ";;<RefMult>",1,0.5,1.5);
    hMeanRefMult3 = new TProfile("hMeanRefMult3", ";;<RefMult3>",1,0.5,1.5);
    hMeanVz       = new TProfile("hMeanVz", ";;<Vz>",1,0.5,1.5);
    hMeanVr       = new TProfile("hMeanVr", ";;<Vr>",1,0.5,1.5);
    hMeanZDCx     = new TProfile("hMeanZDCx", ";;<ZDCx>",1,0.5,1.5);
    hMeanEta      = new TProfile("hMeanEta", ";;<Eta>",1,0.5,1.5);
    hMeanPhi      = new TProfile("hMeanPhi", ";;<Phi>",1,0.5,1.5);
    hMeanDca      = new TProfile("hMeanDca", ";;<Dca>",1,0.5,1.5);
    hMeanPt       = new TProfile("hMeanPt", ";;<Pt>",1,0.5,1.5);
    hNumOfEvent   = new TH1D("hNumOfEvent",";;Event counts",15,0.5,15.5);
    hMeanProton   = new TProfile("hMeanProton",";run index;<#Proton>",1,0.5,1.5);
    hMeanAProton  = new TProfile("hMeanAProton",";run index;<#Anti-Proton>",1,0.5,1.5);
    //--Tof efficiency---
    noTofFlag_eta = new TH1F("noTofFlag_eta","",2000,-10,10);
    TofFlag_eta   = new TH1F("TofFlag_eta","",2000,-10,10);

    noTofFlag_phi = new TH1F("noTofFlag_phi","",628,-3.14,3.14);
    TofFlag_phi   = new TH1F("TofFlag_phi","",628,-3.14,3.14);
    
    noTofFlag_pt = new TH1F("noTofFlag_pt","",1500,0,15);
    TofFlag_pt   = new TH1F("TofFlag_pt","",1500,0,15);

    noTofFlag_eta_pt_pos = new TH2F("noTofFlag_eta_pt_pos","",1500,0,15,2000,-10,10);
    noTofFlag_eta_pt_neg = new TH2F("noTofFlag_eta_pt_neg","",1500,0,15,2000,-10,10);
    TofFlag_eta_pt_pos   = new TH2F("TofFlag_eta_pt_pos","",1500,0,15,2000,-10,10);
    TofFlag_eta_pt_neg   = new TH2F("TofFlag_eta_pt_neg","",1500,0,15,2000,-10,10);
    
    noTofFlag_eta_phi_pos = new TH2F("noTofFlag_eta_phi_pos","",628,-3.14,3.14,2000,-10,10);
    noTofFlag_eta_phi_neg = new TH2F("noTofFlag_eta_phi_neg","",628,-3.14,3.14,2000,-10,10);
    TofFlag_eta_phi_pos   = new TH2F("TofFlag_eta_phi_pos","",628,-3.14,3.14,2000,-10,10);
    TofFlag_eta_phi_neg   = new TH2F("TofFlag_eta_phi_neg","",628,-3.14,3.14,2000,-10,10);
    
    noTofFlag_pro  = new TH1F("noTofFlag_pro","",300,0,15);
    TofFlag_pro    = new TH1F("TofFlag_pro","",300,0,15);
    noTofFlag_apro = new TH1F("noTofFlag_apro","",300,0,15);
    TofFlag_apro   = new TH1F("TofFlag_apro","",300,0,15);

    noTofFlag_pro_lowpt   = new TH1F("noTofFlag_pro_lowpt","",300,0,15);
    TofFlag_pro_lowpt     = new TH1F("TofFlag_pro_lowpt","",300,0,15);
    noTofFlag_apro_lowpt  = new TH1F("noTofFlag_apro_lowpt","",300,0,15);
    TofFlag_apro_lowpt    = new TH1F("TofFlag_apro_lowpt","",300,0,15);

    noTofFlag_pro_highpt  = new TH1F("noTofFlag_pro_highpt","",300,0,15);
    TofFlag_pro_highpt    = new TH1F("TofFlag_pro_highpt","",300,0,15);
    noTofFlag_apro_highpt = new TH1F("noTofFlag_apro_highpt","",300,0,15);
    TofFlag_apro_highpt   = new TH1F("TofFlag_apro_highpt","",300,0,15);
  


    return kStOK;
}

Int_t StMyAnalysisMaker::Finish() {

    outfile->cd();
    outfile->Write();
    outfile->Close();
    return kStOK;
}

void StMyAnalysisMaker::Clear(Option_t *opt) {}

Int_t StMyAnalysisMaker::Make() {
    if(!mPicoDstMaker) {
        LOG_WARN << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }

    mPicoDst = mPicoDstMaker->picoDst();
    if(!mPicoDst) {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
    }
    getEventInfo(mPicoDst);

    return kStOK;
}
void StMyAnalysisMaker::savehistograms(){

        
}
