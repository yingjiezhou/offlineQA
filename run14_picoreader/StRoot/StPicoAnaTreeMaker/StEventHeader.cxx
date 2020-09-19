#include "StEventHeader.h"
#include "StEventTypes.h"
#include "StTree.h"
#include "StuRefMult.hh"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "TVector2.h"
#include "TRandom3.h"
#include "StMuDSTMaker/COMMON/StMuDst.h" 
#include "StMuDSTMaker/COMMON/StMuTrack.h" 
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuMtdHeader.h"
#include "StBTofHeader.h"
#include "StMessMgr.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StEventHeader)

StEventHeader::StEventHeader()
{}

StEventHeader::StEventHeader(const StPicoDst& picoDst, const Float_t *recenterCor, Bool_t doEvtPlane)
{
  prepareEventInfo(picoDst,recenterCor,doEvtPlane);
  //do the refmult correction
  StRefMultCorr* grefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr() ;
  grefmultCorrUtil->init(mRunId);
  grefmultCorrUtil->initEvent(mGRefMult, mVz, mZDCx) ;

  const Int_t cent16_grefmult = grefmultCorrUtil->getCentralityBin16() ;
  const Int_t cent9_grefmult  = grefmultCorrUtil->getCentralityBin9() ;
  const Double_t reweight = grefmultCorrUtil->getWeight() ;
  // NOTE: type should be double or float, not integer
  const Double_t grefmultCor = grefmultCorrUtil->getRefMultCorr() ;

  mGRefMultCorr = (UShort_t)grefmultCor;
  mReweight_GRefMult = (float)reweight;
  mCentrality9 = (UShort_t) cent9_grefmult ;
  mCentrality16  = (UShort_t) cent16_grefmult;

}

StEventHeader::StEventHeader(const StPicoDst& picoDst)
{
  prepareEventInfo(picoDst,0,false);
  StPicoEvent* ev = picoDst.event() ;

  //do the refmult correction
  StRefMultCorr* grefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr() ;
  grefmultCorrUtil->init(mRunId);
  grefmultCorrUtil->initEvent(mGRefMult, mVz, mZDCx) ;

  const Int_t cent16_grefmult = grefmultCorrUtil->getCentralityBin16() ;
  const Int_t cent9_grefmult  = grefmultCorrUtil->getCentralityBin9() ;
  const Double_t reweight = grefmultCorrUtil->getWeight() ;
  // NOTE: type should be double or float, not integer
  const Double_t grefmultCor = grefmultCorrUtil->getRefMultCorr() ;

  mGRefMultCorr = (UShort_t)grefmultCor;
  mReweight_GRefMult = (float)reweight;
  mCentrality9 = (UShort_t) cent9_grefmult ;
  mCentrality16  = (UShort_t) cent16_grefmult;

}

StEventHeader::~StEventHeader()
{ }

void StEventHeader::prepareEventInfo(const StPicoDst& picoDst, const Float_t *recenterCor=0, Bool_t doEvtPlane=false)
{
  TRandom3 *mRandom = new TRandom3();
  StPicoEvent* ev = picoDst.event() ;

  mRunId = ev->runId();
  mEventId = ev->eventId();
  mBField = ev->bField();

  StThreeVectorF mPrimaryVertex = ev->primaryVertex();
  if( mPrimaryVertex.x()==mPrimaryVertex.y()&&mPrimaryVertex.y()==mPrimaryVertex.z() ){
    mVx = -999;
    mVy = -999;
    mVz = -999;
  }else{
    mVx = mPrimaryVertex.x();
    mVy = mPrimaryVertex.y();
    mVz = mPrimaryVertex.z();
  }
  mVzVpd = ev->vzVpd();

  mTriggerIds = ev->triggerIds();
  mRefMultNeg = (UShort_t)(ev->refMultNeg());
  mRefMult = (UShort_t)(ev->refMultNeg()) + (UShort_t)(ev->refMultPos());
  mGRefMult = (UShort_t)(ev->grefMult());

  mRanking = ev->ranking() ;
  mNBEMCMatch = ev->nBEMCMatch() ;
  mNBTOFMatch = ev->nBTOFMatch() ;

  mZDCx = (UInt_t)(ev->ZDCx());
  mBBCx = (UInt_t)(ev->BBCx());
  mBackgroundRate = ev->backgroundRate();
  mbTofTrayMultiplicity = ev->btofTrayMultiplicity() ;
  mNumberOfGlobalTracks = ev->numberOfGlobalTracks() ;
  mNHitsHFT[0] = ev->numberOfPxlInnerHits();
  mNHitsHFT[1] = ev->numberOfPxlOuterHits();
  mNHitsHFT[2] = ev->numberOfIstHits();
  mNHitsHFT[3] = ev->numberOfSsdHits();

  mRandom->SetSeed(mEventId*100000+mZDCx);
  //--------------------------
  for(int i=0;i<4;i++) {
    setHT_Th(i,0);
  }

  StThreeVectorF vtxPos = picoDst.event()->primaryVertex();
  Int_t nTracks = picoDst.numberOfTracks();

  Double_t cossum = 0;
  Double_t sinsum = 0;
  Double_t cossum0 = 0;
  Double_t sinsum0 = 0;
  Double_t cossumQ1 = 0.;
  Double_t sinsumQ1 = 0.;
  Double_t sinsumQ2 = 0.;
  Double_t cossumQ2 = 0.;
  Double_t sinsumQplus = 0.;
  Double_t cossumQplus = 0.;
  Double_t sinsumQminus = 0.;
  Double_t cossumQminus = 0.;
  Double_t sinsumQetaplus = 0.;
  Double_t cossumQetaplus = 0.;
  Double_t sinsumQetaminus = 0.;
  Double_t cossumQetaminus = 0.;

  if(doEvtPlane){
    for(int n=0;n<nTracks;n++){
      StPicoTrack *t = (StPicoTrack*)picoDst.track(n);
      Int_t q = t->charge();
      Float_t dca = t->helix().distance(vtxPos,kFALSE);;
      double ratio = (double)t->nHitsFit()*1./t->nHitsMax();
      Double_t vertZ = vtxPos.z();

      StThreeVectorF gmom = t->gMom(vtxPos, picoDst.event()->bField());
      Float_t pt   = gmom.perp();
      Float_t eta  = gmom.pseudoRapidity();
      Float_t phi  = gmom.phi();

      Float_t nSigmaElectron= t->nSigmaElectron();
      if(phi<0.0) phi += (2.*TMath::Pi());
      Int_t nHitsFit = t->nHitsFit();

      if(fabs(q)!=1) continue;
      if(fabs(nHitsFit)<16) continue;
      if((pt>=2.0)||(pt<0.15)) continue;
      if(fabs(eta)>=1.0) continue;
      if(ratio < 0.52) continue;
      if(ratio >= 1.05) continue;
      if(dca >= 3.0) continue;

      Double_t cos_part =0.;
      Double_t sin_part =0.;
      Double_t coscorrect =0;
      Double_t sincorrect =0;

      if(eta > 0.0 && vertZ > 0.0){
        coscorrect = recenterCor[0];
        sincorrect = recenterCor[1];
      }
      else if(eta > 0.0 && vertZ < 0.0){
        coscorrect = recenterCor[2];
        sincorrect = recenterCor[3];
      }
      else if(eta < 0.0 && vertZ > 0.0){
        coscorrect = recenterCor[4];
        sincorrect = recenterCor[5];
      }
      else{
        coscorrect = recenterCor[6];
        sincorrect = recenterCor[7];
      }

      cos_part = pt*cos(2.*phi)-coscorrect;
      sin_part = pt*sin(2.*phi)-sincorrect;

      //Event Plane Vector
      cossum0 += pt*cos(2.*phi);
      sinsum0 += pt*sin(2.*phi);
      cossum += cos_part;
      sinsum += sin_part;
      if(q >0){
        cossumQplus += cos_part;
        sinsumQplus += sin_part;
      }
      else{
        cossumQminus += cos_part;
        sinsumQminus += sin_part;
      }
      if(eta>0){
        cossumQetaplus += cos_part;
        sinsumQetaplus += sin_part;
      }

      else{
        cossumQetaminus += cos_part;
        sinsumQetaminus += sin_part;
      }

      //randomly divide tracks in this event into two subevent1 and subevent2
      double ran = mRandom->Uniform(0,1);
      if(ran > 0.5){
        cossumQ1 += cos_part;
        sinsumQ1 += sin_part;
      } else {//subevent 2
        cossumQ2 += cos_part;
        sinsumQ2 += sin_part;
      }
    }
  }

  mRandom->Delete();

  mQ0x 			= cossum0*100;
  mQ0y 			= sinsum0*100;

  mQx 			= cossum*100;
  mQy 			= sinsum*100;

  mQ1x 			= cossumQ1*100;
  mQ1y 			= sinsumQ1*100;
  mQ2x 			= cossumQ2*100;
  mQ2y 			= sinsumQ2*100;

  mQplusx 		= cossumQplus*100;
  mQplusy 		= sinsumQplus*100;
  mQminusx 	= cossumQminus*100;
  mQminusy 	= sinsumQminus*100;

  mQetaplusx 	= cossumQetaplus*100;
  mQetaplusy 	= sinsumQetaplus*100;
  mQetaminusx = cossumQetaminus*100;
  mQetaminusy = sinsumQetaminus*100;


  //calculate the original eventplane without any correction
  TVector2 *Q0 = new TVector2(cossum0, sinsum0);
  if((Q0->Mod())!=0.0){
    double  eventPlane0 = 0.5*Q0->Phi();
    if(eventPlane0<0.0)eventPlane0+=TMath::Pi();
    mEventplane0 = eventPlane0*100;
  }


  //calculate the eventplant with the recenter factor applied
  TVector2 *Q = new TVector2(cossum, sinsum);
  if((Q->Mod())!=0.0){
    double  eventPlane = 0.5*Q->Phi();
    if(eventPlane<0.0)eventPlane+=TMath::Pi();
    mEventplane = eventPlane*100;
  }

  //subevent1:
  TVector2 *Q1 = new TVector2(cossumQ1, sinsumQ1);
  if((Q1->Mod())!=0.0){
    double  eventPlane1 = 0.5*Q1->Phi();
    if(eventPlane1<0.0) eventPlane1 += TMath::Pi();
    mEventplane1 = eventPlane1*100;
  }

  //subevent2:
  TVector2 *Q2 = new TVector2(cossumQ2, sinsumQ2);
  if((Q2->Mod())!=0.0){
    double  eventPlane2 = 0.5*Q2->Phi();
    if(eventPlane2<0.0) eventPlane2 += TMath::Pi();
    mEventplane2 = eventPlane2*100;
  }
}

bool StEventHeader::isTrigger(unsigned int id) const
{
    return std::find(mTriggerIds.begin(), mTriggerIds.end(), id) != mTriggerIds.end();
}
