#include "StMuonTrack.h"
#include "StMessMgr.h"
#include "TVector2.h"
#include "TMath.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoMtdPidTraits.h"
#include "StPicoDstMaker/StPicoMtdHit.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StBTofUtil/tofPathLength.hh"
#include "PhysicalConstants.h"
#include <climits>

ClassImp(StMuonTrack)

//----------------------------------------------------------------------------------
StMuonTrack::StMuonTrack() : mId(-1), mPMom(0., 0., 0.), mGMom(0.,0.,0.), 
   mDca(32768), mDcaXY(0), mDcaZ(0), mDcaZLine(32768),
   mNHitsFit(0), mNHitsDedx(0), 
   mNSigmaPion(32768),
   mOrigin(0,0,0),
   mBeta(0), mLocalY(32768), 
   mMatchFlag(0), mChannel(0), mdT(32768),mdY(32768), mdZ(32768),mTriggerFlag(-1)
{

}

//----------------------------------------------------------------------------------
StMuonTrack::StMuonTrack(StPicoDst *picoDst, StPicoTrack* t, Int_t idx)
{
   mId        = idx;
   //mChi2      = (t->chi2() * 1000. > 65536) ? 65536 : (UShort_t)(TMath::Nint(t->chi2() * 1000.));
   mPMom = t->pMom();
   int q      = t->charge();
   //mDedx      = (t->dEdx()*1000. > 65536) ? 65536 : (UShort_t)(TMath::Nint(t->dEdx()*1000.));
   mNHitsFit  = t->nHitsFit()*q;
   //mNHitsMax  = t->nHitsMax();
   mNHitsDedx = (UChar_t)(t->nHitsDedx());
   mNSigmaPion = (fabs(t->nSigmaPion() * 1000.) > 32768) ? 32768 : (Short_t)(TMath::Nint(t->nSigmaPion() * 1000.));

   //mMap0 = t->map0(); // see hitMap definition in StTrackTopologyMap
   //mMap1 = t->map1();
   //const float* params = t->params();
   //const float* errMatrix = t->errMatrix();
   //for (int i = 0; i < 6; i++) mPar[i] = params[i];
   //for (int i = 0; i < 15; i++) mErrMatrix[i] = errMatrix[i];

   StThreeVectorF vertexPos = picoDst->event()->primaryVertex();
   StPhysicalHelixD helix = t->helix();

   StThreeVectorF dcaPoint = helix.at(helix.pathLength(vertexPos.x(), vertexPos.y()));
   float dcaZ = (dcaPoint.z() - vertexPos.z())*10000.;
   float dcaXY = (helix.geometricSignedDistance(vertexPos.x(),vertexPos.y()))*10000.;
   mDcaZ = dcaZ>32768?32768:(Short_t)dcaZ;
   mDcaXY = dcaXY>32768?32768:(Short_t)dcaXY;

   double thePath = helix.pathLength(vertexPos);
   StThreeVectorF dcaPos = helix.at(thePath);
   mGMom = helix.momentumAt(thePath,picoDst->event()->bField()*kilogauss);
   mOrigin = dcaPos;

   double dca = (dcaPos-vertexPos).mag();
   bool isHft = t->isHFTTrack();
   if(isHft) dca *= -1;
   mDca = fabs(dca*10000.)>32768? 32768: (Short_t)(dca*10000.);
   
   StThreeVectorF posDiff = dcaPos - vertexPos;
   if(sin(mGMom.theta())==0) mDcaZLine = 32768;
   else mDcaZLine = (-posDiff.x()*cos(mGMom.phi())-posDiff.y()*sin(mGMom.phi()))*cos(mGMom.theta())/sin(mGMom.theta())+posDiff.z();

   int index2TofPid = t->bTofPidTraitsIndex();
   if (index2TofPid>=0){
      StPicoBTofPidTraits *tofPid = picoDst->btofPidTraits(index2TofPid);
      //mTofMatchFlag = tofPid->btofMatchFlag();
      mLocalY = tofPid->btofYLocal();
      //mLocalZ = tofPid->btofZLocal();
      Float_t beta = tofPid->btofBeta();
      if(beta<1e-4||beta>=(USHRT_MAX-1)/20000){
         Float_t tof = tofPid->btof();
         StThreeVectorF btofHitPos = tofPid->btofHitPos();
         float L = tofPathLength(&vertexPos, &btofHitPos, helix.curvature()); 
         beta = L/(tof*(c_light/1.0e9));
      }
      mBeta = (UShort_t)(beta*20000);
   }

   int index2MtdPid = t->mtdPidTraitsIndex();
   if (index2MtdPid>=0){
      StPicoMtdPidTraits *mtdPid = picoDst->mtdPidTraits(index2MtdPid);
      mMatchFlag = mtdPid->matchFlag();
      mdT = mtdPid->deltaTimeOfFlight()*1000;
      mdZ = mtdPid->deltaZ()*100;
      mdY = mtdPid->deltaY()*100;
      mChannel = (mtdPid->backleg()-1)*60+(mtdPid->module()-1)*12+mtdPid->cell();
      mTriggerFlag = 0;
      for(int i=0;i<picoDst->numberOfMtdHits();i++){
         StPicoMtdHit *hit = (StPicoMtdHit*)picoDst->mtdHit(i);
         if(mChannel == hit->gChannel()){
            mTriggerFlag = hit->triggerFlag();
         }
      }
   }
}
//----------------------------------------------------------------------------------
StMuonTrack::~StMuonTrack()
{
   /* noop */
}
//----------------------------------------------------------------------------------
void StMuonTrack::Print(const Char_t *option) const
{
   if (strcmp(option, "tpc") == 0 || strcmp(option, "") == 0)
   {
      LOG_INFO << "id=" << id() <<endm;
      LOG_INFO << "pMom=" << pMom() << endm;
      LOG_INFO << " nHitsFit = " << nHitsFit() << " nHitsdEdx = " << nHitsDedx() << endm;
   }
}
