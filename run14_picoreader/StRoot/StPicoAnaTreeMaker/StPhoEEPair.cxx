#include "StPhoEEPair.h"
#include "StElectronTrack.h"
#include "TVector2.h"
#include "TRandom.h"
#include "StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#define ElectronMass 0.000510999

ClassImp(StPhoEEPair)

    //----------------------------------------------------------------------------------
StPhoEEPair::StPhoEEPair()
    : mType(0), mPrimEIndex(-1), mPartEIndex(-1),
	mDauDcaDist(32768), mPairPhiV(32768),
mPairPt(0),mPairEta(0),mPairPhi(0),mPairMass(0), 
   mPairPMass(0),
	mPairOx(0),mPairOy(0),mPairOz(0)
{

}

/////////////////////////////////////////////////////////////////////////////////////////
// t - the global track.  p - the associated primary track from the first primary vertex
/////////////////////////////////////////////////////////////////////////////////////////
//----------------------------------------------------------------------------------
StPhoEEPair::StPhoEEPair(
      Char_t   type, Short_t primEIndex, Short_t partEIndex,
		Float_t  dauDcaDist, Float_t  pairPhiV, 
		Float_t  pairPt, Float_t  pairEta, Float_t  pairPhi, Float_t  pairMass, 
		Float_t  pairPMass, Float_t  pairOx, Float_t pairOy, Float_t pairOz
        )
   : mType(type), mPrimEIndex(primEIndex), mPartEIndex(partEIndex),
	mDauDcaDist(dauDcaDist*10000), mPairPhiV(pairPhiV*10000),
   mPairPt(pairPt*1000),mPairEta(pairEta*10000),mPairPhi(pairPhi*10000),mPairMass(pairMass*10000), 
   mPairPMass(pairPMass*10000), mPairOx(pairOx*100), mPairOy(pairOy*100), mPairOz(pairOz*100)
{
}

//----------------------------------------------------------------------------------
StPhoEEPair::~StPhoEEPair()
{
    /* noop */
}
//----------------------------------------------------------------------------------
void StPhoEEPair::Print(const Char_t *option) const
{
   LOG_INFO << "mom=" << pairMom() << " mass = "<<pairMass()<<endm;
}
