#include "StEMuPair.h"
#include "StElectronTrack.h"
#include "TVector2.h"
#include "TRandom.h"
#include "StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#define ElectronMass 0.000510999

ClassImp(StEMuPair)

    //----------------------------------------------------------------------------------

StEMuPair::StEMuPair()
    : mType(0), mDauIndex1(-1), mDauIndex2(-1), mDauDcaDist(0),
   mPairPt(0),mPairEta(0),mPairPhi(0),mPairMass(0),
   mPairPPt(0),mPairPEta(0),mPairPPhi(0),
   mPairPMass(0),mPairOx(0),mPairOy(0),mPairOz(0)
{

}

//----------------------------------------------------------------------------------
StEMuPair::StEMuPair(Char_t type, UShort_t dauIndex1, UShort_t dauIndex2, Float_t dauDcaDist,
		Float_t  pairPt, Float_t  pairEta, Float_t  pairPhi, Float_t  pairMass,
		Float_t  pairPPt, Float_t  pairPEta, Float_t  pairPPhi, 
		Float_t  pairPMass, Float_t pairOx, Float_t pairOy, Float_t pairOz)
    : mType(type), mDauIndex1(dauIndex1), mDauIndex2(dauIndex2), mDauDcaDist(dauDcaDist*1000),
   mPairPt(pairPt), mPairEta(pairEta*10000), mPairPhi(pairPhi*10000), mPairMass(pairMass*1000),
   mPairPPt(pairPt), mPairPEta(pairPEta*10000), mPairPPhi(pairPPhi*10000),
   mPairPMass(pairPMass*1000), mPairOx(pairOx*100), mPairOy(pairOy*100), mPairOz(pairOz*100)
{
}

//----------------------------------------------------------------------------------
StEMuPair::~StEMuPair()
{
    /* noop */
}
//----------------------------------------------------------------------------------
void StEMuPair::Print(const Char_t *option) const
{
   LOG_INFO << "mom=" << pairMom() << " mass = "<<pairMass()<<endm;
}
