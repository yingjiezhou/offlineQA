#ifndef StMuonTrack_hh
#define StMuonTrack_hh

#include <cmath>

class StPicoTrack;
class StPicoDst;
class StPicoMtdPidTraits;
class StDcaGeometry;

#include "TObject.h"
#include "StThreeVectorF.hh"
#include "TVector2.h"
#include <stdio.h>
#include <math.h>
#include "StEvent/StDcaGeometry.h"
#include "PhysicalConstants.h"

// Macro to control EMC variables
#define EMCON 1

class StPicoEvent;

class StMuonTrack : public TObject {
 public:
  StMuonTrack();
  ~StMuonTrack();
  StMuonTrack(StPicoDst *picoDst, StPicoTrack *t, Int_t idx);
  virtual void Print(const Char_t *option = "") const;  ///< Print track info
            
  Int_t   id() const             { return (Int_t)mId; }
  //Float_t chi2() const           { return (Float_t)mChi2/1000.; }
  Float_t gPt() const;
  Float_t gEta() const;
  Float_t gPhi() const;
  StThreeVectorF gMom() const;
  StThreeVectorF const & pMom() const    { return mPMom; }
  Short_t charge() const         { return (mNHitsFit>0) ? +1 : -1; }
  Int_t   nHitsFit() const       { return (mNHitsFit>0) ? (Int_t)mNHitsFit : (Int_t)(-1*mNHitsFit); }
  Int_t   nHitsDedx() const      { return (Int_t)mNHitsDedx; }
  //Int_t   nHitsMax() const       { return (Int_t)mNHitsMax; }
  //Int_t   nHitsMapHFT() const    { return (Int_t)(mMap0 >> 1 & 0x7F); }
  //UInt_t  map0() const { return (UInt_t)mMap0; }
  //UInt_t  map1() const { return (UInt_t)mMap1; }
  //Float_t dEdx() const           { return (Float_t)mDedx/1000.; }
  //Float_t tofMatchFlag() const   { return mTofMatchFlag; }
  Float_t beta() const           { return (Float_t)mBeta/20000.; }
  Float_t localY() const          {return (Float_t)mLocalY/1000.;}
  //Float_t localZ() const          {return (Float_t)mLocalZ/1000.;}
  Float_t nSigmaPion() const     { return (Float_t)mNSigmaPion/1000.; }
  Float_t dca() const           { return (Float_t)abs(mDca)/10000.; }
  Float_t dcaXY() const           { return (Float_t)mDcaXY/10000.; }
  Float_t dcaZ() const           { return (Float_t)mDcaZ/10000.; }
  Float_t dcaZLine() const           { return (Float_t)mDcaZLine/10000.; }
  Int_t matchFlag() const      { return (Int_t)mMatchFlag;}
  Float_t deltaTimeOfFlight() const { return (Float_t)mdT/1000.;}
  Float_t deltaY()  const       { return (Float_t)mdY/100.;}
  Float_t deltaZ()  const       { return (Float_t)mdZ/100.;}
  Int_t  triggerFlag()  const       { return (Int_t)mTriggerFlag;}
  Int_t  channel()  const       { return mChannel;} //
  Int_t  backleg()  const       { return mChannel/60+1;} //1-30
  Int_t  module()  const       { return (mChannel%60)/12+1;} //1-5
  Int_t  cell()  const       { return (mChannel%60)%12;} //0-11


//  const Float_t* params() const     { return mPar; }
//  const Float_t* errMatrix() const  { return mErrMatrix; }

  StPhysicalHelixD helix(float bField) const;
  Bool_t isHFTTrack() const { return mDca<0?true:false; }
          
 protected:
  Short_t mId;               // track Id
  StThreeVectorF mPMom;  // primary momentum, (0.,0.,0.) if none
  StThreeVectorF mGMom;
  //UShort_t mDedx;             // dEdx*1000
  Short_t  mDca;              // dca * 10000 * (isHFT?-1:1)
  Short_t  mDcaXY;            // dcaXY * 10000
  Short_t  mDcaZ;             // dcaZ * 10000
  Short_t  mDcaZLine;             // dcaZLine * 10000
  Char_t   mNHitsFit;         // q*nHitsFit
  //Char_t   mNHitsMax;         // nHitsMax
  UChar_t  mNHitsDedx;        // nHitsDedx
  Short_t  mNSigmaPion;       // nsigmaPi * 1000
  //UInt_t   mMap0;             // TopologyMap data0 HFT + TPC
  //UInt_t   mMap1;             // TopologyMap data1 TPC + Others
  StThreeVectorF mOrigin;
 
  // pidTraits
  //Char_t   mTofMatchFlag;
  UShort_t mBeta;  // *20000
  Short_t  mLocalY; // *1000
 
  // pidTraits
  Char_t   mMatchFlag;
  Short_t  mChannel;
  Short_t  mdT; // *1000
  Short_t  mdY; // *100
  Short_t  mdZ; // *100
  Char_t   mTriggerFlag;

  friend class StAnaTree;

  ClassDef(StMuonTrack, 1)
};
inline Float_t StMuonTrack::gPt() const
{
  //return 1./fabs(mPar[3]);
  return mGMom.perp();
}
inline Float_t StMuonTrack::gEta() const
{

  //float ptt = gPt();
  //StThreeVectorF gmom(ptt*std::cos(mPar[2]),ptt*std::sin(mPar[2]),ptt*mPar[4]);
  //return gmom.pseudoRapidity();
  return mGMom.pseudoRapidity();
}
inline Float_t StMuonTrack::gPhi() const
{
  //float ptt = gPt();
  //StThreeVectorF gmom(ptt*std::cos(mPar[2]),ptt*std::sin(mPar[2]),ptt*mPar[4]);
  //return gmom.phi();
  return mGMom.phi();;
}


inline StThreeVectorF StMuonTrack::gMom() const
{
  //float ptt = gPt();
  //return StThreeVectorF(ptt*std::cos(mPar[2]),ptt*std::sin(mPar[2]),ptt*mPar[4]);
  return mGMom;
}

//inline StDcaGeometry StMuonTrack::dcaGeometry() const
//{
//  StDcaGeometry a;
//  a.set(mPar, mErrMatrix);
//  return a;
//}
      
//inline StPhysicalHelixD StMuonTrack::helix() const
//{
//  return dcaGeometry().helix();
//}        

inline StPhysicalHelixD StMuonTrack::helix(float bField) const
{
  	//return StPhysicalHelixD(mCurv,mDip,mPhase,mOrigin,mH);
  	return StPhysicalHelixD(mGMom,mOrigin,bField*kilogauss,charge());
}        


#endif
