#ifndef StPhoEEPair_hh
#define StPhoEEPair_hh

#include <cmath>

class StPicoTrack;
class StPicoDst;
class StDcaGeometry;
class StElectronTrack;

#include "TObject.h"
#include "TLorentzVector.h"
#include "StThreeVectorF.hh"
#include "TVector3.h"
#include <stdio.h>
#include <math.h>
#include "StEvent/StDcaGeometry.h"

class StPhoEEPair : public TObject {
 public:
  StPhoEEPair();
  ~StPhoEEPair();
	StPhoEEPair(
      Char_t   type, Short_t primEIndex, Short_t partEIndex,
		Float_t  dauDcaDist, Float_t  pairPhiV, 
		Float_t  pairPt, Float_t  pairEta, Float_t  pairPhi, Float_t  pairMass,  
		Float_t  pairPMass, Float_t pairOx, Float_t pairOy, Float_t pairOz
        );

  virtual void Print(const Char_t *option = "") const;  ///< Print track info
            
  Int_t    type() const             { return (Int_t)mType; }
  Int_t    primEIndex() const             { return (Int_t)mPrimEIndex; }
  Int_t    partEIndex() const             { return (Int_t)mPartEIndex; }
  Float_t  pairDca()      const          { return (Float_t)mDauDcaDist/10000.; }

  Float_t  pairPt() const			{return mPairPt/1000.;}
  Float_t  pairY() const;
  Float_t  pairEta() const              { return mPairEta/10000.;}
  Float_t  pairPhi() const              { return mPairPhi/10000.;}
  Float_t  pairMass() const              { return (Float_t)mPairMass/10000.;}
  Float_t  pairPMass() const              { return (Float_t)mPairPMass/10000.;}

  StThreeVectorF pairMom() const;
  StThreeVectorF pairPMom() const;
  StThreeVectorF pairOrigin() const;
  Float_t  pairPhiV() const              { return mPairPhiV/10000.;}

 protected:
  Char_t   mType; // 1 = +-, 2 = ++, 3 = -- 
  Short_t  mPrimEIndex; // primary electron index in StElectronTrack branch
  Short_t  mPartEIndex; // parterner electron index in StPartElectronTrack branch

  Short_t  mDauDcaDist;// *10000 dca between two daughters at pair decay point 
  Short_t  mPairPhiV; // *10000

  // global 
  UShort_t  mPairPt; // *1000
  Short_t  mPairEta;// *10000
  Short_t  mPairPhi;// *10000
  UShort_t mPairMass;// *10000

  // primary
  UShort_t mPairPMass; // *10000

  Short_t  mPairOx; //*100
  Short_t  mPairOy; //*100
  Short_t  mPairOz; //*100

  friend class StPicoDst;

  ClassDef(StPhoEEPair, 1)
};

inline StThreeVectorF StPhoEEPair::pairMom() const
{
    TVector3 m(0,0,0);
    m.SetPtEtaPhi(pairPt(),pairEta(),pairPhi());
  return StThreeVectorF(m.X(),m.Y(),m.Z());
}

inline StThreeVectorF StPhoEEPair::pairOrigin() const
{
  return StThreeVectorF(mPairOx/100.,mPairOy/100.,mPairOz/100.);
}

inline Float_t StPhoEEPair::pairY() const
{
	TLorentzVector ee(0,0,0,0);
	ee.SetPtEtaPhiM(pairPt(),pairEta(),pairPhi(),pairMass());
	return ee.Rapidity();
}
#endif

