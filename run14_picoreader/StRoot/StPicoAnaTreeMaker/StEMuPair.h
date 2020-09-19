#ifndef StEMuPair_hh
#define StEMuPair_hh

#include <cmath>

class StPicoTrack;
class StPicoDst;
class StDcaGeometry;
class StElectronTrack;

#include "TObject.h"
#include "StThreeVectorF.hh"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <stdio.h>
#include <math.h>
#include "StEvent/StDcaGeometry.h"

class StEMuPair : public TObject {
 public:
  StEMuPair();
  ~StEMuPair();
	StEMuPair(Char_t type, UShort_t dauIndex1, UShort_t dauIndex2, Float_t dauDcaDist,
		Float_t  pairPt, Float_t  pairEta, Float_t  pairPhi, Float_t  pairMass,
		Float_t  pairPPt, Float_t  pairPEta, Float_t  pairPPhi, 
      Float_t  pairPMass, Float_t pairOx, Float_t pairOy, Float_t pairOz
        );

  virtual void Print(const Char_t *option = "") const;  ///< Print track info
            
  Int_t   type() const             { return (Int_t)mType; }
  Int_t   dauIndex1() const             { return (Int_t)mDauIndex1; }
  Int_t   dauIndex2() const             { return (Int_t)mDauIndex2; }
  Float_t  pairDca()      const          { return (Float_t)mDauDcaDist/1000.; }  /// dca between two daughters at pair decay point
  Float_t  pairPt() const			{return mPairPt;}
  Float_t  pairEta() const              { return mPairEta/10000.;}
  Float_t  pairPhi() const              { return mPairPhi/10000.;}
  Float_t  pairMass() const              { return (Float_t)mPairMass/1000.;}
  Float_t  pairY() const;
  
  Float_t  pairPPt() const			{return mPairPPt;}
  Float_t  pairPEta() const              { return mPairPEta/10000.;}
  Float_t  pairPPhi() const              { return mPairPPhi/10000.;}
  Float_t  pairPMass() const              { return (Float_t)mPairPMass/1000.;}
  Float_t  pairPY() const;
  StThreeVectorF pairMom() const;
  StThreeVectorF pairPMom() const;
  StThreeVectorF pairOrigin() const;

 protected:
  Char_t   mType;  // 1: e+ mu- 2: e- mu+ 3: e+ mu+ 4: e- mu-
  UShort_t mDauIndex1; // StElectronTrack
  UShort_t mDauIndex2; // StMuonTrack

  Short_t  mDauDcaDist;// *1000

  // global 
  Float_t  mPairPt;
  Short_t  mPairEta;// *10000
  Short_t  mPairPhi;// *10000
  UShort_t mPairMass;// *1000

  // primary
  Float_t  mPairPPt;
  Short_t  mPairPEta;// *10000
  Short_t  mPairPPhi;// *10000
  UShort_t mPairPMass; // *1000
  
  Short_t  mPairOx; //*100
  Short_t  mPairOy; //*100
  Short_t  mPairOz; //*100

  friend class StPicoDst;

  ClassDef(StEMuPair, 1)
};

inline Float_t StEMuPair::pairY() const
{
	TLorentzVector ee(0,0,0,0);
	ee.SetPtEtaPhiM(pairPt(),pairEta(),pairPhi(),pairMass());
	return ee.Rapidity();
}

inline Float_t StEMuPair::pairPY() const
{
	TLorentzVector ee(0,0,0,0);
	ee.SetPtEtaPhiM(pairPPt(),pairPEta(),pairPPhi(),pairPMass());
	return ee.Rapidity();
}

inline StThreeVectorF StEMuPair::pairMom() const
{
   TVector3 m(0,0,0);
   m.SetPtEtaPhi(mPairPt,pairEta(),pairPhi());
   return StThreeVectorF(m.X(),m.Y(),m.Z());
}

inline StThreeVectorF StEMuPair::pairPMom() const
{
    TVector3 m(0,0,0);
    m.SetPtEtaPhi(mPairPPt,pairPEta(),pairPPhi());
  return StThreeVectorF(m.X(),m.Y(),m.Z());
}


inline StThreeVectorF StEMuPair::pairOrigin() const
{
  return StThreeVectorF(mPairOx/100.,mPairOy/100.,mPairOz/100.);
}

#endif
