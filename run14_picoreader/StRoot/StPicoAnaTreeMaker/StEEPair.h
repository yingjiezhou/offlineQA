#ifndef StEEPair_hh
#define StEEPair_hh

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

class StEEPair : public TObject {
 public:
  StEEPair();
  ~StEEPair();
	StEEPair(Char_t mType, Short_t dauIndex1, Short_t dauIndex2, Float_t  dauDcaDist,
        Float_t  pairDcaToVtx, Float_t  pointingAngle, Float_t  pairPhiV, 
		Float_t  pairPt, Float_t  pairEta, Float_t  pairPhi, Float_t  pairMass,  
		Float_t  pairPMass,  
		Float_t pairCtau, Float_t pairOx, Float_t pairOy, Float_t pairOz
        );

  virtual void Print(const Char_t *option = "") const;  ///< Print track info
            
  Int_t   type() const             { return (Int_t)mType; }
  Int_t   dauIndex1() const             { return (Int_t)mDauIndex1; }
  Int_t   dauIndex2() const             { return (Int_t)mDauIndex2; }
  Float_t  pairDca()      const          { return (Float_t)mDauDcaDist/10000.; }  /// dca between two daughters at pair decay point

  Float_t  pairPt() const			{return mPairPt/1000.;}
  Float_t  pairEta() const              { return mPairEta/10000.;}
  Float_t  pairPhi() const              { return mPairPhi/10000.;}
  Float_t  pairMass() const              { return (Float_t)mPairMass/1000.;}
  Float_t  pairY() const;

  //Float_t  pairPPt() const			{return mPairPPt/1000;}
  //Float_t  pairPEta() const              { return mPairPEta/10000.;}
  //Float_t  pairPPhi() const              { return mPairPPhi/10000.;}
  Float_t  pairPMass() const              { return (Float_t)mPairPMass/1000.;}
  //Float_t  pairPY() const;

  StThreeVectorF pairMom() const;
  StThreeVectorF pairOrigin() const;

  //StThreeVectorF pairPMom() const;
  //Float_t  cosThetaStar() const              { return mCosThetaStar/10000.;}
  Float_t  pointingAngle() const              { return mPointingAngle/10000.;}
  Float_t  pairPhiV() const              { return mPairPhiV/10000.;}
  Float_t  pairDcaToVtx()      const       { return (Float_t)mPairDcaToVtx; }
  Float_t  pairCtau()      const       { return (Float_t)mPairCtau; }

 protected:
  Char_t   mType; // 1 = +-, 2 = ++, 3 = -- 
  Short_t  mDauIndex1;
  Short_t  mDauIndex2;

  //Short_t  mDauDcaToVtx1; // *1000
  //Short_t  mDauDcaToVtx2;// *1000
  Short_t  mDauDcaDist;// *10000

  Float_t  mPairDcaToVtx;
  //Short_t  mCosThetaStar; // *10000
  Short_t  mPointingAngle; // *10000
  Short_t  mPairPhiV; // *10000

  // global 
  UShort_t mPairPt; // *1000
  Short_t  mPairEta;// *10000
  Short_t  mPairPhi;// *10000
  UShort_t mPairMass;// *1000

  // primary
  //UShort_t  mPairPPt;  //*1000
  //Short_t  mPairPEta; //*10000
  //Short_t  mPairPPhi; //*10000
  UShort_t mPairPMass; // *1000
  Float_t  mPairCtau; 

  Short_t  mPairOx; //*100
  Short_t  mPairOy; //*100
  Short_t  mPairOz; //*100

  friend class StPicoDst;

  ClassDef(StEEPair, 1)
};

inline StThreeVectorF StEEPair::pairMom() const
{
    TVector3 m(0,0,0);
    m.SetPtEtaPhi(pairPt(),pairEta(),pairPhi());
  return StThreeVectorF(m.X(),m.Y(),m.Z());
}
//inline StThreeVectorF StEEPair::pairPMom() const
//{
//    TVector3 m(0,0,0);
//    m.SetPtEtaPhi(pairPPt(),pairPEta(),pairPPhi());
//  return StThreeVectorF(m.X(),m.Y(),m.Z());
//}

inline StThreeVectorF StEEPair::pairOrigin() const
{
  return StThreeVectorF(mPairOx/100.,mPairOy/100.,mPairOz/100.);
}

inline Float_t StEEPair::pairY() const
{
	TLorentzVector ee(0,0,0,0);
	ee.SetPtEtaPhiM(pairPt(),pairEta(),pairPhi(),pairMass());
	return ee.Rapidity();
}

//inline Float_t StEEPair::pairPY() const
//{
//	TLorentzVector ee(0,0,0,0);
//	ee.SetPtEtaPhiM(pairPPt(),pairPEta(),pairPPhi(),pairPMass());
//	return ee.Rapidity();
//}
#endif
