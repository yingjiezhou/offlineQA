#ifndef STMTDPRODQAMAKER_HH
#define STMTDPRODQAMAKER_HH

/***************************************************************************
 *
 * StMtdProdQAMaker - class to make MTD related QA plots during production
 * Author: Rongrong Ma
 *--------------------------------------------------------------------------
 *
 ***************************************************************************/

#include "StMaker.h"
#include <vector>
#ifndef ST_NO_NAMESPACES
using std::vector;
#endif

class TH1F;
class TH2F;
class TProfile;
class TTree;
class TFile;

class StMuDst;
class StMuTrack;
class StMuMtdHit;
class StTriggerData;
class StEmcPosition;
class StEmcGeom;

class StPicoDst;
class StPicoTrack;
class StPicoMtdHit;

#include "StPhysicalHelixD.hh"
#include "StThreeVector.hh"
#include "StThreeVectorF.hh"
#include "StThreeVectorD.hh"
#include "StMtdUtil/StMtdConstants.h"

#if !defined(ST_NO_TEMPLATE_DEF_ARGS) || defined(__CINT__)
typedef vector<Int_t> IntVec;
#else
typedef vector<Int_t, allocator<Int_t>> IntVec;
#endif

const Double_t muMass = 0.10566;

class StMtdProdQAMaker : public StMaker {
 public:
  StMtdProdQAMaker(const Char_t *name = "MtdQAMaker");
  ~StMtdProdQAMaker();

  Int_t    Init();
  Int_t    InitRun(const Int_t runNumber);
  Int_t    Make();
  Int_t    Finish();
  void     setYear(const int year)             { mYear = year;              }
  void     setDataType(const int type)         { mDataType = type;          }
  void     setOutputFileName(const char *name) { mOutputFileName = name;    }
  void     setPrintConfig(const Bool_t print)  { mPrintConfig = print;      }
  void     setSeparateTrig(const bool sep)     { mSeparateTrig = sep;       }

  void     setMaxTpcVz(const Double_t max)     { mMaxTpcVz = max;           }
  void     setMaxDiffVz(const Double_t max)    { mMaxDiffVz = max;          }
  void     setMinNHitsFit(const Int_t min)     { mMinNHitsFit = min;        }
  void     setMinNHitsDedx(const Int_t min)    { mMinNHitsDedx = min;       } 
  void     setMinNHitsFrac(const Double_t min) { mMinNHitsFrac = min;       }
  void     setMaxDca(const Double_t max)       { mMaxDca = max;             }

  void     setTrackPtLimits(const Double_t min, const Double_t max)
  { mMinTrkPt = min; mMaxTrkPt = max; }
  void     setTrackPhiLimits(const Double_t min, const Double_t max)
  { mMinTrkPhi = min; mMaxTrkPhi = max; }
  void     setTrackEtaLimits(const Double_t min, const Double_t max)
  { mMinTrkEta = min; mMaxTrkPhi = max;}
  void     setNsigmaPiCut(const Double_t min, const Double_t max)
  { mMinNsigmaPi = min; mMaxNsigmaPi = max; }
  void     setMuonDeltaZ(const Double_t min, const Double_t max)
  { mMinMuonDeltaZ = min; mMaxMuonDeltaZ = max; }
  void setMuonDeltaY(const Double_t min, const Double_t max)
  { mMinMuonDeltaY = min; mMaxMuonDeltaY = max; }
  void setMuonDeltaTof(const Double_t min, const Double_t max)
  { mMinMuonDeltaTof = min; mMaxMuonDeltaTof = max; }
  void setMuonDca(const double max)
  { mMaxMuonDca = max; }


 protected:
  void     initTrigIDs();
  void     printConfig();
  void     bookHistos();
  void     getTrigIndex();
  void     triggerData();
  Int_t    processMuDst();
  Int_t    processPicoDst();

  void     fillHisto(TH1F **h, float x);
  void     fillTProfile(TProfile **h, float y);
  void     fillHisto(TH2F **h, float x, float y);

  bool     getBemcInfo(StMuTrack *track, int &maxadc, float &energy, float &zdist, float &phidist, int &neta, int &nphi);

  Bool_t   isValidTrack(StMuTrack *t);
  Bool_t   isValidTrack(StPicoTrack *t, const StThreeVectorF mom, StThreeVectorD vtxPos);
  Int_t    getMtdHitTHUB(const Int_t backleg);

  Bool_t   isMtdHitFiredTrigger(const StMuMtdHit *hit);
  Bool_t   isMtdHitFiredTrigger(const StPicoMtdHit *hit);
  Bool_t   isQTFiredTrigger(const Int_t qt, const Int_t pos);

  Int_t    getMtdHitIndex(const StPicoTrack *track);
  Int_t    getMtdPidTraitsIndex(const StPicoMtdHit *hit);
  Bool_t   isMuonCandidate(const StMuTrack *track);
  Bool_t   isMuonCandidate(const StPicoTrack *track, const StThreeVectorD vtxPos);
  Bool_t   isMuonCandidate(const Double_t nSigmaPi, const Int_t tofindex, const Double_t dz, const Double_t dy, const Double_t dtof, const Double_t dca, const Bool_t isTrig);


  float    rotatePhi(float phi) const;
  void     addCutToHisto(TH1 *h, const Int_t bin, const char *label, const Float_t value = -9999999);

  enum {kNQTboard = 8};
  Int_t    mModuleToQT[gMtdNBacklegs][gMtdNModules];
  Int_t    mModuleToQTPos[gMtdNBacklegs][gMtdNModules]; 
  Int_t    mQTtoModule[kNQTboard][8];
  Int_t    mQTSlewBinEdge[kNQTboard][16][8]; 
  Int_t    mQTSlewCorr[kNQTboard][16][8];
  Int_t    mTrigQTpos[kNQTboard][2]; 
 

  enum     { kNtrig = 4 };

 private:
  StMuDst          *mMuDst;                                    // Pointer to MuDst event
  StPicoDst        *mPicoDst;                                  // 
  Int_t            mDataType;                                  // 0 - MuDst; 1 - PicoDst
  Int_t            mYear;                                      // Year
  Int_t            mRunId;                                     // Run number
  Int_t            mFirstRun;
  Int_t            mLastRun;
  IntVec           mTriggerIDs[kNtrig];   
  Bool_t           mSeparateTrig;                              // Separate different triggers
  Bool_t           mTrigIndex[kNtrig];                         // Valid trigger indices
  Bool_t           mPrintConfig;                               // Print configuration
  TString          mOutputFileName;                            //
  TFile            *mOutputFile;                               //
  Double_t         mTrigTime[2];
  StEmcPosition    *mEmcPosition;
  StEmcGeom        *mEmcGeom[4];
  
  // track quality cuts
  Double_t         mMaxTpcVz;       
  Double_t         mMaxDiffVz;
  Double_t         mMinTrkPt;                                  
  Double_t         mMaxTrkPt;                                  
  Double_t         mMinTrkPhi;                                 
  Double_t         mMaxTrkPhi;                                 
  Double_t         mMinTrkEta;                                 
  Double_t         mMaxTrkEta;                                 
  Int_t            mMinNHitsFit;                               
  Int_t            mMinNHitsDedx;                             
  Double_t         mMinNHitsFrac;                      
  Double_t         mMaxDca;                                 

  // PID cuts
  Double_t         mMinNsigmaPi;                              
  Double_t         mMaxNsigmaPi;                              
  Double_t         mMinMuonDeltaZ;                             
  Double_t         mMaxMuonDeltaZ;
  Double_t         mMinMuonDeltaY;
  Double_t         mMaxMuonDeltaY;
  Double_t         mMinMuonDeltaTof;
  Double_t         mMaxMuonDeltaTof;
  Double_t         mMaxMuonDca;
  Bool_t           mBTofMatch;
  Bool_t           mMtdHitTrigger;

  // Histograms
  TH1F             *mhEventStat[kNtrig];                              
  TH1F             *mhEventCuts;                               
  
  /// vertex
  TH2F             *mhVertexXYRanking0[kNtrig];                       
  TH1F             *mhTpcVzRanking0[kNtrig];                          
  TH1F             *mhDiffVzRanking0[kNtrig];
  TH1F             *mhTpcVtxIndex[kNtrig];
  TH2F             *mhVtxIndexVsClosest[kNtrig];
  TH2F             *mhVertexXY[kNtrig];                       
  TH1F             *mhTpcVz[kNtrig];
  TH1F             *mhVpdVz[kNtrig];
  TH1F             *mhDiffVz[kNtrig];
  TH2F             *mhDiffVzVsTpcVz[kNtrig];
  TH2F             *mhDiffVzVsVpdVz[kNtrig];


  /// reference multiplicity
  TH1F             *mhRefMult[kNtrig];
  TH1F             *mhgRefMult[kNtrig];
  TH2F             *mhgRefMultVsRefMult[kNtrig];
  TH2F             *mhTpcVzVsgRef[kNtrig];
  TH2F             *mhDiffVzVsgRef[kNtrig];
  TH2F             *mhZdcRateVsgRef[kNtrig];
  TH2F             *mhBbcRateVsgRef[kNtrig];

  /// primary tracks
  TH2F             *mhpTrkNHitsVsPt[kNtrig];
  TH1F             *mhpTrkNHits[kNtrig];
  TH2F             *mhpTrkNDedxVsPt[kNtrig];
  TH1F             *mhpTrkNDedx[kNtrig];
  TH2F             *mhpTrkNHitsFracVsPt[kNtrig];
  TH1F             *mhpTrkNHitsFrac[kNtrig];
  TH2F             *mhpTrkDcaVsPt[kNtrig];
  TH1F             *mhpTrkN[kNtrig];
  TH1F             *mhpTrkNMthMtd[kNtrig];
  TH1F             *mhpTrkPt[kNtrig];
  TH2F             *mhpTrkEtaVsPt[kNtrig];
  TH2F             *mhpTrkPhiVsPt[kNtrig];
  TH2F             *mhpTrkPhiVsEta[kNtrig];
  TH2F             *mhpTrkDedxVsMom[kNtrig];
  TH2F             *mhNsigmaEVsMom[kNtrig];
  TH2F             *mhNsigmaPiVsMom[kNtrig];
  TH2F             *mhNsigmaKVsMom[kNtrig];
  TH2F             *mhNsigmaPVsMom[kNtrig];
  TH2F             *mhBetaVsMom[kNtrig];
  TH2F             *mhM2VsMom[kNtrig];

  /// electron PID
  TH1F             *mhAdc0[kNtrig];
  TH2F             *mhEoverPVsPt[kNtrig];
  TH2F             *mhNetaVsPt[kNtrig];
  TH2F             *mhNphiVsPt[kNtrig];
  TH2F             *mhDistZVsPt[kNtrig];
  TH2F             *mhDistPhiVspt[kNtrig];
  TH1F             *mhNElecPos[kNtrig];
  TH1F             *mhNElecNeg[kNtrig];
  TH1F             *mhElecPt[kNtrig];
  TH2F             *mhElecPhiVsEta[kNtrig];
  
  /// trigger performance
  TH1F             *mhNQtSignal[kNtrig];
  TH1F             *mhNMT101Signal[kNtrig];
  TH1F             *mhNTF201Signal[kNtrig];
  TH2F             *mhMtdTacSumMixvsMxq[kNQTboard][2][kNtrig];
  TH2F             *mhMtdVpdTacDiffMT001[kNtrig];
  TH2F             *mhMtdVpdMthTacDiffMT001[kNtrig];
  TH2F             *mhMtdVpdTacDiffMT101[kNtrig];
  TH2F             *mhMtdVpdMthTacDiffMT101[kNtrig];

  /// MTD hits
  TH1F             *mhMtdNRawHits[kNtrig];
  TH2F             *mhMtdRawHitMap[kNtrig];
  TH1F             *mhMtdNHits[kNtrig];
  TH2F             *mhMtdHitMap[kNtrig];
  TH2F             *mhMtdHitTrigTime[kNtrig];
  TH2F             *mhMtdHitLeTimeDiff[kNtrig];
  TH1F             *mhMtdNTrigHits[kNtrig];
  TH2F             *mhMtdTrigHitMap[kNtrig];
  TH1F             *mhMtdNMthHits[kNtrig];
  TH2F             *mhMtdMthHitMap[kNtrig];
  TH1F             *mhMtdNMthTrigHits[kNtrig];
  TH2F             *mhMtdMthTrigHitMap[kNtrig];

  // Matching information
  TH2F             *mhLocalYVsgChan[kNtrig];
  TH2F             *mhLocalZVsgChan[kNtrig];
  TH2F             *mhMtdTofVsgChan[kNtrig];
  TH2F             *mhExpTofVsgChan[kNtrig];
  TH1F             *mhDeltaZ[kNtrig];
  TH2F             *mhDzVsPt[kNtrig];
  TH2F             *mhDzVsgChan[kNtrig];
  TH1F             *mhDeltaY[kNtrig];
  TH2F             *mhDyVsPt[kNtrig];
  TH2F             *mhDyVsgChan[kNtrig];
  TH1F             *mhDeltaTof[kNtrig];
  TH2F             *mhDTofVsPt[kNtrig];
  TH2F             *mhDTofVsgChan[kNtrig];
  
  // Muon analysis
  TH1F             *mhNMuonPos[kNtrig];
  TH1F             *mhNMuonNeg[kNtrig];
  TH1F             *mhMuonPt[kNtrig];
  TH2F             *mhMuonPhiVsEta[kNtrig];
  TH2F             *mhMuonMap[kNtrig];
  TH1F             *mhNULpair[kNtrig];
  TH1F             *mhNLSpairPos[kNtrig];
  TH1F             *mhNLSpairNeg[kNtrig];
  TH2F             *mhInvMvsPtUL[kNtrig];
  TH2F             *mhInvMvsPtLSpos[kNtrig];
  TH2F             *mhInvMvsPtLSneg[kNtrig];
  TH1F             *mhInvMUL[kNtrig];
  TH1F             *mhInvMLSpos[kNtrig];
  TH1F             *mhInvMLSneg[kNtrig];

  
  // Run Dependence
  TH1F             *mhRunStat[kNtrig];
  TProfile         *mhBBCrateVsRun[kNtrig];
  TProfile         *mhZDCrateVsRun[kNtrig];
  TProfile         *mhRefMultVsRun[kNtrig];
  TProfile         *mhgRefMultVsRun[kNtrig];
  TProfile         *mhTpcVxVsRun[kNtrig];
  TProfile         *mhTpcVyVsRun[kNtrig];
  TProfile         *mhTpcVzVsRun[kNtrig];
  TProfile         *mhVpdVzVsRun[kNtrig];
  TProfile         *mhDiffVzVsRun[kNtrig];
  TProfile         *mhpTrkPtVsRun[kNtrig];
  TProfile         *mhpTrkEtaVsRun[kNtrig];
  TProfile         *mhpTrkPhiVsRun[kNtrig];
  TProfile         *mhpTrkDcaVsRun[kNtrig];
  TProfile         *mhNHitsFitVsRun[kNtrig];
  TProfile         *mhNHitsPossVsRun[kNtrig];
  TProfile         *mhNHitsDedxVsRun[kNtrig];
  TProfile         *mhDedxVsRun[kNtrig];
  TProfile         *mhNsigmaPiVsRun[kNtrig];
  TProfile         *mhNsigmaEVsRun[kNtrig];
  TProfile         *mhNsigmaKVsRun[kNtrig];
  TProfile         *mhNsigmaPVsRun[kNtrig];
  TProfile         *mhBetaVsRun[kNtrig];
  TProfile         *mhNElectron[kNtrig];
  TProfile         *mhNPositron[kNtrig];
  TProfile         *mhNMtdRawHitsVsRun[kNtrig];
  TProfile         *mhNMtdHitsVsRun[kNtrig];
  TProfile         *mhNMtdTrigHitsVsRun[kNtrig];
  TProfile         *mhNMtdMthHitsVsRun[kNtrig];
  TProfile         *mhNMuonPosVsRun[kNtrig];
  TProfile         *mhNMuonNegVsRun[kNtrig];
  TProfile         *mhNMuonPairULVsRun[kNtrig];
  TProfile         *mhNMuonPairLSPosVsRun[kNtrig];
  TProfile         *mhNMuonPairLSNegVsRun[kNtrig];
  TProfile         *mhNJpsiVsRun[kNtrig];


  virtual const char *GetCVS() const {
    static const char cvs[]="Tag $Name:  $Id: built " __DATE__ " " __TIME__ ; 
    return cvs;
  }
  
  ClassDef(StMtdProdQAMaker, 0)
};

#endif
