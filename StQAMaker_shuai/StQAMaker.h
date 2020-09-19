#ifndef STQAMAKER_HH
#define STQAMAKER_HH

/***************************************************************************
 *
 * $Id: StQAMaker.h 2018/03/26  Exp $
 * StQAMaker - class to produce run by run QA histograms for ht2 related analysis
 * Author: Shuai Yang
 *--------------------------------------------------------------------------
 *
 ***************************************************************************/

#include "StMaker.h"

#include "StThreeVectorF.hh"
#include "TLorentzVector.h"

#include <vector>
#include <map>
#ifndef ST_NO_NAMESPACES
using std::vector;
#endif

class TProfile;
class TH1D;
class TH2D;
class TString;
class TTree;
class TFile;

class StMuDstMaker;
class StMuDst;
class StEmcCollection;
class StEmcPosition;
class StEmcGeom;
class StMuTrack;

#if !defined(ST_NO_TEMPLATE_DEF_ARGS) || defined(__CINT__)
typedef vector<Int_t> IntVec;
typedef vector<Double_t> DoubleVec;
typedef vector<TLorentzVector> LorentzVec;
#else
typedef vector<Int_t, allocator<Int_t>> IntVec;
typedef vector<Double_t, allocator<Double_t>> DoubleVec;
typedef vector<TLorentzVector, allocator<TLorentzVector>> LorentzVec;
#endif

class StQAMaker : public StMaker {
public:
  StQAMaker(const Char_t *name = "StQAMaker");
  ~StQAMaker();
  
  Int_t    Init();
  Int_t    Make();
  Int_t    Finish();
  
  void	 setTriggerIDs(const IntVec triggerids);
  void     setUseDefaultVtx(const Bool_t flag);
  void     setSelectVtxRanking(const Bool_t flag);
  void     setMaxVtxR(const Double_t max);
  void     setMaxVtxZ(const Double_t max);
  void     setMaxVzDiff(const Double_t max);
  void     setMinTrackPt(const Double_t min);
  void     setMaxTrackEta(const Double_t max);
  void     setMinNHitsFit(const Int_t min);
  void     setMinNHitsFitRatio(const Double_t min);
  void     setMinNHitsDedx(const Int_t min);
  void     setMaxDca(const Double_t max);
  void     setMaxnSigmaE(const Double_t max);
  void     setMaxBeta2TOF(const Double_t max);
  void     setMinBemcPt(const Double_t min);
  void     setMinAdc0(const Int_t min);
  void     setMinMaxPoverE(const Double_t min, const Double_t max);
  void     setMaxZDist(const Double_t max);
  void     setMaxPhiDist(const Double_t max);
  void     setMinNEta(const Int_t min);
  void     setMinNPhi(const Int_t min);
  void     setOutFileName(const TString name);
  void     setStreamName(const TString name);
  void     setPrintMemory(const Bool_t pMem);
  void     setPrintCpu(const Bool_t pCpu);
  void     setPrintConfig(const Bool_t print);
  
protected:
  void     printConfig();
  void     bookHistos();
  Bool_t   processMuDstEvent();
  Bool_t   isValidTrack(StMuTrack *pMuTrack) const;
  Bool_t   getBemcInfo(StMuTrack *pMuTrack, const Int_t runIdx, Int_t &nMthTrks, Int_t &nTrigTrks, Int_t &nBEMCeCans);
    
private:
  static const Int_t    mTotalRuns =  1503;
  
  StMuDstMaker    *mMuDstMaker;        // Pointer to StMuDstMaker
  StMuDst         *mMuDst;             // Pointer to MuDst event
  StEmcCollection *mEmcCollection;
  StEmcPosition   *mEmcPosition;
  StEmcGeom       *mEmcGeom[4];
  
  map<Int_t, Int_t> mTotalRunId;
  
  Bool_t         mPrintMemory;         // Flag to print out memory usage
  Bool_t         mPrintCpu;            // Flag to print out CPU usage
  Bool_t         mPrintConfig;         // Flag to print out task configuration
  TString        mStreamName;          // Data stream name
  Bool_t         mDefaultVtx;          // Use Default Vertex
  Bool_t         mSelectVtxRank;       // Vertex ranking > 0
  Double_t       mMaxVtxR;             // Maximum vertex r
  Double_t       mMaxVtxZ;             // Maximum vertex z
  Double_t       mMaxVzDiff;           // Maximum VpdVz-TpcVz
  Double_t       mMinTrkPt;            // Minimum track pt
  Double_t       mMaxTrkEta;           // Maximum track eta
  Int_t          mMinNHitsFit;         // Minimum number of hits used for track fit
  Double_t       mMinNHitsFitRatio;    // Minimum ratio of hits used for track fit
  Int_t          mMinNHitsDedx;        // Minimum number of hits used for de/dx
  Double_t       mMaxDca;              // Maximum track dca
  Double_t       mMaxnSigmaE;          // Maximum nSigmaE cut
  Double_t       mMaxBeta2TOF;         // Maximum |1-1./beta| for TpcE
  Double_t       mMinBemcPt;           // Minimum BEMC triggered track pt
  Int_t          mMinAdc0;             // Minimum BECM triggered track Adc0
  Double_t       mMinPoverE;           // Minimum BMEC electron p/E
  Double_t       mMaxPoverE;           // Maximum BECM electron p/E
  Double_t       mMaxZDist;            // Maximum BEMC electron #DeltaZ
  Double_t       mMaxPhiDist;          // Maximum BEMC electron #Delta#phi
  Int_t          mMinNEta;             // Minimum BEMC electron n_{#eta}
  Int_t          mMinNPhi;             // Minimum BEMC electron n_{#phi}
  
  TFile          *fOutFile;            // Output file
  TString        mOutFileName;         // Name of the output file
  TTree          *mEvtTree;            // Pointer to the event tree
  
  IntVec         mTriggerIDs;
  
  //define histograms
  
  // default vertex
  TH1D *hEvent_DefVtx;
  TH2D *hVyvsVx_DefVtx;
  TH2D *hVpdVzvsTpcVz_DefVtx;
  TH2D *hVzDiffvsTpcVz_DefVtx;
  TH2D *hVzDiffvsRefMult_DefVtx;
  TH1D *hRefMult_VzVrCut_DefVtx;
  TH1D *hRefMult_EvtCut_DefVtx;
  
  // selected vertex
  TH1D *hEvent;
  TH2D *hVyvsVx;
  TH2D *hVpdVzvsTpcVz;
  TH2D *hTpcVzvsRefMult;
  TH2D *hRawVpdVzvsRefMult;
  TH2D *hVpdVzvsRefMult;
  TH2D *hVzDiffvsTpcVz;
  TH2D *hVzDiffvsRefMult;
  TH1D *hRefMult_VzVrCut;
  TH1D *hRefMult_EvtCut;
  TH2D *hVtxIdxvsRefMult;
  TH2D *hVtxIdxvsRefMult_VzDiffCut;
  TH2D *hRefMultvsZdcX;
  TH2D *hRefMultvsBbcX;
  
  // inclusive QA
  TH1D *hNHitsFit;
  TH1D *hNHitsPoss;
  TH1D *hNHitsDedx;
  TH2D *hEtavsPt;
  TH2D *hFFPhivsPt;
  TH2D *hRFFPhivsPt;
  TH2D *hPosTrkEtavsPhi;
  TH2D *hNegTrkEtavsPhi;
  TH2D *hDcavsPt;
  TH2D *hdEdxvsP;
  TH2D *hdEdxvsPhi;
  TH2D *hdEdxvsEta;
  TH2D *hBetavsP;
  TH2D *hBEMCeEtavsPt;
  TH2D *hBEMCePhivsPt;
  TH2D *hBEMCeEtavsPhi;
  
  // run by run QA
  // event level QA
  TH1D     *hnHT2EvtsvsRun;
  TProfile *hBFieldvsRun;
  TProfile *hZdcXvsRun;
  TProfile *hBbcXvsRun;
  TProfile *hZdcXoverBbcXvsRun;
  TProfile *hTpcVxvsRun;
  TProfile *hTpcVyvsRun;
  TProfile *hTpcVzvsRun;
  TProfile *hRawVpdVzvsRun;
  TProfile *hVpdVzvsRun;
  TProfile *hVzDiffvsRun;
  TProfile *hRefMultvsRun;
  
  // primary tracks
  TProfile *hNTrksvsRun;
  TProfile *hPtvsRun;
  TProfile *hEtavsRun;
  TProfile *hPhivsRun;
  TProfile *hDcavsRun;
  TProfile *hNHitsFitvsRun;
  TProfile *hNHitsPossvsRun;
  TProfile *hNHitsDedxvsRun;
  TProfile *hDedxvsRun;
  TProfile *hNSigmaEvsRun;
  TProfile *hNSigmaPivsRun;
  TProfile *hNSigmaKvsRun;
  TProfile *hNSigmaPvsRun;
  TProfile *hBetavsRun;
  
  // BEMC match tracks
  TProfile *hNMthTrksvsRun;
  TProfile *hMthTrkPtvsRun;
  TProfile *hMthTrkEtavsRun;
  TProfile *hMthTrkPhivsRun;
  TProfile *hMthTrkNSigmaEvsRun;
  TProfile *hMthTrkBetavsRun;
  TProfile *hMthTrkAdc0vsRun;
  TProfile *hMthTrkE0vsRun;
  TProfile *hMthTrkEvsRun;
  TProfile *hMthTrkZDistvsRun;
  TProfile *hMthTrkPhiDistvsRun;
  TProfile *hMthTrkNEtavsRun;
  TProfile *hMthTrkNPhivsRun;
  
  // BEMC trigger tracks
  TProfile *hNTrigTrksvsRun;
  TProfile *hTrigTrkPtvsRun;
  TProfile *hTrigTrkEtavsRun;
  TProfile *hTrigTrkPhivsRun;
  TProfile *hTrigTrkNSigmaEvsRun;
  TProfile *hTrigTrkAdc0vsRun;
  TProfile *hTrigTrkE0vsRun;
  TProfile *hTrigTrkEvsRun;
  TProfile *hTrigTrkZDistvsRun;
  TProfile *hTrigTrkPhiDistvsRun;
  TProfile *hTrigTrkNEtavsRun;
  TProfile *hTrigTrkNPhivsRun;
  
  // BEMC electron candidates
  TProfile *hNBemcEsvsRun;
  
  ClassDef(StQAMaker, 1)
};

inline void	StQAMaker::setTriggerIDs(const IntVec triggerids) { mTriggerIDs.clear(); mTriggerIDs = triggerids; }
inline void StQAMaker::setUseDefaultVtx(const Bool_t flag) { mDefaultVtx = flag; }
inline void StQAMaker::setSelectVtxRanking(const Bool_t flag) { mSelectVtxRank = flag; }
inline void StQAMaker::setMaxVtxR(const Double_t max) { mMaxVtxR = max; }
inline void StQAMaker::setMaxVtxZ(const Double_t max) { mMaxVtxZ = max; }
inline void StQAMaker::setMaxVzDiff(const Double_t max) { mMaxVzDiff = max; }
inline void StQAMaker::setMinTrackPt(const Double_t min){ mMinTrkPt = min;}
inline void StQAMaker::setMaxTrackEta(const Double_t max){ mMaxTrkEta = max; }
inline void StQAMaker::setMinNHitsFit(const Int_t min) { mMinNHitsFit = min; }
inline void StQAMaker::setMinNHitsFitRatio(const Double_t min) { mMinNHitsFitRatio = min; }
inline void StQAMaker::setMinNHitsDedx(const Int_t min) { mMinNHitsDedx = min; }
inline void StQAMaker::setMaxDca(const Double_t max) { mMaxDca = max; }
inline void StQAMaker::setMaxnSigmaE(const Double_t max) { mMaxnSigmaE = max; }
inline void StQAMaker::setMaxBeta2TOF(const Double_t max) { mMaxBeta2TOF = max; }
inline void StQAMaker::setMinBemcPt(const Double_t min) { mMinBemcPt = min; }
inline void StQAMaker::setMinAdc0(const Int_t min) { mMinAdc0 = min; }
inline void StQAMaker::setMinMaxPoverE(const Double_t min, const Double_t max) { mMinPoverE = min; mMaxPoverE = max; }
inline void StQAMaker::setMaxZDist(const Double_t max) { mMaxZDist = max; }
inline void StQAMaker::setMaxPhiDist(const Double_t max) { mMaxPhiDist = max; }
inline void StQAMaker::setMinNEta(const Int_t min) { mMinNEta = min; }
inline void StQAMaker::setMinNPhi(const Int_t min) { mMinNPhi = min; }
inline void StQAMaker::setOutFileName(const TString name) { mOutFileName = name; }
inline void StQAMaker::setStreamName(const TString name) { mStreamName = name; }
inline void StQAMaker::setPrintMemory(const Bool_t pMem) { mPrintMemory = pMem; }
inline void StQAMaker::setPrintCpu(const Bool_t pCpu) { mPrintCpu = pCpu; }
inline void StQAMaker::setPrintConfig(const Bool_t print) { mPrintConfig = print; }
#endif
