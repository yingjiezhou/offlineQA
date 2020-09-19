#ifndef StPicoDstarMixedMaker_h
#define StPicoDstarMixedMaker_h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoDstarEvent
 *  simultaneously and do analysis.
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include "TChain.h"
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
class TString;
class TFile;
class TNtuple;
class StPicoTrack;
class StPicoDstMaker;
class StPicoEvent;

struct ParticleInfo_Electron
{
    Int_t charge;
    Float_t pt;
    Float_t eta;
    Float_t phi;
    Float_t p;
    Float_t nSigmaPi;
    Float_t beta;
   // Float_t betaElectron;
   // Float_t dBetaElectron;
    Float_t energy;
    Float_t p1;
    Float_t p2;
    Float_t p3;
};

//vector<ParticleInfo> ElectronInfo;

class StPicoDstarMixedMaker : public StMaker
{
  public:
    StPicoDstarMixedMaker(char const * name, TString const inputFilesList,
        TString const outBaseName, StPicoDstMaker* picoDstMaker);
    virtual ~StPicoDstarMixedMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();
  private:

    StPicoDstarMixedMaker() {}
    void initHists();
    bool isGoodTrigger(StPicoEvent const*) const;
    bool isGoodQaEvent(StPicoEvent const* const picoEvent) const;
    bool isGoodEvent(StPicoEvent const* const picoEvent) const;
    bool isGoodQaTrack(StPicoTrack const* const trk) const;
    bool isGoodTrack(StPicoTrack const* trk,float dca) const;
    float getTofBeta(StPicoTrack const* const trk) const;
    StPicoDstMaker* mPicoDstMaker;
    TString mInputFilesList;
    TString mOutFileBaseName;
    bool isBadrun(Int_t runId);
    
    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate 

  public:
    void setRunNumList(string list){
      mRunNumList = list;
    }
    void setRunbyRunQA(bool b){
      mRunbyRunQA = b;
    }

    void getBadruns(string inputFileName);
  private: 
    TFile* mFile;

    TFile* mFile_RunID;//check ToF problems
    TH2F* hinvBetavsP_RunID[78];//check ToF problems

    std::map<int,int> runnum;
    string mRunNumList;
    vector<int> mBadRun;

    //Event level
    int  mRunId;
    float mVpdVz;
    float mRefmult;
    float mVpdHitEast;
    float mVpdHitWest;
    //primaryVertex vertex
    float mVx;
    float mVy;
    float mVz;
    float mVr;
    //event QA
    TH3F* hVxVyVz;
    TH1F* hVz;
    TH1F* hVpdVz;
    TH1F* hVr;
    TH1F* hVzVpdVz;
    TH2F* hnEvsPhivsVz;
    TH2F* hnEvsEtavsVz;
    TH2F* hnTofMulvsRef; 
    TH2F* hnTofMatvsRef;
    TH2F* hnTofHitvsRef;
    TH2F* h_Vx_Vy;
    TH1F* h_p_E0;
    TH1F* h_bemcdz;
    TH1F* h_bemcDphi;
    TH1F* h_nSMDphi;
    TH1F* h_nSMDeta;
    TH1F* h_btowPhiDz;
    TH1F* h_btowEtaDz;

    TH1F* h_mtdDeltaY;
    TH1F* h_mtdDeltaZ;
    TH1F* h_mtdDeltaTOF;
    TProfile* ph_p_E0;
    TProfile* ph_bemcdz;
    TProfile* ph_bemcDphi;
    TProfile* ph_nSMDphi;
    TProfile* ph_nSMDeta;
    TProfile* ph_btowPhiDz;
    TProfile* ph_btowEtaDz;

    TProfile* ph_mtdDeltaY;
    TProfile* ph_mtdDeltaZ;
    TProfile* ph_mtdDeltaTOF;

    TH1D* hevt;
    TH1D* hevtcut;
    TH1D* hevtbadcut;
    TH1D* hpassevtcut;
    TH1F* hrefmult;
     
   //tracl level QA
    TH3F* hNsigEvsinvBeta;
    TH1F* hnHitsFit;
    TH1F* hpDca;
    TH1F* hgDca;
    TH2F* hinvBetavsP;
    TH2F* hinvBetavsY;
    TH2F* hdEdx;
    TH2F* h_mTpc;
    TH1F* hpt;
    TH1F* hGpt;
    TH1F* hEta;
    TH1F* hPhi;
    TH3F* hBadTofId;
    TH2F* h_nSigmaElectron_P;//tof problem cell calibration
    TH2F* h_nSigmaElectron_P_tpc;//tof problem cell calibration

    //invariant mass electron
    //TH1F* hMeeCount;
    TH1F* hMeeCount;//unlike sign
    TH1F* hMeeCount_like1;//like sign electron
    TH1F* hMeeCount_like2;//like sign positron
    TH2F* hMeeCountPt;//unlike sign
    TH2F* hMeeCountPt_like1;//like sign electron
    TH2F* hMeeCountPt_like2;//like sign positron
    float M_electron=0.000511;//GeV

    //MRPC ToF 
    TH1F* ModuleId_1;//1/Beta 0.8-0.9, P 0.4-3 GeV
    TH1F* TofId_1;
    TH1F* TrayId_1;
    TH1F* ModuleId_2;//1/Beta 0.82-0.9, P 0.4-3 GeV
    TH1F* TofId_2;
    TH1F* ModuleId_3;//1/Beta 0.82-0.88, P 0.4-3 GeV
    TH1F* TofId_3;
    TH1F* ModuleId_4;//1/Beta 0.84-0.88, P 0.4-3 GeV
    TH1F* TofId_4;
    TH1F* ModuleId_5;//1/Beta 0.84-0.9, P 0.4-3 GeV
    TH1F* TofId_5;


    //Run by run QA
    bool mRunbyRunQA;
    TProfile* pVpdVz;
    TProfile* pVzVpdVz;
    TProfile* pRefmult;
    TProfile* pVpdHitEast;
    TProfile* pVpdHitWest;
    TProfile* pVx;
    TProfile* pVy;
    TProfile* pVz;
    TProfile* pVr;
    //Run by run track level
    TProfile* pTof;
    TProfile* pDedx;
    TProfile* pRMSDedx;
    TProfile* pgDCA;
    TProfile* pgPt;
    TProfile* pgPhi;
    TProfile* pgEta;
    TProfile* pNFits; 
    TProfile* ppPt;
    TProfile* ppEta;
    TProfile* ppPhi;
    ClassDef(StPicoDstarMixedMaker, 1)
};

inline void StPicoDstarMixedMaker::getBadruns(string inputFileName){
    ifstream fin(inputFileName.c_str());
    if(!fin){
      cout <<"no Bad runs list" << endl;
      return;
    }
    cout << "  " << inputFileName << flush;

    Int_t runId = 0 ;
    while( fin >> runId ) {
      mBadRun.push_back(runId);
    }
    cout << "get Bad runs list [OK]" << endl;
}
inline  bool StPicoDstarMixedMaker::isBadrun(Int_t runId){
    vector<Int_t>::iterator iter = std::find(mBadRun.begin(), mBadRun.end(), runId);
    return ( iter != mBadRun.end() ) ;
}

#endif
