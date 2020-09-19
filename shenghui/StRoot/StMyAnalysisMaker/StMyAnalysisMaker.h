#ifndef STAR_picoMaker
#define STAR_picoMaker
#include "StMaker.h"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StThreeVectorD.hh"
#include "TString.h"
#include "TVector3.h"

class StPicoBTofPidTraits;
class StPicoDst;
class StPicoTrack;
class StPicoEvent;
class StPicoDstMaker;
class TH1D;
class TH1F;
class TH2D;
class TH2F;
class TProfile;

class StMyAnalysisMaker : public StMaker {
  public:
    StMyAnalysisMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName);
    virtual ~StMyAnalysisMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();

    void savehistograms();
    Float_t getEventInfo(StPicoDst* mPicoDst);
    Char_t *    outDir;
    void   setOutputName( Char_t* dir, Char_t* name="eFinder");
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    StPicoEvent    *mPicoEvent;
    StPicoTrack    *mPicoTrack;

    TString mOutName;
    TString mOutputName;
    int     mNEvent;
    char    mFileName[50];
    TFile *outfile;

    TH2D* hVzVsVpdVz_noCut;
    TH2D* hVxVy_noCut;
    TH2D* hVzVsVpdVz;
    TH2D* hVxVy;

    TH1D* hRef;
    TH1D* hRef3;
    TH1D* hRef3_vz1;
    TH2D* hRefVsTofMult;
    TH2D* hRef3VsTofMult;

    TH1D* hPt;
    TH1D* hEta;
    TH1D* hPhi;
    TH1D* hDca;
    TH1D* hCharge;

    TH1D* hnHits;
    TH1D* hnHitsdEdx;
    TH2D* hRagiVsMass2;

    TH2D* hPVsTOF;
    TH2D* hPVsBeta;
    TH2D* hPVsInvBeta;
    TH2D* hdEdxVsRagi;
    TProfile* mZDCx;
    TH2D* hnHitsPt;

    TProfile* hMeanRefMult;
    TProfile* hMeanRefMult3;
    TProfile* hMeanVz;
    TProfile* hMeanVr;
    TProfile* hMeanZDCx;
    TProfile* hMeanEta;
    TProfile* hMeanPhi;
    TProfile* hMeanDca;
    TProfile* hMeanPt;
    TProfile* hMeanProton;
    TProfile* hMeanAProton;


    TH2D *hTofMatchHisto;
    TH2D* hBeta_eta1Histo;

    TH2D*  hnTrksVsnTofMatch;
    TH1D *hNumOfEvent;
    //track level 
	  TH1F* noTofFlag_pt;
	  TH1F* TofFlag_pt;
	  
	  TH1F* noTofFlag_eta;
	  TH1F* TofFlag_eta;
	  
	  TH1F* noTofFlag_phi;
	  TH1F* TofFlag_phi;
	  
	  TH2F* noTofFlag_eta_pt_pos;
	  TH2F* TofFlag_eta_pt_pos;
	  TH2F* noTofFlag_eta_pt_neg;
	  TH2F* TofFlag_eta_pt_neg;
	  
	  TH2F* noTofFlag_eta_phi_pos;
	  TH2F* TofFlag_eta_phi_pos;
	  TH2F* noTofFlag_eta_phi_neg;
	  TH2F* TofFlag_eta_phi_neg;
	  
	  TH1F* noTofFlag_pro;
	  TH1F* TofFlag_pro;
	  TH1F* noTofFlag_apro;
	  TH1F* TofFlag_apro;

	  TH1F* noTofFlag_pro_lowpt;
	  TH1F* TofFlag_pro_lowpt;
	  TH1F* noTofFlag_apro_lowpt;
	  TH1F* TofFlag_apro_lowpt;
	  
	  TH1F* noTofFlag_pro_highpt;
	  TH1F* TofFlag_pro_highpt;
	  TH1F* noTofFlag_apro_highpt;
	  TH1F* TofFlag_apro_highpt;

    ClassDef(StMyAnalysisMaker,1)
};

#endif
