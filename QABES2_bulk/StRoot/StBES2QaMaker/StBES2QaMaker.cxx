#include "StBES2QaMaker.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StThreeVectorF.hh"
//#include "StarClassLibrary/PhysicalConstants.h"  
#include "StPicoEvent/StPicoEpdHit.h"
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StEpdEpFinder.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"
//#include "StV0TofCorrection/StV0TofCorrection.h"

// tof
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "phys_constants.h"
// mtd
#include "StPicoEvent/StPicoMtdHit.h"
// btow
#include "StPicoEvent/StPicoBTowHit.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
// etof
#include "StPicoEvent/StPicoETofPidTraits.h"


const float Mass[3] = { //Gev
	0.13957, // pion 
	0.49368, // kaon
	0.93827  // proton
};
const Float_t massLamMin = 1.1105;
const Float_t massLamMax = 1.1205;

Double_t StBES2QaMaker::mZDCSMDCenterex = 0;
Double_t StBES2QaMaker::mZDCSMDCenterey = 0;
Double_t StBES2QaMaker::mZDCSMDCenterwx = 0;
Double_t StBES2QaMaker::mZDCSMDCenterwy = 0;

ClassImp(StBES2QaMaker)

//-----------------------------------------------------------------------------
StBES2QaMaker::StBES2QaMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName)
  : StMaker(name)
{
	mPicoDstMaker = picoMaker;
	mPicoDst = 0;
	mOutName = outName;

	gRandom->SetSeed(0);
	mCalibMode        = 3;
	mTrgSelection     = 0;
	mDcaMax           = 3.0;
	mNHitsFitMax      = 15;
	mVzAbsMax         = 6.0;
	mTrgEffCorrection = 0;
	mTrkEffCorrection = 0;
	mTrkEP            = 0; // 0=primary, 1=global
	mEventCounter     = 0;
	mTopologicalCut   = 0; // 0=fixed cut
	//mV0TofCorrection  = 0; // 0=not apply
	mEvMixVtxShift    = 0; // 0=off
	mFixedMode        = 0; // 0=collider mode,  1=fixed target

	creatingWeights = false;


	mDcaP[0] = 0.5;  mDcaPi[0] = 1.6;  mDcaV0[0] = 0.9;  mDL[0] = 5.0;  mDcaDh[0] = 1.0; //  0-10%
	mDcaP[1] = 0.3;  mDcaPi[1] = 1.6;  mDcaV0[1] = 0.9;  mDL[1] = 5.0;  mDcaDh[1] = 1.0; // 10-20%
	mDcaP[2] = 0.3;  mDcaPi[2] = 1.3;  mDcaV0[2] = 1.0;  mDL[2] = 4.0;  mDcaDh[2] = 1.0; // 20-30%
	mDcaP[3] = 0.3;  mDcaPi[3] = 1.2;  mDcaV0[3] = 1.0;  mDL[3] = 3.0;  mDcaDh[3] = 1.0; // 30-40%
	mDcaP[4] = 0.2;  mDcaPi[4] = 1.0;  mDcaV0[4] = 1.0;  mDL[4] = 3.0;  mDcaDh[4] = 1.0; // 40-50%
	mDcaP[5] = 0.1;  mDcaPi[5] = 0.9;  mDcaV0[5] = 1.0;  mDL[5] = 3.0;  mDcaDh[5] = 1.0; // 50-60%
	mDcaP[6] = 0.1;  mDcaPi[6] = 0.8;  mDcaV0[6] = 1.0;  mDL[6] = 3.0;  mDcaDh[6] = 1.0; // 60-70%
	mDcaP[7] = 0.1;  mDcaPi[7] = 0.7;  mDcaV0[7] = 1.0;  mDL[7] = 3.0;  mDcaDh[7] = 1.0; // 70-80%

	refMult_pre    = -999.9;
	vz_pre         = -999.9;
	zdcAdcSumE_pre = -999.9;

	mData  = 5; // 1=run11, 2=run14, 3=run16, 4=run18-isobar, 5=run19-19GeV, 6=run19-14GeV, 7=run19-7GeV, 8=run19-9GeV, 9=run19-FXT
}

//----------------------------------------------------------------------------- 
StBES2QaMaker::~StBES2QaMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StBES2QaMaker::Init() {

	DeclareHistograms();

	// For QV Maker
	pQVMaker = new QVMaker(mCalibMode);

	//string infile3_ = "OutEpdFinder";
	//const Char_t *Input = "InputEpdFinder";

	// ---------------------------------------------------------------------
	//Isaac-Mike's block to define EPD, copied from Isaac's Analysis code for 27 GeV
	// ---------------------------------------------------------------------
	// set up EPD Event Plane finder - malisa
	//if(!creatingWeights) mEpFinder = new StEpdEpFinder(9, infile3_.c_str(), Input); // 9 centrality bins.
	//else                 mEpFinder = new StEpdEpFinder(9, infile3_.c_str()); // 9 centrality bins.
	//mEpFinder->SetnMipThreshold(0.3);
	//mEpFinder->SetMaxTileWeight(2.0);
	//mEpFinder->SetEpdHitFormat(2);     // 2=pico

	//// Here, I set eta-based weights
	//// --->   This is only to be used for Au+Au 27 GeV <---
	//// https://drupal.star.bnl.gov/STAR/blog/lisa/optimizing-ep1-resolution-au27au-ring-dependent-weights-eta-dependent-weights
	//double lin[9] = {-1.950, -1.900, -1.850, -1.706, -1.438, -1.340, -1.045, -0.717, -0.700};
	//double cub[9] = {0.1608, 0.1600, 0.1600, 0.1595, 0.1457, 0.1369, 0.1092, 0.0772, 0.0700};
	//TH2D wt("Order1etaWeight","Order1etaWeight",100,1.5,6.5,9,0,9);
	//for (int ix=1; ix<101; ix++){
	//	for (int iy=1; iy<10; iy++){
	//		double eta = wt.GetXaxis()->GetBinCenter(ix);
	//		wt.SetBinContent(ix,iy,lin[iy-1]*eta+cub[iy-1]*pow(eta,3));
	//	}
	//}
	//mEpFinder->SetEtaWeights(1,wt);
	// ---------------------------------------------------------------------
	// Isaac-Mike's block ends here
	// ---------------------------------------------------------------------
	

	//pRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
	// run14
	if( mData==2 ){ // run14
		pRefMultCorr = CentralityMaker::instance()->getgRefMultCorr_P16id();
		pRefMultCorr->setVzForWeight(6,-6.0,6.0);
		pRefMultCorr->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P16id.txt");
	}else if( mData==3 ){ // run16
		pRefMultCorr = CentralityMaker::instance()->getgRefMultCorr_P16id(); // <- somehow same as run14, but it's correct
		pRefMultCorr->setVzForWeight(6,-6.0,6.0);
		pRefMultCorr->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_VpdnoVtx_Vpd5_Run16.txt");
	}

	// Load caribration parameters
	// -------------------------------------------------------
	CalibParamInit();
	SetBBCTilesPhi();
	if( LoadBBCalibParam( GetRunNum() )==-1 ) return kStOK;
	//if( mCalibMode>1 && LoadEPCalibParam( GetRunNum() )==-1 ) return kStOK;

	// Cut parameters
	string sData[] = { "Run10AuAu200", "Run11AuAu200", "Run14AuAu200", "Run16AuAu200", "Run18Isobar200", "Run19AuAu19", "Run19AuAu14", "Run19AuAu7", "Run19AuAu9", "Run19FXT", "Run19FXT31", "Run19AuAu200",
   						"Run20AuAu11p5", "Run20AuAu9p2", "Run20AuAuXX" };
	cout << "=============================================================" << endl;
	cout << " Data setup = " << sData[mData] << endl; 
	cout << " Calibration mode = " << mCalibMode << endl;
	cout << " Selected trigger ID = " << mTrgSelection << " ";
	if     ( mTrgSelection==0 ) cout << " (all MB triggers)" << endl; 
	else if( mTrgSelection==1 ) cout << " (MB triggers except period-2)" << endl; 
	else                        cout << " trigger selection is strange!!" << endl; 
	cout << " Trigger efficiency correction = " << mTrgEffCorrection << endl;
	cout << " Track efficiency correction = " << mTrkEffCorrection << endl;
	cout << " z-vertex cut:  |vz|<=" << mVzAbsMax << endl;
	cout << " radial vertex cut:  vr<=" << mVrAbsMax << endl;
	cout << " DCA cut     :  dca<=" << mDcaMax << " (not for run16)" << endl;
	cout << " nHitsFit cut:  nHitsFit>=" << mNHitsFitMax << endl;
	cout << " Track type (0=primary, 1=global) : "<< mTrkEP << endl;
	cout << "=============================================================" << endl;

	for( int ich=0; ich<2; ich++ ){
		pEvt [ich].reserveTrkCap(200);
		piEvt[ich].reserveTrkCap(300);
	}

	//// only for FXT31GeV
	////TFile *fcal = TFile::Open("/direct/star+u/taknn/ana/QABES2/jobs/submit/Out/31GeV_qtest_step1.root");
	//TFile *fcal = TFile::Open("/direct/star+u/taknn/ana/QABES2/jobs/submit/Out/31GeV_qtest_step1_w.root");
	//if( !fcal ) return kStOK;
	//for( int ih=0; ih<2; ih++ ){
	//	for( int iep=0; iep<2; iep++ ){
	//		for( int ict=0; ict<10; ict++ ){
	//			for( int ixy=0; ixy<2; ixy++ ){
	//				hInRunvsQvPar[ih][iep][ict][ixy] = (TProfile*)fcal->Get( Form("hRunvsQvPar_%d_%d_%d_%d",ih,iep,ict,ixy) );
	//				if( !hInRunvsQvPar[ih][iep][ict][ixy] ) return kStOK;
	//			}
	//		}
	//	}
	//}
	

	return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StBES2QaMaker::Finish() {

	if(mOutName!="") {
		TFile *fout = new TFile(mOutName.Data(),"RECREATE");
		fout->cd();
		WriteHistograms();
		fout->Close();
	}
	cout << "StBES2QaMaker::Finish()" << endl;
	cout << "# of processed events = " << mEventCounter << endl;
	return kStOK;
}

//-----------------------------------------------------------------------------
void StBES2QaMaker::DeclareHistograms() {

	// Event QA histograms
	// __________________________________________________________________
	// Vertex 
	if( mFixedMode==0 ){ 
		hVz = new TH1F("hVz","hVz;V_{z}^{TPC} [cm] ",300,-150,150);
		hVz_VrCut = new TH1F("hVz_VrCut","hVz_VrCut;V_{z}^{TPC} [cm] ",300,-150,150);
		hVzVpdz = new TH2F("hVzVpdz","hVzVpdz;V_{z}^{TPC} [cm] ;V_{z}^{vpd} [cm] ",140,-35,35,140,-35,35);
		hVzVr = new TH2F("hVzVr","hVzVr;V_{z} [cm] ;V_{r} [cm] ",200,-150,150,80,0,8);
		hVxVy   = new TH2F("hVxVy","hVxVy;V_{x} [cm] ;V_{y} [cm] ",100,-10,10,100,-10,10);
	}else{                
		hVz = new TH1F("hVz","hVz;V_{z}^{TPC} [cm] ",200,150,250);
		hVz_VrCut = new TH1F("hVz_VrCut","hVz_VrCut;V_{z}^{TPC} [cm] ",200,150,250);
		hVzVpdz = new TH2F("hVzVpdz","hVzVpdz;V_{z}^{TPC} [cm] ;V_{z}^{vpd} [cm] ",120,170,230,120,170,230);
		hVzVr = new TH2F("hVzVr","hVzVr;V_{z} [cm] ;V_{r} [cm] ",100,150,250,80,0,8);
		hVxVy   = new TH2F("hVxVy","hVxVy;V_{x} [cm] ;V_{y} [cm] ",200,-10,10,200,-10,10); // from 9p8 GeV
	}
	hdVz    = new TH1F("hdVz","hdVz; V_{z}^{TPC}-V_{z}^{vpd} [cm] ",120,-30,30);
	hVxVyCut= new TH2F("hVxVyCut","hVxVyCut;V_{x} [cm] ;V_{y} [cm] ",100,-10,10,100,-10,10);


	// multiplicity
	if( mData<9 ) hgrefMult = new TH1F("hgrefMult","hgrefMult;grefMult ",600,0,600);
	else          hgrefMult = new TH1F("hgrefMult","hgrefMult;grefMult ",900,0,900);
	htofMult  = new TH1F("htofMult","htofMult;tofMult ",600,0,2400);
	if( mData==11 ){
		hrefMult  = new TH1F("hrefMult","hrefMult;refMult ",900,0,900);
		hrefMultvsTofMult  = new TH2F( "hrefMultvsTofMult", "hrefMultvsTofMult;refMult ;tofMult ", 200, 0, 900, 200, 0, 2600 );
		hrefMultvsTofMatch = new TH2F( "hrefMultvsTofMatch", "hrefMultvsTofMatch;refMult ;tofMatch ", 200, 0, 900, 200, 0, 900 );
		hrefMultvsETofMult  = new TH2F( "hrefMultvsETofMult", "hrefMultvsETofMult;refMult ;etofMult ", 200, 0, 900, 200, 0, 2400 );
	}else{
		hrefMult  = new TH1F("hrefMult","hrefMult;refMult ",600,0,600);
		hrefMultvsTofMult  = new TH2F( "hrefMultvsTofMult", "hrefMultvsTofMult;refMult ;tofMult ", 200, 0, 600, 200, 0, 2400 );
		hrefMultvsTofMatch = new TH2F( "hrefMultvsTofMatch", "hrefMultvsTofMatch;refMult ;tofMatch ", 200, 0, 600, 200, 0, 600 );
		hrefMultvsETofMult  = new TH2F( "hrefMultvsETofMult", "hrefMultvsETofMult;refMult ;etofMult ", 200, 0, 600, 200, 0, 1600 );
	}
	
	hTrg = new TH1F( "hTrg", "hTrg", ntrg, 0, ntrg );
	for( int k=0; k<2; k++ ) hVtxRank[k] = new TH1F( Form("hVtxRank_%d",k), Form("hVtxRank_%d",k), 50, -5, 5 );
	hnTofT0 = new TH1F( "hnTofT0", "hnTofT0", 100, 0, 200 );

	// EPD 
	hNEpdHits  = new TH1F("hNEpdHits","EPD nhits", 200, 0, 1000);
	hEpdHitAdc = new TH1F("hEpdHitAdc","EPD hits Adc", 120, -0.5, 119.5);
	hEpdEtaPhiEast = new TH2F("hEpdEtaPhiEast", "#eta - #phi;#eta ;#phi ", 100, -5.1,-2.1, 300, -TMath::Pi(), TMath::Pi());
	hEpdEtaPhiWest = new TH2F("hEpdEtaPhiWest", "#eta - #phi;#eta ;#phi ", 100, 2.1, 5.1, 300, -TMath::Pi(), TMath::Pi());
	for( int iring=0; iring<16; iring++ ){
		if( mData==11 ){
			hrefMultvsEpdRingE[iring] = new TH2F( Form("hrefMultvsEpdRingE_%d",iring), Form("hrefMultvsEpdRingE_%d;refMult ;EpdMipE per ring ",iring), 200, 0, 900, 200, 0, 200 );
			hrefMultvsEpdRingW[iring] = new TH2F( Form("hrefMultvsEpdRingW_%d",iring), Form("hrefMultvsEpdRingW_%d;refMult ;EpdMipW per ring ",iring), 200, 0, 900, 200, 0, 200 );

			hrefMultvsEpdRingEraw[iring] = new TH2F( Form("hrefMultvsEpdRingEraw_%d",iring), Form("hrefMultvsEpdRingEraw_%d;refMult ;EpdMipE per ring ",iring), 200, 0, 900, 200, 0, 500 );
			hrefMultvsEpdRingWraw[iring] = new TH2F( Form("hrefMultvsEpdRingWraw_%d",iring), Form("hrefMultvsEpdRingWraw_%d;refMult ;EpdMipW per ring ",iring), 200, 0, 900, 200, 0, 500 );
		}else{
			hrefMultvsEpdRingE[iring] = new TH2F( Form("hrefMultvsEpdRingE_%d",iring), Form("hrefMultvsEpdRingE_%d;refMult ;EpdMipE per ring ",iring), 200, 0, 600, 200, 0, 200 );
			hrefMultvsEpdRingW[iring] = new TH2F( Form("hrefMultvsEpdRingW_%d",iring), Form("hrefMultvsEpdRingW_%d;refMult ;EpdMipW per ring ",iring), 200, 0, 600, 200, 0, 200 );

			hrefMultvsEpdRingEraw[iring] = new TH2F( Form("hrefMultvsEpdRingEraw_%d",iring), Form("hrefMultvsEpdRingEraw_%d;refMult ;EpdMipE per ring ",iring), 200, 0, 600, 200, 0, 200 );
			hrefMultvsEpdRingWraw[iring] = new TH2F( Form("hrefMultvsEpdRingWraw_%d",iring), Form("hrefMultvsEpdRingWraw_%d;refMult ;EpdMipW per ring ",iring), 200, 0, 600, 200, 0, 200 );
		}
	}
	if( mData==11 ) hEpdEvsW = new TH2F( "hEpdEvsW", "hEpdEvsW;EpdMipSumE ;EpdMipSumW ", 200, 0, 1800, 200, 0, 1800 );
	else            hEpdEvsW = new TH2F( "hEpdEvsW", "hEpdEvsW;EpdMipSumE ;EpdMipSumW ", 200, 0, 1000, 200, 0, 1000 );


	// BBC, ZDC
	hBBCEvsW= new TH2F("hBBCEvsW","hBBCEvsW;BbcAdcSumE ;BbcAdcSumW ",100,0,70000,100,0,70000);
	if( mData==11 ) hZDCEvsW= new TH2F("hZDCEvsW","hZDCEvsW;ZdcAdcSumE ;ZdcAdcSumW ",100,0,4200,100,0,4200);
	else            hZDCEvsW= new TH2F("hZDCEvsW","hZDCEvsW;ZdcAdcSumE ;ZdcAdcSumW ",100,0,2200,100,0,2200);

	for( int iep=0; iep<2; iep++ ){
		hBBCvsMult[iep] = new TH2F( Form("hBBCvsMult_%d",iep), Form("hBBCvsMult_%d",iep), 100, 0, 60000, 100, 0, 600 );
		hZDCvsMult[iep] = new TH2F( Form("hZDCvsMult_%d",iep), Form("hZDCvsMult_%d",iep), 200, 0, 4200, 200, 0, 600 );
		hBBCvsGMult[iep] = new TH2F( Form("hBBCvsGMult_%d",iep), Form("hBBCvsGMult_%d",iep), 100, 0, 60000, 100, 0, 600 );
		hZDCvsGMult[iep] = new TH2F( Form("hZDCvsGMult_%d",iep), Form("hZDCvsGMult_%d",iep), 200, 0, 4200, 200, 0, 600 );
		hBBCvsTMult[iep] = new TH2F( Form("hBBCvsTMult_%d",iep), Form("hBBCvsTMult_%d",iep), 100, 0, 60000, 100, 0, 4000 );
	}
	hBBCAdc2D = new TH2F("hBBCAdc2D","hBBCAdc2D",48,0,48,100,0,4200);
	hZDCAdc2D = new TH2F("hZDCAdc2D","hZDCAdc2D",32,0,32,120,0,2400);
	hBBCAdc2DHi = new TH2F("hBBCAdc2DHi","hBBCAdc2DHi",48,0,48,200,3900,4100);
	for( int iring=0; iring<4; iring++ ){
		if( mData==11 ){
			hrefMultvsBBCRingE[iring] = new TH2F( Form("hrefMultvsBBCRingE_%d",iring), Form("hrefMultvsBBCRingE_%d;refMult ;BbcAdcSumE per ring ",iring), 200, 0, 800, 200, 0, 24000);
			hrefMultvsBBCRingW[iring] = new TH2F( Form("hrefMultvsBBCRingW_%d",iring), Form("hrefMultvsBBCRingW_%d;refMult ;BbcAdcSumW per ring ",iring), 200, 0, 800, 200, 0, 24000);
			//hrefMultvsBBCRingE[iring] = new TH2F( Form("hrefMultvsBBCRingE_%d",iring), Form("hrefMultvsBBCRingE_%d;refMult ;BbcAdcSumE per ring ",iring), 200, 0, 600, 200, 0, 54000); // <- 31GeVFXT only for day189
			//hrefMultvsBBCRingW[iring] = new TH2F( Form("hrefMultvsBBCRingW_%d",iring), Form("hrefMultvsBBCRingW_%d;refMult ;BbcAdcSumW per ring ",iring), 200, 0, 600, 200, 0, 54000);
		}else{
			hrefMultvsBBCRingE[iring] = new TH2F( Form("hrefMultvsBBCRingE_%d",iring), Form("hrefMultvsBBCRingE_%d;refMult ;BbcAdcSumE per ring ",iring), 200, 0, 600, 200, 0, 24000);
			hrefMultvsBBCRingW[iring] = new TH2F( Form("hrefMultvsBBCRingW_%d",iring), Form("hrefMultvsBBCRingW_%d;refMult ;BbcAdcSumW per ring ",iring), 200, 0, 600, 200, 0, 24000);
		}
	}

	// Rate
	for( int k=0; k<3; k++ ){
		double     bbc_max = 1000;
		if( k==2 ) bbc_max = 200;
		hBbcRate[k] = new TH1F( Form("hBbcRate_%d",k), Form("hBbcRate_%d;BBC rate [kHz] ",k), 200, 0, bbc_max );
		hZdcRate[k] = new TH1F( Form("hZdcRate_%d",k), Form("hZdcRate_%d;ZDC rate [kHz] ",k), 200, 0, 100 );
	}
	hBGRate = new TH1F("hBGRate","BGRate",200,0,200);


	// ZDC-SMD
	//=====================================
	for( int iep=0; iep<2; iep++ ){
		hZDCAdc[iep]    = new TH2F( Form("hZDCAdc_%d",iep), Form("hZDCAdc_%d",iep), 20, 0, 20, 400, 0, 2000 );
		hZDCAdcCor[iep] = new TH2F( Form("hZDCAdcCor_%d",iep), Form("hZDCAdcCor_%d",iep), 20, 0, 20, 400, 0, 2000 );
	}
	hBbcAdc   = new TProfile( Form("hBbcAdc"),    Form("hBbcAdc"),    48, 0, 48, 0, 4100 );
	hBbcAdcCor= new TProfile( Form("hBbcAdcCor"), Form("hBbcAdcCor"), 48, 0, 48, 0, 4100 );

	// Vz vs Multiplicity
	//=====================================
	if( mFixedMode==0 ){ 
		hVzvsrefMult = new TH2F("hVzvsrefMult","hVzvsrefMult;V_{z} [cm] ;refMult ",200,-100,100,200,0,600);
		hVzvsEpdMipE = new TH2F("hVzvsEpdMipE","hVzvsEpdMipE;V_{z} [cm] ;EpdMipE ",200,-100,100,200,0,1000);
		hVzvsEpdMipW = new TH2F("hVzvsEpdMipW","hVzvsEpdMipW;V_{z} [cm] ;EpdMipW ",200,-100,100,200,0,1000);
		hprVzvsrefMult = new TProfile("hprVzvsrefMult","hprVzvsrefMult;V_{z} [cm] ;<refMult> ",200,-100,100);
		hprVzvsEpdMipE = new TProfile("hprVzvsEpdMipE","hprVzvsEpdMipE;V_{z} [cm] ;<EpdMipE> ",200,-100,100);
		hprVzvsEpdMipW = new TProfile("hprVzvsEpdMipW","hprVzvsEpdMipW;V_{z} [cm] ;<EpdMipW> ",200,-100,100);
	}else{
		hVzvsrefMult = new TH2F("hVzvsrefMult","hVzvsrefMult;V_{z} [cm] ;refMult ",200,160,240,200,0,600);
		hVzvsEpdMipE = new TH2F("hVzvsEpdMipE","hVzvsEpdMipE;V_{z} [cm] ;EpdMipE ",200,160,240,200,0,1000);
		hVzvsEpdMipW = new TH2F("hVzvsEpdMipW","hVzvsEpdMipW;V_{z} [cm] ;EpdMipW ",200,160,240,200,0,1000);
		hprVzvsrefMult = new TProfile("hprVzvsrefMult","hprVzvsrefMult;V_{z} [cm] ;<refMult> ",200,160,240);
		hprVzvsEpdMipE = new TProfile("hprVzvsEpdMipE","hprVzvsEpdMipE;V_{z} [cm] ;<EpdMipE> ",200,160,240);
		hprVzvsEpdMipW = new TProfile("hprVzvsEpdMipW","hprVzvsEpdMipW;V_{z} [cm] ;<EpdMipW> ",200,160,240);
	}

	// Mtd and BEMC
	//=====================================
	hMtdHit = new TH2F("hMtdHit","hMtdHit;backleg-id ;(module-id)*ncell + cell-id ",30,0,30,60,0,60);
	hBTowHit = new TH2F("hBTowHit","hBTowHit;phi ;#eta ",100,-TMath::Pi(),TMath::Pi(),100,-1,1);

	// track info
	//=====================================
	float etaMin = -1.8; 
	float etaMax =  1.8; 
	if( mFixedMode==1 ){
		etaMin = -2.6;
		etaMax =  1.0;
	}
	for( int ich=0; ich<2; ich++ ){ 
		hPPhivsEta[ich]  = new TH2F( Form("hPPhivsEta_%d",ich), Form("hPPhivsEta_%d;#phi [rad] ;#eta ",ich), 200, -TMath::Pi(), TMath::Pi(), 200, etaMin, etaMax );
		hGPhivsEta[ich]  = new TH2F( Form("hGPhivsEta_%d",ich), Form("hGPhivsEta_%d;#phi [rad] ;#eta ",ich), 200, -TMath::Pi(), TMath::Pi(), 200, etaMin, etaMax );
		hGPhivsEta2[ich] = new TH2F( Form("hGPhivsEta2_%d",ich), Form("hGPhivsEta2_%d;#phi [rad] ;#eta ",ich), 200, -TMath::Pi(), TMath::Pi(), 200, etaMin, etaMax );
		hPPhivsEtaTof[ich]  = new TH2F( Form("hPPhivsEtaTof_%d",ich), Form("hPPhivsEtaTof_%d;#phi [rad] ;#eta ",ich), 200, -TMath::Pi(), TMath::Pi(), 200, etaMin, etaMax );
	}
	hPPhivsEtaMtd   = new TH2F( Form("hPPhivsEtaMtd"), Form("hPPhivsEtaMtd;#phi [rad] ;#eta "), 200, -TMath::Pi(), TMath::Pi(), 120, -1.2, 1.2 );
	hPPhivsEtaBemc  = new TH2F( Form("hPPhivsEtaBemc"), Form("hPPhivsEtaBemc;#phi [rad] ;#eta "), 200, -TMath::Pi(), TMath::Pi(), 120, -1.2, 1.2 );

	hPdEdx    = new TH1F( "hPdEdx", "hPdEdx;dE/dx [keV/cm] ", 100, 0, 10 );
	hGdEdx    = new TH1F( "hGdEdx", "hGdEdx;dE/dx [keV/cm] ", 100, 0, 10 );
	hPdEdxvsP = new TH2F( "hPdEdxvsP",  "hPdEdxvsP;p*charge [GeV/c] ;dE/dx [keV/cm] ", 200, -5, 5, 200, 1, 15 );
	hGdEdxvsP = new TH2F( "hGdEdxvsP",  "hGdEdxvsP;p*charge [GeV/c] ;dE/dx [keV/cm] ", 200, -5, 5, 200, 1, 15 );

	hNhitsFit = new TH1F("hNhitsFit","hNhitsFit;nHitsFit ",80,0,80);
	hNhitsMax = new TH1F("hNhitsMax","hNhitsMax;nHitsMax ",80,0,80);
	hgDca  = new TH1F( Form("hgDca"), Form("hgDca;dca(global) [cm] "), 150, 0, 15);
	hDca   = new TH1F( Form("hDca"), Form("hDca;dca(primary) [cm] "), 100, 0, 5);
	hDcaXY = new TH1F( Form("hDcaXY"), Form("hDcaXY;dca_{XY} [cm] "), 100, 0, 5);
	hDcaZ  = new TH1F( Form("hDcaZ"), Form("hDcaZ;dcaZ [cm] "), 100, -5, 5);
	for( int ich=0; ich<2; ich++ ){ 
		hPPt[ich] = new TH1F( Form("hPPt_%d",ich), Form("hP:Pt_%d;p_{T} [GeV/c] ",ich), 100, 0, 10 );
		hGPt[ich] = new TH1F( Form("hGPt_%d",ich), Form("hGPt_%d;p_{T} [GeV/c] ",ich), 100, 0, 10 );
	}
	hRatioPrGlTrk = new TH1F("hRatioPrGlTrk","hRatioPrGlTrk",100,0,1);

	hM2vsP     = new TH2F( Form("hM2vsP"), Form("hM2vsP;p*charge [GeV/c] ;mass^{2} [GeV^{2}/c^{4}] "), 400, -5, 5, 200, -0.5, 1.5 );
	hM2vsPdca  = new TH2F( Form("hM2vsPdca"), Form("hM2vsPdca;p*charge [GeV/c] ;mass^{2} [GeV^{2}/c^{4}] "), 400, -5, 5, 200, -0.5, 1.5 );
	hBetaIvsP  = new TH2F( "hBetaIvsP", "hBetaIvsP", 200, -5, 5, 200, 0.6, 2.8 );
	hBetaIvsPe = new TH2F( "hBetaIvsPe", "hBetaIvsPe", 200, -5, 5, 200, 0.6, 2.8 );
	hETOFxy    = new TH2F( "hETOFxy", "hETOFxy", 200, -250, 250, 200, -250, 250 );


	for( int k=0; k<6; k++ ){
		hDcavsP[k]  = new TH2F( Form("hDcavsP_%d",k), Form("hDcavsP_%d;p*charge [GeV/c] ;dca [cm] ",k), 100, -2.5, 2.5, 200, 0, 10 );
	}

	string spid[4] = { "#pi", "K", "p", "e" };
	for( int k=0; k<4; k++ ){
		hNsigmaP[k] = new TH2F( Form("hNsigmaP_%d",k), Form("hNsigmaP_%d;p [GeV/c] ;n#sigma_{%s} ",k,spid[k].c_str()), 100, 0, 5, 280, -35, 35  );
	}

	for( int ich=0; ich<2; ich++ ){
		hLambdaMass[ich] = new TH1F( Form("hLambdaMass_%d",ich), Form("hLambdaMass_%d;M_{inv} [GeV/c^{2}] ",ich), 200, 1.07, 1.17 );
		hLambdaMassTopoCut[ich] = new TH1F( Form("hLambdaMassTopoCut_%d",ich), Form("hLambdaMassTopoCut_%d;M_{inv} [GeV/c^{2}] ",ich), 200, 1.07, 1.17 );
	}


	// EP calibration
	//=====================================
	if( mCalibMode>0 && mCalibMode<3 ){

		// SMD EP parameters
		for( int iep=0; iep<NsubSMD; iep++ ){ 
			hSMDrecX[iep] = new TProfile( Form("hSMDrecX_%d",iep), Form("hSMDrecX_%d",iep), NqvCent*NqvZvtx, 0.5, NqvCent*NqvZvtx+0.5, -9000, 9000, "s" );
			hSMDrecY[iep] = new TProfile( Form("hSMDrecY_%d",iep), Form("hSMDrecY_%d",iep), NqvCent*NqvZvtx, 0.5, NqvCent*NqvZvtx+0.5, -9000, 9000, "s" );

			if( mCalibMode==2 ){
				hSMDfltC[iep] = new TProfile( Form("hSMDfltC_%d",iep), Form("hSMDfltC_%d",iep), Nflt*NqvCent*NqvZvtx, 0.5, Nflt*NqvCent*NqvZvtx+0.5, -10, 10 );
				hSMDfltS[iep] = new TProfile( Form("hSMDfltS_%d",iep), Form("hSMDfltS_%d",iep), Nflt*NqvCent*NqvZvtx, 0.5, Nflt*NqvCent*NqvZvtx+0.5, -10, 10 );
		}	}

		// BBC EP parameters
		for( int ith=0; ith<NordBBC; ith++ ){
			for( int iep=0; iep<NsubBBC; iep++ ){
				hBBCrecX[ith][iep] = new TProfile( Form("hBBCrecX_%d_%d",ith,iep), Form("hBBCrecX_%d_%d",ith,iep), NqvCent*NqvZvtx, 0.5, NqvCent*NqvZvtx+0.5, -9000, 9000, "s" );
				hBBCrecY[ith][iep] = new TProfile( Form("hBBCrecY_%d_%d",ith,iep), Form("hBBCrecY_%d_%d",ith,iep), NqvCent*NqvZvtx, 0.5, NqvCent*NqvZvtx+0.5, -9000, 9000, "s" );

				if( mCalibMode==2 ){
					hBBCfltC[ith][iep] = new TProfile( Form("hBBCfltC_%d_%d",ith,iep), Form("hBBCfltC_%d_%d",ith,iep), Nflt*NqvCent*NqvZvtx, 0.5, Nflt*NqvCent*NqvZvtx+0.5, -10, 10 );
					hBBCfltS[ith][iep] = new TProfile( Form("hBBCfltS_%d_%d",ith,iep), Form("hBBCfltS_%d_%d",ith,iep), Nflt*NqvCent*NqvZvtx, 0.5, Nflt*NqvCent*NqvZvtx+0.5, -10, 10 );
		}	}	}

		// TPC EP parameters
		for( int ith=0; ith<NordTPC; ith++ ){
			for( int iep=0; iep<NsubTPC; iep++ ){
				hTPCrecX[ith][iep] = new TProfile( Form("hTPCrecX_%d_%d",ith,iep), Form("hTPCrecX_%d_%d",ith,iep), NqvCent*NqvZvtx, 0.5, NqvCent*NqvZvtx+0.5, -9000, 9000, "s" );
				hTPCrecY[ith][iep] = new TProfile( Form("hTPCrecY_%d_%d",ith,iep), Form("hTPCrecY_%d_%d",ith,iep), NqvCent*NqvZvtx, 0.5, NqvCent*NqvZvtx+0.5, -9000, 9000, "s" );

				if( mCalibMode==2 ){
					hTPCfltC[ith][iep] = new TProfile( Form("hTPCfltC_%d_%d",ith,iep), Form("hTPCfltC_%d_%d",ith,iep), Nflt*NqvCent*NqvZvtx, 0.5, Nflt*NqvCent*NqvZvtx+0.5, -10, 10 );
					hTPCfltS[ith][iep] = new TProfile( Form("hTPCfltS_%d_%d",ith,iep), Form("hTPCfltS_%d_%d",ith,iep), Nflt*NqvCent*NqvZvtx, 0.5, Nflt*NqvCent*NqvZvtx+0.5, -10, 10 );
		}	}	}

		// EPD EP parameters
		for( int ith=0; ith<NordEPD; ith++ ){
			for( int iep=0; iep<NsubEPD; iep++ ){
				hEPDrecX[ith][iep] = new TProfile( Form("hEPDrecX_%d_%d",ith,iep), Form("hEPDrecX_%d_%d",ith,iep), NqvCent*NqvZvtx, 0.5, NqvCent*NqvZvtx+0.5, -9000, 9000, "s" );
				hEPDrecY[ith][iep] = new TProfile( Form("hEPDrecY_%d_%d",ith,iep), Form("hEPDrecY_%d_%d",ith,iep), NqvCent*NqvZvtx, 0.5, NqvCent*NqvZvtx+0.5, -9000, 9000, "s" );

				if( mCalibMode==2 ){
					hEPDfltC[ith][iep] = new TProfile( Form("hEPDfltC_%d_%d",ith,iep), Form("hEPDfltC_%d_%d",ith,iep), Nflt*NqvCent*NqvZvtx, 0.5, Nflt*NqvCent*NqvZvtx+0.5, -10, 10 );
					hEPDfltS[ith][iep] = new TProfile( Form("hEPDfltS_%d_%d",ith,iep), Form("hEPDfltS_%d_%d",ith,iep), Nflt*NqvCent*NqvZvtx, 0.5, Nflt*NqvCent*NqvZvtx+0.5, -10, 10 );
		}	}	}

	}//mCalibMode

	// EP distributions
	//=====================================
	if( mCalibMode==3 ){
		//SMD
		for( int iep=0; iep<NsubSMD; iep++ ){ 
			for( int ict=0; ict<NqvCent; ict++ ){
				for( int ic=0; ic<3; ic++ ){
					hSMDEP[iep][ict][ic] = new TH1F( Form("hSMDEP_%d_%d_%d",iep,ict,ic), Form("hSMDEP_%d_%d_%d",iep,ict,ic), 100, -TMath::Pi(), TMath::Pi() );
		}	}	}

		//BBC
		for( int ith=0; ith<NordBBC; ith++ ){
			for( int iep=0; iep<NsubBBC; iep++ ){
				for( int ict=0; ict<NqvCent; ict++ ){
					for( int icc=0; icc<3; icc++ ){
						hBBCEP[ith][iep][ict][icc] = new TH1F( Form("hBBCEP_%d_%d_%d_%d",ith,iep,ict,icc), Form("hBBCEP_%d_%d_%d_%d",ith,iep,ict,icc), 100, -TMath::Pi()/(ith+1), TMath::Pi()/(ith+1) );
		}	}	}	}

		//TPC
		for( int ith=0; ith<NordTPC; ith++ ){
			for( int iep=0; iep<NsubTPC; iep++ ){
				for( int ict=0; ict<NqvCent; ict++ ){
					for( int ic=0; ic<3; ic++ ){
						hTPCEP[ith][iep][ict][ic] = new TH1F( Form("hTPCEP_%d_%d_%d_%d",ith,iep,ict,ic), Form("hTPCEP_%d_%d_%d_%d",ith,iep,ict,ic), 100, -TMath::Pi()/(ith+1.), TMath::Pi()/(ith+1.) );
		}	}	}	}

		//EPD
		for( int ith=0; ith<NordEPD; ith++ ){
			for( int iep=0; iep<NsubEPD; iep++ ){
				for( int ict=0; ict<NqvCent; ict++ ){
					for( int ic=0; ic<3; ic++ ){
						hEPDEP[ith][iep][ict][ic] = new TH1F( Form("hEPDEP_%d_%d_%d_%d",ith,iep,ict,ic), Form("hEPDEP_%d_%d_%d_%d",ith,iep,ict,ic), 100, -TMath::Pi()/(ith+1.), TMath::Pi()/(ith+1.) );
		}	}	}	}

	}//mCalibMode

/*
	// Q-vector check
	//=====================================
	// ZDC-SMD
	for( int iep=0; iep<NsubSMD; iep++ ){ 
		for( int ic=0; ic<2; ic++ ){
			hQSMD2D[iep][ic] = new TH2F( Form("hQSMD2D_%d_%d",iep,ic), Form("hQSMD2D_%d_%d",iep,ic), 200, -20, 20, 200, -20, 20 );
		}
	}

	// BBC Q-vector
	for( int ith=0; ith<2; ith++ ){
		for( int iep=0; iep<NsubBBC; iep++ ){
			for( int ict=0; ict<NCent; ict++ ){
				hQBBC[ith][iep][ict] = new TH1F( Form("hQBBC_%d_%d_%d",ith,iep,ict), Form("hQBBC_%d_%d_%d",ith,iep,ict), 100, 0, 10. );
	}	 }	}

	// TPC Q-vector
	for( int ith=0; ith<3; ith++ ){
		for( int iep=0; iep<NsubTPC; iep++ ){
			for( int ict=0; ict<NCent; ict++ ){
				hQTPC[ith][iep][ict] = new TH1F( Form("hQTPC_%d_%d_%d",ith,iep,ict), Form("hQTPC_%d_%d_%d",ith,iep,ict), 100, 0, 10 );
	}	 }	}


	// EP correlation
	for( int ith=0; ith<NordSMD; ith++ ){
		for( int ics=0; ics<2; ics++ ){
			hSMDEPCrrvsMult[ith][ics] = new TProfile( Form("hSMDEPCrrvsMult_%d_%d",ith,ics), Form("hSMDEPCrrvsMult_%d_%d",ith,ics), NCent, 0, NCent, -10, 10 );
			if( mTrgEffCorrection ) hSMDEPCrrvsMult[ith][ics]->Sumw2();

	}	}

	for( int ith=0; ith<NordBBC; ith++ ){
		for( int ics=0; ics<4; ics++ ){
			hBBCEPCrrvsMult[ith][ics] = new TProfile( Form("hBBCEPCrrvsMult_%d_%d",ith,ics), Form("hBBCEPCrrvsMult_%d_%d",ith,ics), NCent, 0, NCent, -10, 10 );
			if( mTrgEffCorrection ) hBBCEPCrrvsMult[ith][ics]->Sumw2();
	}	}

*/


	for( int ih=0; ih<2; ih++ ){
		for( int iep=0; iep<4; iep++ ){
			hQvCor[ih][iep] = new TProfile( Form("hQvCor_%d_%d",ih,iep), Form("hQvCor_%d_%d",ih,iep), 11, 0, 11 );
		}
	}

	for( int ih=0; ih<2; ih++ ){
		for( int iep=0; iep<2; iep++ ){
			for( int ict=0; ict<10; ict++ ){
				for( int ixy=0; ixy<2; ixy++ ){
					hRunvsQvPar[ih][iep][ict][ixy] = new TProfile( Form("hRunvsQvPar_%d_%d_%d_%d",ih,iep,ict,ixy), Form("hRunvsQvPar_%d_%d_%d_%d",ih,iep,ict,ixy), mRunbins, mRunmin, mRunmax, -10, 10, "s" );
				}
			}
		}
	}


	// run dependece
	// event info.
	hRunidvsrefMult = new TProfile( "hRunidvsrefMult", "hRunidvsrefMult;RunId ;<refMult> ", mRunbins, mRunmin, mRunmax );
	hRunidvsgrefMult = new TProfile( "hRunidvsgrefMult", "hRunidvsgrefMult;RunId ;<grefMult> ", mRunbins, mRunmin, mRunmax );
	hRunidvstofMult = new TProfile( "hRunidvstofMult", "hRunidvstofMult;RunId ;<tofMult> ", mRunbins, mRunmin, mRunmax );
	hRunidvstofMatch = new TProfile( "hRunidvstofMatch", "hRunidvstofMatch;RunId ;<nBTOFMatch> ", mRunbins, mRunmin, mRunmax );
	hRunidvsBemcMatch = new TProfile( "hRunidvsBemcMatch", "hRunidvsBemcMatch;RunId ;<nBMECMatch> ", mRunbins, mRunmin, mRunmax );
	hRunidvsEpdMipE = new TProfile( "hRunidvsEpdMipE", "hRunidvsEpdMipE;RunId ;<EpdMipE> ", mRunbins, mRunmin, mRunmax );
	hRunidvsEpdMipW = new TProfile( "hRunidvsEpdMipW", "hRunidvsEpdMipW;RunId ;<EpdMipW> ", mRunbins, mRunmin, mRunmax );
	hRunidvsEtofMult = new TProfile( "hRunidvsEtofMult", "hRunidvsEtofMult;RunId ;<nHitETof> ", mRunbins, mRunmin, mRunmax );

	hRunidvsVx = new TProfile( "hRunidvsVx", "hRunidvsVx;RunId ;<Vx> ", mRunbins, mRunmin, mRunmax );
	hRunidvsVy = new TProfile( "hRunidvsVy", "hRunidvsVy;RunId ;<Vy> ", mRunbins, mRunmin, mRunmax );
	hRunidvsVz = new TProfile( "hRunidvsVz", "hRunidvsVz;RunId ;<Vz> ", mRunbins, mRunmin, mRunmax );
	hRunidvsRank = new TProfile( "hRunidvsRank", "hRunidvsRank;RunId ;<Rank> ", mRunbins, mRunmin, mRunmax );

	hRunidvsZdcx = new TProfile( "hRunidvsZdcx", "hRunidvsZdcx;RunId ;<ZDC rate> ", mRunbins, mRunmin, mRunmax );
	hRunidvsBbcx = new TProfile( "hRunidvsBbcx", "hRunidvsBbcx;RunId ;<BBC rate> ", mRunbins, mRunmin, mRunmax );
	hRunidvsBgrate = new TProfile( "hRunidvsBgrate", "hRunidvsBgrate;RunId ;<BG rate> ", mRunbins, mRunmin, mRunmax );

	for( int ih=0; ih<2; ih++ ){
		hRunidvsTpcQx[ih] = new TProfile( Form("hRunidvsTpcQx_%d",ih), Form("hRunidvsTpcQx_%d;RunId ;TPC <Q_{x,%d}> ",ih,ih+1), mRunbins, mRunmin, mRunmax );
		hRunidvsTpcQy[ih] = new TProfile( Form("hRunidvsTpcQy_%d",ih), Form("hRunidvsTpcQy_%d;RunId ;TPC <Q_{y,%d}> ",ih,ih+1), mRunbins, mRunmin, mRunmax );
	}

	for( int ih=0; ih<2; ih++ ){
		for( int isub=0; isub<2; isub++ ){
			if( isub==0 ){
				hRunidvsEpdQx[ih][isub] = new TProfile( Form("hRunidvsEpdQx_%d_%d",ih,isub), Form("hRunidvsEpdQx_%d_%d;RunId ;EPD.E <Q_{x,%d}> ",ih,isub,ih+1), mRunbins, mRunmin, mRunmax );
				hRunidvsEpdQy[ih][isub] = new TProfile( Form("hRunidvsEpdQy_%d_%d",ih,isub), Form("hRunidvsEpdQy_%d_%d;RunId ;EPD.E <Q_{y,%d}> ",ih,isub,ih+1), mRunbins, mRunmin, mRunmax );
			}else{
				hRunidvsEpdQx[ih][isub] = new TProfile( Form("hRunidvsEpdQx_%d_%d",ih,isub), Form("hRunidvsEpdQx_%d_%d;RunId ;EPD.W <Q_{x,%d}> ",ih,isub,ih+1), mRunbins, mRunmin, mRunmax );
				hRunidvsEpdQy[ih][isub] = new TProfile( Form("hRunidvsEpdQy_%d_%d",ih,isub), Form("hRunidvsEpdQy_%d_%d;RunId ;EPD.W <Q_{y,%d}> ",ih,isub,ih+1), mRunbins, mRunmin, mRunmax );
			}
		}
	}

	hRunidvsAllEvt = new TH1F( "hRunidvsAllEvt", "hRunidvsAllEvt;RunId ; # of events ", mRunbins, mRunmin, mRunmax );
	for( int k=0; k<4; k++ ) hRunidvsGoodEvt[k] = new TH1F( Form("hRunidvsGoodEvt_%d",k), Form("hRunidvsGoodEvt_%d;RunId ; # of good events ",k), mRunbins, mRunmin, mRunmax );
	
	// track info.
	hRunidvsPt         = new TProfile( "hRunidvsPt", "hRunidvsPt;RunId ;<p_{T}> [GeV/c] ", mRunbins, mRunmin, mRunmax );
	hRunidvsEta        = new TProfile( "hRunidvsEta", "hRunidvsEta;RunId ;<#eta> ", mRunbins, mRunmin, mRunmax );
	hRunidvsPhi        = new TProfile( "hRunidvsPhi", "hRunidvsPhi;RunId ;<#phi> ", mRunbins, mRunmin, mRunmax );
	hRunidvsNchp       = new TProfile( "hRunidvsNchp", "hRunidvsNchp;RunId ;<N_{trk}>, q>0 ", mRunbins, mRunmin, mRunmax );
	hRunidvsNchm       = new TProfile( "hRunidvsNchm", "hRunidvsNchm;RunId ;<N_{trk}>, q<0 ", mRunbins, mRunmin, mRunmax );
	hRunidvsNproton    = new TProfile( "hRunidvsNproton", "hRunidvsNproton;RunId ;<N_{p}> ", mRunbins, mRunmin, mRunmax );
	hRunidvsNprotonbar = new TProfile( "hRunidvsNprotonbar", "hRunidvsNprotonbar;RunId ;<N_{pbar}> ", mRunbins, mRunmin, mRunmax );
	hRunidvsNhits      = new TProfile( "hRunidvsNhits", "hRunidvsNhits;RunId ;<NhitsFit> ", mRunbins, mRunmin, mRunmax );
	hRunidvsDedx       = new TProfile( "hRunidvsDedx", "hRunidvsDedx;RunId ;<dEdx> ", mRunbins, mRunmin, mRunmax );
	hRunidvsDca        = new TProfile( "hRunidvsDca", "hRunidvsDca;RunId ;<dca> ", mRunbins, mRunmin, mRunmax );
	if( mFixedMode==1 ) hRunidvsEtaDist = new TH2F( "hRunidvsEtaDist", "hRunidvsEtaDist;RunId ;#eta ", mRunbins, mRunmin, mRunmax, 90, -2.6, 1.0 );
	else                hRunidvsEtaDist = new TH2F( "hRunidvsEtaDist", "hRunidvsEtaDist;RunId ;#eta ", mRunbins, mRunmin, mRunmax, 90, -1.8, 1.8 );

	hPt    = new TH1F("hPt","hPt",200,0,10);
	hPhi   = new TH1F("hPhi","hPhi",200,-3.2,3.2);
	hNhits = new TH1F("hNhits","hNhits",80,0,80);
	hDedx  = new TH1F("hDedx","hDedx",200,0,30);
	hDcaTmp= new TH1F("hDcaTmp","hDcaTmp",200,0,5);
	if( mFixedMode==1 ) hEta = new TH1F("hEta","hEta",200,-3,1);
	else                hEta = new TH1F("hEta","hEta",200,-2,2);

	hVzAll   = new TH1F("hVzAll","hVzAll",400,-200,200);

	// only for test
	for( int iq=0; iq<2; iq++ ){
		hAveDedx2D [iq] = new TProfile2D( Form("hAveDedx2D_%d",iq), Form("hAveDedx2D_%d",iq), 200, -3.2, 3.2, 200, -2, 2 );
		hAveDedxPhi[iq] = new TProfile( Form("hAveDedx_%d",iq), Form("hAveDedx_%d",iq), 200, -3.2, 3.2 );
		hDedxPhi   [iq] = new TH2F( Form("hDedxPhi_%d",iq), Form("hDedxPhi_%d",iq), 100, -3.2, 3.2, 200, 0, 30 );
	}
}

//-----------------------------------------------------------------------------
void StBES2QaMaker::WriteHistograms() {

	// QA plots


	// event info
	// __________________________________
	// Vertex
	hVz      ->Write();
	hdVz     ->Write();
	hVxVy    ->Write();
	hVxVyCut ->Write();
	hVzVpdz  ->Write();
	hVz_VrCut->Write();
	hVzVr    ->Write();

	// multiplicity
	hrefMult          ->Write();
	hgrefMult         ->Write();
	htofMult          ->Write();
	hrefMultvsTofMult ->Write();
	hrefMultvsTofMatch->Write();
	hrefMultvsETofMult->Write();

	hTrg->Write();
	for( int k=0; k<2; k++ ) hVtxRank[k]->Write();
	hnTofT0->Write();

	// BBC/ZDC
	hBBCEvsW->Write();
	hZDCEvsW->Write();
	for( int iep=0; iep<2; iep++ ){
		hBBCvsMult[iep]->Write();
		hZDCvsMult[iep]->Write();
		hBBCvsGMult[iep]->Write();
		hZDCvsGMult[iep]->Write();
		hBBCvsTMult[iep]->Write();
	}
	hBBCAdc2D->Write();
	hZDCAdc2D->Write();
	hBBCAdc2DHi->Write();
	for( int iring=0; iring<4; iring++ ){
		hrefMultvsBBCRingE[iring]->Write(); 
		hrefMultvsBBCRingW[iring]->Write(); 
	}

	for( int k=0; k<3; k++ ){
		hZdcRate[k]->Write();
		hBbcRate[k]->Write();
	}
	hBGRate->Write();

	for( int iep=0; iep<2; iep++ ){
		hZDCAdc   [iep]->Write();
		hZDCAdcCor[iep]->Write();
	}
	hBbcAdc   ->Write();
	hBbcAdcCor->Write();

	// EPD
	hNEpdHits ->Write();
	hEpdHitAdc->Write();
	hEpdEtaPhiEast->Write();
	hEpdEtaPhiWest->Write(); 
	for( int iring=0; iring<16; iring++ ){
		hrefMultvsEpdRingE[iring]->Write();
		hrefMultvsEpdRingW[iring]->Write();

		hrefMultvsEpdRingEraw[iring]->Write();
		hrefMultvsEpdRingWraw[iring]->Write();
	}
	hEpdEvsW->Write();

	hVzvsrefMult->Write(); 
	hVzvsEpdMipE->Write(); 
	hVzvsEpdMipW->Write(); 
	hprVzvsrefMult->Write(); 
	hprVzvsEpdMipE->Write(); 
	hprVzvsEpdMipW->Write(); 

	hMtdHit->Write();
	hBTowHit->Write();

	// track info
	for( int ich=0; ich<2; ich++ ){ 
		hPPt[ich]->Write();
		hGPt[ich]->Write();
		hPPhivsEta[ich]->Write();
		hGPhivsEta[ich]->Write();
		hGPhivsEta2[ich]->Write();
		hPPhivsEtaTof[ich]->Write();
	}
	hPPhivsEtaMtd->Write();
	hPPhivsEtaBemc->Write();

	hPdEdx   ->Write();
	hGdEdx   ->Write();
	hPdEdxvsP->Write();
	hGdEdxvsP->Write();
	hNhitsFit->Write();
	hNhitsMax->Write();
	hgDca    ->Write();
	hDca     ->Write();
	hDcaXY   ->Write();
	hDcaZ    ->Write();
	hRatioPrGlTrk->Write();


	hM2vsP->Write();
	//hM2vsPdca->Write();
	for( int k=0; k<6; k++ ) hDcavsP[k]->Write();
	for( int k=0; k<4; k++ ) hNsigmaP[k]->Write();
	hBetaIvsP->Write();
	hBetaIvsPe->Write();
	hETOFxy->Write();
	
	for( int ich=0; ich<2; ich++ ){
		hLambdaMass       [ich]->Write();
		hLambdaMassTopoCut[ich]->Write();
	}



/*
	// EP calibration
	//=====================================
	if( mCalibMode>0 && mCalibMode<3 ){
		// SMD
		for( int iep=0; iep<NsubSMD; iep++ ){ 
			hSMDrecX[iep]->Write();
			hSMDrecY[iep]->Write();
			if( mCalibMode==2 ){
				hSMDfltC[iep]->Write();
				hSMDfltS[iep]->Write();
		}	}

		// BBC
		for( int ith=0; ith<NordBBC; ith++ ){
			for( int iep=0; iep<NsubBBC; iep++ ){
				hBBCrecX[ith][iep]->Write();
				hBBCrecY[ith][iep]->Write();
				if( mCalibMode==2 ){
					hBBCfltC[ith][iep]->Write();
					hBBCfltS[ith][iep]->Write();
		}	}	}

		// TPC
		for( int ith=0; ith<NordTPC; ith++ ){
			for( int iep=0; iep<NsubTPC; iep++ ){
				hTPCrecX[ith][iep]->Write();
				hTPCrecY[ith][iep]->Write();
				if( mCalibMode==2 ){
					hTPCfltC[ith][iep]->Write();
					hTPCfltS[ith][iep]->Write();
		}	}	}
	}

	// EP distributions
	//=====================================
	if( mCalibMode==3 ){
		//SMD
		for( int iep=0; iep<NsubSMD; iep++ ){ 
			for( int ict=0; ict<NqvCent; ict++ ){
				for( int ic=0; ic<3; ic++ ){
					hSMDEP[iep][ict][ic]->Write();
		}	}	}

		//BBC
		for( int ith=0; ith<NordBBC; ith++ ){
			for( int iep=0; iep<NsubBBC; iep++ ){
				for( int ict=0; ict<NqvCent; ict++ ){
					for( int icc=0; icc<3; icc++ ){
						hBBCEP[ith][iep][ict][icc]->Write();
		}	}	}	}
		//TPC
		for( int ith=0; ith<NordTPC; ith++ ){
			for( int iep=0; iep<NsubTPC; iep++ ){
				for( int ict=0; ict<NqvCent; ict++ ){
					for( int ic=0; ic<3; ic++ ){
						hTPCEP[ith][iep][ict][ic]->Write();
		}	}	}	}



	}//mCalibMode
*/

	// Q-vector check
	//=====================================
	/*
	if( mCalibMode>0 ){
		// ZDC-SMD
		for( int iep=0; iep<NsubSMD; iep++ ){ 
			for( int ic=0; ic<2; ic++ ){
				hQSMD2D[iep][ic]->Write();
			}
		}

		// BBC Q-vector
		for( int ith=0; ith<2; ith++ ){
			for( int iep=0; iep<NsubBBC; iep++ ){
				for( int ict=0; ict<NCent; ict++ ){
					hQBBC[ith][iep][ict]->Write();
		}   }	}

		// TPC Q-vector
		for( int ith=0; ith<3; ith++ ){
			for( int iep=0; iep<NsubTPC; iep++ ){
				for( int ict=0; ict<NCent; ict++ ){
					hQTPC[ith][iep][ict]->Write();
		}	}	}
	}
	//=====================================

	if( mCalibMode!=3 ) return;

	for( int ith=0; ith<NordSMD; ith++ ){
		for( int ics=0; ics<2; ics++ ) hSMDEPCrrvsMult[ith][ics]->Write();
	}

	for( int ith=0; ith<NordBBC; ith++ ){
		for( int ics=0; ics<4; ics++ ){
			hBBCEPCrrvsMult[ith][ics]->Write();
		}
	}

	*/
	for( int ih=0; ih<2; ih++ ){
		for( int iep=0; iep<4; iep++ ){
			hQvCor[ih][iep]->Write();
		}
	}
	for( int ih=0; ih<2; ih++ ){
		for( int iep=0; iep<2; iep++ ){
			for( int ict=0; ict<10; ict++ ){
				for( int ixy=0; ixy<2; ixy++ ){
					hRunvsQvPar[ih][iep][ict][ixy]->Write();
				}
			}
		}
	}

	// run dependence
	// event info.
	hRunidvsrefMult ->Write();
	hRunidvsgrefMult->Write();
	hRunidvstofMult ->Write();
	hRunidvstofMatch->Write();
	hRunidvsBemcMatch->Write();
	hRunidvsEpdMipE ->Write();
	hRunidvsEpdMipW ->Write();
	hRunidvsEtofMult->Write();
	hRunidvsVx      ->Write();
	hRunidvsVy      ->Write();
	hRunidvsVz      ->Write();
	hRunidvsRank    ->Write();
	hRunidvsZdcx    ->Write();
	hRunidvsBbcx    ->Write();
	hRunidvsBgrate  ->Write();
	for( int ih=0; ih<2; ih++ ){
		hRunidvsTpcQx[ih]->Write();
		hRunidvsTpcQy[ih]->Write();
	}
	for( int ih=0; ih<2; ih++ ){
		for( int isub=0; isub<2; isub++ ){
			hRunidvsEpdQx[ih][isub]->Write();
			hRunidvsEpdQy[ih][isub]->Write();
		}
	}
	hRunidvsAllEvt->Write();
	for( int k=0; k<4; k++ ) hRunidvsGoodEvt[k]->Write(); 

	// track info.
	hRunidvsPt        ->Write();	
	hRunidvsEta       ->Write();
	hRunidvsPhi       ->Write();
	hRunidvsNhits     ->Write();
	hRunidvsDedx      ->Write();
	hRunidvsDca       ->Write();
	hRunidvsNchp      ->Write(); 
	hRunidvsNchm      ->Write();
	hRunidvsNproton   ->Write();
	hRunidvsNprotonbar->Write();
	hRunidvsEtaDist   ->Write();

	// only for test
	for( int iq=0; iq<2; iq++ ){
		hAveDedx2D [iq]->Write();
		hAveDedxPhi[iq]->Write();
		hDedxPhi   [iq]->Write();
	}
}

//----------------------------------------------------------------------------- 
void StBES2QaMaker::Clear(Option_t *opt) {
	//pEve->InitValue();
}

//----------------------------------------------------------------------------- 
Int_t StBES2QaMaker::Make() {

	// check dst pointers
	if(!mPicoDstMaker) {
		LOG_WARN << " No PicoDstMaker! Skip! " << endm;
		return kStWarn;
	}

	mPicoDst = mPicoDstMaker->picoDst();
	if(!mPicoDst) {
		LOG_WARN << " No PicoDst! Skip! " << endm;
		return kStWarn;
	}

	StPicoEvent *picoEvent = (StPicoEvent*)mPicoDst->event();
	if( !picoEvent ){
		LOG_WARN << " No PicoEvent! Skip! " << endm;
		return kStWarn;
	}


	// Bad run rejection
	Int_t runId   = picoEvent->runId();
	mRunNumber = runId;
	//pRefMultCorr->init(mRunNumber);
	//if( pRefMultCorr->isBadRun(mRunNumber) ) return kStOK; 
	mEventCounter++;

	// event information
	Double_t refMult  = picoEvent->refMult(); 
	Double_t grefMult = picoEvent->grefMult(); 
	UShort_t tofMult  = picoEvent->btofTrayMultiplicity();
	if( mFixedMode ){
		int myrefM  = 0;
		int mygrefM = 0;
		TrackLoop(myrefM,mygrefM);
		refMult  = myrefM;
		grefMult = mygrefM;
	}

	Float_t vx     = picoEvent->primaryVertex().x();
	Float_t vy     = picoEvent->primaryVertex().y();
	Float_t vz     = picoEvent->primaryVertex().z();
	Float_t vpd_vz = picoEvent->vzVpd();
	//cout <<"vx="<<vx<<" vy="<<vy<<" vz="<<vz<< endl;

	Double_t zdcRate = picoEvent->ZDCx(); // zdcCoincidenceRate()
	Double_t bbcRate = picoEvent->BBCx(); // bbcCoincidenceRate()

	Float_t zdcAdcSumE = picoEvent->ZdcSumAdcEast();
	Float_t zdcAdcSumW = picoEvent->ZdcSumAdcWest();


	// EPD info.
	// __________________________________________________
	// based on StEpdEpFinder by Mike 
	Int_t nepdHits = mPicoDst->numberOfEpdHits();
	hNEpdHits->Fill(nepdHits);

	// Hit loop
	float nepdMIPsE = 0;
	float nepdMIPsW = 0;

	float nepdRingE[16] = {0};
	float nepdRingW[16] = {0};
	float nepdRingEraw[16] = {0};
	float nepdRingWraw[16] = {0};

	StEpdGeom * mEpdGeom = new StEpdGeom();
	for(Int_t iHit=0; iHit<nepdHits; iHit++) {

		StPicoEpdHit *epdHit = mPicoDst->epdHit(iHit);
		if( !epdHit ) continue;
		hEpdHitAdc->Fill(epdHit->nMIP());
		Short_t side_EW = epdHit->side(); //+1 for West and -1 for East
		//	if (side_EW == -1 && epdHit->nMIP()>0.3 && epdHit->nMIP()<3.) nepdMIPsE++;
		//	if (side_EW == 1 && epdHit->nMIP()>0.3 && epdHit->nMIP()<3.) nepdMIPsW++;
		//if (side_EW == -1 && epdHit->nMIP()>0.3) nepdMIPsE += (epdHit->nMIP()>5 ? epdHit->nMIP() : 1);
		//if (side_EW ==  1 && epdHit->nMIP()>0.3) nepdMIPsW += (epdHit->nMIP()>5 ? epdHit->nMIP() : 1);
		if( side_EW==-1 && epdHit->nMIP()>0.3 ) nepdMIPsE += (epdHit->nMIP()>5) ? 5 : epdHit->nMIP();
		if( side_EW== 1 && epdHit->nMIP()>0.3 ) nepdMIPsW += (epdHit->nMIP()>5) ? 5 : epdHit->nMIP();


		//TVector3 StraightLine = mEpdGeom->TileCenter(epdHit->id()) - pRcVx; // This give eta-phi plot with gap that I don't like
		TVector3 StraightLine = mEpdGeom->RandomPointOnTile(epdHit->id()) - picoEvent->primaryVertex(); //Becasue I like this more than TileCenter
		double phi_epd = StraightLine.Phi();
		double eta_epd = StraightLine.Eta();
		if(side_EW == -1) hEpdEtaPhiEast->Fill(eta_epd,phi_epd,epdHit->nMIP());
		if(side_EW ==  1) hEpdEtaPhiWest->Fill(eta_epd,phi_epd,epdHit->nMIP());

		float mip = epdHit->nMIP(); // gain calibrated energy loss in tile, in units of Landau MPV for one MIP
		int iring = epdHit->row() - 1; // (1-16)-1 -> 0-15
		if( mip>0.3 ){
			if( side_EW==-1 ) nepdRingE[iring] += (mip>5) ? 5 : mip;
			if( side_EW== 1 ) nepdRingW[iring] += (mip>5) ? 5 : mip;

			if( side_EW==-1 ) nepdRingEraw[iring] += epdHit->nMIP();
			if( side_EW== 1 ) nepdRingWraw[iring] += epdHit->nMIP();
		}

		//	cout<<" iHit= "<<iHit<<" id= "<<epdHit->id()<<" ADC= "<<epdHit->adc()<<" nMIP= "<<epdHit->nMIP()<<endl;
	} // end of iHit loop 
	// End of Epd block __________________________________________________


	// check data
	if( mEventCounter>1 ){ 
		if( refMult_pre==refMult && vz_pre==vz && zdcAdcSumE_pre==zdcAdcSumE ){
			cout << "skip this event... possible duplicate...  RM=" << refMult_pre <<" "<< refMult <<" vz="<< vz_pre <<" "<< vz << endl;
			return kStOK;
		} 
	}
	refMult_pre    = refMult;
	vz_pre         = vz;
	zdcAdcSumE_pre = zdcAdcSumE;


	Int_t trgWord = -1;
	//Int_t trgWord = picoEvent->triggerWord();
	// Run11 AuAu200, vpd-zdc-mb-protected
	// <high> hgfe dcba <low>
	// c=350003, d=350013, e=350023, f=350033, g=350043
	
	
	Int_t trgID = -1;
	Int_t trgMatch[ntrg] = {0};
	if( mData==1 ){ //run11
		if     ( trgWord & (Int_t)pow(2.,2.) ) trgID = 0;
		else if( trgWord & (Int_t)pow(2.,3.) ) trgID = 1;
		else if( trgWord & (Int_t)pow(2.,4.) ) trgID = 2;
		else if( trgWord & (Int_t)pow(2.,5.) ) trgID = 3;
		else if( trgWord & (Int_t)pow(2.,6.) ) trgID = 4;

	}else if( mData==2 ){ //run14
		// MB trgID 
		// 0=450050,  1=450060,  2=450005,  3=450015,  4=450025
		if( picoEvent->isTrigger(450050) ){ trgMatch[0] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(450060) ){ trgMatch[1] = 1;  trgID = 1; }
		if( picoEvent->isTrigger(450005) ){ trgMatch[2] = 1;  trgID = 2; }
		if( picoEvent->isTrigger(450015) ){ trgMatch[3] = 1;  trgID = 3; }
		if( picoEvent->isTrigger(450025) ){ trgMatch[4] = 1;  trgID = 4; }

	}else if( mData==3 ){ //run16
		// VPDMB-5-p-sst
		// 520001, 520011, 520021, 520031, 520041, 520051
		if( picoEvent->isTrigger(520001) ){ trgMatch[0] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(520011) ){ trgMatch[1] = 1;  trgID = 1; }
		if( picoEvent->isTrigger(520021) ){ trgMatch[2] = 1;  trgID = 2; }
		if( picoEvent->isTrigger(520031) ){ trgMatch[3] = 1;  trgID = 3; }
		if( picoEvent->isTrigger(520041) ){ trgMatch[4] = 1;  trgID = 4; }
		if( picoEvent->isTrigger(520051) ){ trgMatch[5] = 1;  trgID = 5; }

	}else if( mData==4 ){ // run18 
		// isobar
		// vpdmb-30     = 600001, 600011, 600021, 600031
		// vpdmb-30-hlt = 600002, 600012, 600022, 600032, 600042
		//
		//if( picoEvent->isTrigger(600001) ){ trgMatch[0] = 1;  trgID = 0; }
		//if( picoEvent->isTrigger(600011) ){ trgMatch[1] = 1;  trgID = 1; }
		//if( picoEvent->isTrigger(600021) ){ trgMatch[2] = 1;  trgID = 2; }
		//if( picoEvent->isTrigger(600031) ){ trgMatch[3] = 1;  trgID = 3; }
		//
		// AuAu 27GeV
		// mb= 610001, 610011, 610021, 610031, 610041, 610051, 
		if( picoEvent->isTrigger(610001) ){ trgMatch[0] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(610011) ){ trgMatch[1] = 1;  trgID = 1; }
		if( picoEvent->isTrigger(610021) ){ trgMatch[2] = 1;  trgID = 2; }
		if( picoEvent->isTrigger(610031) ){ trgMatch[3] = 1;  trgID = 3; }
		if( picoEvent->isTrigger(610041) ){ trgMatch[4] = 1;  trgID = 4; }
		if( picoEvent->isTrigger(610051) ){ trgMatch[5] = 1;  trgID = 5; }

	}else if( mData==5 ){ // run19 19GeV
		// AuAu 19GeV MB
		if( picoEvent->isTrigger(640001) ){ trgMatch[0] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(640011) ){ trgMatch[1] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(640021) ){ trgMatch[2] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(640031) ){ trgMatch[3] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(640041) ){ trgMatch[4] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(640051) ){ trgMatch[5] = 1;  trgID = 0; }

	}else if( mData==6 ){ // run19 14GeV
		if( picoEvent->isTrigger(650000) ){ trgMatch[0] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(650010) ){ trgMatch[1] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(650020) ){ trgMatch[2] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(650030) ){ trgMatch[3] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(650040) ){ trgMatch[4] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(650050) ){ trgMatch[5] = 1;  trgID = 0; }

	}else if( mData==7 ){ // run19 7GeV, minbias
		if( picoEvent->isTrigger(660000) ){ trgMatch[0] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(660010) ){ trgMatch[1] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(660020) ){ trgMatch[2] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(660030) ){ trgMatch[3] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(660040) ){ trgMatch[4] = 1;  trgID = 0; }

	}else if( mData==8 ){ // run19 9GeV, minbias
		if( picoEvent->isTrigger(670000) ){ trgMatch[0] = 1;  trgID = 0; }

	}else if( mData==9 || mData==10 ){ // run19 FXT

		// run20 5p75 FXT, epde-or-bbce-or-vpde-tof1-etof
		if( picoEvent->isTrigger(720007) ){ trgMatch[0] = 1;  trgID = 0; }

	}else if( mData==11 ){ // run19 AuAu200
		if( picoEvent->isTrigger(700001) ){ trgMatch[0] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(700011) ){ trgMatch[0] = 1;  trgID = 0; }

	}else if( mData==12 ){ // run20 AuAu11p5
		if( picoEvent->isTrigger(710000) ){ trgMatch[0] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(710010) ){ trgMatch[0] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(710020) ){ trgMatch[0] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(710030) ){ trgMatch[0] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(710040) ){ trgMatch[0] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(710050) ){ trgMatch[0] = 1;  trgID = 0; }
	}else if( mData==13 ){ // run20 9p2
		if( picoEvent->isTrigger(780010) ){ trgMatch[0] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(780020) ){ trgMatch[0] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(780030) ){ trgMatch[0] = 1;  trgID = 0; }
		if( picoEvent->isTrigger(780040) ){ trgMatch[0] = 1;  trgID = 0; }
	}

	// mData==9 for FXT, no trogger selection
	if( trgID<0 ) return kStOK;
	hTrg->Fill(trgID);

	bool IsETOF; 
	if( mFixedMode==0 ){ // minbias_withetof
		// 710018 for 11p2, 780008 and 780018 for 9p2
		IsETOF = (picoEvent->isTrigger(710018) || picoEvent->isTrigger(780008) || picoEvent->isTrigger(780018)); 
	}else IsETOF = true;

	// Fill event info
	hVz    ->Fill(vz);
	hVxVy  ->Fill(vx,vy);
	hdVz   ->Fill(vz-vpd_vz);
	hVzVpdz->Fill(vz,vpd_vz);
	hVzAll ->Fill(vz);
	hVzVr  ->Fill(vz,sqrt(vx*vx+vy*vy));
	hRunidvsAllEvt->Fill( picoEvent->runId() );
	hVtxRank[0]->Fill( picoEvent->ranking() );


	Float_t vx_ave = 0.0; // 0.056; // 0.0;
	Float_t vy_ave = 0.0; //-0.326; // 0.0;
	if( mFixedMode==1 ){
		vx_ave = -0.4;
		vy_ave = -2.0;
	}

	Float_t vxc = vx - vx_ave;
	Float_t vyc = vy - vy_ave;
	Float_t vr  = sqrt( vxc*vxc+vyc*vyc );
	if( vr<mVrAbsMax ){ 
		hVz_VrCut->Fill(vz);
		hVxVyCut ->Fill(vx,vy);
		hVzvsrefMult->Fill(vz,refMult);
		hVzvsEpdMipE->Fill(vz,nepdMIPsE);
		hVzvsEpdMipW->Fill(vz,nepdMIPsW);
		hprVzvsrefMult->Fill(vz,refMult);
		hprVzvsEpdMipE->Fill(vz,nepdMIPsE);
		hprVzvsEpdMipW->Fill(vz,nepdMIPsW);

		// collider
		if( mFixedMode==0 ){
			if( fabs(vz)<mVzAbsMax ){ 
				hRunidvsGoodEvt[0]->Fill( picoEvent->runId() );
				if( picoEvent->nBTOFMatch()>0 ){ 
					hRunidvsGoodEvt[1]->Fill( picoEvent->runId() );
					if( fabs(vz-vpd_vz)<mVpdVzDif ) hRunidvsGoodEvt[2]->Fill( picoEvent->runId() );
					if( IsETOF ) hRunidvsGoodEvt[3]->Fill( picoEvent->runId() );
				}
			}

		// fixed
		}else{
			if( vz>195 && vz<205 ){ 
				hRunidvsGoodEvt[0]->Fill( picoEvent->runId() );
				if( picoEvent->nBTOFMatch()>0 ){
					hRunidvsGoodEvt[1]->Fill( picoEvent->runId() );
					if( fabs(vz-vpd_vz)<mVpdVzDif ) hRunidvsGoodEvt[2]->Fill( picoEvent->runId() );
					if( IsETOF ) hRunidvsGoodEvt[3]->Fill( picoEvent->runId() );
				}
			}
		}
	}

	if( mFixedMode==0 ){
		if( fabs(vz)>mVzAbsMax )        return kStOK;
		if( picoEvent->nBTOFMatch()<1 ) return kStOK;
	}else{
		if( vz<195 || vz>205 )      return kStOK;
	}
	if( vr>mVrAbsMax )              return kStOK;
	if( mFixedMode==0 && fabs(vz-vpd_vz)>mVpdVzDif ) return kStOK;


	hrefMult ->Fill(refMult);
	hgrefMult->Fill(grefMult);
	htofMult ->Fill(tofMult);
	hrefMultvsTofMult->Fill(refMult,tofMult);
	hrefMultvsTofMatch->Fill(refMult,picoEvent->nBTOFMatch());
	hrefMultvsETofMult->Fill(refMult,picoEvent->etofHitMultiplicity());
	hVtxRank[1]->Fill( picoEvent->ranking() );

	// nTofT0 info.
	Int_t  nTofT0 = picoEvent->nTofT0();
	hnTofT0->Fill(nTofT0);



	// TofMult cut
	//bool fTofMatch = IsTofMatch( refMult, tofMult, 0 );
	//if( !fTofMatch ) return kStOK;


	// Set Centrality
	////pRefMultCorr->initEvent( refMult, vz, zdcRate );
	//pRefMultCorr->initEvent(grefMult,vz,zdcRate); 
	//grefMult = pRefMultCorr->getRefMultCorr(); // return grefMult
	//mTrgEff = pRefMultCorr->getWeight();
	mTrgEff = 1.0;


	// Centrality
	int fQzvtx = static_cast<int>( ( vz + mVzAbsMax ) / (mVzAbsMax*2.0) * NqvZvtx );
	//int fQcent = pRefMultCorr->getCentralityBin16(); // 0=75-80%, 15=0-5%
	//int fCent  = pRefMultCorr->getCentralityBin9();  // 0=70-80%, 8=0-5%
	// temporary setting
	int fQcent = 0; //pRefMultCorr->getCentralityBin16(); // 0=75-80%, 15=0-5%
	int fCent  = 0; //pRefMultCorr->getCentralityBin9();  // 0=70-80%, 8=0-5%
	if( fQcent<0 || fQcent>=NqvCent ) return kStOK;
	if( fCent<0 || fCent>=NCent )     return kStOK;
	if( mFixedMode==0 ){
		if( fQzvtx<0 || fQzvtx>=NqvZvtx ) return kStOK;
	}else fQzvtx = 0;

	// when using RefMultCor
	int fCent_tmp = NCent - fCent - 1;
	fCent = fCent_tmp;


	// check ZDC and BBC 
	hZDCEvsW->Fill( zdcAdcSumE, zdcAdcSumW );
	hZDCvsMult[0]->Fill( zdcAdcSumE, refMult );
	hZDCvsMult[1]->Fill( zdcAdcSumW, refMult );
	hZDCvsGMult[0]->Fill( zdcAdcSumE, picoEvent->grefMult() );
	hZDCvsGMult[1]->Fill( zdcAdcSumW, picoEvent->grefMult() );

	Double_t bbcAdcSumE = 0.0;
	Double_t bbcAdcSumW = 0.0;
	Double_t bbcAdcRingE[4] = {0}; // inner only
	Double_t bbcAdcRingW[4] = {0};
	for( int ipmt=0; ipmt<24; ipmt++ ){
		Double_t bbce = picoEvent->bbcAdcEast(ipmt);
		Double_t bbcw = picoEvent->bbcAdcWest(ipmt);
		hBBCAdc2D->Fill( ipmt, bbce );
		hBBCAdc2D->Fill( ipmt+24, bbcw );
		if( ipmt<16 ){
			bbcAdcSumE += bbce;
			bbcAdcSumW += bbcw;
		}
		hBBCAdc2DHi->Fill( ipmt, bbce );
		hBBCAdc2DHi->Fill( ipmt+24, bbcw );

		// ring-by-ring
		if( ipmt<6 ){ // most inner
			bbcAdcRingE[0] += bbce;
			bbcAdcRingW[0] += bbcw;
		}else if( ipmt<16 ){ // second inner
			bbcAdcRingE[1] += bbce;
			bbcAdcRingW[1] += bbcw;
		}else if( ipmt<20 ){
			bbcAdcRingE[2] += bbce;
			bbcAdcRingW[2] += bbcw;
		}else{
			bbcAdcRingE[3] += bbce;
			bbcAdcRingW[3] += bbcw;
		}
	}
	hBBCEvsW->Fill( bbcAdcSumE, bbcAdcSumW );
	hBBCvsMult [0]->Fill( bbcAdcSumE, refMult );
	hBBCvsMult [1]->Fill( bbcAdcSumW, refMult );
	hBBCvsGMult[0]->Fill( bbcAdcSumE, picoEvent->grefMult() );
	hBBCvsGMult[1]->Fill( bbcAdcSumW, picoEvent->grefMult() );
	hBBCvsTMult[0]->Fill( bbcAdcSumE, tofMult );
	hBBCvsTMult[1]->Fill( bbcAdcSumW, tofMult );

	for( int iring=0; iring<4; iring++ ){
		hrefMultvsBBCRingE[iring]->Fill( refMult, bbcAdcRingE[iring] );
		hrefMultvsBBCRingW[iring]->Fill( refMult, bbcAdcRingW[iring] );
	}

	for( int ist=0; ist<8; ist++ ){
		hZDCAdc2D->Fill( ist,    picoEvent->ZdcSmdEastHorizontal(ist) );
		hZDCAdc2D->Fill( ist+8,  picoEvent->ZdcSmdEastVertical  (ist) );
		hZDCAdc2D->Fill( ist+16, picoEvent->ZdcSmdWestHorizontal(ist) );
		hZDCAdc2D->Fill( ist+24, picoEvent->ZdcSmdWestVertical  (ist) );

		//cout << picoEvent->ZdcSmdEastHorizontal(ist) <<" ";
	}
	//cout << endl;

	// rate
	hZdcRate[0]->Fill( picoEvent->zdcEastRate()/1000. );
	hBbcRate[0]->Fill( picoEvent->bbcEastRate()/1000. );
	hZdcRate[1]->Fill( picoEvent->zdcWestRate()/1000. );
	hBbcRate[1]->Fill( picoEvent->bbcWestRate()/1000. );
	hZdcRate[2]->Fill( zdcRate/1000. );
	hBbcRate[2]->Fill( bbcRate/1000. );
	hBGRate    ->Fill( picoEvent->backgroundRate()/1000. );



	for( int iring=0; iring<16; iring++ ){
		hrefMultvsEpdRingE[iring]->Fill( refMult, nepdRingE[iring] );
		hrefMultvsEpdRingW[iring]->Fill( refMult, nepdRingW[iring] );
		hrefMultvsEpdRingEraw[iring]->Fill( refMult, nepdRingEraw[iring] );
		hrefMultvsEpdRingWraw[iring]->Fill( refMult, nepdRingWraw[iring] );
	}
	hEpdEvsW->Fill( nepdMIPsE, nepdMIPsW );


	// Mtd hit info.
	// ==================================================
	unsigned int nHitsMtd = mPicoDst->numberOfMtdHits(); 
	for( unsigned int l=0; l<nHitsMtd; l++ ){
		StPicoMtdHit* hit = mPicoDst->mtdHit(l);
		if( !hit ) continue;
		
		int backleg = hit->backleg(); // 1-30
		int module  = hit->module();  // 1-5
		int cell    = hit->cell();    // 0-11
		hMtdHit->Fill( backleg-1, cell+(module-1)*12 );
		//cout <<"bl="<< backleg <<" mod="<< module <<" cell="<< cell << endl;
	}
	
	// Bemc hit info.
	// ==================================================
	//unsigned int nHitsBTow = mPicoDst->numberOfBTowHits(); 
	//StEmcGeom *mBTowGeom = new StEmcGeom(1); // bemc=1 from StEmcGeom.cxx
	//for( unsigned int l=0; l<nHitsBTow; l++ ){
	//	StPicoBTowHit* hit = mPicoDst->btowHit(l);
	//	if( !hit ) continue;
	//	
	//	int softId = hit->numericIndex2SoftId(l);
	//	int adc    = hit->adc();
	//	float energy = hit->energy();

	//	float btow_eta, btow_phi;
	//	mBTowGeom->getEtaPhi( softId, btow_eta, btow_phi );
	//	hBTowHit->Fill( btow_phi, btow_eta );
	//	//cout <<"bl="<< backleg <<" mod="<< module <<" cell="<< cell << endl;
	//}
	


	// Event plane reconstruction
	// ==================================================
	// ZDC-SMD event plane
	// ------------------------------------------------------
	//QV Q_SMD[NsubSMD];
	//CalcQvSMD( picoEvent, Q_SMD, fQcent, fQzvtx );

	// BBC event plane
	// ------------------------------------------------------
	//QV Q_BBC[NordBBC][NsubBBC];
	//CalcQvBBC( picoEvent, Q_BBC, fQcent, fQzvtx );

	// EPD event plane
	// ------------------------------------------------------
	QV Q_EPD[NordEPD][NsubEPD];
	CalcQvEPD( mPicoDst, Q_EPD, fQcent, fQzvtx );

	if( mCalibMode==0 ) return kStOK;

	// TPC event plane
	// ------------------------------------------------------
	//QV Q_TPC[NordTPC][NsubTPC];
	//CalcQvTPC( mPicoDst, Q_TPC, fQcent, fQzvtx );
	// ------------------------------------------------------

	// EPD 
	//StEpdEpInfo result = mEpFinder->Results(pico->picoArray(StPicoArrays::EpdHit),pRcVx,myCentral);

	// BBC Q-vector
	//for( int ith=0; ith<2; ith++ ){
	//	for( int iep=0; iep<NsubBBC; iep++ ){
	//		if( Q_BBC[ith][iep].get_Psi()>-9000 ){
	//			//float fq = sqrt( pow(Q_BBC[ith][iep].get_Qx(), 2.0) + pow(Q_BBC[ith][iep].get_Qy(), 2.0) );
	//			hQBBC[ith][iep][fCent]->Fill( Q_BBC[ith][iep].get_Qabs() );
	//		}
	//}	}

	// TPC Q-vector
	//for( int ith=0; ith<3; ith++ ){
	//	for( int iep=0; iep<NsubTPC; iep++ ){
	//		if( Q_TPC[ith][iep].get_Psi()>-9000 ){
	//			float fq = sqrt( pow(Q_TPC[ith][iep].get_Qx(), 2.0) + pow(Q_TPC[ith][iep].get_Qy(), 2.0) );
	//			hQTPC[ith][iep][fCent]->Fill( fq );
	//		}
	//}	}



	// EP correlation
	// ===========================================
	//for( int ith=0; ith<NordSMD; ith++ ){
	//	if( Q_SMD[0].get_Psi()>-9000 && Q_SMD[1].get_Psi()>-9000 ){
	//		float dpsiSMD = ( ith + 1.0 ) * ( Q_SMD[0].get_Psi() - Q_SMD[1].get_Psi() );
	//		hSMDEPCrrvsMult[ith][0]->Fill( fCent, cos(dpsiSMD), mTrgEff );
	//		hSMDEPCrrvsMult[ith][1]->Fill( fCent, sin(dpsiSMD), mTrgEff );
	//	}
	//}

	//for( int ith=0; ith<NordBBC; ith++ ){
	//	if( Q_BBC[ith][0].get_Psi()>-9000 && Q_BBC[ith][1].get_Psi()>-9000 ){
	//		float ndpsi = ( ith + 1.0 ) * ( Q_BBC[ith][0].get_Psi() - Q_BBC[ith][1].get_Psi() );
	//		hBBCEPCrrvsMult[ith][0]->Fill( fCent, cos( ndpsi ), mTrgEff );
	//		hBBCEPCrrvsMult[ith][1]->Fill( fCent, sin( ndpsi ), mTrgEff );
	//	}

	//	if( Q_BBC[ith][3].get_Psi()>-9000 && Q_BBC[ith][4].get_Psi()>-9000 ){
	//		float ndpsi = ( ith + 1.0 ) * ( Q_BBC[ith][3].get_Psi() - Q_BBC[ith][4].get_Psi() );
	//		hBBCEPCrrvsMult[ith][2]->Fill( fCent, cos( ndpsi ), mTrgEff );
	//		hBBCEPCrrvsMult[ith][3]->Fill( fCent, sin( ndpsi ), mTrgEff );
	//	}
	//}

	// ===========================================


	TVector3 pvtx_3d = picoEvent->primaryVertex();
	StThreeVectorF pvtx( pvtx_3d.X(), pvtx_3d.Y(), pvtx_3d.Z() );
	//float ZDCEP[3], BBCEP[3], TPCEP[3];
	//for( int iep=0; iep<3; iep++ ){ 
	//	ZDCEP[iep] = Q_SMD[iep].get_Psi();
	//	BBCEP[iep] = Q_BBC[1][iep+3].get_Psi();
	//	TPCEP[iep] = Q_TPC[1][iep+3].get_Psi(); // |eta|>0.3
	//	TPCEP[iep] = Q_TPC[1][iep+9].get_Psi(); // |eta|>0.5
	//}
	for(int ich=0; ich<2; ich++){
		pEvt[ich].set_pVtx(pvtx);
		pEvt[ich].set_BField(picoEvent->bField());
		pEvt[ich].set_Cent(fCent);
		pEvt[ich].set_Pid(2);
		//pEvt[ich].set_Psi( ZDCEP );
		//pEvt[ich].set_Psi2B( BBCEP );
		//pEvt[ich].set_Psi2 ( TPCEP );
		//pEvt[ich].set_Q1B( Q_BBC[0][5].get_Qabs() );
		pEvt[ich].set_Trgeff(mTrgEff);

		piEvt[ich].set_pVtx(pvtx);
		piEvt[ich].set_BField(picoEvent->bField());
		piEvt[ich].set_Cent(fCent);
		piEvt[ich].set_Pid(0);
		//piEvt[ich].set_Psi( ZDCEP );
		//piEvt[ich].set_Psi2B( BBCEP );
		//piEvt[ich].set_Psi2 ( TPCEP );
		//piEvt[ich].set_Q1B( Q_BBC[0][5].get_Qabs() );
		piEvt[ich].set_Trgeff(mTrgEff);
	}

	int Nch[2] = {0};
	int Np [2] = {0};

	float fqvCh[4][2][2] = {{{0}}};
	float fqwCh   [2][2] = {{0}};

	float fqvFXT[4][2][2] = {{{0}}};
	float fqwFXT   [2][2] = {{0}};

	unsigned int nprtrks = 0;
	unsigned int ngltrks = 0;

	// track loop
	for( unsigned int i=0; i<mPicoDst->numberOfTracks(); i++ ){

		StPicoTrack *ptrk = (StPicoTrack*)mPicoDst->track(i);
		if( !ptrk ) continue;  

		// track info.
		float dca = ptrk->gDCA( pvtx_3d ).Mag();
		int nHitsFit = ptrk->nHitsFit();
		short ch  = ptrk->charge();
		short iCH = (ch>0) ? 0 : 1;
		hNhitsFit->Fill( nHitsFit );
		hNhitsMax->Fill( ptrk->nHitsMax() );
		if( nHitsFit<15 )  continue; 
		if( ptrk->nHitsFit()<0.52*ptrk->nHitsMax() ) continue;

		// global tracks
		float gpt  = ptrk->gMom().Perp(); 
		float gphi = ptrk->gMom().Phi();  
		float geta = ptrk->gMom().PseudoRapidity();
		float gmom = ptrk->gMom().Mag();

		if( gpt>0.15 ) hgDca->Fill( dca );
		if( fabs(dca)>10 ) continue;

		// primary track info.
		float pt  = ptrk->pMom().Perp(); 
		float phi = ptrk->pMom().Phi();  
		float eta = ptrk->pMom().PseudoRapidity();
		float mom = ptrk->pMom().Mag();
		if( pt>0.15 ){ 
			hDca->Fill( dca );
			hDcaXY->Fill( ptrk->gDCAxy(pvtx_3d.X(), pvtx_3d.Y()) );
			hDcaZ ->Fill( ptrk->gDCAz(pvtx_3d.Z()) );
		}


		// good primary
		bool isGoodPrimary = false;
		if( pt>0.15 && pt<10 && fabs(dca)<3.0 ){ 
			//hPDca          ->Fill( dca );
			//hPNhitsFit     ->Fill( ptrk->nHitsFit() );
			hPPhivsEta[iCH]->Fill( phi, eta );

			if( mFixedMode==1 ){
				hPPt[iCH]->Fill( pt );
				hPdEdx   ->Fill( ptrk->dEdx() );
				hPdEdxvsP->Fill( ptrk->pMom().Mag()*ch, ptrk->dEdx() );
				hPt      ->Fill(pt);
				hEta     ->Fill(eta);
				hPhi     ->Fill(phi);
				hNhits   ->Fill(nHitsFit);
				hDedx    ->Fill(ptrk->dEdx());
				hDcaTmp  ->Fill(dca);
				Nch[iCH]++;
				nprtrks++;
				isGoodPrimary = true;
			}else if( mFixedMode==0 && fabs(eta)<1.8 ){
				hPPt[iCH]->Fill( pt );
				hPdEdx   ->Fill( ptrk->dEdx() );
				hPdEdxvsP->Fill( ptrk->pMom().Mag()*ch, ptrk->dEdx() );
				hPt      ->Fill(pt);
				hEta     ->Fill(eta);
				hPhi     ->Fill(phi);
				hNhits   ->Fill(nHitsFit);
				hDedx    ->Fill(ptrk->dEdx());
				hDcaTmp  ->Fill(dca);
				Nch[iCH]++;
				nprtrks++;
				isGoodPrimary = true;
			}

			// only for test
			hAveDedx2D [iCH]->Fill( phi, eta, ptrk->dEdx() );
			hAveDedxPhi[iCH]->Fill( phi, ptrk->dEdx() );
			hDedxPhi   [iCH]->Fill( phi, ptrk->dEdx() );
		}
		if( isGoodPrimary && ptrk->isMtdTrack() ) hPPhivsEtaMtd->Fill( phi, eta );
		if( isGoodPrimary && ptrk->isBemcTrack() ) hPPhivsEtaBemc->Fill( phi, eta );

		// good global
		if( gpt>0.15 && gpt<10 ){ 
			//hGDca->Fill( dca );
			//hGNhitsFit     ->Fill( ptrk->nHitsFit() );
			hGPhivsEta[iCH]->Fill( gphi, geta );
			if( fabs(dca)<3 ) hGPhivsEta2[iCH]->Fill( gphi, geta );
			if( mFixedMode==1 ){
				hGPt[iCH]      ->Fill( gpt );
				hGdEdx         ->Fill( ptrk->dEdx() );
				hGdEdxvsP      ->Fill( ptrk->gMom().Mag()*ch, ptrk->dEdx() );
				ngltrks++;
			}else if( mFixedMode==0 && fabs(geta)<1.8 ){
				hGPt[iCH]      ->Fill( gpt );
				hGdEdx         ->Fill( ptrk->dEdx() );
				hGdEdxvsP      ->Fill( ptrk->gMom().Mag()*ch, ptrk->dEdx() );
				ngltrks++;
			}
		}


		// PID
		//===========================================
		// TPC dE/dx
		float nsPi = ptrk->nSigmaPion();
		float nsP  = ptrk->nSigmaProton();
		float nsK  = ptrk->nSigmaKaon();
		bool isPionTpc   = ( fabs(nsPi)<2 && fabs(nsK)>2 && fabs(nsP)>2 ) ? true : false;
		//bool isKaonTpc   = ( fabs(nsPi)>2 && fabs(nsK)<2 && fabs(nsP)>2 ) ? true : false;
		bool isProtonTpc = ( fabs(nsPi)>2 && fabs(nsK)>2 && fabs(nsP)<2 ) ? true : false;
		if( isProtonTpc && pt<0.4 ) isProtonTpc = false;

		// TOF m2
		float m2   = -9999;
		float beta = getTofBeta( ptrk, picoEvent->primaryVertex(), picoEvent->bField() ); //ptrk->btofBeta();
		if( nTofT0<3 ) beta = -1;
		if( beta>0 ) m2 = mom*mom*( 1.0 / beta / beta - 1.0 ); 

		if( isGoodPrimary ){
			hM2vsP->Fill( mom*ch, m2 );			
			//if( dca<3 ) hM2vsPdca->Fill( mom*ch, m2 );			
			hNsigmaP[0]->Fill(mom,nsPi);
			hNsigmaP[1]->Fill(mom,nsK);
			hNsigmaP[2]->Fill(mom,nsP);
			hNsigmaP[3]->Fill(mom,ptrk->nSigmaElectron());
			if( beta>0 && pt>0.15 ) hBetaIvsP->Fill( mom*ch, 1.0/beta );
		}
		if( beta>0 && pt>0.15 && pt<10 ) hPPhivsEtaTof[iCH]->Fill( phi, eta );

		// dca vs p for global tracks
		if( fabs(nsPi)<2 ){
			if( beta>0 ) hDcavsP[0]->Fill( gmom*ch, dca ); 
			else         hDcavsP[1]->Fill( gmom*ch, dca ); 
		}
		if( fabs(nsK)<2 ){
			if( beta>0 ) hDcavsP[2]->Fill( gmom*ch, dca ); 
			else         hDcavsP[3]->Fill( gmom*ch, dca ); 
		}
		if( fabs(nsP)<2 ){
			if( beta>0 ) hDcavsP[4]->Fill( gmom*ch, dca ); 
			else         hDcavsP[5]->Fill( gmom*ch, dca ); 
		}

		// count (anti)protons
		if( isGoodPrimary && isProtonTpc && pt<2 ) Np[iCH]++;


		// ETOF
		int index2etof = ptrk->eTofPidTraitsIndex();
		float betaETOF = -1.0;
		if( index2etof>=0 ){
			StPicoETofPidTraits const* const etofPid = mPicoDstMaker->picoDst()->etofPidTraits(index2etof);
			if( etofPid ){
				betaETOF = etofPid->beta();
				float crsX = etofPid->crossingX();
				float crsY = etofPid->crossingY();
				if( betaETOF>0 && pt>0.15 ){ 
					hBetaIvsPe->Fill( mom*ch, 1.0/betaETOF );
					hETOFxy->Fill( crsX, crsY );
				}
				//cout <<"beta="<< betaETOF <<" x="<< crsX <<" y="<< crsY << endl;
			}
		}

		// Set buffer
		if( isPionTpc || isProtonTpc ){

			SgTrack v0trk(ptrk,picoEvent,mPicoDst);
			//if( i_pid==0 && fabs(dca)>(mDcaPi[(15-fQcent)/2][0]-0.1) ) piEvt[iCH].AddTrk( v0trk );
			//if( i_pid==1 && fabs(dca)>(mDcaP [(15-fQcent)/2][0]-0.1) ) pEvt [iCH].AddTrk( v0trk );

			// rough DCA cuts for Lambdas
			if( isPionTpc   && fabs(dca)>mDcaPi[(15-fQcent)/2] ) piEvt[iCH].AddTrk( v0trk );
			if( isProtonTpc && fabs(dca)>mDcaP [(15-fQcent)/2] ) pEvt [iCH].AddTrk( v0trk );
		}


		// Calculate TPC Q-vectors
		//________________________________________________________
		// track cuts first
		if( !ptrk->isPrimary() ) continue;
		if( pt<0.15 || pt>2.0 )  continue;
		if( fabs(dca)>3.0 )      continue;
		if( mFixedMode==0 && fabs(eta)>1.5 )     continue;
		if( mFixedMode==1 && (eta<-2 || eta>0) ) continue;

		// raw flow vectors, charge separately 
		for( int ih=0; ih<4; ih++ ){
			fqvCh[ih][iCH][0] += cos( (ih+1.0)*phi );
			fqvCh[ih][iCH][1] += sin( (ih+1.0)*phi );
			if( ih==0 ){
				fqwCh[iCH][0] += 1;
				fqwCh[iCH][1] += pt;
			}
		}

		if( mFixedMode==1 ){
			int iTrk = (int)(gRandom->Rndm()*2.0);
			if( iTrk==0 || iTrk==1 ){
				for( int ih=0; ih<4; ih++ ){
					//fqvFXT[ih][iTrk][0] += cos( (ih+1.0)*phi );
					//fqvFXT[ih][iTrk][1] += sin( (ih+1.0)*phi );

					float wgt = 1.0;
					if( ih==0 ) wgt = pt*(eta+2.1);
					else        wgt = pt;
					fqvFXT[ih][iTrk][0] += wgt*cos( (ih+1.0)*phi );
					fqvFXT[ih][iTrk][1] += wgt*sin( (ih+1.0)*phi );

					if( ih==0 ){
						fqwFXT[iTrk][0] += 1;
						fqwFXT[iTrk][1] += pt;
					}
				}
			}
		}

	}// track loop

	// skip pair calculation 
	//if( mCalibMode<3 ) return kStOK;
	  
	// make V0 particles in real events
	//MakeV0Pair( pEvt[0], piEvt[1], 0, REAL ); // p(2)    & pim(0) for lambda 
	//MakeV0Pair( pEvt[1], piEvt[0], 1, REAL ); // pbar(2) & pip(0) for anti-lambda 

	if( ngltrks!=0 ) hRatioPrGlTrk->Fill( (float)nprtrks/(float)ngltrks );


	// Run dependence
	// ________________________________________________________________________
	// event info.
	hRunidvsrefMult  ->Fill( picoEvent->runId(), refMult );
	hRunidvsgrefMult ->Fill( picoEvent->runId(), grefMult );
	hRunidvstofMult  ->Fill( picoEvent->runId(), tofMult );
	hRunidvstofMatch ->Fill( picoEvent->runId(), picoEvent->nBTOFMatch() );
	hRunidvsBemcMatch->Fill( picoEvent->runId(), picoEvent->nBEMCMatch() );
	hRunidvsEpdMipE  ->Fill( picoEvent->runId(), nepdMIPsE );
	hRunidvsEpdMipW  ->Fill( picoEvent->runId(), nepdMIPsW );
	if( IsETOF || picoEvent->etofHitMultiplicity()>0 ){
		hRunidvsEtofMult ->Fill( picoEvent->runId(), picoEvent->etofHitMultiplicity() );
	}

	hRunidvsVx  ->Fill( picoEvent->runId(), vx );
	hRunidvsVy  ->Fill( picoEvent->runId(), vy );
	hRunidvsVz  ->Fill( picoEvent->runId(), vz );
	hRunidvsRank->Fill( picoEvent->runId(), picoEvent->ranking() );

	hRunidvsZdcx  ->Fill( picoEvent->runId(), picoEvent->ZDCx() );
	hRunidvsBbcx  ->Fill( picoEvent->runId(), picoEvent->BBCx() );
	hRunidvsBgrate->Fill( picoEvent->runId(), picoEvent->backgroundRate() );

	for( int ih=0; ih<2; ih++ ){
		if( (fqwCh[0][0]+fqwCh[1][0])>0 ){
			hRunidvsTpcQx[ih]->Fill( picoEvent->runId(), (fqvCh[ih][0][0]+fqvCh[ih][1][0])/(fqwCh[0][0]+fqwCh[1][0]) );
			hRunidvsTpcQy[ih]->Fill( picoEvent->runId(), (fqvCh[ih][0][1]+fqvCh[ih][1][1])/(fqwCh[0][0]+fqwCh[1][0]) );
		}
	}
	for( int ih=0; ih<2; ih++ ){
		for( int isub=0; isub<2; isub++ ){
			if( Q_EPD[ih][isub].get_Psi()>-100 ){
				hRunidvsEpdQx[ih][isub]->Fill( picoEvent->runId(), Q_EPD[ih][isub].get_Qx() );
				hRunidvsEpdQy[ih][isub]->Fill( picoEvent->runId(), Q_EPD[ih][isub].get_Qy() );
			}
		}
	}

	// primary track info.
	hRunidvsPt        ->Fill( picoEvent->runId(), hPt    ->GetMean() );
	hRunidvsEta       ->Fill( picoEvent->runId(), hEta   ->GetMean() );
	hRunidvsPhi       ->Fill( picoEvent->runId(), hPhi   ->GetMean() );
	hRunidvsNhits     ->Fill( picoEvent->runId(), hNhits ->GetMean() );
	hRunidvsDedx      ->Fill( picoEvent->runId(), hDedx  ->GetMean() );
	hRunidvsDca       ->Fill( picoEvent->runId(), hDcaTmp->GetMean() );
	hRunidvsNchp      ->Fill( picoEvent->runId(), Nch[0] );
	hRunidvsNchm      ->Fill( picoEvent->runId(), Nch[1] );
	hRunidvsNproton   ->Fill( picoEvent->runId(), Np[0] );
	hRunidvsNprotonbar->Fill( picoEvent->runId(), Np[1] );

	if( hEpdEvsW->GetEntries()>0 ){
		for( int ix=0; ix<hEta->GetNbinsX(); ix++ ){
			hRunidvsEtaDist->Fill( picoEvent->runId(), hEta->GetBinCenter(ix+1), hEta->GetBinContent(ix+1) );
		}
	}


	// test only for 31GeV (FXT7GeV)
	// ==========================================
	/*
	int centFXT;
	if     ( refMult<= 9 ) centFXT = 10;
	else if( refMult<=17 ) centFXT =  9;
	else if( refMult<=29 ) centFXT =  8;
	else if( refMult<=43 ) centFXT =  7;
	else if( refMult<=63 ) centFXT =  6;
	else if( refMult<=87 ) centFXT =  5;
	else if( refMult<=114) centFXT =  4;
	else if( refMult<=148) centFXT =  3;
	else if( refMult<=167) centFXT =  2;
	else if( refMult<=188) centFXT =  1;
	else                   centFXT =  0;

	if( centFXT<10 ){
		for( int ih=0; ih<2; ih++ ){
			for( int iep=0; iep<2; iep++ ){
				hRunvsQvPar[ih][iep][centFXT][0]->Fill(picoEvent->runId(),fqvFXT[ih][iep][0]);
				hRunvsQvPar[ih][iep][centFXT][1]->Fill(picoEvent->runId(),fqvFXT[ih][iep][1]);
			}
		}

		// recentering
		for( int ih=0; ih<2; ih++ ){
			for( int iep=0; iep<2; iep++ ){
				for( int ixy=0; ixy<2; ixy++ ){
					int ibin = hInRunvsQvPar[ih][iep][centFXT][ixy]->FindBin(picoEvent->runId());
					float mn = hInRunvsQvPar[ih][iep][centFXT][ixy]->GetBinContent(ibin);
					float sg = hInRunvsQvPar[ih][iep][centFXT][ixy]->GetBinError  (ibin);
					if( sg==0 ) continue;
					fqvFXT[ih][iep][ixy] -= mn;
					fqvFXT[ih][iep][ixy] /= sg;
				}
			}
		}

		float psiFXT[2][2];
		for( int ih=0; ih<2; ih++ ){
			for( int iep=0; iep<2; iep++ ){
				if( fqwFXT[iep][0]>0 ) psiFXT[ih][iep] = atan2( fqvFXT[ih][iep][1], fqvFXT[ih][iep][0] ) / (ih+1.0);
				else                   psiFXT[ih][iep] = -9999;
			}
		}

		for( int ih=0; ih<2; ih++ ){
			for( int iep=0; iep<2; iep++ ){
				if( psiFXT[ih][iep]>-100 ){
					hQvCor[ih][iep  ]->Fill( centFXT, cos( (ih+1.0)*(psiFXT[ih][0]-psiFXT[ih][1]) ) );
					hQvCor[ih][iep+2]->Fill( centFXT, sin( (ih+1.0)*(psiFXT[ih][0]-psiFXT[ih][1]) ) );
				}
			}
		}
	}
	// ==========================================
	*/



	// reset buffers
	// ________________________________________________________________________
	for( int ich=0; ich<2; ich++ ){
		pEvt [ich].Clear();
		piEvt[ich].Clear();
	}
	hPt    ->Reset();
	hEta   ->Reset();
	hPhi   ->Reset();
	hNhits ->Reset();
	hDedx  ->Reset();
	hDcaTmp->Reset();

	return kStOK;
}

//---------------------------------------------------------------------------------------
Float_t StBES2QaMaker::ZDCSMD( StPicoEvent *pEv, int eastwest, int verthori, int strip ) const {

	float val = 0;
	//val = mEv->zdcTriggerDetector().zdcSmd( (StBeamDirection)eastwest, verthori, strip );
	if     ( eastwest==0 && verthori==0 ) val = pEv->ZdcSmdEastVertical  (strip-1);
	else if( eastwest==0 && verthori==1 ) val = pEv->ZdcSmdEastHorizontal(strip-1);
	else if( eastwest==1 && verthori==0 ) val = pEv->ZdcSmdWestVertical  (strip-1);
	else if( eastwest==1 && verthori==1 ) val = pEv->ZdcSmdWestHorizontal(strip-1);
	return val;
}

Float_t StBES2QaMaker::ZDCSMD_GetPosition( int eastwest, int verthori, int strip ) {
// Get position of each slat;strip starts from 1

	Float_t zdcsmd_x[7] = {0.5,2,3.5,5,6.5,8,9.5};
	Float_t zdcsmd_y[8] = {1.25,3.25,5.25,7.25,9.25,11.25,13.25,15.25};

	if(eastwest==0 && verthori==0) return zdcsmd_x[strip-1]-mZDCSMDCenterex;
	if(eastwest==1 && verthori==0) return mZDCSMDCenterwx-zdcsmd_x[strip-1];
	if(eastwest==0 && verthori==1) return zdcsmd_y[strip-1]/sqrt(2.)-mZDCSMDCenterey;
	if(eastwest==1 && verthori==1) return zdcsmd_y[strip-1]/sqrt(2.)-mZDCSMDCenterwy;

  return 0;
}

bool StBES2QaMaker::IsTofMatch( double& refM, unsigned short& tofM, int runnum ){

	bool fmatch = true;

	// Run10 Au+Au 200
	if( mData==0 ){

		// TofMult cut parameters 
		double parUp[5][2];
		// run10 RF
		parUp[0][0] = 140;  parUp[0][1] = 11;
		parUp[1][0] = 130;  parUp[1][1] = 9.5;
		parUp[2][0] = 130;  parUp[2][1] = 7.5;
		// run10 FF
		parUp[3][0] = 130;  parUp[3][1] = 9.5;
		parUp[4][0] = 130;  parUp[4][1] = 9.5;

		double parLw[5][2];
		// run10 RF
		parLw[0][0] = -70;  parLw[0][1] = 4.5;
		parLw[1][0] = -75;  parLw[1][1] = 4.5;
		parLw[2][0] = -75;  parLw[2][1] = 3.0;
		// run10 FF
		parLw[3][0] = -150;  parLw[3][1] = 3.5;
		parLw[4][0] = -180;  parLw[4][1] = 4.2;

		int id = -1;
		if     ( runnum>=11044006 && runnum<=11054067 ) id = 0; // RF-1
		else if( runnum>=11055095 && runnum<=11064043 ) id = 1; // RF-2
		else if( runnum>=11064057 && runnum<=11077018 ) id = 2; // RF-3

		if     ( runnum>=11002122 && runnum<=11028088 ) id = 3; // FF trg=260001
		else if( runnum>=11028094 && runnum<=11035026 ) id = 4; // FF trg=260011

		if( id>=0 ){ // mark bad events in RF runs
			if( tofM < (parLw[id][0]+parLw[id][1]*refM) ) fmatch = false;
			if( tofM > (parUp[id][0]+parUp[id][1]*refM) ) fmatch = false;
		}	

	// Run11 Au+Au 200
	}else if( mData==1 ){

		// test
		//if( tofM < (-160+3.1*refM) ) fmatch = false;
		//if( tofM > ( 130+5.5*refM) ) fmatch = false;

		// adjust
		if( tofM < (-160+3.5*refM) ) fmatch = false;
		if( tofM > ( 140+5.0*refM) ) fmatch = false;

	// Run14 Au+Au 200
	}else if( mData==2 ){

		// rough cut
		if( mCalibMode==0 ){
			if( tofM < (-160+3.0*refM) ) fmatch = false;
			if( tofM > ( 140+10.*refM) ) fmatch = false;

		// adjusted cut
		}else{
			if( tofM < (-240+4.5*refM) ) fmatch = false;
			if( tofM > ( 140+10.*refM) ) fmatch = false;
		}

	// Run16 Au+Au 200, grefMult cut
	}else if( mData==3 ){

		// rough cut
		if( mCalibMode==0 ){
			if( tofM < (-220+3.0*refM) ) fmatch = false;
			if( tofM > ( 180+7.0*refM) ) fmatch = false;
		}else{
			if( tofM < (-200+3.5*refM) ) fmatch = false;
			if( tofM > ( 180+5.8*refM) ) fmatch = false;
		}

	// Run18 isobar
	}else if( mData==4 ){
	}

	return fmatch;
}

void StBES2QaMaker::CalcQvSMD( StPicoEvent *pEv, QV fQv[NsubSMD], int fCent, int fVz ){

	if( fCent<0 || fCent>=NqvCent || fVz<0 || fVz>=NqvZvtx ) return;
	

	// Calculate raw flow vectors
	// -----------------------------------------------
	float fqv[NsubSMD][4] = {{0}}; //east-west, xy+weight
	for( int iep=0; iep<2; iep++ ){ // east-west

		int nstrip;
		for( int ixy=0; ixy<2; ixy++ ){
			if( ixy==0 ) nstrip = 8; // vertical strips (7 in x-direction)
			else         nstrip = 9; // horizontal strips (8 in y-direction)
			for( int is=1; is<nstrip; is++ ){
				Float_t zdc_adc = ZDCSMD( pEv, iep, ixy, is );
				hZDCAdc[iep]->Fill( ixy*9+is, zdc_adc );
				cout << "Zdc_Adc = " << zdc_adc << endl;

				// pedestal subtraction for Au+Au 2010
				//int id = ixy*7 + is;
				//if( ixy==1 ) id += 1;
				if( mData==0 ){ // only for Run10 and/or older runs
					Float_t ped = hZDCSMDped->GetBinContent( mZDCSMDmap[iep][ixy][is-1]+1, mZDCPEDRun-11001000 );
					zdc_adc -= ped;
					if( zdc_adc<0 ) zdc_adc = 0;
				}

				fqv[iep][ixy]   += ZDCSMD_GetPosition( iep, ixy, is ) * zdc_adc;
				fqv[iep][ixy+2] += zdc_adc;
				hZDCAdcCor[iep]->Fill( ixy*9+is, zdc_adc );
			}

			if( fqv[iep][ixy+2]>0 ) fqv[iep][ixy] /= fqv[iep][ixy+2];
			else                    fqv[iep][ixy] = -9999;

			// only for blind data
			//if( iep==0 && pEv->ZdcSumAdcEast()<1 ) fqv[iep][ixy] = -9999;
			//if( iep==1 && pEv->ZdcSumAdcWest()<1 ) fqv[iep][ixy] = -9999;
		}
		hQSMD2D[iep][0]->Fill( fqv[iep][0], fqv[iep][1] );
	}
	// combined
	for( int iq=0; iq<2; iq++ ){ 
		fqv[2][iq] = fqv[0][iq] - fqv[1][iq];
		if( fqv[0][iq]<-9999 || fqv[1][iq]<-9999 ) fqv[2][iq] = -9999;
	}
	hQSMD2D[2][0]->Fill( fqv[2][0], fqv[2][1] );


	// Set raw Q-vectors and Calibration //
	// =======================================================
	int fth   = 0;
	float fQw = 1.0;
	for( int isub=0; isub<NsubSMD; isub++ ){
		if( fqv[isub][0]==-9999 || fqv[isub][1]==-9999 ) fQw = -1;
		fQv[isub].set_Qv( fqv[isub][0], fqv[isub][1], fQw, dZDC, isub, fth, fCent, fVz );
		pQVMaker->setProfiles( hSMDrecX[isub], hSMDrecY[isub], hSMDfltC[isub], hSMDfltS[isub], hSMDEP[isub][fCent] );
		pQVMaker->doCalibration( &(fQv[isub]) );
		hQSMD2D[isub][1]->Fill( fQv[isub].get_Qx(), fQv[isub].get_Qy() );
	}
	// =======================================================
}

void StBES2QaMaker::CalcQvTPC( StPicoDst *pD, QV fQv[NordTPC][NsubTPC], int fCent, int fVz ){

	StPicoEvent *picoEvent = (StPicoEvent*)pD->event();

	float fqvCh[4][4][2] = {{{0}}};
	float fqwCh[4][2] = {{0}};
	float fqvTPC[NordTPC][NsubTPC][4] = {{{0}}};
	for( unsigned int l=0; l<pD->numberOfTracks(); l++ ){

    	StPicoTrack *ptrk =  (StPicoTrack*)pD->track(l); 
    	if( !ptrk ) continue;
		if( !isGoodTrackEP(ptrk) ) continue;

	    //StPhysicalHelix helix = ptrk->dcaGeometry().helix();
	    //float dca = helix.geometricSignedDistance(picoEvent->primaryVertex());

		// run16
		//float dca = ptrk->dca().mag();
		
		float dca = ptrk->gDCA( picoEvent->primaryVertex() ).Mag();
		if( fabs(dca)>3.0 ) continue;

		float pt, eta, phi;
		if( mTrkEP==0 ){
			pt  = ptrk->pMom().Perp();
			eta = ptrk->pMom().PseudoRapidity();
			phi = ptrk->pMom().Phi();
		}else{
			pt  = ptrk->gMom().Perp();
			eta = ptrk->gMom().PseudoRapidity();
			phi = ptrk->gMom().Phi();
		}
		short ch  = ptrk->charge();

		int i_EP01 = -1;
		if     ( eta<-0.1 ) i_EP01 = 0;
		else if( eta> 0.1 ) i_EP01 = 1;

		int i_EP03 = -1;
		if     ( eta<-0.3 ) i_EP03 = 0;
		else if( eta> 0.3 ) i_EP03 = 1;

		int i_EP04 = -1;
		if     ( eta<-0.4 ) i_EP04 = 0;
		else if( eta> 0.4 ) i_EP04 = 1;

		int i_EP05 = -1;
		if     ( eta<-0.5 ) i_EP05 = 0;
		else if( eta> 0.5 ) i_EP05 = 1;

		int i_EP02 = -1;
		if( fabs(eta)<0.2 ) i_EP02 = 0;


		int iCH = 0;
		if( ch<0 ) iCH = 1;

		for( int ith=0; ith<NordTPC; ith++ ){

			int sgn = 1.0;
			//if( ith==0 && i_EP==1 ) sgn = -1.0;


			// eta<-0.1 or eta>0.1 or combined
			if( i_EP01>=0 ){
				float wgt  = pt;
				float wgtC = pt;

				fqvTPC[ith][i_EP01][0] += wgt * cos( (ith+1.0)*phi );
				fqvTPC[ith][i_EP01][1] += wgt * sin( (ith+1.0)*phi );
				fqvTPC[ith][i_EP01][2] += wgt;
				if( wgt!=0 ) fqvTPC[ith][i_EP01][3] += 1;

				fqvTPC[ith][2][0] += sgn * wgtC * cos( (ith+1.0)*phi );
				fqvTPC[ith][2][1] += sgn * wgtC * sin( (ith+1.0)*phi );
				fqvTPC[ith][2][2] += wgtC;
				if( wgtC!=0 ) fqvTPC[ith][2][3] += 1;
			}

			// eta<-0.3 or eta>0.3 or combined
			if( i_EP03>=0 ){
				float wgt  = pt;
				float wgtC = pt;

				fqvTPC[ith][i_EP03+3][0] += wgt * cos( (ith+1.0)*phi );
				fqvTPC[ith][i_EP03+3][1] += wgt * sin( (ith+1.0)*phi );
				fqvTPC[ith][i_EP03+3][2] += wgt;
				if( wgt!=0 ) fqvTPC[ith][i_EP03+3][3] += 1;

				fqvTPC[ith][5][0] += sgn * wgtC * cos( (ith+1.0)*phi );
				fqvTPC[ith][5][1] += sgn * wgtC * sin( (ith+1.0)*phi );
				fqvTPC[ith][5][2] += wgtC;
				if( wgtC!=0 ) fqvTPC[ith][5][3] += 1;
			}

			// |eta|>0.4 
			if( i_EP04>=0 ){
				float wgt = pt;
				float wgtC = pt;

				fqvTPC[ith][i_EP04+6][0] += wgt * cos( (ith+1.0)*phi );
				fqvTPC[ith][i_EP04+6][1] += wgt * sin( (ith+1.0)*phi );
				fqvTPC[ith][i_EP04+6][2] += wgt;
				if( wgt!=0 ) fqvTPC[ith][i_EP04+6][3] += 1;

				fqvTPC[ith][8][0] += sgn * wgtC * cos( (ith+1.0)*phi );
				fqvTPC[ith][8][1] += sgn * wgtC * sin( (ith+1.0)*phi );
				fqvTPC[ith][8][2] += wgtC;
				if( wgtC!=0 ) fqvTPC[ith][8][3] += 1;
			}

			// |eta|>0.5
			if( i_EP05>=0 ){
				float wgt  = pt;
				float wgtC = pt;

				fqvTPC[ith][i_EP05+9][0] += wgt * cos( (ith+1.0)*phi );
				fqvTPC[ith][i_EP05+9][1] += wgt * sin( (ith+1.0)*phi );
				fqvTPC[ith][i_EP05+9][2] += wgt;
				if( wgt!=0 ) fqvTPC[ith][i_EP05+9][3] += 1;

				fqvTPC[ith][11][0] += sgn * wgtC * cos( (ith+1.0)*phi );
				fqvTPC[ith][11][1] += sgn * wgtC * sin( (ith+1.0)*phi );
				fqvTPC[ith][11][2] += wgtC;
				if( wgtC!=0 ) fqvTPC[ith][11][3] += 1;
			}

			// h+ or h-
			//if( i_EP03>=0 ){
			//	float wgt  = pt;
			//	fqvTPC[ith][i_EP03+12+iCH*2][0] += wgt * cos( (ith+1.0)*phi );
			//	fqvTPC[ith][i_EP03+12+iCH*2][1] += wgt * sin( (ith+1.0)*phi );
			//	fqvTPC[ith][i_EP03+12+iCH*2][2] += wgt;
			//	if( wgt!=0 ) fqvTPC[ith][i_EP03+12+iCH*2][3] += 1;
			//}

			// |eta|<0.2
			if( i_EP02==0 ){
				float wgt = pt;
				//fqvTPC[ith][16][0] += wgt * cos( (ith+1.0)*phi );
				//fqvTPC[ith][16][1] += wgt * sin( (ith+1.0)*phi );
				//fqvTPC[ith][16][2] += wgt;
				//if( wgt!=0 ) fqvTPC[ith][16][3] += 1;
				fqvTPC[ith][12][0] += wgt * cos( (ith+1.0)*phi );
				fqvTPC[ith][12][1] += wgt * sin( (ith+1.0)*phi );
				fqvTPC[ith][12][2] += wgt;
				if( wgt!=0 ) fqvTPC[ith][12][3] += 1;
			}

			// raw flow vectors
			// eta<-0.1 or eta>0.1 charge separately 
			if( i_EP01>=0 ){
				fqvCh[ith][i_EP01+iCH*2][0] += cos( (ith+1.0)*phi );
				fqvCh[ith][i_EP01+iCH*2][1] += sin( (ith+1.0)*phi );
				if( ith==1 ){
					fqwCh[i_EP01+iCH*2][0] += 1;
					fqwCh[i_EP01+iCH*2][1] += pt;
				}
			}//i_EP01

		}// ith

	}// track loop

	for( int ith=0; ith<4; ith++ ){
		for( int isub=0; isub<4; isub++ ){
			for( int ixy=0; ixy<2; ixy++ ){
				mQvCh[ith][isub][ixy] = fqvCh[ith][isub][ixy];
			}
		}
	}
	for( int isub=0; isub<4; isub++ ){
		for( int ixy=0; ixy<2; ixy++ ){
			mQwCh[isub][ixy] = fqwCh[isub][ixy];
		}
	}


	// for n=1, weight = mult
	for( int isub=0; isub<NsubTPC; isub++ ) fqvTPC[0][isub][2] = fqvTPC[0][isub][3]; 

	// Set raw Q-vectors and Calibration //
	// =======================================================
	for( int ith=0; ith<NordTPC; ith++ ){
		for( int isub=0; isub<NsubTPC; isub++ ){

			fQv[ith][isub].set_Qv( fqvTPC[ith][isub][0], fqvTPC[ith][isub][1], fqvTPC[ith][isub][2], fqvTPC[ith][isub][3], dTPC, isub, ith, fCent, fVz );
			pQVMaker->setProfiles( hTPCrecX[ith][isub], hTPCrecY[ith][isub], hTPCfltC[ith][isub], hTPCfltS[ith][isub], hTPCEP[ith][isub][fCent] );
			pQVMaker->doCalibration( &(fQv[ith][isub]) );
		}
	}
	// =======================================================
}

void StBES2QaMaker::CalcQvBBC( StPicoEvent* pEv, QV fQv[NordBBC][NsubBBC], int fCent, int fVz ){

	if( fCent<0 || fCent>=NqvCent || fVz<0 || fVz>=NqvZvtx ) return;

	// run14             0     1     2     3     4     5     6     7     8     9     10    11    12    13    14    15
	float adcMax[48] = { 4030, 4031, 4018, 4022, 4019, 4043, 4018, 4015, 4016, 4013, 4017, 4030, 4007, 4046, 4030, 4042, 
		4200, 4200, 4200, 4200, 4200, 4200, 4200, 4200,
   						 4022, 4022, 4027, 4013, 4012, 4022, 4027, 4019, 4021, 4013, 4013, 4014, 4023, 4015, 4021, 4013,
		4200, 4200, 4200, 4200, 4200, 4200, 4200, 4200 };
	int Saturate[2] = {0};

	float fqv[NordBBC][NsubBBC][3] = {{{0.0}}}; //nth, EP, (qx,qy,qw)
	for( int iep=0; iep<2; iep++ ){  // iep=0->east, iep=1->west

		// Random selected
		//-------------------------------------------------
		float adc0[24] = {0};
		for( int ipmt=0; ipmt<24; ipmt++ ){

			//float adc = (float)mEv->bbcTriggerDetector().adc(ipmt+iep*24);
			//float tdc = (float)mEv->bbcTriggerDetector().tdc(ipmt+iep*24);
			float adc; 
			if( iep==0 ) adc = pEv->bbcAdcEast(ipmt);
			else         adc = pEv->bbcAdcWest(ipmt);
			//if( tdc<10 || adc<1 || tdc>3500 ) continue; 
			if( ipmt>15 ) continue;
			if( adc<1 )   continue; 
			if( adc>adcMax[ipmt+iep*24] ) Saturate[iep] = 1;
			hBbcAdc->Fill( ipmt+iep*24, adc );

			if( BBCgc ) adc *= BBCgc->GetBinContent(ipmt+1+24*iep);
			adc0[ipmt] = adc;
			hBbcAdcCor->Fill( ipmt+iep*24, adc );

			float phi = GetBBCPhi( iep, ipmt );
			for( int ith=0; ith<NordBBC; ith++ ){
				// east or west
				fqv[ith][iep][0] += adc * cos( (ith+1.0)*phi );
				fqv[ith][iep][1] += adc * sin( (ith+1.0)*phi );
				fqv[ith][iep][2] += adc;
				// both	
				float sgn = 1.0;
				if( ith==0 && iep==1 ) sgn = -1.0; 
				fqv[ith][2][0] += sgn * adc * cos( (ith+1.0)*phi );
				fqv[ith][2][1] += sgn * adc * sin( (ith+1.0)*phi );
				fqv[ith][2][2] += adc;

			}//ith
		}//ipmt

		// Weighted average
		//-------------------------------------------------
		float AdcTile[18] = {0};
		SetBBCAdc( adc0, AdcTile, iep );

		for( int itl=0; itl<18; itl++ ){

			float phi = GetBBCTilePhi( iep, itl );
			for( int ith=0; ith<NordBBC; ith++ ){
				// east or west
				fqv[ith][iep+3][0] += AdcTile[itl] * cos( (ith+1.0)*phi );
				fqv[ith][iep+3][1] += AdcTile[itl] * sin( (ith+1.0)*phi );
				fqv[ith][iep+3][2] += AdcTile[itl];
				// both
				float sgn = 1.0;
				if( ith==0 && iep==1 ) sgn = -1.0;
				fqv[ith][5][0] += sgn * AdcTile[itl] * cos( (ith+1.0)*phi );
				fqv[ith][5][1] += sgn * AdcTile[itl] * sin( (ith+1.0)*phi );
				fqv[ith][5][2] += AdcTile[itl];
			}//ith
			//hBbcTileAdc[fCent]->Fill( itl+iep*18, AdcTile[itl] );

		}//itile	

	}//iep

	for( int iep=0; iep<2; iep++ ){  // iep=0->east, iep=1->west
		if( Saturate[iep] ){
			for( int ith=0; ith<NordBBC; ith++ ){
				fqv[ith][iep  ][2] = -1;
				fqv[ith][iep+3][2] = -1;
				fqv[ith][  2  ][2] = -1;
				fqv[ith][  5  ][2] = -1;
			}
		}
	}

	// Set raw Q-vectors and Calibration //
	// =======================================================
	for( int ith=0; ith<NordBBC; ith++ ){
		for( int isub=0; isub<NsubBBC; isub++ ){

			fQv[ith][isub].set_Qv( fqv[ith][isub][0], fqv[ith][isub][1], fqv[ith][isub][2], dBBC, isub, ith, fCent, fVz );
			pQVMaker->setProfiles( hBBCrecX[ith][isub], hBBCrecY[ith][isub], hBBCfltC[ith][isub], hBBCfltS[ith][isub], hBBCEP[ith][isub][fCent] );
			pQVMaker->doCalibration( &(fQv[ith][isub]) );
		}
	}
	// =======================================================
}

void StBES2QaMaker::SetBBCAdc( float adcP[24], float adcT[18], int iep ){

		for( int ipmt=0; ipmt<6; ipmt++ ) adcT[ipmt] = adcP[ipmt];

		for( int ipmt=6; ipmt<16; ipmt++ ){ 
			if( ipmt==6 || ipmt==11 ) continue;
			int itl = ipmt;
			if( ipmt>=8  ) itl += 1;
			if( ipmt>=13 ) itl += 1;
			adcT[itl] = adcP[ipmt];
		}

		float adc1 = 0; 
		float adc2 = 0;
		float sum1 = adcP[7] + adcP[8]; 
		float sum2 = adcP[7] + adcP[15];
		if( (sum1+sum2)!=0 ){
			adc1 = sum1 * adcP[6] / ( sum1 + sum2 ) * BBCgc->GetBinContent(50+iep); // for tile-9
			adc2 = sum2 * adcP[6] / ( sum1 + sum2 ) * BBCgc->GetBinContent(50+iep); // for tile-7
		}else{
			adc1 = (gRandom->Rndm()>0.5) ? 0.5 * adcP[6] * BBCgc->GetBinContent(50+iep) : 0;
			adc2 = 0.0;
			if( adc1==0 ) adc2 = 0.5 * adcP[6] * BBCgc->GetBinContent(50+iep);
		}
		adcT[6] = adc1;
		adcT[8] = adc2;
		
		sum1 = adcP[12] + adcP[10]; 
		sum2 = adcP[12] + adcP[13]; 
		if( (sum1+sum2)!=0 ){
			adc1 = sum1 * adcP[11] / ( sum1 + sum2 ) * BBCgc->GetBinContent(50+iep); // for tile-13
			adc2 = sum2 * adcP[11] / ( sum1 + sum2 ) * BBCgc->GetBinContent(50+iep); // for tile-15
		}else{
			adc1 = (gRandom->Rndm()>0.5) ? 0.5 * adcP[11] * BBCgc->GetBinContent(50+iep) : 0;
			adc2 = 0.0;
			if( adc1==0 ) adc2 = 0.5 * adcP[11] * BBCgc->GetBinContent(50+iep);
		}
		adcT[12] = adc1;
		adcT[14] = adc2;
}

Int_t StBES2QaMaker::LoadBBCalibParam( int runnum ){

	FileStat_t buf;
	int check = gSystem->GetPathInfo(Form("%s/%d.root",gSystem->Getenv("BBCALIB"),runnum),buf);
	if( check==1 ){ cout << "LoadBBCalibParam(): File does not exist : " << Form("%s/%d.root",gSystem->Getenv("BBCALIB"),runnum) << endl; }

	TFile *in; 
	if( mCalibMode==0 ){
		BBCgc = new TH1F("hBbcGainCrr","hBbcGainCrr",52,0,52);
		for( int ix=0; ix<BBCgc->GetNbinsX(); ix++ ) BBCgc->SetBinContent( ix+1, 1.0 );
		cout << "LoadBBCalibParam(): Gain correction factors were set to 1. continue..." << endl;
	}else{
		if( check==0 ) in = TFile::Open( Form("%s/%d.root",gSystem->Getenv("BBCALIB"),runnum) );
		else           in = NULL;
		if( !in ){ 
			cout << "Could not open a calibration file for BBC gain correction !! Gain correction factors were set to 1." << endl;
		
			BBCgc = new TH1F("hBbcGainCrr","hBbcGainCrr",52,0,52);
			for( int ix=0; ix<BBCgc->GetNbinsX(); ix++ ) BBCgc->SetBinContent( ix+1, 1.0 );
		}else { 
			
			cout << endl << Form("BBC gain correction -----> %s/BBCgain_%d.root was loaded",gSystem->Getenv("BBCALIB"),runnum) << endl; 

			TH1F *hist = (TH1F*)in->Get("hBbcGainCrr");
			if( hist==NULL || hist->IsZombie() ){ 
				cout <<"Could not get correction histogram!"<< endl; 
				BBCgc = new TH1F("hBbcGainCrr","hBbcGainCrr",52,0,52);
				for( int ix=0; ix<BBCgc->GetNbinsX(); ix++ ) BBCgc->SetBinContent( ix+1, 1.0 );
				cout << "LoadBBCalibParam(): Gain correction factors were set to 1." << endl;
				return -1; 
			}
			BBCgc = (TH1F*)hist->Clone();
			BBCgc->SetDirectory(0);
			in->Close();
		}
	}

	return 0;
}

Float_t StBES2QaMaker::GetBBCPhi(const Int_t eastWest, const Int_t tileId) const
{
  //float GetPhiInBBC(int eastWest, int bbcN) { //tileId=0 to 23
  const float Pi = TMath::Pi() ;
  const float phi_div=Pi/6;
  float bbc_phi=phi_div;
  switch(tileId) {
  case 0: bbc_phi=3*phi_div;
          break;
  case 1: bbc_phi=phi_div;
          break;
  case 2: bbc_phi=-1*phi_div;
          break;
  case 3: bbc_phi=-3*phi_div;
          break;
  case 4: bbc_phi=-5*phi_div;
          break;
  case 5: bbc_phi=5*phi_div;
          break;
  case 6: bbc_phi= (gRandom->Rndm()>0.5) ? 2*phi_div:4*phi_div;
          break;
  case 7: bbc_phi=3*phi_div;
          break;
  case 8: bbc_phi=phi_div;
          break;
  case 9: bbc_phi=0.;
          break;
  case 10: bbc_phi=-phi_div;
           break;
  case 11: bbc_phi=(gRandom->Rndm()>0.5) ? -2*phi_div:-4*phi_div;
           break;
  case 12: bbc_phi=-3*phi_div;
           break;
  case 13: bbc_phi=-5*phi_div;
           break;
  case 14: bbc_phi=Pi;
           break;
  case 15: bbc_phi=5*phi_div;
           break;
  case 16: bbc_phi=3*phi_div;
           break;
  case 17: bbc_phi=0.;
           break;
  case 18: bbc_phi=-3*phi_div;
           break;
  case 19: bbc_phi= Pi;
           break;
  case 20: bbc_phi=3*phi_div;
           break;
  case 21: bbc_phi=0.;
           break;
  case 22: bbc_phi=-3*phi_div;
           break;
  case 23: bbc_phi= Pi;
           break;
  }

  if(eastWest==0) {
    if(bbc_phi > -0.001) {bbc_phi = Pi-bbc_phi;}
    else {bbc_phi = -Pi-bbc_phi;}
  }
  if (bbc_phi < 0.) { bbc_phi += 2.0*Pi;}

  return bbc_phi;
}

float StBES2QaMaker::GetBBCTilePhi(const Int_t eastWest, const Int_t tileId){

	if( tileId<0 || tileId>17 ){ cout<<"GetBBCTilePhi(): Invalid tile-ID !!!"<<endl; return -9999; }
	return BbcTilePhi[eastWest][tileId];
}

void StBES2QaMaker::SetBBCTilesPhi(){

	const float Pi = TMath::Pi() ;
	
	// Most inner tiles
	BbcTilePhi[0][0] =  3.0*Pi/6.0;
	BbcTilePhi[0][1] =  1.0*Pi/6.0;
	BbcTilePhi[0][2] = -1.0*Pi/6.0;
	BbcTilePhi[0][3] = -3.0*Pi/6.0;
	BbcTilePhi[0][4] = -5.0*Pi/6.0;
	BbcTilePhi[0][5] =  5.0*Pi/6.0;
	
	// Middle inner tiles
	BbcTilePhi[0][6]  = 4.0*Pi/6.0; 
	BbcTilePhi[0][7]  = 3.0*Pi/6.0; 
	BbcTilePhi[0][8]  = 2.0*Pi/6.0; 
	BbcTilePhi[0][9]  = 1.0*Pi/6.0; 
	BbcTilePhi[0][10] = 0.0; 
	BbcTilePhi[0][11] = -1.0*Pi/6.0; 
	BbcTilePhi[0][12] = -2.0*Pi/6.0; 
	BbcTilePhi[0][13] = -3.0*Pi/6.0; 
	BbcTilePhi[0][14] = -4.0*Pi/6.0; 
	BbcTilePhi[0][15] = -5.0*Pi/6.0; 
	BbcTilePhi[0][16] = Pi; 
	BbcTilePhi[0][17] = 5.0*Pi/6.0; 

	for( int itl=0; itl<18; itl++ ) BbcTilePhi[1][itl] = BbcTilePhi[0][itl];

	for( int itl=0; itl<18; itl++ ){ 
		if( BbcTilePhi[0][itl]>-0.001 ) BbcTilePhi[0][itl] =  Pi - BbcTilePhi[0][itl];
		else                            BbcTilePhi[0][itl] = -Pi - BbcTilePhi[0][itl];

		if( BbcTilePhi[0][itl]< 0 ) BbcTilePhi[0][itl] += 2.0*Pi;
		if( BbcTilePhi[1][itl]< 0 ) BbcTilePhi[1][itl] += 2.0*Pi;
	}
}

bool StBES2QaMaker::isGoodPTrack( StPicoTrack *trk ){

	bool flg = true;
	float pt = trk->pMom().Perp();
	if( pt<0.15 || pt>10.0 )                     flg = false;
	if( trk->nHitsFit()<mNHitsFitMax )           flg = false;
	if( fabs(trk->pMom().PseudoRapidity())>1.0 ) flg = false;
	if( trk->nHitsFit()<0.52*trk->nHitsMax() )   flg = false;
	//if( fabs(trk->dca().mag())>mDcaMax )         flg = false; // applied in the original loop
	
	return flg;
}

bool StBES2QaMaker::isGoodGTrack( StPicoTrack *trk ){
	// No DCA cut

	bool flg = true;
	float pt = trk->gMom().Perp();
	if( pt<0.15 || pt>10.0 )                     flg = false;
	if( trk->nHitsFit()<mNHitsFitMax )           flg = false;
	if( fabs(trk->gMom().PseudoRapidity())>1.0 ) flg = false;
	if( trk->nHitsFit()<0.52*trk->nHitsMax() )   flg = false;
	
	return flg;
}

bool StBES2QaMaker::isGoodTrackEP( StPicoTrack *trk ){

	bool flg = true;
	float pt; 
	if( mTrkEP==0 ) pt = trk->pMom().Perp();
	else            pt = trk->gMom().Perp();
	if( pt<0.15 || pt>2.0 )                      flg = false;
	if( trk->nHitsFit()<15 )                     flg = false;
	if( fabs(trk->pMom().PseudoRapidity())>1.0 ) flg = false;
	if( trk->nHitsFit()<0.52*trk->nHitsMax() )   flg = false;
	//if( trk->dca()>3.0 )                         flg = false;
	
	return flg;
}

//Int_t StBES2QaMaker::LoadEPCalibParam( int runnum ){
//
//	TFile *in = TFile::Open( Form("%s/%d.root",gSystem->Getenv("EPCALIB"),runnum) );
//	if( !in ){ cout << "Could not open a calibration file!!  run=" << runnum << endl; return -1; }
//	else { cout << Form("EP calibration -----> %s/OutMuAna_%d.root was loaded",gSystem->Getenv("EPCALIB"),runnum) << endl << endl; }
//	
//	pQVMaker->LoadParameters(in);
//
//	in->Close();
//	return 0;
//}

void StBES2QaMaker::CalibParamInit(){

	// Initialization
	for( int ith=0; ith<NordBBC; ith++ ){
		for( int ict=0; ict<NqvCent; ict++ ){
			for( int iz=0; iz<NqvZvtx; iz++ ){
				for( int iep=0; iep<NsubBBC; iep++ ){
					BBCXYmn[ith][ict][iz][iep][0] = 0.0;
					BBCXYmn[ith][ict][iz][iep][1] = 0.0;
					BBCXYsg[ith][ict][iz][iep][0] = 1.0;
					BBCXYsg[ith][ict][iz][iep][1] = 1.0;

					for( int ik=0; ik<Nflt; ik++ ){
						BBCfltC[ith][ict][iz][iep][ik] = 0.0;
						BBCfltS[ith][ict][iz][iep][ik] = 0.0;
					}
	}	}	}	}

	for( int ict=0; ict<NqvCent; ict++ ){
		for( int iz=0; iz<NqvZvtx; iz++ ){
			for( int iep=0; iep<NsubSMD; iep++ ){
				for( int ixy=0; ixy<2; ixy++ ){
					SMDXYmn[ict][iz][iep][ixy] = 0.0;
					SMDXYsg[ict][iz][iep][ixy] = 1.0;
				}
				for( int ik=0; ik<Nflt; ik++ ){
					SMDfltC[ict][iz][iep][ik] = 0.0;
					SMDfltS[ict][iz][iep][ik] = 0.0;
				}
	}	}	}	

	for( int ith=0; ith<NordTPC; ith++ ){
		for( int ict=0; ict<NqvCent; ict++ ){
			for( int iz=0; iz<NqvZvtx; iz++ ){
				for( int iep=0; iep<NsubTPC; iep++ ){
					TPCXYmn[ith][ict][iz][iep][0] = 0.0;
					TPCXYmn[ith][ict][iz][iep][1] = 0.0;
					TPCXYsg[ith][ict][iz][iep][0] = 1.0;
					TPCXYsg[ith][ict][iz][iep][1] = 1.0;
					for( int ik=0; ik<Nflt; ik++ ){
						TPCfltC[ith][ict][iz][iep][ik] = 0.0;
						TPCfltS[ith][ict][iz][iep][ik] = 0.0;
					}
	}	}	}	}

}

//int StBES2QaMaker::GetZvtxBin( float zvtx ){
//
//	//int fZvtx = static_cast<int>( ( zvtx + 30.0 ) / 60.0 * NqvZvtx );
//	int fZvtx = static_cast<int>( ( zvtx + 6.0 ) / 12.0 * NqvZvtx );
//	return fZvtx;
//}

void StBES2QaMaker::SetCalibrationMode( int mode ){ 
	mCalibMode = mode; 
	cout << "StBES2QaMaker: CalibrationMode = " << mode << endl;
}

void StBES2QaMaker::SetDCAMax( float fDcaMax ){
	mDcaMax = fDcaMax;
	cout << "StBES2QaMaker: Set DCA cut parameter to be  dca<" << mDcaMax << endl;
}

void StBES2QaMaker::SetNHitsFitMax( int fNHitsFitMax ){
	mNHitsFitMax = fNHitsFitMax;
	cout << "StBES2QaMaker: Set nHitsFit cut parameter to be  nHitsFit>=" << mNHitsFitMax << endl;
}

void StBES2QaMaker::SetZVertexMax( float fVzMax ){
	mVzAbsMax = fVzMax;
	cout << "StBES2QaMaker: Set z-vertex cut parameter to be  vz<=" << mVzAbsMax << endl;
}

void StBES2QaMaker::SetTrgEffCorrection( int fTrgEff ){
	mTrgEffCorrection = fTrgEff;
	cout << "StBES2QaMaker: Set TrgEffCorrection = " << mTrgEffCorrection << endl;
}

void StBES2QaMaker::SetData( int flg ){ 
	mData = flg; 
	cout << "StBES2QaMaker: Data set = " << flg << endl;
}

void StBES2QaMaker::SetTopologicalCut( int fCut ){
	mTopologicalCut = fCut;
	cout << "StBES2QaMaker: Set mTopologicalCut = " << fCut << endl;
}

//void StBES2QaMaker::SetV0TofCorrection( int fCor ){
//	mV0TofCorrection = fCor;
//	cout << "StBES2QaMaker: Set mV0TofCorrection = " << fCor << endl;
//}

//void StBES2QaMaker::SetEvMixVtxShift( int flg ){
//	mEvMixVtxShift = flg;
//	cout << "StBES2QaMaker: Set mEvMixVtxShift = " << flg << endl;
//}

bool StBES2QaMaker::isProtonTof(float fm2){

	bool isOK = true;
	if( !(fm2>0.5 && fm2<1.5) ) isOK = false;
	return isOK;
}

bool StBES2QaMaker::isPionTof(float fm2, float fp){

	bool isOK = true;
	if( !(fm2>(0.017-0.013*fp) && fm2<0.04) ) isOK = false;
	return isOK;
}

//bool StBES2QaMaker::TofCorrM2Cut(StThreeVectorF& pvtx, StThreeVectorF& v0vtx, StLorentzVectorD& v0p, SgTrack& trk, StPhysicalHelixD& helix, int ipart){
//
//	Float_t tof = trk.tof(); 
//	Float_t betaCor = -1.0;
//	pTofCorr->setVectors3D(pvtx)(v0vtx)(trk.tofPos());
//	pTofCorr->setMotherTracks(v0p);
//	pTofCorr->correctBeta( helix, tof, betaCor );
//
//	float mom = trk.mom();
//	float m2 = mom * mom * ( 1.0 / betaCor / betaCor - 1.0 );
//
//	bool isOK = true;
//	// pion
//	if( ipart==0 && !isPionTof(m2,mom) ) isOK = false;
//	// proton
//	if( ipart==1 && !isProtonTof(m2) ) isOK = false; 
//
//	pTofCorr->clearContainers();
//	return isOK;
//}


// From the following code:
// https://github.com/MustafaMustafa/auau200GeVRun14Ana/blob/master/StRoot/StPicoD0AnaMaker/StPicoD0AnaMaker.cxx
// modified for run16
//float StBES2QaMaker::getTofBeta(StPicoTrack const* const trk, StThreeVectorF const& vtx, float BField) const {
float StBES2QaMaker::getTofBeta(StPicoTrack const* const trk, TVector3 const& vtx, float BField) const {

	int index2tof = trk->bTofPidTraitsIndex();
	float beta = -1.0; // std::numeric_limits<float>::quiet_NaN();

	if (index2tof >= 0)
	{
		StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

		if (tofPid)
		{
			beta = tofPid->btofBeta();

			if (beta < 1e-4)
			{
				TVector3 tofPos = tofPid->btofHitPos();
				StThreeVectorF const btofHitPos( tofPos.X(), tofPos.Y(), tofPos.Z() );
				StThreeVectorF vtxSt( vtx.X(), vtx.Y(), vtx.Z() );
				StThreeVectorF gmom( trk->gMom().X(), trk->gMom().Y(), trk->gMom().Z() );
				StThreeVectorF org( trk->origin().X(), trk->origin().Y(), trk->origin().Z() );

				//StPhysicalHelixD helix = trk->helix(BField);
				StPhysicalHelixD helix( gmom, org, BField*kilogauss, static_cast<float>(trk->charge()) );
				float L   = tofPathLength(&vtxSt, &btofHitPos, helix.curvature());
				float tof = tofPid->btof();
				if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
				else         beta = -1.0; //std::numeric_limits<float>::quiet_NaN();
			}
		}
	}

	return beta;
}


void StBES2QaMaker::MakeV0Pair(SgPEvent& eveA, SgPEvent& eveB, int fCH, PAIRTYPE fPair ){

	Trks trkA = eveA.get_trks();
	Trks trkB = eveB.get_trks();
	if( trkA.size()<1 || trkB.size()<1 ) return;

	
	int iCent = eveA.get_Cent() - 1;
	if( iCent==-1 ) iCent = 0; // for 0-5%
	if( iCent<0 ) return;

	double BField = eveA.get_BField();
	StThreeVectorF vtx_shift(0.0,0.0,0.0);

	StThreeVectorF pvtx = eveA.get_pVtx();
	if( fPair==MIX ){ 
		if( mEvMixVtxShift ) vtx_shift = eveA.get_pVtx() - eveB.get_pVtx();
	}

	for( unsigned int itrkA=0; itrkA<trkA.size(); itrkA++ ){

		StPhysicalHelixD helixA0 = trkA[itrkA].helix();
		StThreeVector<double> newOrigin = helixA0.origin() - vtx_shift;

		// StPhysicalHelix.hh
    	// Requires: momentum, origin, signed Magnetic Field and Charge of particle (+/- 1)
    	//StPhysicalHelix(const StThreeVector<double>&, const StThreeVector<double>&, double, double);
		StPhysicalHelixD helixA( trkA[itrkA].momentum(), newOrigin, eveA.get_BField()*kilogauss, (Double_t)trkA[itrkA].charge() );
		//if( fPair!=MIX ) cout << "mom_org=" << helixA0.momentum(BField*kilogauss).mag() <<" mom_new="<< helixA.momentum(BField*kilogauss).mag() << endl;
		
		for( unsigned int itrkB=0; itrkB<trkB.size(); itrkB++ ){

			StPhysicalHelixD helixB = trkB[itrkB].helix();

			// daughter particle information
			pair<double,double> s = helixA.pathLengths(helixB);
			StThreeVectorF pA = helixA.momentumAt( s.first, BField*kilogauss );   // momentum of daughter A (proton for lambda)
			StThreeVectorF pB = helixB.momentumAt( s.second, BField*kilogauss );  // momentum of daughter B
			float eA = sqrt( pA.mag()*pA.mag() + Mass[eveA.get_Pid()]*Mass[eveA.get_Pid()] ); // energy of A
			float eB = sqrt( pB.mag()*pB.mag() + Mass[eveB.get_Pid()]*Mass[eveB.get_Pid()] ); // energy of B

			
			// info. on topological cuts  
			StThreeVectorF pV0 = pA + pB; // momentum of parent
			float eV0 = eA + eB;          // energy of parent
			float invM = sqrt( eV0*eV0 - pV0.mag()*pV0.mag() ); // invariant mass
			StThreeVectorF dcaA = helixA.at(s.first);
			StThreeVectorF dcaB = helixB.at(s.second);
			StThreeVectorF v0  = (dcaA+dcaB)*0.5;
			if( invM>1.17 ) continue; // rough Minv cut
			hLambdaMass[fCH]->Fill( invM );


			float fdcaA = trkA[itrkA].dca();
			float fdcaB = trkB[itrkB].dca();

			TLorentzVector mom4d_v0( TVector3( pV0.x(), pV0.y(), pV0.z() ), eV0 );
			float etav0 = mom4d_v0.Eta();
			float ptv0  = mom4d_v0.Pt();
			//float phiv0 = mom4d_v0.Phi();

			// dca between two daughters
			float dcaDaughters = (dcaA-dcaB).mag();

			// v0 decay length
			StThreeVectorF v0toPV = v0 - pvtx;
			float decayLength = v0toPV.mag(); 
			float v0dir = v0toPV.x()*pV0.x() + v0toPV.y()*pV0.y() + v0toPV.z()*pV0.z();

			// v0 dca to PV
			float angle = v0toPV.angle(pV0);  // angle uses acos() which retuns a value within 0-pi
			float dcaV0 = v0toPV.mag()*TMath::Sin(angle); // always positive value
			if( v0dir<0 ) continue;
			if( ptv0<0.5 ) continue;
			if( fabs(etav0)>=1 ) continue;

			if( !(dcaDaughters<(mDcaDh[iCent]) && decayLength>(mDL[iCent]) && dcaV0<(mDcaV0[iCent])) ) continue;
			if( !(fdcaA>(mDcaP[iCent]) && fdcaB>(mDcaPi[iCent])) ) continue;
			hLambdaMassTopoCut[fCH]->Fill( invM );
			

		} // itrkB
	}// itrkA

}

void StBES2QaMaker::TrackLoop(int& myrefMult, int& mygrefMult){

	int countPTrk = 0;
	int countGTrk = 0;
	TVector3 pvtx_3d = mPicoDst->event()->primaryVertex();

	// track loop
	for( unsigned int i=0; i<mPicoDst->numberOfTracks(); i++ ){

		StPicoTrack *ptrk = (StPicoTrack*)mPicoDst->track(i);
		if( !ptrk ) continue;  

		// track info.
		float dca = ptrk->gDCA( pvtx_3d ).Mag();
		int nHitsFit = ptrk->nHitsFit();
		//short ch  = ptrk->charge();
		//short iCH = (ch>0) ? 0 : 1;
		if( nHitsFit<10 )  continue; 
		if( ptrk->nHitsFit()<0.52*ptrk->nHitsMax() ) continue;

		// global tracks
		//float gpt  = ptrk->gMom().Perp(); 
		//float gphi = ptrk->gMom().Phi();  
		float geta = ptrk->gMom().PseudoRapidity();
		float gmom = ptrk->gMom().Mag();
		if( gmom>1.e-10 && geta<0 && geta>-2 ) countGTrk++;
		

		// primary track info.
		//float pt  = ptrk->pMom().Perp(); 
		//float phi = ptrk->pMom().Phi();  
		float eta = ptrk->pMom().PseudoRapidity();
		float mom = ptrk->pMom().Mag();
		if( mom>1.e-10 && fabs(dca)<3 && eta<0 && eta>-2 ) countPTrk++;
	}
	myrefMult  = countPTrk;
	mygrefMult = countGTrk;
}

void StBES2QaMaker::CalcQvEPD( StPicoDst* pD, QV fQv[NordEPD][NsubEPD], int fCent, int fVz ){

	if( fCent<0 || fCent>=NqvCent || fVz<0 || fVz>=NqvZvtx ) return;

	Int_t nepdHits = pD->numberOfEpdHits();
	//hNEpdHits->Fill(nepdHits);

	float nepdMIPsE = 0;
	float nepdMIPsW = 0;
	float nepdRingE[16] = {0};
	float nepdRingW[16] = {0};
	float fqv[NordEPD][3][3] = {{{0.0}}}; //nth-order, 3(E,W,E+W) x nSub, (qx,qy,qw)
	if( NsubEPD!=3 ){ cout <<"CalcQvEPD(): Recheck NsubEPD!!! Skip EP calc. for EPD"<< endl; return; }

	StEpdGeom * mEpdGeom = new StEpdGeom();
	for(Int_t iHit=0; iHit<nepdHits; iHit++) {

		StPicoEpdHit *epdHit = pD->epdHit(iHit);
		if( !epdHit ) continue;
		//hEpdHitAdc->Fill(epdHit->nMIP());

		Short_t side_EW = epdHit->side(); //+1 for West and -1 for East
		float mip = epdHit->nMIP(); // gain calibrated energy loss in tile, in units of Landau MPV for one MIP
		float weight = (mip>3) ? 3 : mip;
		//float weight = (mip>6) ? 6 : mip;

		int iEW = (side_EW + 1)/2;
		//int itile = epdHit->tile() - 1; // 0-30, id in each supersector, 0=most inner
		//int isect = epdHit->position() - 1; // 0-11, supersector id on a wheel

		if( mip<0.3 ) continue;
		if( side_EW==-1 ) nepdMIPsE += weight;
		if( side_EW== 1 ) nepdMIPsW += weight;


		TVector3 StraightLine = mEpdGeom->TileCenter(epdHit->id()) - pD->event()->primaryVertex(); // This give eta-phi plot with gap 
		//TVector3 StraightLine = mEpdGeom->RandomPointOnTile(epdHit->id()) - picoEvent->primaryVertex();  // I guess this is just for plot-visualization
		double phi_epd = StraightLine.Phi();
		//double eta_epd = StraightLine.Eta();
		//if(side_EW == -1) hEpdEtaPhiEast->Fill(eta_epd,phi_epd,mip);
		//if(side_EW ==  1) hEpdEtaPhiWest->Fill(eta_epd,phi_epd,mip);

		int iring = epdHit->row() - 1; // (1-16)-1 -> 0-15

		if( side_EW==-1 ) nepdRingE[iring] += weight;
		if( side_EW== 1 ) nepdRingW[iring] += weight;
		

		// flow vectors
		//if( iring<2 ) continue; // only for Isobar because of saturation
		for( int ih=0; ih<NordEPD; ih++ ){
			// iEW=(0,1)=(E,W)
			fqv[ih][iEW][0] += weight * cos( (ih+1.0)*phi_epd );
			fqv[ih][iEW][1] += weight * sin( (ih+1.0)*phi_epd );
			fqv[ih][iEW][2] += weight;

			// combined
			float sgn = 1.0;
			if( ih==0 && iEW==1 ) sgn = -1.0; // just to be constent with other detectors
			fqv[ih][2][0] += sgn * weight * cos( (ih+1.0)*phi_epd );
			fqv[ih][2][1] += sgn * weight * sin( (ih+1.0)*phi_epd );
			fqv[ih][2][2] += weight;
		}

		//	cout<<" iHit= "<<iHit<<" id= "<<epdHit->id()<<" ADC= "<<epdHit->adc()<<" nMIP= "<<epdHit->nMIP()<<endl;
	} // end of iHit loop 
	// End of Epd block __________________________________________________

	//for( int iring=0; iring<16; iring++ ){
	//	hrefMultvsEpdRingE[iring]->Fill( pD->event()->refMult(), nepdRingE[iring] );
	//	hrefMultvsEpdRingW[iring]->Fill( pD->event()->refMult(), nepdRingW[iring] );
	//}
	//hEpdEvsW->Fill( nepdMIPsE, nepdMIPsW );
	//if( mCalibMode==0 ){
	//	hRunidvsEpdMipE->Fill( pD->event()->runId(), nepdMIPsE );
	//	hRunidvsEpdMipW->Fill( pD->event()->runId(), nepdMIPsW );
	//}
	if( nepdHits<1 ) return;

	// Set raw Q-vectors and Calibration //
	// =======================================================
	for( int ih=0; ih<NordEPD; ih++ ){
		for( int isub=0; isub<3; isub++ ){

			fQv[ih][isub].set_Qv( fqv[ih][isub][0], fqv[ih][isub][1], fqv[ih][isub][2], dEPD, isub, ih, fCent, fVz );
			pQVMaker->setProfiles( hEPDrecX[ih][isub], hEPDrecY[ih][isub], hEPDfltC[ih][isub], hEPDfltS[ih][isub], hEPDEP[ih][isub][fCent] );
			pQVMaker->doCalibration( &(fQv[ih][isub]) );
		}
	}
	// =======================================================
}
