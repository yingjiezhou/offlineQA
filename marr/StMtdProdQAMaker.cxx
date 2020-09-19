#include <iostream>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <iterator>
#include <bitset>

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TLorentzVector.h"

#include "StEventTypes.h"
#include "StThreeVectorF.hh"
#include "PhysicalConstants.h"
#include "StMemoryInfo.hh"
#include "StMessMgr.h"
#include "StTimer.hh"
#include "StEnumerations.h"

#include "StEvent.h"
#include "StVertex.h"
#include "StTriggerData.h"
#include "StTrack.h"
#include "StDcaGeometry.h"
#include "StDedxPidTraits.h"
#include "StTrackPidTraits.h"
#include "StBTofPidTraits.h"
#include "StBTofCollection.h"
#include "StBTofHit.h"
#include "StBTofRawHit.h"
#include "StBTofHeader.h"
#include "StMtdCollection.h"
#include "StMtdHeader.h"
#include "StMtdRawHit.h"
#include "StMtdHit.h"
#include "StMtdPidTraits.h"
#include "StTpcDedxPidAlgorithm.h"
#include "StMuDSTMaker/COMMON/StMuBTofPidTraits.h"
#include "StMuDSTMaker/COMMON/StMuBTofHit.h"
#include "StarClassLibrary/StParticleDefinition.hh"

#include "StMuDSTMaker/COMMON/StMuEmcCollection.h"
#include "StMuDSTMaker/COMMON/StMuEmcPoint.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StEmcCollection.h"
#include "StEmcCluster.h"
#include "StEmcDetector.h"
#include "StEmcModule.h"
#include "StEmcClusterCollection.h"
#include "StEmcPoint.h"
#include "StEmcRawHit.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/others/emcDetectorName.h"
#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "StEmcRawMaker/StEmcRawMaker.h"
#include "StEmcRawMaker/defines.h"

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuMtdCollection.h"
#include "StMuDSTMaker/COMMON/StMuMtdHeader.h"
#include "StMuDSTMaker/COMMON/StMuMtdRawHit.h"
#include "StMuDSTMaker/COMMON/StMuMtdHit.h"
#include "StMuDSTMaker/COMMON/StMuMtdPidTraits.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoMtdHit.h"
#include "StPicoEvent/StPicoMtdTrigger.h"
#include "StPicoEvent/StPicoMtdPidTraits.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"

// #include "StPicoDstMaker/StPicoDst.h"
// #include "StPicoDstMaker/StPicoEvent.h"
// #include "StPicoDstMaker/StPicoTrack.h"
// #include "StPicoDstMaker/StPicoMtdHit.h"
// #include "StPicoDstMaker/StPicoMtdTrigger.h"
// #include "StPicoDstMaker/StPicoMtdPidTraits.h"
// #include "StPicoDstMaker/StPicoBTofPidTraits.h"
// #include "StPicoDstMaker/StPicoEmcPidTraits.h"

#include "StMtdUtil/StMtdGeometry.h"
#include "StMtdProdQAMaker.h"
#include "tables/St_mtdModuleToQTmap_Table.h"
#include "tables/St_mtdQTSlewingCorr_Table.h"
#include "tables/St_mtdQTSlewingCorrPart2_Table.h"

ClassImp(StMtdProdQAMaker)

//_____________________________________________________________________________
StMtdProdQAMaker::StMtdProdQAMaker(const Char_t *name) : 
  StMaker(name),
  mMuDst(NULL), mPicoDst(NULL), mDataType(0),
  mYear(-1),
  mRunId(-1), mFirstRun(0), mLastRun(0),
  mSeparateTrig(kTRUE),
  mPrintConfig(kTRUE),
  mOutputFileName(""), mOutputFile(NULL),
  mEmcPosition(NULL),
  mMaxTpcVz(1e4), mMaxDiffVz(3),
  mMinTrkPt(1.), mMaxTrkPt(1e4), 
  mMinTrkPhi(0.), mMaxTrkPhi(2*pi), 
  mMinTrkEta(-1), mMaxTrkEta(1),
  mMinNHitsFit(15), mMinNHitsDedx(10), mMinNHitsFrac(0.52), 
  mMaxDca(3),
  mMinNsigmaPi(-2.), mMaxNsigmaPi(3.),
  mMinMuonDeltaZ(-20.), mMaxMuonDeltaZ(20.),
  mMinMuonDeltaY(-20), mMaxMuonDeltaY(20),
  mMinMuonDeltaTof(-1e4), mMaxMuonDeltaTof(1.0),
  mMaxMuonDca(3),
  mBTofMatch(kFALSE), mMtdHitTrigger(kTRUE)
{
  // default constructor
  mTrigTime[0] = 0;
  mTrigTime[1] = 0;
  for(int i=0;i<4;i++) mEmcGeom[i] = 0;
  for(int i=0; i<kNtrig; i++)
    {      
      mhEventStat[i] = NULL;
      mhVertexXYRanking0[i] = NULL;
      mhTpcVzRanking0[i] = NULL;
      mhDiffVzRanking0[i] = NULL;
      mhTpcVtxIndex[i] = NULL;
      mhVtxIndexVsClosest[i] = NULL;
      mhVertexXY[i] = NULL;
      mhTpcVz[i] = NULL;
      mhVpdVz[i] = NULL;
      mhDiffVz[i] = NULL;
      mhDiffVzVsTpcVz[i] = NULL;
      mhDiffVzVsVpdVz[i] = NULL;
      mhRefMult[i] = NULL;
      mhgRefMult[i] = NULL;
      mhgRefMultVsRefMult[i] = NULL;
      mhTpcVzVsgRef[i] = NULL;
      mhDiffVzVsgRef[i] = NULL;
      mhZdcRateVsgRef[i] = NULL;
      mhBbcRateVsgRef[i] = NULL;
      mhpTrkNHitsVsPt[i] = NULL;
      mhpTrkNHits[i] = NULL;
      mhpTrkNDedxVsPt[i] = NULL;
      mhpTrkNDedx[i] = NULL;
      mhpTrkNHitsFracVsPt[i] = NULL;
      mhpTrkNHitsFrac[i] = NULL;
      mhpTrkDcaVsPt[i] = NULL;
      mhpTrkN[i] = NULL;
      mhpTrkNMthMtd[i] = NULL;
      mhpTrkPt[i] = NULL;
      mhpTrkEtaVsPt[i] = NULL;
      mhpTrkPhiVsPt[i] = NULL;
      mhpTrkPhiVsEta[i] = NULL;
      mhpTrkDedxVsMom[i] = NULL;
      mhNsigmaEVsMom[i] = NULL;
      mhNsigmaPiVsMom[i] = NULL;
      mhNsigmaKVsMom[i] = NULL;
      mhNsigmaPVsMom[i] = NULL;
      mhBetaVsMom[i] = NULL;
      mhM2VsMom[i] = NULL;
      mhpTrkN[i] = NULL;
      mhpTrkPt[i] = NULL;
      mhpTrkNMthMtd[i] = NULL;
      mhAdc0[i]          = NULL;
      mhEoverPVsPt[i]    = NULL;
      mhNetaVsPt[i]      = NULL;
      mhNphiVsPt[i]      = NULL;
      mhDistZVsPt[i]     = NULL;
      mhDistPhiVspt[i]   = NULL;
      mhNElecPos[i]      = NULL;
      mhNElecNeg[i]      = NULL;
      mhElecPt[i]        = NULL;
      mhElecPhiVsEta[i]  = NULL;
      mhNQtSignal[i] = NULL;
      mhNMT101Signal[i] = NULL;
      mhNTF201Signal[i] = NULL;
      mhMtdVpdTacDiffMT001[i] = NULL;
      mhMtdVpdMthTacDiffMT001[i] = NULL;
      mhMtdVpdTacDiffMT101[i] = NULL;
      mhMtdVpdMthTacDiffMT101[i] = NULL;
      for(Int_t j=0; j<kNQTboard; j++)
	{
	  for(Int_t k=0; k<2; k++)
	    {
	      mhMtdTacSumMixvsMxq[j][k][i] = NULL;
	    }
	}
      mhMtdNRawHits[i] = NULL;
      mhMtdRawHitMap[i] = NULL;
      mhMtdNHits[i] = NULL;
      mhMtdHitMap[i] = NULL;
      mhMtdHitTrigTime[i]= NULL;
      mhMtdHitLeTimeDiff[i]= NULL;
      mhMtdNTrigHits[i] = NULL;
      mhMtdTrigHitMap[i] = NULL;
      mhMtdNMthHits[i] = NULL;
      mhMtdMthHitMap[i] = NULL;
      mhMtdNMthTrigHits[i] = NULL;
      mhMtdMthTrigHitMap[i] = NULL;
      mhLocalYVsgChan[i] = NULL;
      mhLocalZVsgChan[i] = NULL;
      mhMtdTofVsgChan[i] = NULL;
      mhExpTofVsgChan[i] = NULL;
      mhDeltaZ[i] = NULL;
      mhDzVsPt[i] = NULL;
      mhDzVsgChan[i] = NULL;
      mhDeltaY[i] = NULL;
      mhDyVsPt[i] = NULL;
      mhDyVsgChan[i] = NULL;
      mhDeltaTof[i] = NULL;
      mhDTofVsPt[i] = NULL;
      mhDTofVsgChan[i] = NULL;
      mhNMuonPos[i] = NULL;
      mhNMuonNeg[i] = NULL;
      mhMuonPt[i] = NULL;
      mhMuonPhiVsEta[i] = NULL;
      mhMuonMap[i] = NULL;
      mhNULpair[i] = NULL;
      mhNLSpairPos[i] = NULL;
      mhNLSpairNeg[i] = NULL;
      mhInvMvsPtUL[i] = NULL;
      mhInvMvsPtLSpos[i] = NULL;
      mhInvMvsPtLSneg[i] = NULL;
      mhInvMUL[i] = NULL;
      mhInvMLSpos[i] = NULL;
      mhInvMLSneg[i] = NULL;
      mhRunStat[i] = NULL;
      mhBBCrateVsRun[i] = NULL;
      mhZDCrateVsRun[i] = NULL;
      mhRefMultVsRun[i] = NULL;
      mhgRefMultVsRun[i] = NULL;
      mhTpcVxVsRun[i] = NULL;
      mhTpcVyVsRun[i] = NULL;
      mhTpcVzVsRun[i] = NULL;
      mhVpdVzVsRun[i] = NULL;
      mhDiffVzVsRun[i] = NULL;
      mhpTrkPtVsRun[i] = NULL;
      mhpTrkEtaVsRun[i] = NULL;
      mhpTrkPhiVsRun[i] = NULL;
      mhpTrkDcaVsRun[i] = NULL;
      mhNHitsFitVsRun[i] = NULL;
      mhNHitsPossVsRun[i] = NULL;
      mhNHitsDedxVsRun[i] = NULL;
      mhDedxVsRun[i] = NULL;
      mhNsigmaPiVsRun[i] = NULL;
      mhNsigmaEVsRun[i] = NULL;
      mhNsigmaKVsRun[i] = NULL;
      mhNsigmaPVsRun[i] = NULL;
      mhBetaVsRun[i] = NULL;
      mhNElectron[i] = NULL;
      mhNPositron[i] = NULL;
      mhNMtdRawHitsVsRun[i] = NULL;
      mhNMtdHitsVsRun[i] = NULL;
      mhNMtdTrigHitsVsRun[i] = NULL;
      mhNMtdMthHitsVsRun[i] = NULL;
      mhNMuonPosVsRun[i] = NULL;
      mhNMuonNegVsRun[i] = NULL;
      mhNMuonPairULVsRun[i] = NULL;
      mhNMuonPairLSPosVsRun[i] = NULL;
      mhNMuonPairLSNegVsRun[i] = NULL;
      mhNJpsiVsRun[i] = NULL;
    }
}
 
//_____________________________________________________________________________
StMtdProdQAMaker::~StMtdProdQAMaker()
{
  // default destructor
  if(mOutputFile)  delete mOutputFile;
  if(mEmcPosition) delete mEmcPosition;
  for(int i=0;i<4;i++) {
    mEmcGeom[i] = 0;
  }
}

//_____________________________________________________________________________
Int_t StMtdProdQAMaker::Init()
{
  mEmcPosition = new StEmcPosition();
  for(int i=0;i<4;i++)
    {
      mEmcGeom[i] = StEmcGeom::getEmcGeom(detname[i].Data());
    }

  if(mPrintConfig) printConfig();

  if(mYear==2013) { mFirstRun = 14130000; mLastRun = 14162000; }
  if(mYear==2014) { mFirstRun = 15074000; mLastRun = 15168000; }
  if(mYear==2015) { mFirstRun = 16040000; mLastRun = 16170000; }
  if(mYear==2016) { mFirstRun = 17030000; mLastRun = 17180000; }
  if(mYear==2017) { mFirstRun = 18050000; mLastRun = 18180000; }
  initTrigIDs();
  bookHistos();

  return kStOK;
}


//_____________________________________________________________________________
Int_t StMtdProdQAMaker::InitRun(const Int_t runNumber)
{
  // initialize maps
  memset(mModuleToQT,-1,sizeof(mModuleToQT));
  memset(mModuleToQTPos,-1,sizeof(mModuleToQTPos));
  memset(mQTtoModule,-1,sizeof(mQTtoModule));

  // obtain maps from DB
  LOG_INFO << "Retrieving mtdModuleToQTmap table from database ..." << endm;
  TDataSet *dataset = GetDataBase("Geometry/mtd/mtdModuleToQTmap");
  St_mtdModuleToQTmap *mtdModuleToQTmap = static_cast<St_mtdModuleToQTmap*>(dataset->Find("mtdModuleToQTmap"));
  if(!mtdModuleToQTmap)
    {
      LOG_ERROR << "No mtdModuleToQTmap table found in database" << endm;
      return kStErr;
    }
  mtdModuleToQTmap_st *mtdModuleToQTtable = static_cast<mtdModuleToQTmap_st*>(mtdModuleToQTmap->GetTable());

  for(Int_t i=0; i<gMtdNBacklegs; i++)
    {
      for(Int_t j=0; j<gMtdNModules; j++)
	{
	  Int_t index = i*5 + j;
	  Int_t qt = mtdModuleToQTtable->qtBoardId[index];
	  Int_t channel = mtdModuleToQTtable->qtChannelId[index];
	  mModuleToQT[i][j]    = qt;
	  if(channel<0)
	    {
	      mModuleToQTPos[i][j] = channel;
	    }
	  else
	    {
	      if(channel%8==1) mModuleToQTPos[i][j] = 1 + channel/8 * 2;
	      else             mModuleToQTPos[i][j] = 2 + channel/8 * 2;
	    }
	  if(mModuleToQT[i][j]>0 && mModuleToQTPos[i][j]>0)
	    mQTtoModule[mModuleToQT[i][j]-1][mModuleToQTPos[i][j]-1] = j + 1;
	}
    }

  // online slewing correction for QT board
  memset(mQTSlewBinEdge,-1,sizeof(mQTSlewBinEdge));
  memset(mQTSlewCorr,-1,sizeof(mQTSlewCorr));
  LOG_INFO << "Retrieving mtdQTSlewingCorr table from database ..." << endm;
  dataset = GetDataBase("Calibrations/mtd/mtdQTSlewingCorr");
  St_mtdQTSlewingCorr *mtdQTSlewingCorr = static_cast<St_mtdQTSlewingCorr*>(dataset->Find("mtdQTSlewingCorr"));
  if(!mtdQTSlewingCorr)
    {
      LOG_ERROR << "No mtdQTSlewingCorr table found in database" << endm;
      return kStErr;
    }
  mtdQTSlewingCorr_st *mtdQTSlewingCorrtable = static_cast<mtdQTSlewingCorr_st*>(mtdQTSlewingCorr->GetTable());
  for(int j=0; j<4; j++)
    {
      for(int i=0; i<16; i++)
	{
	  for(Int_t k=0; k<8; k++)
	    {
	      Int_t index = j*16*8 + i*8 + k;
	      mQTSlewBinEdge[j][i][k] = (int) mtdQTSlewingCorrtable->slewingBinEdge[index];
	      mQTSlewCorr[j][i][k] = (int) mtdQTSlewingCorrtable->slewingCorr[index];
	    }
	}
    }
  if(mYear==2016)
    {
      LOG_INFO << "Retrieving mtdQTSlewingCorrPart2 table from database ..." << endm;
      dataset = GetDataBase("Calibrations/mtd/mtdQTSlewingCorrPart2");
      if(dataset)
	{
	  St_mtdQTSlewingCorrPart2 *mtdQTSlewingCorr2 = static_cast<St_mtdQTSlewingCorrPart2*>(dataset->Find("mtdQTSlewingCorrPart2"));
	  mtdQTSlewingCorrPart2_st *mtdQTSlewingCorrtable2 = static_cast<mtdQTSlewingCorrPart2_st*>(mtdQTSlewingCorr2->GetTable());
	  for(int j=0; j<4; j++)
	    {
	      for(int i=0; i<16; i++)
		{
		  for(Int_t k=0; k<8; k++)
		    {
		      Int_t index = j*16*8 + i*8 + k;
		      mQTSlewBinEdge[j+4][i][k] = (int) mtdQTSlewingCorrtable2->slewingBinEdge[index];
		      mQTSlewCorr[j+4][i][k] = (int) mtdQTSlewingCorrtable2->slewingCorr[index];
		    }
		}
	    }
	}
    }
  LOG_INFO << "===== End retrieving mtdQTSlewingCorr =====" << endm;

  return kStOK;
}

//_____________________________________________________________________________
void StMtdProdQAMaker::initTrigIDs()
{
  printf("==================================\n");
  LOG_INFO << "Initialize triggers for year: " << mYear << endm;

  for(int i=0; i<kNtrig; i++) mTriggerIDs[i].clear();

  if (mYear == 2013)
    {
      const Int_t ntrig = 2; 
      Int_t di_muon[ntrig] = {430103, 430113};
      Int_t sg_muon[ntrig] = {430101, 430111};
      Int_t el_muon[ntrig] = {430102, 430112};
      for(Int_t i=0; i<ntrig; i++)
	{
	  mTriggerIDs[0].push_back(di_muon[i]);
	  mTriggerIDs[1].push_back(sg_muon[i]);
	  mTriggerIDs[2].push_back(el_muon[i]);
	}
      mTriggerIDs[2].push_back(430122);
    }
  else if(mYear == 2014)
    {
      const Int_t ntrig = 5; 
      Int_t di_muon[ntrig] = {450601, 450611, 450621, 450631, 450641};
      Int_t sg_muon[ntrig] = {450600, 450610, 450620, 450630, 450640};
      Int_t el_muon[ntrig] = {450602, 450612, 450622, 450632, 450642};
      for(Int_t i=0; i<ntrig; i++)
	{
	  mTriggerIDs[0].push_back(di_muon[i]);
	  mTriggerIDs[1].push_back(sg_muon[i]);
	  mTriggerIDs[2].push_back(el_muon[i]);
	}
      mTriggerIDs[3].push_back(450604);
      mTriggerIDs[3].push_back(450605);
      mTriggerIDs[3].push_back(450606);
    }
  else if(mYear == 2015)
    {
      // pp 200 GeV
      const Int_t ntrig = 5; 
      Int_t di_muon[ntrig] = {470602, 480602, 490602, 500602, 510602};
      Int_t sg_muon[ntrig] = {470600, 480600, 490600, 500600, 510600};
      Int_t el_muon[ntrig] = {470601, 480601, 490601, 500601, 510601};
      for(Int_t i=0; i<ntrig; i++)
	{
	  mTriggerIDs[0].push_back(di_muon[i]);
	  mTriggerIDs[1].push_back(sg_muon[i]);
	  mTriggerIDs[2].push_back(el_muon[i]);
	}
    }
  else if (mYear == 2016)
    {
      const Int_t ntrig = 3; 
      Int_t di_muon[ntrig] = {520602, 520612, 520622};
      Int_t sg_muon[ntrig] = {520604, 520614, 520624};
      Int_t el_muon[ntrig] = {520606, 520616, 520626};
      for(Int_t i=0; i<ntrig; i++)
	{
	  mTriggerIDs[0].push_back(di_muon[i]);
	  mTriggerIDs[1].push_back(sg_muon[i]);
	  mTriggerIDs[2].push_back(el_muon[i]);
	}
      mTriggerIDs[0].push_back(46);
      mTriggerIDs[0].push_back(47);
      mTriggerIDs[0].push_back(520803);

      mTriggerIDs[3].push_back(520005);  //VPD-ZDC-novtx
      mTriggerIDs[3].push_back(520015);  //VPD-ZDC-novtx
      mTriggerIDs[3].push_back(570005);  //VPD-ZDC-novtx
    }
  else if (mYear == 2017)
    {
      // pp 510 GeV
      const Int_t ntrig = 2; 
      Int_t di_muon[ntrig] = {570602, 580602};
      Int_t sg_muon[ntrig] = {0, 0};
      Int_t el_muon[ntrig] = {0, 0};
      for(Int_t i=0; i<ntrig; i++)
	{
	  mTriggerIDs[0].push_back(di_muon[i]);
	  mTriggerIDs[1].push_back(sg_muon[i]);
	  mTriggerIDs[2].push_back(el_muon[i]);
	}
      mTriggerIDs[3].push_back(570004);  //VPD-novtx
    }
}


//_____________________________________________________________________________
Int_t StMtdProdQAMaker::Finish()
{  
  if(mOutputFile)
    {
      mOutputFile->cd();
      mOutputFile->Write();
      mOutputFile->Close();
      LOG_INFO << "StMtdProdQAMaker::Finish() -> write out histograms to " << mOutputFileName.Data() << endm;
    }

  return kStOK;
}

//_____________________________________________________________________________
Int_t StMtdProdQAMaker::Make()
{

  // Check the availability of input data
  Int_t iret;
  if(mDataType==0)
    {
      StMuDstMaker *muDstMaker = (StMuDstMaker*) GetMaker("MuDst");
      if(muDstMaker) 
	{
	  mMuDst = muDstMaker->muDst();
	  if(mMuDst)  iret = processMuDst();
	}
    }

  else if(mDataType==1)
    {
      LOG_DEBUG << "Running on PicoDst ..." << endm;
      StPicoDstMaker *picoDstMaker = (StPicoDstMaker*) GetMaker("picoDst");
      if(picoDstMaker)
	{
	  mPicoDst = picoDstMaker->picoDst();
	  if(mPicoDst) iret = processPicoDst();
	}
    }

  return iret;
}

//_____________________________________________________________________________
void StMtdProdQAMaker::getTrigIndex()
{
  for(int i=0; i<kNtrig; i++) mTrigIndex[i] = false;
  if(!mSeparateTrig)
    {
      for(Int_t j=0; j<kNtrig; j++)
	{
	  if(j==0) mTrigIndex[j] = true;
	  else     mTrigIndex[j] = false;
	}
      return;
    }

  if(mDataType==0)
    {
      for(Int_t j=0; j<kNtrig; j++)
	{
	  bool found = false;
	  for(UInt_t i=0; i<mTriggerIDs[j].size(); i++)
	    {
	      if(mMuDst->event()->triggerIdCollection().nominal().isTrigger(mTriggerIDs[j][i]))
		{
		  //cout << mTriggerIDs[j][i] << endl;
		  found = true;
		  break;
		}
	    }
	  mTrigIndex[j] = found;
	}
      // reject di-muon event overlap with other triggers
      // but should have been rejected
      StMuMtdHeader *muMtdHeader = mMuDst->mtdHeader();
      if(muMtdHeader && muMtdHeader->shouldHaveRejectEvent()==1)
	mTrigIndex[0] = false;
    }
  else if(mDataType==1)
    {
      for(Int_t j=0; j<kNtrig; j++)
	{
	  bool found = false;
	  for(UInt_t i=0; i<mTriggerIDs[j].size(); i++)
	    {
	      if(mPicoDst->event()->isTrigger(mTriggerIDs[j][i]))
		{
		  found = true;
		  break;
		}
	    }
	  mTrigIndex[j] = found;
	}
      StPicoMtdTrigger *trigger = (StPicoMtdTrigger*)mPicoDst->mtdTrigger(0);
      if(mYear==2014 && trigger && trigger->shouldHaveRejectEvent()==1)
	mTrigIndex[0] = false;
    }
}

//_____________________________________________________________________________
void StMtdProdQAMaker::fillHisto(TH1F **h, float x)
{
  if(!mSeparateTrig)
    {
      if(h[0]) h[0]->Fill(x);
    }
  else
    {
      for(int i=0; i<kNtrig; i++)
	{
	  if(mTrigIndex[i] && h[i]) h[i]->Fill(x);
	}
    }
}

//_____________________________________________________________________________
void StMtdProdQAMaker::fillHisto(TH2F **h, float x, float y)
{
  if(!mSeparateTrig)
    {
      if(h[0]) h[0]->Fill(x,y);
    }
  else
    {
      for(int i=0; i<kNtrig; i++)
	{
	  if(mTrigIndex[i] && h[i]) h[i]->Fill(x,y);
	}
    }
}

//_____________________________________________________________________________
void StMtdProdQAMaker::fillTProfile(TProfile **h, float y)
{
  if(mRunId<0) return;

  if(!mSeparateTrig)
    {
      if(h[0]) h[0]->Fill(mRunId,y);
    }
  else
    {
      for(int i=0; i<kNtrig; i++)
	{
	  if(mTrigIndex[i] && h[i]) h[i]->Fill(mRunId,y);
	}
    }
}

//_____________________________________________________________________________
Int_t StMtdProdQAMaker::processMuDst()
{
  //cout << "Work on MuDst: " << mMuDst->event()->eventInfo().id() << endl;
  mRunId   = mMuDst->event()->runId();
  getTrigIndex();

  // Event statistics
  for(Int_t j=0; j<kNtrig; j++)
    {
      if(mhEventStat[j]) mhEventStat[j]->Fill(0.5);
    }

  // count triggers
  fillHisto(mhEventStat,1.5);
  fillHisto(mhRunStat,mRunId);

  // Cut on vertex z
  /// vpd vz
  float vpdVz = -1e4;
  StBTofHeader *tofHeader = mMuDst->btofHeader();
  if(tofHeader) 
    vpdVz = tofHeader->vpdVz();
  fillHisto(mhVpdVz,vpdVz);

  if(abs(vpdVz)>500)
    {
      LOG_DEBUG << "Bad vpd vertext: " << vpdVz << endm;
      return kStOK;
    }
  fillHisto(mhEventStat,2.5);

  /// default primary vertex
  StMuPrimaryVertex* priVertex = mMuDst->primaryVertex();
  if(priVertex)
    {
      StThreeVectorF verPos = priVertex->position();
      fillHisto(mhTpcVzRanking0, verPos.z());
      fillHisto(mhVertexXYRanking0, verPos.x(), verPos.y());
      fillHisto(mhDiffVzRanking0, verPos.z() - vpdVz);
    }

  //////////////////////////////////////
  // select the right vertex using VPD
  int index = 0;
  float dzmax = 999;
  int closest = -1;
  for(unsigned int i=0;i<mMuDst->numberOfPrimaryVertices();i++)
    {
      StMuPrimaryVertex *vtx = mMuDst->primaryVertex(i);
      if(!vtx) continue;
      float vz = vtx->position().z();
      if(abs(vpdVz-vz)<dzmax)
	{
	  dzmax = abs(vpdVz-vz);
	  closest = i;
	}
      if(fabs(vpdVz-vz)<mMaxDiffVz) 
	{
	  index = i;
	  break;
	}
    }
  fillHisto(mhTpcVtxIndex, index);
  fillHisto(mhVtxIndexVsClosest, closest, index);
  mMuDst->setVertexIndex(index);
  
  priVertex = mMuDst->primaryVertex();
  if(!priVertex) return kStWarn;
  StThreeVectorF verPos = priVertex->position();
  float tpcVz = verPos.z();
  fillHisto(mhTpcVz, tpcVz);
  fillHisto(mhVertexXY, verPos.x(), verPos.y());
  fillHisto(mhDiffVz, tpcVz - vpdVz);
  fillHisto(mhDiffVzVsTpcVz, tpcVz, tpcVz-vpdVz);
  fillHisto(mhDiffVzVsVpdVz, vpdVz, tpcVz-vpdVz);
  fillTProfile(mhTpcVxVsRun, verPos.x());
  fillTProfile(mhTpcVyVsRun, verPos.y());
  fillTProfile(mhTpcVzVsRun, tpcVz);
  fillTProfile(mhVpdVzVsRun, vpdVz);
  fillTProfile(mhDiffVzVsRun, tpcVz-vpdVz);  

  /// vertex cut
  if(mMaxTpcVz<1e4 && abs(tpcVz)>mMaxTpcVz) return kStOK;
  if(mMaxDiffVz<1e4 && abs(tpcVz-vpdVz)>mMaxDiffVz) return kStOK;
  fillHisto(mhEventStat,3.5);

  /// reference multiplicity
  StMuEvent *muEvent = mMuDst->event();
  StRunInfo runInfo = muEvent->runInfo();
  int refMult = muEvent->refMult();
  int gRefMult = muEvent->grefmult();
  float bbcRate = runInfo.bbcCoincidenceRate() * 1e-3; // kHz
  float zdcRate = runInfo.zdcCoincidenceRate() * 1e-3; // kHz
  fillHisto(mhRefMult, refMult);
  fillHisto(mhgRefMult, gRefMult);
  fillHisto(mhgRefMultVsRefMult, refMult, gRefMult);
  fillHisto(mhTpcVzVsgRef, gRefMult, tpcVz);
  fillHisto(mhDiffVzVsgRef, gRefMult, tpcVz-vpdVz);
  fillHisto(mhZdcRateVsgRef, gRefMult, zdcRate);
  fillHisto(mhBbcRateVsgRef, gRefMult, bbcRate);
  fillTProfile(mhBBCrateVsRun, bbcRate);
  fillTProfile(mhZDCrateVsRun, zdcRate);
  fillTProfile(mhRefMultVsRun, refMult);
  fillTProfile(mhgRefMultVsRun, gRefMult);
  
   
  /// primary tracks
  IntVec trkId;
  trkId.clear();
  int maxadc, neta, nphi;
  float energy, zdist, phidist;
  int nElec = 0, nPosi = 0;

  float nGoodTrack = 0, nMthGoodTrack = 0;
  int nTracks = mMuDst->numberOfPrimaryTracks();
  for(int i=0; i<nTracks; i++)
    {
      StMuTrack* pTrack = mMuDst->primaryTracks(i);
      if(!pTrack) continue;
      StThreeVectorF mom = pTrack->momentum();
      float p   = mom.mag();
      float pt  = mom.perp();
      float eta = mom.pseudoRapidity();
      float phi = rotatePhi(mom.phi());
      int charge = pTrack->charge();  
      int nHitsFit  = pTrack->nHitsFit(kTpcId);
      int nHitsDedx = pTrack->nHitsDedx();
      int nHitsMax  = pTrack->nHitsPoss(kTpcId);
      float nHitsFrac = nHitsFit*1./nHitsMax;
      float dca = pTrack->dca().mag();
      fillHisto(mhpTrkNHitsVsPt, pt, nHitsFit);
      fillHisto(mhpTrkNHits, nHitsFit*charge);
      fillHisto(mhpTrkNDedxVsPt, pt, nHitsDedx);
      fillHisto(mhpTrkNDedx, nHitsDedx*charge);
      fillHisto(mhpTrkNHitsFracVsPt, pt, nHitsFrac);
      fillHisto(mhpTrkNHitsFrac, nHitsFrac);
      fillHisto(mhpTrkDcaVsPt, pt, dca);

      if(!isValidTrack(pTrack)) continue;
      nGoodTrack++;
      float dedx = pTrack->dEdx()*1e6;
      float nSigmaE  = pTrack->nSigmaElectron();
      float nSigmaPi = pTrack->nSigmaPion();
      float nSigmaK  = pTrack->nSigmaKaon();
      float nSigmaP  = pTrack->nSigmaProton();
      fillHisto(mhpTrkPt, pt);
      fillHisto(mhpTrkEtaVsPt, pt, eta);
      fillHisto(mhpTrkPhiVsPt, pt, phi);
      fillHisto(mhpTrkPhiVsEta, eta, phi);
      fillHisto(mhpTrkDedxVsMom, p, dedx);
      fillHisto(mhNsigmaEVsMom, p, nSigmaE);
      fillHisto(mhNsigmaPiVsMom, p, nSigmaPi);
      fillHisto(mhNsigmaKVsMom, p, nSigmaK);
      fillHisto(mhNsigmaPVsMom, p, nSigmaP);
      fillTProfile(mhpTrkPtVsRun, pt);
      fillTProfile(mhpTrkEtaVsRun, eta);
      fillTProfile(mhpTrkPhiVsRun, phi);
      fillTProfile(mhpTrkDcaVsRun, dca);
      fillTProfile(mhNHitsFitVsRun, nHitsFit);
      fillTProfile(mhNHitsPossVsRun, nHitsMax);
      fillTProfile(mhNHitsDedxVsRun, nHitsDedx);
      fillTProfile(mhDedxVsRun, dedx);
      fillTProfile(mhNsigmaPiVsRun, nSigmaPi);
      fillTProfile(mhNsigmaEVsRun, nSigmaE);
      fillTProfile(mhNsigmaKVsRun, nSigmaK);
      fillTProfile(mhNsigmaPVsRun, nSigmaP);

      /// TOF
      if(pTrack->index2BTofHit()>-1)
	{
	  StMuBTofPidTraits btofPid = pTrack->btofPidTraits(); 
	  float beta = btofPid.beta();
	  if(beta!=0)
	    {
	      float m2 = pow(p,2) * (1/pow(beta,2)-1);
	      fillHisto(mhBetaVsMom, p, 1./beta);
	      fillHisto(mhM2VsMom, p, m2);
	      fillTProfile(mhBetaVsRun, 1./beta);
	    }
	}

      /// electrons
      if(getBemcInfo(pTrack,maxadc,energy,zdist,phidist,neta,nphi))
	{
	  fillHisto(mhAdc0, maxadc);
	  fillHisto(mhEoverPVsPt, pt, p/energy);
	  fillHisto(mhNetaVsPt, pt, neta);
	  fillHisto(mhNphiVsPt, pt, nphi);
	  fillHisto(mhDistZVsPt, pt, zdist);
	  fillHisto(mhDistPhiVspt, pt, phidist);
	  if(abs(pTrack->nSigmaElectron())<2 &&
	     p/energy>0.3 && p/energy<1.5)
	    {
	      fillHisto(mhElecPt, pt);
	      fillHisto(mhElecPhiVsEta, eta, phi);
	      if(charge>0) nPosi++;
	      else         nElec++;
	    }
	}

      /// muons
      int index = pTrack->index2MtdHit();
      if(index>-1)  
	{
	  nMthGoodTrack++;
	  const StMuMtdPidTraits mtdPid = pTrack->mtdPidTraits();
	  float dy     = mtdPid.deltaY();
	  float dz     = mtdPid.deltaZ();
	  float localy = mtdPid.yLocal();
	  float localz = mtdPid.zLocal();
	  float mtdtof = mtdPid.timeOfFlight();
	  float exptof = mtdPid.expTimeOfFlight();
	  float dtof   = mtdtof - exptof;

	  StMuMtdHit *hit = mMuDst->mtdHit(index);
	  float gChan   = (hit->backleg()-1)*60 + (hit->module()-1)*12 + hit->cell();
	  fillHisto(mhLocalYVsgChan, gChan, localy);
	  fillHisto(mhLocalZVsgChan, gChan, localz);
	  fillHisto(mhMtdTofVsgChan, gChan, mtdtof);
	  fillHisto(mhExpTofVsgChan, gChan, exptof);
	  fillHisto(mhDeltaZ, dz);
	  fillHisto(mhDzVsPt, pt, dz);
	  fillHisto(mhDzVsgChan, gChan, dz);
	  fillHisto(mhDeltaY, dy);
	  fillHisto(mhDyVsPt, pt, dy);
	  fillHisto(mhDyVsgChan, gChan, dy);
	  fillHisto(mhDeltaTof, dtof);
	  fillHisto(mhDTofVsPt, pt, dtof);
	  fillHisto(mhDTofVsgChan, gChan, dtof);

	  if(isMuonCandidate(pTrack))
	    {
	      trkId.push_back(i);
	      fillHisto(mhMuonMap, hit->backleg(), (hit->module()-1)*12+hit->cell()+1);
	    }
	}
    }
  fillHisto(mhpTrkN, nGoodTrack);
  fillHisto(mhpTrkNMthMtd, nMthGoodTrack);
  fillHisto(mhNElecPos, nPosi);
  fillHisto(mhNElecNeg, nElec);

  fillTProfile(mhNPositron, nPosi);
  fillTProfile(mhNElectron, nElec);

  /// trigger performance
  triggerData();

  /// MTD hits
  int nMtdRawHits = mMuDst->numberOfBMTDRawHit();
  for(int i=0; i<nMtdRawHits; i++)
    {
      StMuMtdRawHit *rawHit = (StMuMtdRawHit*)mMuDst->mtdRawHit(i);
      if(!rawHit) continue;
      fillHisto(mhMtdRawHitMap, rawHit->backleg(), rawHit->channel());
    }
  fillHisto(mhMtdNRawHits, nMtdRawHits);
  fillTProfile(mhNMtdRawHitsVsRun, nMtdRawHits);

  StMuMtdHeader *muMtdHeader = mMuDst->mtdHeader();
  if(muMtdHeader)
    {
      mTrigTime[0] = 25.*(muMtdHeader->triggerTime(1)&0xfff);
      mTrigTime[1] = 25.*(muMtdHeader->triggerTime(2)&0xfff);
    }

  int nMtdHits = mMuDst->numberOfMTDHit();
  int nTrigHit = 0, nMthHit = 0, nMthTrigHit = 0;
  for(int i=0; i<nMtdHits; i++)
    {
      StMuMtdHit *hit = mMuDst->mtdHit(i);
      if(!hit) continue;
      int backleg = hit->backleg();
      int module  = hit->module();
      int cell    = hit->cell();
      int gChan   = (backleg-1)*60 + (module-1)*12 + cell;
      fillHisto(mhMtdHitMap, backleg, (module-1)*12+cell+1);

      int tHub = getMtdHitTHUB(backleg);
      float tDiff = (hit->leadingEdgeTime().first+hit->leadingEdgeTime().second)/2 - mTrigTime[tHub-1];
      while(tDiff<0) tDiff += 51200;
      fillHisto(mhMtdHitTrigTime, gChan, tDiff);
      fillHisto(mhMtdHitLeTimeDiff, gChan, hit->leadingEdgeTime().second-hit->leadingEdgeTime().first);

      bool isTrig = isMtdHitFiredTrigger(hit);
      int index2track = hit->index2Global();
      bool isMth = false;
      if(index2track>-1)
	{
	  StMuTrack* gTrack = mMuDst->globalTracks(index2track);
	  if(isValidTrack(gTrack)) isMth = true;
	}

      if(isTrig)
	{
	  nTrigHit ++;
	  fillHisto(mhMtdTrigHitMap, backleg, (module-1)*12+cell+1);
	}

      if(isMth)
	{
	  nMthHit ++;
	  fillHisto(mhMtdMthHitMap, backleg, (module-1)*12+cell+1);
	}

      if(isMth && isTrig)
	{
	  nMthTrigHit ++;
	  fillHisto(mhMtdMthTrigHitMap, backleg, (module-1)*12+cell+1);
	}
    }
  fillHisto(mhMtdNHits, nMtdHits);
  fillHisto(mhMtdNTrigHits, nTrigHit);
  fillHisto(mhMtdNMthHits, nMthHit);
  fillHisto(mhMtdNMthTrigHits, nMthTrigHit);

  fillTProfile(mhNMtdHitsVsRun, nMtdHits);
  fillTProfile(mhNMtdTrigHitsVsRun, nTrigHit);
  fillTProfile(mhNMtdMthHitsVsRun, nMthHit);


  /// muon analysis
  int nPosMuon = 0, nNegMuon = 0;
  int nULpair = 0, nLSpairPos = 0, nLSpairNeg = 0;
  int nJpsi = 0;
  unsigned int nMuon = trkId.size();
  for(unsigned int i=0; i<nMuon; i++)
    {
      StMuTrack *ipTrack = mMuDst->primaryTracks(trkId[i]); 
      int iq = ipTrack->charge();
      const StThreeVectorF imom = ipTrack->momentum();
      if(iq>0) nPosMuon++;
      else     nNegMuon++;
      fillHisto(mhMuonPt, imom.perp());
      fillHisto(mhMuonPhiVsEta, imom.pseudoRapidity(),rotatePhi(imom.phi()));

      TLorentzVector imuon;
      imuon.SetXYZM(imom.x(),imom.y(),imom.z(),muMass);
      for(UInt_t j=i+1; j<nMuon; j++)
	{
	  StMuTrack *jpTrack = mMuDst->primaryTracks(trkId[j]);
	  int jq = jpTrack->charge();
	  const StThreeVectorF jmom = jpTrack->momentum();
	  TLorentzVector jmuon;
	  jmuon.SetXYZM(jmom.x(),jmom.y(),jmom.z(),muMass);

	  float pt1 = imom.perp(), pt2 = jmom.perp();
	  if(pt1<pt2) 
	    {
	      pt1 = jmom.perp();
	      pt2 = imom.perp();
	    }

	  TLorentzVector jpsi = imuon + jmuon;
	  float invmass = jpsi.M();
	  if(iq*jq<0)
	    {
	      nULpair++;
	      if(pt1>1.5) 
		{
		  fillHisto(mhInvMvsPtUL, invmass, jpsi.Pt());
		  fillHisto(mhInvMUL, invmass);
		  if(invmass>3.0 && invmass<3.2) nJpsi++;
		}
	    }
	  else 
	    {
	      if(pt1>1.5 && invmass>3.0 && invmass<3.2) nJpsi--;
	      if(iq>0)
		{
		  nLSpairPos++;
		  if(pt1>1.5) 
		    {
		      fillHisto(mhInvMvsPtLSpos, invmass, jpsi.Pt());
		      fillHisto(mhInvMLSpos, invmass);
		    }
		}
	      else
		{
		  nLSpairNeg++;
		  if(pt1>1.5) 
		    {
		      fillHisto(mhInvMvsPtLSneg, invmass, jpsi.Pt());
		      fillHisto(mhInvMLSneg, invmass);
		    }
		}
	    }
	}
    }

  fillHisto(mhNMuonPos, nPosMuon); 
  fillHisto(mhNMuonNeg, nNegMuon); 
  fillHisto(mhNULpair, nULpair); 
  fillHisto(mhNLSpairPos, nLSpairPos); 
  fillHisto(mhNLSpairNeg, nLSpairNeg);

  fillTProfile(mhNMuonPosVsRun, nPosMuon);
  fillTProfile(mhNMuonNegVsRun, nNegMuon);
  fillTProfile(mhNMuonPairULVsRun, nULpair);
  fillTProfile(mhNMuonPairLSPosVsRun, nLSpairPos);
  fillTProfile(mhNMuonPairLSNegVsRun, nLSpairNeg);
  fillTProfile(mhNJpsiVsRun, nJpsi);

  return kStOK;
}

//_____________________________________________________________________________
Int_t StMtdProdQAMaker::processPicoDst()
{
  StPicoEvent *picoEvent = mPicoDst->event();
  if(!picoEvent) return kStOK;

  getTrigIndex();
  mRunId = picoEvent->runId();

  // Event statistics
  for(Int_t j=0; j<kNtrig; j++)
    {
      if(mhEventStat[j]) mhEventStat[j]->Fill(0.5);
    }

  // count triggers
  fillHisto(mhEventStat,1.5);
  fillHisto(mhRunStat,mRunId);

  // Cut on vertex z
  /// vpd vz
  float vpdVz = picoEvent->vzVpd();
  fillHisto(mhVpdVz,vpdVz);
  if(abs(vpdVz)>300)
    {
      LOG_DEBUG << "Bad vpd vertext: " << vpdVz << endm;
      return kStOK;
    }
  fillHisto(mhEventStat,2.5);

  /// default primary vertex
  //StThreeVectorF verPos = picoEvent->primaryVertex();
  TVector3 priVtxPos = picoEvent->primaryVertex();
  StThreeVectorF verPos(priVtxPos.x(), priVtxPos.y(), priVtxPos.z());
  float tpcVz = verPos.z();
  fillHisto(mhTpcVz, tpcVz);
  fillHisto(mhVertexXY, verPos.x(), verPos.y());
  fillHisto(mhDiffVz, tpcVz - vpdVz);
  fillHisto(mhDiffVzVsTpcVz, tpcVz, tpcVz-vpdVz);
  fillHisto(mhDiffVzVsVpdVz, vpdVz, tpcVz-vpdVz);
  fillTProfile(mhTpcVxVsRun, verPos.x());
  fillTProfile(mhTpcVyVsRun, verPos.y());
  fillTProfile(mhTpcVzVsRun, tpcVz);
  fillTProfile(mhVpdVzVsRun, vpdVz);
  fillTProfile(mhDiffVzVsRun, tpcVz-vpdVz);  

  //printf("PiDst: tpc = %3.2f, vpd = %3.2f, dz= %3.2f\n",tpcVz,vpdVz,tpcVz-vpdVz);

  /// vertex cut
  if(mMaxTpcVz<1e4 && abs(tpcVz)>mMaxTpcVz) return kStOK;
  if(mMaxDiffVz<1e4 && abs(tpcVz-vpdVz)>mMaxDiffVz) return kStOK;
  fillHisto(mhEventStat,3.5);

  /// reference multiplicity
  int refMult = picoEvent->refMult();
  int gRefMult = picoEvent->grefMult();
  float bbcRate = picoEvent->BBCx() * 1e-3; // kHz
  float zdcRate = picoEvent->ZDCx() * 1e-3; // kHz
  fillHisto(mhRefMult, refMult);
  fillHisto(mhgRefMult, gRefMult);
  fillHisto(mhgRefMultVsRefMult, refMult, gRefMult);
  fillHisto(mhTpcVzVsgRef, gRefMult, tpcVz);
  fillHisto(mhDiffVzVsgRef, gRefMult, tpcVz-vpdVz);
  fillHisto(mhZdcRateVsgRef, gRefMult, zdcRate);
  fillHisto(mhBbcRateVsgRef, gRefMult, bbcRate);
  fillTProfile(mhBBCrateVsRun, bbcRate);
  fillTProfile(mhZDCrateVsRun, zdcRate);
  fillTProfile(mhRefMultVsRun, refMult);
  fillTProfile(mhgRefMultVsRun, gRefMult);
  
   
  /// primary tracks
  IntVec trkId;
  trkId.clear();
  int nPosi = 0, nElec = 0;
  float nGoodTrack = 0, nMthGoodTrack = 0;
  int nGlobal = mPicoDst->numberOfTracks();
  for(int i=0; i<nGlobal; i++)
    {
      StPicoTrack* gTrack = mPicoDst->track(i);
      if(!gTrack) continue;
      //StThreeVectorF mom = gTrack->pMom();
      StThreeVectorF mom(gTrack->pMom().x(), gTrack->pMom().y(), gTrack->pMom().z());
      if(mom.mag()<=0) continue;
      float p   = mom.mag();
      float pt  = mom.perp();
      float eta = mom.pseudoRapidity();
      float phi = rotatePhi(mom.phi());
      int charge = gTrack->charge();  
      int nHitsFit  = gTrack->nHitsFit();
      int nHitsDedx = gTrack->nHitsDedx();
      int nHitsMax  = gTrack->nHitsMax();
      float nHitsFrac = nHitsFit*1./nHitsMax;
      //double dca = gTrack->helix().distance(verPos);
      //double dca = (gTrack->dca()-verPos).mag();
      //double dca = (gTrack->dcaPoint()-verPos).mag();
      double dca = gTrack->gDCA(picoEvent->primaryVertex()).Mag();
      fillHisto(mhpTrkNHitsVsPt, pt, nHitsFit);
      fillHisto(mhpTrkNHits, nHitsFit*charge);
      fillHisto(mhpTrkNDedxVsPt, pt, nHitsDedx);
      fillHisto(mhpTrkNDedx, nHitsDedx*charge);
      fillHisto(mhpTrkNHitsFracVsPt, pt, nHitsFrac);
      fillHisto(mhpTrkNHitsFrac, nHitsFrac);
      fillHisto(mhpTrkDcaVsPt, pt, dca);

      if(!isValidTrack(gTrack,mom,verPos)) continue;
      nGoodTrack++;
      float dedx = gTrack->dEdx();
      float nSigmaE  = gTrack->nSigmaElectron();
      float nSigmaPi = gTrack->nSigmaPion();
      float nSigmaK  = gTrack->nSigmaKaon();
      float nSigmaP  = gTrack->nSigmaProton();
      fillHisto(mhpTrkPt, pt);
      fillHisto(mhpTrkEtaVsPt, pt, eta);
      fillHisto(mhpTrkPhiVsPt, pt, phi);
      fillHisto(mhpTrkPhiVsEta, eta, phi);
      fillHisto(mhpTrkDedxVsMom, p, dedx);
      fillHisto(mhNsigmaEVsMom, p, nSigmaE);
      fillHisto(mhNsigmaPiVsMom, p, nSigmaPi);
      fillHisto(mhNsigmaKVsMom, p, nSigmaK);
      fillHisto(mhNsigmaPVsMom, p, nSigmaP);
      fillTProfile(mhpTrkPtVsRun, pt);
      fillTProfile(mhpTrkEtaVsRun, eta);
      fillTProfile(mhpTrkPhiVsRun, phi);
      fillTProfile(mhpTrkDcaVsRun, dca);
      fillTProfile(mhNHitsFitVsRun, nHitsFit);
      fillTProfile(mhNHitsPossVsRun, nHitsMax);
      fillTProfile(mhNHitsDedxVsRun, nHitsDedx);
      fillTProfile(mhDedxVsRun, dedx);
      fillTProfile(mhNsigmaPiVsRun, nSigmaPi);
      fillTProfile(mhNsigmaEVsRun, nSigmaE);
      fillTProfile(mhNsigmaKVsRun, nSigmaK);
      fillTProfile(mhNsigmaPVsRun, nSigmaP);

      /// TOF
      if(gTrack->bTofPidTraitsIndex()>-1)
	{
	  StPicoBTofPidTraits *btofPid = mPicoDst->btofPidTraits(gTrack->bTofPidTraitsIndex());
	  float beta = btofPid->btofBeta();
	  if(beta!=0)
	    {
	      float m2 = pow(p,2) * (1/pow(beta,2)-1);
	      fillHisto(mhBetaVsMom, p, 1./beta);
	      fillHisto(mhM2VsMom, p, m2);
	      fillTProfile(mhBetaVsRun, 1./beta);
	    }
	}

      /// electrons
      int emcIndex = gTrack->bemcPidTraitsIndex();
      //int emcIndex = gTrack->emcPidTraitsIndex();
      if(emcIndex>=0)
	{
	  StPicoBEmcPidTraits *emcPid = mPicoDst->bemcPidTraits(emcIndex);
	  int maxadc = emcPid->bemcAdc0();
	  float energy = emcPid->bemcE();
	  int neta = emcPid->bemcSmdNEta();
	  int nphi = emcPid->bemcSmdNPhi();
	  float zdist = emcPid->bemcZDist();
	  float phidist = emcPid->bemcPhiDist();

	  // StPicoEmcPidTraits *emcPid = mPicoDst->emcPidTraits(emcIndex);
	  // int maxadc = emcPid->adc0();
	  // float energy = emcPid->e();
	  // int neta = emcPid->nEta();
	  // int nphi = emcPid->nPhi();
	  // float zdist = emcPid->zDist();
	  // float phidist = emcPid->phiDist();

	  fillHisto(mhAdc0, maxadc);
	  fillHisto(mhEoverPVsPt, pt, p/energy);
	  fillHisto(mhNetaVsPt, pt, neta);
	  fillHisto(mhNphiVsPt, pt, nphi);
	  fillHisto(mhDistZVsPt, pt, zdist);
	  fillHisto(mhDistPhiVspt, pt, phidist);
	  if(abs(gTrack->nSigmaElectron())<2 &&
	     p/energy>0.3 && p/energy<1.5)
	    {
	      fillHisto(mhElecPt, pt);
	      fillHisto(mhElecPhiVsEta, eta, phi);
	      if(charge>0) nPosi++;
	      else         nElec++;
	    }
	}

      /// muons
      int index = gTrack->mtdPidTraitsIndex();
      if(index>=0)  
	{
	  nMthGoodTrack++;
	  StPicoMtdPidTraits *mtdPid = mPicoDst->mtdPidTraits(index);
	  float dy     = mtdPid->deltaY();
	  float dz     = mtdPid->deltaZ();
	  float dtof   = mtdPid->deltaTimeOfFlight();
	  float gChan   = (mtdPid->backleg()-1)*60 + (mtdPid->module()-1)*12 + mtdPid->cell();
	  fillHisto(mhDeltaZ, dz);
	  fillHisto(mhDzVsPt, pt, dz);
	  fillHisto(mhDzVsgChan, gChan, dz);
	  fillHisto(mhDeltaY, dy);
	  fillHisto(mhDyVsPt, pt, dy);
	  fillHisto(mhDyVsgChan, gChan, dy);
	  fillHisto(mhDeltaTof, dtof);
	  fillHisto(mhDTofVsPt, pt, dtof);
	  fillHisto(mhDTofVsgChan, gChan, dtof);
	  if(isMuonCandidate(gTrack, verPos))
	    {
	      trkId.push_back(i);
	      fillHisto(mhMuonMap, mtdPid->backleg(),(mtdPid->module()-1)*12+mtdPid->cell()+1);
	    }
	}
    }
  fillHisto(mhpTrkN, nGoodTrack);
  fillHisto(mhpTrkNMthMtd, nMthGoodTrack);
  fillHisto(mhNElecPos, nPosi);
  fillHisto(mhNElecNeg, nElec);

  fillTProfile(mhNPositron, nPosi);
  fillTProfile(mhNElectron, nElec);

  /// trigger performance
  triggerData();

  /// MTD hits
  int nMtdHits = mPicoDst->numberOfMtdHits();
  int nTrigHit = 0, nMthHit = 0, nMthTrigHit = 0;
  for(int i=0; i<nMtdHits; i++)
    {
      StPicoMtdHit *hit = mPicoDst->mtdHit(i);
      if(!hit) continue;
      int backleg = hit->backleg();
      int module  = hit->module();
      int cell    = hit->cell();
      int gChan   = (backleg-1)*60 + (module-1)*12 + cell;
      fillHisto(mhMtdHitMap, backleg, (module-1)*12+cell+1);

      fillHisto(mhMtdHitLeTimeDiff, gChan, hit->leadingEdgeTime().second-hit->leadingEdgeTime().first);

      bool isTrig = isMtdHitFiredTrigger(hit);
      bool isMth  = false;
      int index2traits = getMtdPidTraitsIndex(hit);
      if(index2traits>-1)
	{
	  StPicoMtdPidTraits *mtdPid = mPicoDst->mtdPidTraits(index2traits);
	  int index2track = mtdPid->trackIndex();
	  StPicoTrack* gTrack = mPicoDst->track(index2track);
	  if(gTrack)
	    {
	      //StThreeVectorF mom = gTrack->gMom(verPos,picoEvent->bField());
	      StThreeVectorF mom(gTrack->gMom().x(), gTrack->gMom().y(), gTrack->gMom().z());
	      if(isValidTrack(gTrack,mom,verPos)) isMth = true;
	    }
	}

      if(isTrig)
	{
	  nTrigHit ++;
	  fillHisto(mhMtdTrigHitMap, backleg, (module-1)*12+cell+1);
	}

      if(isMth)
	{
	  nMthHit ++;
	  fillHisto(mhMtdMthHitMap, backleg, (module-1)*12+cell+1);
	}

      if(isMth && isTrig)
	{
	  nMthTrigHit ++;
	  fillHisto(mhMtdMthTrigHitMap, backleg, (module-1)*12+cell+1);
	}
    }
  fillHisto(mhMtdNHits, nMtdHits);
  fillHisto(mhMtdNTrigHits, nTrigHit);
  fillHisto(mhMtdNMthHits, nMthHit);
  fillHisto(mhMtdNMthTrigHits, nMthTrigHit);

  fillTProfile(mhNMtdHitsVsRun, nMtdHits);
  fillTProfile(mhNMtdTrigHitsVsRun, nTrigHit);
  fillTProfile(mhNMtdMthHitsVsRun, nMthHit);


  /// muon analysis
  int nPosMuon = 0, nNegMuon = 0;
  int nULpair = 0, nLSpairPos = 0, nLSpairNeg = 0;
  int nJpsi = 0;
  unsigned int nMuon = trkId.size();
  for(unsigned int i=0; i<nMuon; i++)
    {
      StPicoTrack *ipTrack = mPicoDst->track(trkId[i]);
      // StThreeVectorF imom = ipTrack->pMom();
      StThreeVectorF imom(ipTrack->pMom().x(), ipTrack->pMom().y(), ipTrack->pMom().z());
      int iq = ipTrack->charge();
      if(iq>0) nPosMuon++;
      else     nNegMuon++;
      fillHisto(mhMuonPt, imom.perp());
      fillHisto(mhMuonPhiVsEta, imom.pseudoRapidity(),rotatePhi(imom.phi()));

      TLorentzVector imuon;
      imuon.SetXYZM(imom.x(),imom.y(),imom.z(),muMass);
      for(UInt_t j=i+1; j<nMuon; j++)
	{
	  StPicoTrack *jpTrack = mPicoDst->track(trkId[j]);
	  int jq = jpTrack->charge();
	  //StThreeVectorF jmom = jpTrack->pMom();
	  StThreeVectorF jmom(jpTrack->pMom().x(), jpTrack->pMom().y(), jpTrack->pMom().z());
	  TLorentzVector jmuon;
	  jmuon.SetXYZM(jmom.x(),jmom.y(),jmom.z(),muMass);

	  float pt1 = imom.perp(), pt2 = jmom.perp();
	  if(pt1<pt2) 
	    {
	      pt1 = jmom.perp();
	      pt2 = imom.perp();
	    }

	  TLorentzVector jpsi = imuon + jmuon;
	  float invmass = jpsi.M();
	  if(iq*jq<0)
	    {
	      nULpair++;
	      if(pt1>1.5) 
		{
		  fillHisto(mhInvMvsPtUL, invmass, jpsi.Pt());
		  fillHisto(mhInvMUL, invmass);
		}
	      if(invmass>3.0 && invmass<3.2) nJpsi++;
	    }
	  else 
	    {
	      if(invmass>3.0 && invmass<3.2) nJpsi--;
	      if(iq>0)
		{
		  nLSpairPos++;
		  if(pt1>1.5) 
		    {
		      fillHisto(mhInvMvsPtLSpos, invmass, jpsi.Pt());
		      fillHisto(mhInvMLSpos, invmass);
		    }
		}
	      else
		{
		  nLSpairNeg++;
		  if(pt1>1.5) 
		    {
		      fillHisto(mhInvMvsPtLSneg, invmass, jpsi.Pt());
		      fillHisto(mhInvMLSneg, invmass);
		    }
		}
	    }
	}
    }
  fillHisto(mhNMuonPos, nPosMuon); 
  fillHisto(mhNMuonNeg, nNegMuon); 
  fillHisto(mhNULpair, nULpair); 
  fillHisto(mhNLSpairPos, nLSpairPos); 
  fillHisto(mhNLSpairNeg, nLSpairNeg);

  fillTProfile(mhNMuonPosVsRun, nPosMuon);
  fillTProfile(mhNMuonNegVsRun, nNegMuon);
  fillTProfile(mhNMuonPairULVsRun, nULpair);
  fillTProfile(mhNMuonPairLSPosVsRun, nLSpairPos);
  fillTProfile(mhNMuonPairLSNegVsRun, nLSpairNeg);
  fillTProfile(mhNJpsiVsRun, nJpsi);

  return kStOK;
}

//_____________________________________________________________________________
void StMtdProdQAMaker::triggerData()
{
  for(int i = 0; i < kNQTboard; i++)
    {
      for(int j=0; j<2; j++)
	{
	  mTrigQTpos[i][j] = -1;
	}
    }

  // QT cut
  UShort_t mtd_qt_tac_max = 4095;
  UShort_t mtd_qt_tac_min = 100;
  if (mRunId >= 16045067) mtd_qt_tac_min = 80;
  UShort_t mtd_qt_tac_diff_range_abs = 600;
  if (mYear == 2015) mtd_qt_tac_diff_range_abs = 1023;

  if(mDataType==0)
    {
      const Int_t ip = 0;
      StTriggerData *trigData = const_cast<StTriggerData*>(mMuDst->event()->triggerData());

      // VPD tac information
      int vpdTacSum = trigData->vpdEarliestTDCHighThr(east,ip) + trigData->vpdEarliestTDCHighThr(west,ip);

      // QT
      int mxq_tacsum[kNQTboard][2];
      int mxq_tacsum_pos[kNQTboard][2];
      for(int i=0; i<kNQTboard; i++)
	{
	  for(int j=0; j<2; j++)
	    {
	      mxq_tacsum[i][j] = 0;
	      mxq_tacsum_pos[i][j] = -1;
	    }
	}

      // extract tac information for each QT board
      int mtdQTtac[kNQTboard][16];
      int mtdQTadc[kNQTboard][16];
      for(Int_t i=0; i<32; i++)
	{
	  Int_t type = (i/4)%2;
	  if(mYear<=2015)
	    {
	      if(type==1)
		{
		  mtdQTtac[0][i-i/4*2-2] = trigData->mtdAtAddress(i,0);
		  mtdQTtac[1][i-i/4*2-2] = trigData->mtdgemAtAddress(i,0);
		  mtdQTtac[2][i-i/4*2-2] = trigData->mtd3AtAddress(i,0);
		  mtdQTtac[3][i-i/4*2-2] = trigData->mtd4AtAddress(i,0);
		}
	      else
		{
		  mtdQTadc[0][i-i/4*2] = trigData->mtdAtAddress(i,0);
		  mtdQTadc[1][i-i/4*2] = trigData->mtdgemAtAddress(i,0);
		  mtdQTadc[2][i-i/4*2] = trigData->mtd3AtAddress(i,0);
		  mtdQTadc[3][i-i/4*2] = trigData->mtd4AtAddress(i,0);
		}
	    }
	  else
	    {
	      for(int im=0; im<kNQTboard; im++)
		{
		  if(mYear!=2016 && im>=4) continue;
		  if(type==0) mtdQTadc[im][i-i/4*2] = trigData->mtdQtAtCh(im+1,i,0);
		  else        mtdQTtac[im][i-i/4*2-2] = trigData->mtdQtAtCh(im+1,i,0);
		}
	    }

	}

      // In each QT board, find the two signals with
      // largest TacSum values
      int nQtSignal = 0;
      int j[2], a[2];
      int sumTac[kNQTboard][8];
      for(Int_t im=0; im<kNQTboard; im++)
	{
	  for(int i=0; i<8; i++)
	    {
	      if(mYear==2016 && i%2 == 0) continue;
	      if(mYear!=2016 && im>3)     continue;
	      // slewing correction
	      for(int k=0; k<2; k++)
		{
		  j[k] = mtdQTtac[im][i*2+k];
		  a[k] = mtdQTadc[im][i*2+k];

		  int slew_bin = -1;
		  if(a[k]>=0 && a[k]<=mQTSlewBinEdge[im][i*2+k][0]) slew_bin = 0;
		  else
		    {
		      for(int l=1; l<8; l++)
			{
			  if(a[k]>mQTSlewBinEdge[im][i*2+k][l-1] && a[k]<=mQTSlewBinEdge[im][i*2+k][l])
			    {
			      slew_bin = l;
			      break;
			    }
			}
		    }
		  if(slew_bin>=0)
		    j[k] += mQTSlewCorr[im][i*2+k][slew_bin];
		}

	      // position correction
	      int module = mQTtoModule[im][i];
	      sumTac[im][i] = int( j[0] + j[1] + abs(module-3)*1./8 * (j[0]-j[1]) );

	      if(j[0]<mtd_qt_tac_min || j[0]>mtd_qt_tac_max || 
		 j[1]<mtd_qt_tac_min || j[1]>mtd_qt_tac_max ||
		 TMath::Abs(j[0]-j[1])>mtd_qt_tac_diff_range_abs) continue;
	      nQtSignal++;
	  
	      if(mxq_tacsum[im][0] < sumTac[im][i])
		{
		  mxq_tacsum[im][1] = mxq_tacsum[im][0];
		  mxq_tacsum[im][0] = sumTac[im][i];

		  mxq_tacsum_pos[im][1] = mxq_tacsum_pos[im][0];
		  mxq_tacsum_pos[im][0] = i+1;
		}
	      else if (mxq_tacsum[im][1] < sumTac[im][i])
		{
		  mxq_tacsum[im][1]  = sumTac[im][i];
		  mxq_tacsum_pos[im][1] = i+1;
		}
	      int bin = 0;
	      if(mYear!=2016) bin = im*8+i+1;
	      else            bin = im*4+i/2+1;
	      fillHisto(mhMtdVpdTacDiffMT001, bin, (sumTac[im][i]-vpdTacSum)/8+1024);
	    }
	}
      fillHisto(mhNQtSignal, nQtSignal);


      // MT101 informaiton
      int mix_tacsum[kNQTboard][2];
      int mix_id[kNQTboard][2];
      int nMixSignal = 0;
      for(int i = 0; i < kNQTboard; i++)
	{
	  if(mYear!=2016 && i>3)     continue;
	  int idx = 0;
	  if(mYear==2016) idx = i/2*3 + i%2*16;
	  else            idx = i*3;
	  mix_tacsum[i][0] = (trigData->mtdDsmAtCh(idx,0)) + ((trigData->mtdDsmAtCh(idx+1,0)&0x3)<<8);
	  mix_id[i][0]     = (trigData->mtdDsmAtCh(idx+1,0)&0xc)>>2;
	  mix_tacsum[i][1] = (trigData->mtdDsmAtCh(idx+1,0)>>4) + ((trigData->mtdDsmAtCh(idx+2,0)&0x3f)<<4);
	  mix_id[i][1]     = (trigData->mtdDsmAtCh(idx+2,0)&0xc0)>>6;

	  for(int j=0; j<2; j++)
	    {
	      if(mix_tacsum[i][j]>0) 
		{
		  nMixSignal ++;
		  fillHisto(mhMtdTacSumMixvsMxq[i][j],mxq_tacsum[i][j]/8,mix_tacsum[i][j]);
		  int bin = 0;
		  if(mYear!=2016) bin = i*8+mxq_tacsum_pos[i][j];
		  else            bin = i*4+mxq_tacsum_pos[i][j]/2;

		  fillHisto(mhMtdVpdTacDiffMT101,bin,mix_tacsum[i][j]-vpdTacSum/8+1024);
		}
	    }
	}
      fillHisto(mhNMT101Signal, nMixSignal);

      // fill the MtdVpdTacDiff histograms for matched hits
      for(int i=0; i<mMuDst->numberOfMTDHit(); i++)
	{
	  StMuMtdHit *hit = mMuDst->mtdHit(i);
	  if(!hit) continue;
	  if(hit->index2Global()<0) continue;
	  int qt = mModuleToQT[hit->backleg()-1][hit->module()-1];
	  int pos = mModuleToQTPos[hit->backleg()-1][hit->module()-1];
	  if(qt>=1 && qt<=kNQTboard && pos>=1 && pos<=8)
	    {
	      int bin = 0;
	      if(mYear!=2016) bin = (qt-1)*8+pos;
	      else            bin = (qt-1)*4+(pos-1)/2+1;
	      fillHisto(mhMtdVpdMthTacDiffMT001, bin, sumTac[qt-1][pos-1]/8-vpdTacSum/8+1024);

	      if(pos==mxq_tacsum_pos[qt-1][0] || pos==mxq_tacsum_pos[qt-1][1])
		{
		  int index = (pos==mxq_tacsum_pos[qt-1][0]) ? 0 : 1;
		  if(mYear!=2016) bin = (qt-1)*8+mxq_tacsum_pos[qt-1][index];
		  else            bin = (qt-1)*4+mxq_tacsum_pos[qt-1][index]/2;
		  fillHisto(mhMtdVpdMthTacDiffMT101,bin,mix_tacsum[qt-1][index]-vpdTacSum/8+1024);
		}
	    }
	}
  
      // TF201
      UInt_t decision = trigData->dsmTF201Ch(0);
      unsigned short decision2 = 0;
      if(mYear==2016) decision2 = trigData->dsmTF201Ch(6);
      Int_t nMuon = 0;
      for(Int_t i = 0; i < 4; i++)
	{
	  for(Int_t j=0; j<2; j++)
	    {
	      if(mYear==2016)
		{
		  if((decision>>(i*2+j+4))&0x1)
		    {
		      nMuon ++;
		      int qt = i*2;
		      mTrigQTpos[qt][j] = mxq_tacsum_pos[qt][j];
		    }
		  if((decision2>>(i*2+j+4))&0x1)
		    {
		      nMuon ++;
		      int qt = i*2+1;
		      mTrigQTpos[qt][j] = mxq_tacsum_pos[qt][j];
		    }
		}
	      else
		{
		  if((decision>>(i*2+j+4))&0x1)
		    {
		      nMuon ++;
		      mTrigQTpos[i][j] = mxq_tacsum_pos[i][j];
		    }
		}
	    }
	}
      fillHisto(mhNTF201Signal, nMuon);
    }
  else if(mDataType==1)
    {
      StPicoMtdTrigger *mtdTrig = mPicoDst->mtdTrigger(0);
      if(!mtdTrig) return;
      int vpdTacSum = mtdTrig->getVpdTacSum();

      // QT
      int mxq_tacsum[kNQTboard][2];
      int mxq_tacsum_pos[kNQTboard][2];
      for(int i=0; i<kNQTboard; i++)
	{
	  for(int j=0; j<2; j++)
	    {
	      mxq_tacsum[i][j] = 0;
	      mxq_tacsum_pos[i][j] = -1;
	    }
	}

      // In each QT board, find the two signals with
      // largest TacSum values
      int nQtSignal = 0;
      int sumTac[kNQTboard][8];
      for(int im=0; im<kNQTboard; im++)
	{
	  if(mYear!=2016 && im>3)     continue;
	  for(int i=0; i<8; i++)
	    {
	      sumTac[im][i] = mtdTrig->getQTtacSum(im+1,i+1);
	      if(sumTac[im][i]<2*mtd_qt_tac_min) continue;
	      nQtSignal++;

	      if(mxq_tacsum[im][0] < sumTac[im][i])
		{
		  mxq_tacsum[im][1] = mxq_tacsum[im][0];
		  mxq_tacsum[im][0] = sumTac[im][i];

		  mxq_tacsum_pos[im][1] = mxq_tacsum_pos[im][0];
		  mxq_tacsum_pos[im][0] = i+1;
		}
	      else if (mxq_tacsum[im][1] < sumTac[im][i])
		{
		  mxq_tacsum[im][1]  = sumTac[im][i];
		  mxq_tacsum_pos[im][1] = i+1;
		}

	      int bin = 0;
	      if(mYear!=2016) bin = im*8+i+1;
	      else            bin = im*4+i/2+1;
	      fillHisto(mhMtdVpdTacDiffMT001, bin, sumTac[im][i]/8-vpdTacSum/8+1024);
	    }
	}
      fillHisto(mhNQtSignal, nQtSignal);

 
      // MT101 informaiton
      int mix_tacsum[kNQTboard][2];
      int mix_id[kNQTboard][2];
      int nMixSignal = 0;
      for(Int_t i = 0; i < kNQTboard; i++)
	{
	  if(mYear!=2016 && i>3)     continue;
	  for(Int_t j=0; j<2; j++)
	    {
	      mix_tacsum[i][j] = mtdTrig->getMT101Tac(i+1,j);
	      mix_id[i][j]     = mtdTrig->getMT101Id(i+1,j);
	      if(mix_tacsum[i][j]>0) 
		{
		  nMixSignal ++;
		  fillHisto(mhMtdTacSumMixvsMxq[i][j],mxq_tacsum[i][j]/8,mix_tacsum[i][j]);
		  int bin = 0;
		  if(mYear!=2016) bin = i*8+mxq_tacsum_pos[i][j];
		  else            bin = i*4+mxq_tacsum_pos[i][j]/2;
		  fillHisto(mhMtdVpdTacDiffMT101,bin,mix_tacsum[i][j]-vpdTacSum/8+1024);
		}
	    }
	}
      fillHisto(mhNMT101Signal, nMixSignal);

      // fill the MtdVpdTacDiff histograms for matched hits
      for(int i=0; i<mPicoDst->numberOfMtdHits(); i++)
	{
	  StPicoMtdHit *hit = mPicoDst->mtdHit(i);
	  if(!hit) continue;
	  if(getMtdPidTraitsIndex(hit)<0) continue;
	  int qt = mModuleToQT[hit->backleg()-1][hit->module()-1];
	  int pos = mModuleToQTPos[hit->backleg()-1][hit->module()-1];
	  if(qt>=1 && qt<=kNQTboard && pos>=1 && pos<=8)
	    {
	      int bin = 0;
	      if(mYear!=2016) bin = (qt-1)*8+pos;
	      else            bin = (qt-1)*4+(pos-1)/2+1;
	      fillHisto(mhMtdVpdMthTacDiffMT001, bin, sumTac[qt-1][pos-1]/8-vpdTacSum/8+1024);

	      if(pos==mxq_tacsum_pos[qt-1][0] || pos==mxq_tacsum_pos[qt-1][1])
		{
		  int index = (pos==mxq_tacsum_pos[qt-1][0]) ? 0 : 1;
		  if(mYear!=2016) bin = (qt-1)*8+mxq_tacsum_pos[qt-1][index];
		  else            bin = (qt-1)*4+mxq_tacsum_pos[qt-1][index]/2;
		  fillHisto(mhMtdVpdMthTacDiffMT101,bin,mix_tacsum[qt-1][index]-vpdTacSum/8+1024);
		}
	    }
	}
  
      // TF201 decision
      unsigned short decision = mtdTrig->getTF201TriggerBit();
      int nMuon = 0;
      for(int i = 0; i < kNQTboard; i++)
	{
	  for(int j=0; j<2; j++)
	    {
	      if((decision>>(i*2+j))&0x1)
		{
		  nMuon ++;
		  mTrigQTpos[i][j] = mxq_tacsum_pos[i][j];
		}
	    }
	}
      fillHisto(mhNTF201Signal, nMuon);
      /*
      if(nMuon==0)
	{
	  printf("\n+++ New Event +++ \n");
	  std::vector<unsigned int> trigs = mPicoDst->event()->triggerIds();
	  for(int i=0; i<trigs.size(); i++)
	    {
	      printf("[i] Trigger %d: %d\n",i,trigs[i]);
	    }

	  cout << "TF201 = " << std::bitset<16>(decision)  << endl;
	}
      */
    }
}

//_____________________________________________________________________________
bool StMtdProdQAMaker::getBemcInfo(StMuTrack *track, int &maxadc, float &energy, float &zdist, float &phidist, int &neta, int &nphi)
{
  maxadc = 0;
  energy = 0.;
  zdist = -999.;
  phidist = -999.;
  neta = 0;
  nphi = 0;

  Float_t maxtowerE = -999.;
  Float_t mindist = 999.;
  Int_t mod = -1, eta=-1, sub=-1;

  StEmcCollection *mEmcCollection = (StEmcCollection*)mMuDst->emcCollection();
  if(!mEmcCollection) 
    {
      //LOG_WARN << " No Emc Collection for this event " << endm;
      return false;
    }

  StThreeVectorD position, momentum;
  StThreeVectorD positionBSMDE, momentumBSMDE;
  StThreeVectorD positionBSMDP, momentumBSMDP;

  Double_t bField = mMuDst->event()->runInfo().magneticField()/10.; //Tesla
  Bool_t ok       = kFALSE;
  Bool_t okBSMDE  = kFALSE;
  Bool_t okBSMDP  = kFALSE;
  if(mEmcPosition) 
    {
      ok = mEmcPosition->projTrack(&position, &momentum, track, bField, mEmcGeom[0]->Radius());
      okBSMDE = mEmcPosition->projTrack(&positionBSMDE, &momentumBSMDE, track, bField, mEmcGeom[2]->Radius());
      okBSMDP = mEmcPosition->projTrack(&positionBSMDP, &momentumBSMDP, track, bField, mEmcGeom[3]->Radius());
    }
  if(!ok) {
    LOG_WARN << " Projection failed for this track ... " << endm;
    return false;
  }

  Bool_t bemcMatchFlag = kFALSE;
  if(ok && okBSMDE && okBSMDP)
    {
      StSPtrVecEmcPoint& bEmcPoints = mEmcCollection->barrelPoints();
      mindist=1.e9;
      mEmcGeom[0]->getBin(positionBSMDP.phi(), positionBSMDE.pseudoRapidity(), mod, eta, sub); //project on SMD plan
      for(StSPtrVecEmcPointIterator it = bEmcPoints.begin(); it != bEmcPoints.end(); it++) 
	{
	  Bool_t associatedPoint = kFALSE;
	  Bool_t associatedCluster = kFALSE;
	  StPtrVecEmcCluster& bEmcClusters = (*it)->cluster(kBarrelEmcTowerId);
	  if(bEmcClusters.size()==0 ) continue;
	  if(bEmcClusters[0]==NULL) continue;
	  for(StPtrVecEmcClusterIterator cIter = bEmcClusters.begin(); cIter != bEmcClusters.end(); cIter++)
	    {
	      StPtrVecEmcRawHit& bEmcHits = (*cIter)->hit();
	      for(StPtrVecEmcRawHitIterator hIter = bEmcHits.begin(); hIter != bEmcHits.end(); hIter++) 
		{
		  if(mod == (Int_t)(*hIter)->module() && eta == (Int_t)(*hIter)->eta() && sub == (Int_t)(*hIter)->sub())
		    {
		      bemcMatchFlag = kTRUE;
		      associatedPoint = kTRUE;
		      associatedCluster = kTRUE;
		      break;
		    }
		}
	      if(associatedCluster)
		{
		  for(StPtrVecEmcRawHitIterator hitit=bEmcHits.begin(); hitit!=bEmcHits.end();hitit++)
		    {
		      if((*hitit)->energy()>maxtowerE) maxtowerE = (*hitit)->energy();
		      if((int)((*hitit)->adc())>maxadc) maxadc = (*hitit)->adc();
		    }
		}
	      associatedCluster = kFALSE;
	    }
	  StPtrVecEmcCluster& smdeClusters = (*it)->cluster(kBarrelSmdEtaStripId);
	  StPtrVecEmcCluster& smdpClusters = (*it)->cluster(kBarrelSmdPhiStripId);

	  if(associatedPoint)
	    {
	      energy += (*it)->energy(); //use point's energy, not tower cluster's energy

	      float deltaphi=(*it)->position().phi()-positionBSMDP.phi();
	      if(deltaphi>=TMath::Pi()) deltaphi=deltaphi-TMath::TwoPi();
	      if(deltaphi<-TMath::Pi()) deltaphi=deltaphi+TMath::TwoPi();
	      
	      float rsmdp=mEmcGeom[3]->Radius();
	      float pointz=(*it)->position().z();
	      float deltaz=pointz-positionBSMDE.z();
	      if(sqrt(deltaphi*deltaphi*rsmdp*rsmdp+deltaz*deltaz)<mindist) 
		{
		  phidist=deltaphi;
		  zdist  =deltaz;
		  if(smdeClusters.size()>=1) neta=smdeClusters[0]->nHits();
		  if(smdpClusters.size()>=1) nphi=smdpClusters[0]->nHits();
		  mindist=sqrt(deltaphi*deltaphi*rsmdp*rsmdp+deltaz*deltaz);
		}
	    }//associated
	}
    } // end if (ok && okBSMDE && okBSMDP)
  return bemcMatchFlag;
}

//_____________________________________________________________________________
void StMtdProdQAMaker::addCutToHisto(TH1 *h, const Int_t bin, const char *label, const Float_t value)
{
  if(!h) return;
  h->GetXaxis()->SetBinLabel(bin,label);
  if(value!=-9999999)
    h->SetBinContent(bin,value);
}

//_____________________________________________________________________________
void StMtdProdQAMaker::bookHistos()
{
  //this array describe the Channel input from which backleg & position & direction 
  const char qtlabel2[32][100] = {"QT1-1  25-1","QT1-2  25-5","QT1-3  25-2","QT1-4  25-4","QT1-5  25-3","QT1-6  30-3","QT1-7  30-1","QT1-8  30-5","QT2-1  05-1","QT2-2  05-5","QT2-3  05-2","QT2-4  05-4","QT2-5  05-3","QT2-6      ","QT2-7  30-2","QT2-8  30-4","QT3-1  10-1","QT3-2  10-5","QT3-3  10-2","QT3-4  10-4","QT3-5  10-3","QT3-6  15-3","QT3-7      ","QT3-8      ","QT4-1  21-1","QT4-2  21-5","QT4-3  20-2","QT4-4  20-4","QT4-5  20-3","QT4-6      ","QT4-7  15-2","QT4-8  15-4"};

  const Int_t nSpecMBins = 56;//PRL mass bin
  Double_t specM[nSpecMBins+1] = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 
				  0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.125, 0.175, 0.2, 0.31, 0.4, 0.51, 
				  0.63, 0.75, 0.78, 0.785, 0.79, 0.8, 0.89, 0.965, 1, 1.01, 1.015, 1.02, 
				  1.035, 1.13, 1.25, 1.45, 1.65, 1.875, 2.075, 2.25, 2.475, 2.665, 2.85, 
				  2.99, 3.02, 3.035, 3.055, 3.07, 3.075, 3.09, 3.095, 3.1, 3.115, 3.13, 
				  3.225, 3.4, 3.85, 4.4, 5.5};

  if(mOutputFileName.Length()==0)
    {
      LOG_WARN << "Input file name is needed for output. Set it to qa.root!" << endm;
      mOutputFileName = "qa.root";
    }
  mOutputFile = new TFile(mOutputFileName.Data(),"recreate");

  mhEventCuts = new TH1F("hEventCuts","Cuts used for analysis",20,0,20);
  addCutToHisto(mhEventCuts, 1,  "|vtx_z|",            mMaxTpcVz);
  addCutToHisto(mhEventCuts, 2,  "dz",                 mMaxDiffVz);
  addCutToHisto(mhEventCuts, 3,  "trk_pt_min",         mMinTrkPt);
  addCutToHisto(mhEventCuts, 4,  "trk_pt_max",         mMaxTrkPt);
  addCutToHisto(mhEventCuts, 5,  "trk_eta",            mMaxTrkEta);
  addCutToHisto(mhEventCuts, 6,  "MinNHitsFit",        mMinNHitsFit);
  addCutToHisto(mhEventCuts, 7,  "MinNHitsDedx",       mMinNHitsDedx);
  addCutToHisto(mhEventCuts, 8,  "MinNHitsFrac",       mMinNHitsFrac);
  addCutToHisto(mhEventCuts, 9,  "mMaxDca",            mMaxDca);
  addCutToHisto(mhEventCuts, 10,  "mMinNsigmaPi",      mMinNsigmaPi);
  addCutToHisto(mhEventCuts, 11,  "mMaxNsigmaPi",      mMaxNsigmaPi);
  addCutToHisto(mhEventCuts, 12, "mMinMuonDeltaZ",     mMinMuonDeltaZ);
  addCutToHisto(mhEventCuts, 13, "mMaxMuonDeltaZ",     mMaxMuonDeltaZ);
  addCutToHisto(mhEventCuts, 14, "mMinMuonDeltaY",     mMinMuonDeltaY);
  addCutToHisto(mhEventCuts, 15, "mMaxMuonDeltaY",     mMaxMuonDeltaY);
  addCutToHisto(mhEventCuts, 16, "mMinMuonDeltaTof",   mMinMuonDeltaTof);
  addCutToHisto(mhEventCuts, 17, "mMaxMuonDeltaTof",   mMaxMuonDeltaTof);
  addCutToHisto(mhEventCuts, 18, "mMaxMuonDca",        mMaxMuonDca);
  addCutToHisto(mhEventCuts, 19, "mDataType",          mDataType);


  const char *trig_name[kNtrig] = {"dimuon","simuon","emu","dimuon_hft"};
  const int nPtBins = 40;
  const double lowPtBin = 0, hiPtBin = 20;
  TString name;
  for(int i=0; i<kNtrig; i++)
    {
      name = trig_name[i];
      if(!mSeparateTrig)
	{
	  if(i==0) name = "all";
	  else continue;
	}
      
      mhEventStat[i] = new TH1F(Form("mhEventStat_%s",name.Data()),"Event statistics",5,0,5);
      mhEventStat[i]->GetXaxis()->SetBinLabel(1,"All events");
      mhEventStat[i]->GetXaxis()->SetBinLabel(2,"Good trigger");
      mhEventStat[i]->GetXaxis()->SetBinLabel(3,"VPD");
      mhEventStat[i]->GetXaxis()->SetBinLabel(4,"Vtx cuts");

      /// vertex 

      mhVertexXYRanking0[i] = new TH2F(Form("mhVertexXYRanking0_%s",name.Data()),"Default TPC v_{y} vs v_{x} (no vertex cut);TPC v_{x} (cm);TPC v_{y} (cm)",100,-5,5,100,-5,5);

      mhTpcVzRanking0[i] = new TH1F(Form("mhTpcVzRanking0_%s",name.Data()),"Default TPC v_{z} distribution (no vertex cut); TPC v_{z} (cm)",201,-201,201);

      mhDiffVzRanking0[i] = new TH1F(Form("mhDiffVzRanking0_%s",name.Data()),"TPC-VPD v_{z} for Default TPC vertex (no vertex cut); TPC-VPD v_{z} (cm)",400,-100,100);

      mhTpcVtxIndex[i] = new TH1F(Form("mhTpcVtxIndex_%s",name.Data()),"Index of selected TPC vertex;Vertex index",50,0,50);

      mhVtxIndexVsClosest[i] = new TH2F(Form("mhVtxIndexVsClosest_%s",name.Data()),"Index TPC vertex: chosen vs closest;Closest index to VPD;Chosen index",50,0,50,50,0,50);

      mhVertexXY[i] = new TH2F(Form("mhVertexXY_%s",name.Data()),"Selected TPC v_{y} vs v_{x} (no vertex cut);TPC v_{x} (cm);TPC v_{y} (cm)",100,-5,5,100,-5,5);

      mhTpcVz[i] = new TH1F(Form("mhTpcVz_%s",name.Data()),"Selected TPC v_{z} distribution (no vertex cut); TPC v_{z} (cm)",201,-201,201);

      mhVpdVz[i] = new TH1F(Form("mhVpdVz_%s",name.Data()),"VPD v_{z} distribution (no vertex cut); VPD v_{z} (cm)",201,-201,201);


      mhDiffVz[i] = new TH1F(Form("mhDiffVz_%s",name.Data()),"Selected TPC-VPD v_{z} (no vertex cut); TPC-VPD v_{z} (cm)",200,-10,10);

      mhDiffVzVsTpcVz[i] = new TH2F(Form("mhDiffVzVsTpcVz_%s",name.Data()),"Selected TPC-VPD v_{z} vs TPC v_{z} (no vertex cut); TPC v_{z} (cm); TPC-VPD v_{z} (cm)",201,-201,201,200,-10,10);
      
      mhDiffVzVsVpdVz[i] = new TH2F(Form("mhDiffVzVsVpdVz_%s",name.Data()),"Selected TPC-VPD v_{z} vs VPD v_{z} (no vertex cut); VPD v_{z} (cm); TPC-VPD v_{z} (cm)",201,-201,201,200,-10,10);


      /// reference multiplicity
      mhRefMult[i] = new TH1F(Form("mhRefMult_%s",name.Data()),"Reference multiplicity;RefMult",800,0,800);

      mhgRefMult[i] = new TH1F(Form("mhgRefMult_%s",name.Data()),"Global reference multiplicity;gRefMult",800,0,800);

      mhgRefMultVsRefMult[i] = new TH2F(Form("mhgRefMultVsRefMult_%s",name.Data()),"Global reference multiplicity vs reference multiplicity;RefMult;gRefMult",800,0,800,800,0,800);

      mhTpcVzVsgRef[i] = new TH2F(Form("mhTpcVzVsgRef_%s",name.Data()),"TPC v_{z} vs gRefMult;gRefMult;TPC v_{z} (cm)",800,0,800,201,-201,201);

      mhDiffVzVsgRef[i] = new TH2F(Form("mhDiffVzVsgRef_%s",name.Data()),"TPC-VPD v_{z} vs gRefMult;gRefMult;TPC-VPD v_{z} (cm)",800,0,800,200,-10,10);

      mhZdcRateVsgRef[i] = new TH2F(Form("mhZdcRateVsgRef_%s",name.Data()),"ZDC rate vs gRefMult;gRefMult;ZDC (kHz)",800,0,800,500,0,5e3);

      mhBbcRateVsgRef[i] = new TH2F(Form("mhBbcRateVsgRef_%s",name.Data()),"BBC rate vs gRefMult;gRefMult;BBC (kHz)",800,0,800,500,0,5e3);

      /// primary tracks
      mhpTrkNHitsVsPt[i] = new TH2F(Form("mhpTrkNHitsVsPt_%s",name.Data()),"Primary track: NHitsFit vs p_{T} (no track quality cuts);p_{T} (GeV/c);NHitsFit",nPtBins,lowPtBin,hiPtBin,45,0,45);

      mhpTrkNHits[i] = new TH1F(Form("mhpTrkNHits_%s",name.Data()),"Primary track: NHitsFit distribution (no track quality cuts);q*NHitsFit",90,-45,45);

      mhpTrkNDedxVsPt[i] = new TH2F(Form("mhpTrkNDedxVsPt_%s",name.Data()),"Primary track: NHitsDedx vs p_{T} (no track quality cuts);p_{T} (GeV/c);NHitsDedx",nPtBins,lowPtBin,hiPtBin,45,0,45);

      mhpTrkNDedx[i] = new TH1F(Form("mhpTrkNDedx_%s",name.Data()),"Primary track: NHitsDedx distribution (no track quality cuts);q*NHitsDedx",90,-45,45);

      mhpTrkNHitsFracVsPt[i] = new TH2F(Form("mhpTrkNHitsFracVsPt_%s",name.Data()),"Primary track: NHitsFit/NHitsPoss vs p_{T} (no track quality cuts);p_{T} (GeV/c);NHitsFrac",nPtBins,lowPtBin,hiPtBin,100,0,1);

      mhpTrkNHitsFrac[i] = new TH1F(Form("mhpTrkNHitsFrac_%s",name.Data()),"Primary track: NHitsFit/NHitsPoss distribution (no track quality cuts);NHitsFit/NHitsPoss",100,0,1);

      mhpTrkDcaVsPt[i] = new TH2F(Form("mhpTrkDcaVsPt_%s",name.Data()),"Primary track: DCA vs p_{T} (no track quality cuts);p_{T} (GeV/c);DCA (cm)",nPtBins,lowPtBin,hiPtBin,90,0,30);

      mhpTrkN[i] = new TH1F(Form("mhpTrkN_%s",name.Data()),"# of primary tracks per event;N_{trk}",800,0,800);

      mhpTrkNMthMtd[i] = new TH1F(Form("mhpTrkNMthMtd_%s",name.Data()),"# of primary tracks matched MTD hits per event;N_{trk}",15,0,15);

      mhpTrkPt[i] = new TH1F(Form("mhpTrkPt_%s",name.Data()),"p_{T} distribution of primary tracks;p_{T} (GeV/c)",nPtBins,lowPtBin,hiPtBin);

      mhpTrkEtaVsPt[i] = new TH2F(Form("mhpTrkEtaVsPt_%s",name.Data()),"#eta vs p_{T} of primary tracks;p_{T} (GeV/c);#eta",nPtBins,lowPtBin,hiPtBin,120,-1.2,1.2);

      mhpTrkPhiVsPt[i] = new TH2F(Form("mhpTrkPhiVsPt_%s",name.Data()),"#varphi vs p_{T} of primary tracks;p_{T} (GeV/c);#varphi",nPtBins,lowPtBin,hiPtBin,120,0,2*pi);

      mhpTrkPhiVsEta[i] = new TH2F(Form("mhpTrkPhiVsEta_%s",name.Data()),"#varphi vs #eta of primary tracks;#eta;#varphi",120,-1.2,1.2,120,0,2*pi);

      mhpTrkDedxVsMom[i] = new TH2F(Form("mhpTrkDedxVsMom_%s",name.Data()),"dE/dx vs momentum of primary tracks;p (GeV/c);dE/dx (keV/cm)",nPtBins,lowPtBin,hiPtBin,100,0,10);

      mhNsigmaEVsMom[i] = new TH2F(Form("mhNsigmaEVsMom_%s",name.Data()),"n#sigma_{e} vs momentum of primary tracks;p (GeV/c);n#sigma_{e}",nPtBins,lowPtBin,hiPtBin,100,-10,10);

      mhNsigmaPiVsMom[i] = new TH2F(Form("mhNsigmaPiVsMom_%s",name.Data()),"n#sigma_{#pi} vs momentum of primary tracks;p (GeV/c);n#sigma_{#pi}",nPtBins,lowPtBin,hiPtBin,100,-10,10);

      mhNsigmaKVsMom[i] = new TH2F(Form("mhNsigmaKVsMom_%s",name.Data()),"n#sigma_{K} vs momentum of primary tracks;p (GeV/c);n#sigma_{K}",nPtBins,lowPtBin,hiPtBin,100,-10,10);

      mhNsigmaPVsMom[i] = new TH2F(Form("mhNsigmaPVsMom_%s",name.Data()),"n#sigma_{P} vs momentum of primary tracks;p (GeV/c);n#sigma_{P}",nPtBins,lowPtBin,hiPtBin,100,-10,10);

      mhBetaVsMom[i] = new TH2F(Form("mhBetaVsMom_%s",name.Data()),"1/#beta vs momentum of primary tracks;p (GeV/c);1/#beta",nPtBins,lowPtBin,hiPtBin,100,0,5);

      mhM2VsMom[i] = new TH2F(Form("mhM2VsMom_%s",name.Data()),"m^{2} vs momentum of primary tracks;p (GeV/c);m^{2} (GeV/c^{2})^{2}",nPtBins,lowPtBin,hiPtBin,100,0,2);

      /// electron PID
      mhAdc0[i] = new TH1F(Form("mhAdc0_%s",name.Data()),"adc0 distribution for BEMC-matched primary tracks;adc0",500,0,500);

      mhEoverPVsPt[i] = new TH2F(Form("mhEoverPVsPt_%s",name.Data()),"P/E vs p_{T} for BEMC-matched primary tracks;p_{T} (GeV/c); P/E",nPtBins,lowPtBin,hiPtBin,20,0,2);

      mhNetaVsPt[i] = new TH2F(Form("mhNetaVsPt_%s",name.Data()),"nEta vs p_{T} for BEMC-matched primary tracks;p_{T} (GeV/c); nEta",nPtBins,lowPtBin,hiPtBin,5,0,5);

      mhNphiVsPt[i] = new TH2F(Form("mhNphiVsPt_%s",name.Data()),"nPhi vs p_{T} for BEMC-matched primary tracks;p_{T} (GeV/c); nPhi",nPtBins,lowPtBin,hiPtBin,5,0,5);

      mhDistZVsPt[i] = new TH2F(Form("mhDistZVsPt_%s",name.Data()),"#Deltaz vs p_{T} for BEMC-matched primary tracks;p_{T} (GeV/c); #Deltaz (cm)",nPtBins,lowPtBin,hiPtBin,50,-50,50);

      mhDistPhiVspt[i] = new TH2F(Form("mhDistPhiVspt_%s",name.Data()),"#Delta#varphi vs p_{T} for BEMC-matched primary tracks;p_{T} (GeV/c); #Delta#varphi",nPtBins,lowPtBin,hiPtBin,40,-0.2,0.2);

      mhNElecPos[i] = new TH1F(Form("mhNElecPos_%s",name.Data()),"Number of positron candidates per event;N",5,0,5);

      mhNElecNeg[i] = new TH1F(Form("mhNElecNeg_%s",name.Data()),"Number of electron candidates per event;N",5,0,5);

      mhElecPt[i] = new TH1F(Form("mhElecPt_%s",name.Data()),"p_{T} distribution of e^{+}/e^{-} candidates;p_{T} (GeV/c)",nPtBins,lowPtBin,hiPtBin);

      mhElecPhiVsEta[i] = new TH2F(Form("mhElecPhiVsEta_%s",name.Data()),"#varphi vs #eta of e^{+}/e^{-} candidates;#eta;#varphi",120,-1.2,1.2,120,0,2*pi);


      /// trigger performance
      mhNQtSignal[i] = new TH1F(Form("mhNQtSignal_%s",name.Data()),"Number of QT signals per event;N",10,0,10);

      mhNMT101Signal[i] = new TH1F(Form("mhNMT101Signal_%s",name.Data()),"Number of MT101 signals per event;N",10,0,10);

      mhNTF201Signal[i] = new TH1F(Form("mhNTF201Signal_%s",name.Data()),"Number of TF201 signals per event;N",10,0,10);

      mhMtdVpdTacDiffMT001[i] = new TH2F(Form("mhMtdVpdTacDiffMT001_%s",name.Data()),"QT: MTD-VPD tac difference (all);;tac_{MTD}-tac_{VPD}",32,0.5,32.5,1000,500,1500);

      mhMtdVpdMthTacDiffMT001[i] = new TH2F(Form("mhMtdVpdMthTacDiffMT001_%s",name.Data()),"QT: MTD-VPD tac difference (track-matched);;tac_{MTD}-tac_{VPD}",32,0.5,32.5,1000,500,1500);

      mhMtdVpdTacDiffMT101[i] = new TH2F(Form("mhMtdVpdTacDiffMT101_%s",name.Data()),"MT101: MTD-VPD tac difference (all);;tac_{MTD}-tac_{VPD}",32,0.5,32.5,1000,500,1500);

      mhMtdVpdMthTacDiffMT101[i] = new TH2F(Form("mhMtdVpdMthTacDiffMT101_%s",name.Data()),"MT101: MTD-VPD tac difference (track-matched);;tac_{MTD}-tac_{VPD}",32,0.5,32.5,1000,500,1500);
      for(Int_t j=0; j<32; j++)
	{
	  mhMtdVpdTacDiffMT001[i]->GetXaxis()->SetBinLabel(j+1,qtlabel2[j]);
	  mhMtdVpdMthTacDiffMT001[i]->GetXaxis()->SetBinLabel(j+1,qtlabel2[j]);
	  mhMtdVpdTacDiffMT101[i]->GetXaxis()->SetBinLabel(j+1,qtlabel2[j]);
	  mhMtdVpdMthTacDiffMT101[i]->GetXaxis()->SetBinLabel(j+1,qtlabel2[j]);
	}

      for(Int_t j=0; j<kNQTboard; j++)
	{
	  for(Int_t k=0; k<2; k++)
	    {
	      mhMtdTacSumMixvsMxq[j][k][i] = new TH2F(Form("mhMtdTacSumMixvsMxq_QT%d_%d_%s",j+1,k,name.Data()),Form("MTD QT%d: MIX vs MXQ at %d;mxq_mtdtacsum;mix_mtdtacsum",j+1,k),1024,0,1024,1024,0,1024);
	    }
	}

      /// MTD hits
      mhMtdNRawHits[i] = new TH1F(Form("mhMtdNRawHits_%s",name.Data()),"Number of raw MTD hits per event;N",150,0,150);

      mhMtdRawHitMap[i] = new TH2F(Form("mhMtdRawHitMap_%s",name.Data()),"Channel vs backleg of MTD raw hits;backleg;channel",30,0.5,30.5,120,0.5,120.5);

      mhMtdNHits[i] = new TH1F(Form("mhMtdNHits_%s",name.Data()),"Number of MTD hits per event;N",30,0,30);

      mhMtdHitMap[i] = new TH2F(Form("mhMtdHitMap_%s",name.Data()),"Channel vs backleg of MTD hits;backleg;channel",30,0.5,30.5,60,0.5,60.5);

      mhMtdHitTrigTime[i]= new TH2F(Form("mhMtdHitTrigTime_%s",name.Data()),"MTD: trigger time of hit (west+east)/2;channel;tdc-t_{trigger} (ns)",1801,-0.5,1800.5,450,2400,3300);

      mhMtdHitLeTimeDiff[i]= new TH2F(Form("mhMtdHitLeTimeDiff_%s",name.Data()),"MTD: (east-west) leading time of hits;channel;#Deltat_{leading} (ns)",1801,-0.5,1800.5,41,-20.5,20.5);

      mhMtdNTrigHits[i] = new TH1F(Form("mhMtdNTrigHits_%s",name.Data()),"Number of triggering MTD hits per event;N",10,0,10);

      mhMtdTrigHitMap[i] = new TH2F(Form("mhMtdTrigHitMap_%s",name.Data()),"Channel vs backleg of triggering MTD hits;backleg;channel",30,0.5,30.5,60,0.5,60.5);

      mhMtdNMthHits[i] = new TH1F(Form("mhMtdNMthHits_%s",name.Data()),"Number of matched MTD hits per event;N",10,0,10);

      mhMtdMthHitMap[i] = new TH2F(Form("mhMtdMthHitMap_%s",name.Data()),"Channel vs backleg of matched MTD hits;backleg;channel",30,0.5,30.5,60,0.5,60.5);

      mhMtdNMthTrigHits[i] = new TH1F(Form("mhMtdNMthTrigHits_%s",name.Data()),"Number of matched triggering MTD hits per event;N",10,0,10);

      mhMtdMthTrigHitMap[i] = new TH2F(Form("mhMtdMthTrigHitMap_%s",name.Data()),"Channel vs backleg of matched triggering MTD hits;backleg;channel",30,0.5,30.5,60,0.5,60.5);

      
      /// Matching information
      mhLocalYVsgChan[i] = new TH2F(Form("mhLocalYVsgChan_%s",name.Data()),"Local y vs channel for matched track-hit pairs;channel;local y (cm)",1801,-0.5,1800.5,160,-80,80);

      mhLocalZVsgChan[i] = new TH2F(Form("mhLocalZVsgChan_%s",name.Data()),"Local z vs channel for matched track-hit pairs;channel;local z (cm)",1801,-0.5,1800.5,200,-100,100);

      mhMtdTofVsgChan[i] = new TH2F(Form("mhMtdTofVsgChan_%s",name.Data()),"MTD measured tof vs channel for matched track-hit pairs;channel;tof_{measured} (ns)",1801,-0.5,1800.5,300,0,30);

      mhExpTofVsgChan[i] = new TH2F(Form("mhExpTofVsgChan_%s",name.Data()),"Expected tof vs channel for matched track-hit pairs;channel;tof_{expected} (ns)",1801,-0.5,1800.5,300,0,30);

      mhDeltaZ[i] = new TH1F(Form("mhDeltaZ_%s",name.Data()),"#Deltaz distribution for matched track-hit pairs;#Deltaz (cm)",201,-201,201);

      mhDzVsPt[i] = new TH2F(Form("mhDzVsPt_%s",name.Data()),"#Deltaz vs p_{T} for matched track-hit pairs;p_{T} (GeV/c);#Deltaz (cm)",nPtBins,lowPtBin,hiPtBin,201,-201,201);

      mhDzVsgChan[i] = new TH2F(Form("mhDzVsgChan_%s",name.Data()),"#Deltaz vs channel for matched track-hit pairs;channel;#Deltaz (cm)",1801,-0.5,1800.5,201,-201,201);

      mhDeltaY[i] = new TH1F(Form("mhDeltaY_%s",name.Data()),"#Deltay distribution for matched track-hit pairs;#Deltay (cm)",101,-101,101);

      mhDyVsPt[i] = new TH2F(Form("mhDyVsPt_%s",name.Data()),"#Deltay vs p_{T} for matched track-hit pairs;p_{T} (GeV/c);#Deltay (cm)",nPtBins,lowPtBin,hiPtBin,101,-101,101);

      mhDyVsgChan[i] = new TH2F(Form("mhDyVsgChan_%s",name.Data()),"#Deltay vs channel for matched track-hit pairs;channel;#Deltay (cm)",1801,-0.5,1800.5,101,-101,101);

      mhDeltaTof[i] = new TH1F(Form("mhDeltaTof_%s",name.Data()),"#Deltatof distribution for matched track-hit pairs;#Deltatof (ns)",300,-10,20);

      mhDTofVsPt[i] = new TH2F(Form("mhDTofVsPt_%s",name.Data()),"#Deltatof vs p_{T} for matched track-hit pairs;p_{T} (GeV/c);#Deltatof (ns)",nPtBins,lowPtBin,hiPtBin,300,-10,20);

      mhDTofVsgChan[i] = new TH2F(Form("mhDTofVsgChan_%s",name.Data()),"#Deltatof vs channel for matched track-hit pairs;channel;#Deltatof (ns)",1801,-0.5,1800.5,300,-10,20);


      /// Muon analysis
      mhNMuonPos[i] = new TH1F(Form("mhNMuonPos_%s",name.Data()),"Number of positive muon candidates per event;N",5,-0.5,4.5);

      mhNMuonNeg[i] = new TH1F(Form("mhNMuonNeg_%s",name.Data()),"Number of negative muon candidates per event;N",5,-0.5,4.5);

      mhMuonPt[i] = new TH1F(Form("mhMuonPt_%s",name.Data()),"p_{T} distribution of muon candidates;p_{T} (GeV/c)",nPtBins,lowPtBin,hiPtBin);

      mhMuonPhiVsEta[i] = new TH2F(Form("mhMuonPhiVsEta_%s",name.Data()),"#varphi vs #eta of muon candidates;#eta;#varphi",120,-1.2,1.2,120,0,2*pi);

      mhMuonMap[i] = new TH2F(Form("mhMuonMap_%s",name.Data()),"Channel vs backleg of muon candidates;backleg;channel",30,0.5,30.5,60,0.5,60.5);

      mhNULpair[i] = new TH1F(Form("mhNULpair_%s",name.Data()),"Number of unlike-sign pairs per event;N",5,-0.5,4.5);

      mhNLSpairPos[i] = new TH1F(Form("mhNLSpairPos_%s",name.Data()),"Number of positive like-sign pairs per event;N",5,-0.5,4.5);
	  
      mhNLSpairNeg[i] = new TH1F(Form("mhNLSpairNeg_%s",name.Data()),"Number of negative like-sign pairs per event;N",5,-0.5,4.5);
	  
      mhInvMvsPtUL[i] = new TH2F(Form("mhInvMvsPtUL_%s",name.Data())," p_{T} vs M_{#mu#mu} for unlike-sign pairs;M_{#mu#mu} (GeV/c^2);p_{T} (GeV/c)",1500,0,15,10,0,10);
	  
      mhInvMvsPtLSpos[i] = new TH2F(Form("mhInvMvsPtLSpos_%s",name.Data())," p_{T} vs M_{#mu#mu} for positive like-sign pairs;M_{#mu#mu} (GeV/c^2);p_{T} (GeV/c)",1500,0,15,10,0,10);
	  
      mhInvMvsPtLSneg[i] = new TH2F(Form("mhInvMvsPtLSneg_%s",name.Data())," p_{T} vs M_{#mu#mu} for negative like-sign pairs;M_{#mu#mu} (GeV/c^2);p_{T} (GeV/c)",1500,0,15,10,0,10);

      mhInvMUL[i] = new TH1F(Form("mhInvMUL_%s",name.Data()),"M_{#mu#mu} for unlike-sign pairs;M_{#mu#mu} (GeV/c^2)",nSpecMBins,specM);
	  
      mhInvMLSpos[i] = new TH1F(Form("mhInvMLSpos_%s",name.Data()),"M_{#mu#mu} for positive like-sign pairs;M_{#mu#mu} (GeV/c^2)",nSpecMBins,specM);
	  
      mhInvMLSneg[i] = new TH1F(Form("mhInvMLSneg_%s",name.Data()),"M_{#mu#mu} for negative like-sign pairs;M_{#mu#mu} (GeV/c^2)",nSpecMBins,specM);

      int nRun = mLastRun-mFirstRun+1;
      double run_lowBin = mFirstRun - 0.5;
      double run_highBin = mLastRun + 0.5;

      mhRunStat[i] = new TH1F(Form("mhRunStat_%s",name.Data()),"Number of events per Run; run index",nRun, run_lowBin, run_highBin);

      mhBBCrateVsRun[i] = new TProfile(Form("mhBBCrateVsRun_%s",name.Data()),"BBC rate vs run; run index; BBC rate (kHz)",nRun, run_lowBin, run_highBin);

      mhZDCrateVsRun[i] = new TProfile(Form("mhZDCrateVsRun_%s",name.Data()),"ZDC rate vs run; run index; ZDC rate (kHz)",nRun, run_lowBin, run_highBin);

      mhRefMultVsRun[i] = new TProfile(Form("mhRefMultVsRun_%s",name.Data()),"Reference multiplicity vs run; run index; RefMult",nRun, run_lowBin, run_highBin);

      mhgRefMultVsRun[i] = new TProfile(Form("mhgRefMultVsRun_%s",name.Data()),"Global reference multiplicity vs run; run index; gRefMult",nRun, run_lowBin, run_highBin);

      mhTpcVxVsRun[i] = new TProfile(Form("mhTpcVxVsRun_%s",name.Data()),"TPC v_{x} vs run; run index; TPC v_{x} (cm)",nRun, run_lowBin, run_highBin);
      
      mhTpcVyVsRun[i] = new TProfile(Form("mhTpcVyVsRun_%s",name.Data()),"TPC v_{y} vs run; run index; TPC v_{y} (cm)",nRun, run_lowBin, run_highBin);
      
      mhTpcVzVsRun[i] = new TProfile(Form("mhTpcVzVsRun_%s",name.Data()),"TPC v_{z} vs run; run index; TPC v_{z} (cm)",nRun, run_lowBin, run_highBin);

      mhVpdVzVsRun[i] = new TProfile(Form("mhVpdVzVsRun_%s",name.Data()),"VPD v_{z} vs run; run index; VPD v_{z} (cm)",nRun, run_lowBin, run_highBin);

      mhDiffVzVsRun[i] = new TProfile(Form("mhDiffVzVsRun_%s",name.Data()),"TPC-VPD v_{z} vs run; run index;TPC-VPD v_{z} (cm)",nRun, run_lowBin, run_highBin);

      mhpTrkPtVsRun[i] = new TProfile(Form("mhpTrkPtVsRun_%s",name.Data()),"Primary track p_{T} vs run; run index; p_{T,trk} (GeV/c)",nRun, run_lowBin, run_highBin);

      mhpTrkEtaVsRun[i] = new TProfile(Form("mhpTrkEtaVsRun_%s",name.Data()),"Primary track #eta vs run; run index; #eta_{trk}",nRun, run_lowBin, run_highBin);

      mhpTrkPhiVsRun[i] = new TProfile(Form("mhpTrkPhiVsRun_%s",name.Data()),"Primary track #varphi vs run; run index; #varphi_{trk}",nRun, run_lowBin, run_highBin);

      mhpTrkDcaVsRun[i] = new TProfile(Form("mhpTrkDcaVsRun_%s",name.Data()),"Primary track DCA vs run; run index; DCA (cm)",nRun, run_lowBin, run_highBin);

      mhNHitsFitVsRun[i] = new TProfile(Form("mhNHitsFitVsRun_%s",name.Data()),"Primary track NHitsFit vs run; run index; NhitsFit",nRun, run_lowBin, run_highBin);

      mhNHitsPossVsRun[i] = new TProfile(Form("mhNHitsPossVsRun_%s",name.Data()),"Primary track NHitsPoss vs run; run index; NHitsPoss",nRun, run_lowBin, run_highBin);

      mhNHitsDedxVsRun[i] = new TProfile(Form("mhNHitsDedxVsRun_%s",name.Data()),"Primary track NHitsDedx vs run; run index; NhitsDedx",nRun, run_lowBin, run_highBin);

      mhDedxVsRun[i] = new TProfile(Form("mhDedxVsRun_%s",name.Data()),"Primary track dE/dx vs run; run index; dE/dx (keV/cm)",nRun, run_lowBin, run_highBin);

      mhNsigmaPiVsRun[i] = new TProfile(Form("mhNsigmaPiVsRun_%s",name.Data()),"Primary track n#sigma_{#pi} vs run; run index; n#sigma_{#pi}",nRun, run_lowBin, run_highBin);

      mhNsigmaEVsRun[i] = new TProfile(Form("mhNsigmaEVsRun_%s",name.Data()),"Primary track n#sigma_{e} vs run; run index; n#sigma_{e}",nRun, run_lowBin, run_highBin);

      mhNsigmaKVsRun[i] = new TProfile(Form("mhNsigmaKVsRun_%s",name.Data()),"Primary track n#sigma_{K} vs run; run index; n#sigma_{K}",nRun, run_lowBin, run_highBin);

      mhNsigmaPVsRun[i] = new TProfile(Form("mhNsigmaPVsRun_%s",name.Data()),"Primary track n#sigma_{P} vs run; run index; n#sigma_{P}",nRun, run_lowBin, run_highBin);

      mhBetaVsRun[i] = new TProfile(Form("mhBetaVsRun_%s",name.Data()),"Primary track 1/#beta vs run; run index; 1/#beta",nRun, run_lowBin, run_highBin);

      mhNElectron[i] = new TProfile(Form("mhNElectron_%s",name.Data()),"Number of electron candidates vs run; run index; N_{electron}",nRun, run_lowBin, run_highBin);

      mhNPositron[i] = new TProfile(Form("mhNPositron_%s",name.Data()),"Number of positron candidates vs run; run index; N_{positron}",nRun, run_lowBin, run_highBin);

      mhNMtdRawHitsVsRun[i] = new TProfile(Form("mhNMtdRawHitsVsRun_%s",name.Data()),"Number of MTD raw hits per event vs run; run index; N",nRun, run_lowBin, run_highBin);

      mhNMtdHitsVsRun[i] = new TProfile(Form("mhNMtdHitsVsRun_%s",name.Data()),"Number of MTD hits per event vs run; run index; N",nRun, run_lowBin, run_highBin);

      mhNMtdTrigHitsVsRun[i] = new TProfile(Form("mhNMtdTrigHitsVsRun_%s",name.Data()),"Number of triggering MTD hits per event vs run; run index; N",nRun, run_lowBin, run_highBin);

      mhNMtdMthHitsVsRun[i] = new TProfile(Form("mhNMtdMthHitsVsRun_%s",name.Data()),"Number of matched MTD hits per event vs run; run index; N",nRun, run_lowBin, run_highBin);

      mhNMuonPosVsRun[i] = new TProfile(Form("mhNMuonPosVsRun_%s",name.Data()),"Number of positive muon candidates per event vs run; run index; N",nRun, run_lowBin, run_highBin);

      mhNMuonNegVsRun[i] = new TProfile(Form("mhNMuonNegVsRun_%s",name.Data()),"Number of negative muon candidates per event vs run; run index; N",nRun, run_lowBin, run_highBin);

      mhNMuonPairULVsRun[i] = new TProfile(Form("mhNMuonPairULVsRun_%s",name.Data()),"Number of unlike-sign muon pairs per event vs run; run index; N",nRun, run_lowBin, run_highBin);

      mhNMuonPairLSPosVsRun[i] = new TProfile(Form("mhNMuonPairLSPosVsRun_%s",name.Data()),"Number of positve like-sign muon pairs per event vs run; run index; N",nRun, run_lowBin, run_highBin);

      mhNMuonPairLSNegVsRun[i] = new TProfile(Form("mhNMuonPairLSNegVsRun_%s",name.Data()),"Number of negative like-sign muon pairs per event vs run; run index; N",nRun, run_lowBin, run_highBin);

      mhNJpsiVsRun[i] = new TProfile(Form("mhNJpsiVsRun_%s",name.Data()),"Number of J/psi candidates per event (US-LS, 3.0 < M_{#mu#mu} < 3.2 GeV/c^{2}); run index; N",nRun, run_lowBin, run_highBin);
    }
}


//_____________________________________________________________________________
void StMtdProdQAMaker::printConfig()
{
  const char *type[2] = {"MuDst", "PicoDst"};
  const char *decision[2] = {"no","yes"};
  printf("=== Configuration for StMtdProdQAMaker ===\n");
  printf("Data type: %s\n",type[mDataType]);
  printf("Run year: %d\n",mYear);
  printf("Maximum TPC z: %1.0f\n",mMaxTpcVz);
  printf("Maximum vz diff: %1.0f\n",mMaxDiffVz);
  printf("Track pt  range: [%1.2f, %1.2f]\n",mMinTrkPt,mMaxTrkPt);
  printf("Track phi range: [%1.2f, %1.2f]\n",mMinTrkPhi,mMaxTrkPhi);
  printf("Track eta range: [%1.2f, %1.2f]\n",mMinTrkEta,mMaxTrkEta);
  printf("Minimum number of fit hits: %d\n",mMinNHitsFit);
  printf("Minimum number of dedx hits: %d\n",mMinNHitsDedx);
  printf("Minimum fraction of fit hits: %4.2f\n",mMinNHitsFrac);
  printf("Maximum dca: %1.2f\n",mMaxDca);
  printf("Muon PID cuts:\n");
  printf("    %1.1f < NsigmaPi < %1.1f\n",mMinNsigmaPi,mMaxNsigmaPi);
  printf("    TOF match: %s\n",decision[mBTofMatch]);
  printf("    MTD hit trigger: %s\n",decision[mMtdHitTrigger]);
  printf("    %1.0f < dz < %1.0f cm\n",mMinMuonDeltaZ,mMaxMuonDeltaZ);
  printf("    %1.0f < dy < %1.0f cm\n",mMinMuonDeltaY,mMaxMuonDeltaY);
  printf("    %1.1f < dtof < %1.1f ns\n",mMinMuonDeltaTof,mMaxMuonDeltaTof);
  printf("    dca < %1.1f cm \n", mMaxMuonDca);
  printf("=======================================\n");
}

//_____________________________________________________________________________
Bool_t StMtdProdQAMaker::isMtdHitFiredTrigger(const StPicoMtdHit *hit)
{
  return isQTFiredTrigger( mModuleToQT[hit->backleg()-1][hit->module()-1], mModuleToQTPos[hit->backleg()-1][hit->module()-1]);
}

//_____________________________________________________________________________
Bool_t StMtdProdQAMaker::isMtdHitFiredTrigger(const StMuMtdHit *hit)
{
  return isQTFiredTrigger( mModuleToQT[hit->backleg()-1][hit->module()-1], mModuleToQTPos[hit->backleg()-1][hit->module()-1]);
}

//_____________________________________________________________________________
Bool_t StMtdProdQAMaker::isQTFiredTrigger(const Int_t qt, const Int_t pos)
{
  return (pos==mTrigQTpos[qt-1][0] || pos==mTrigQTpos[qt-1][1]);
}

//_____________________________________________________________________________
Bool_t StMtdProdQAMaker::isMuonCandidate(const StMuTrack *track)
{
  if(!track) return kFALSE;
  double nSigmaPi = -999.;
  int    tofIndex = -1;
  double dz = -999., dy = -999, dtof = -999., dca = -999.;
  bool isMtdTrig = kFALSE;

  // nSigmaPi cut
  nSigmaPi = track->nSigmaPion();

  // tof match
  tofIndex = track->index2BTofHit();

  // dca cut
  dca = track->dcaGlobal().mag();

  // dz cut
  int iMtd = track->index2MtdHit();
  if(iMtd>=0)
    {
      const StMuMtdPidTraits mtdPid = track->mtdPidTraits();
      dy     = mtdPid.deltaY();
      dz     = mtdPid.deltaZ();
      dtof   = mtdPid.timeOfFlight() - mtdPid.expTimeOfFlight();

      StMuMtdHit *hit = mMuDst->mtdHit(iMtd);
      isMtdTrig = isMtdHitFiredTrigger(hit);
    }

  return isMuonCandidate(nSigmaPi, tofIndex, dz, dy, dtof, dca, isMtdTrig);
}

//_____________________________________________________________________________
Bool_t StMtdProdQAMaker::isMuonCandidate(const StPicoTrack *track, const StThreeVectorD vtxPos)
{
  if(!track) return kFALSE;
  double nSigmaPi = -999.;
  int    tofIndex = -1;
  double dz = -999., dy = -999, dtof = -999., dca = -999;
  bool isMtdTrig = kFALSE;

  // nSigmaPi cut
  nSigmaPi = track->nSigmaPion();

  // tof match
  tofIndex = track->bTofPidTraitsIndex();

  // dca cut
  // dca = track->helix().distance(vtxPos);
  // dca = (track->dca()-vtxPos).mag();
  // dca = (track->dcaPoint()-vtxPos).mag();
  dca = track->gDCA(mPicoDst->event()->primaryVertex()).Mag();

  // dz cut
  int iMtd = track->mtdPidTraitsIndex();
  if(iMtd>=0)
    {
      StPicoMtdPidTraits *mtdPid = mPicoDst->mtdPidTraits(iMtd);
      dy = mtdPid->deltaY();
      dz = mtdPid->deltaZ();
      dtof = mtdPid->deltaTimeOfFlight();

      int hitIndex = getMtdHitIndex(track);
      StPicoMtdHit *hit = mPicoDst->mtdHit(hitIndex);
      isMtdTrig = hit->triggerFlag();
    }
	  
  return isMuonCandidate(nSigmaPi, tofIndex, dz, dy, dtof, dca, isMtdTrig);
}


//_____________________________________________________________________________
Bool_t StMtdProdQAMaker::isMuonCandidate(const Double_t nSigmaPi, const Int_t tofindex, const Double_t dz, const Double_t dy, const Double_t dtof, const Double_t dca, const Bool_t isTrig)
{
  if(mMaxMuonDeltaZ<1e4 && (dz<mMinMuonDeltaZ || dz>mMaxMuonDeltaZ))            return kFALSE;
  if(mMaxMuonDeltaY<1e4 && (dy<mMinMuonDeltaY || dy>mMaxMuonDeltaY))            return kFALSE;
  if(mMaxMuonDeltaTof<1e4 && (dtof<mMinMuonDeltaTof || dtof>mMaxMuonDeltaTof))  return kFALSE;
  if(mMaxNsigmaPi<1e4 && (nSigmaPi<mMinNsigmaPi || nSigmaPi>mMaxNsigmaPi))      return kFALSE;
  if(mBTofMatch && tofindex<0)                        return kFALSE;
  if(mMtdHitTrigger && !isTrig)                       return kFALSE;
  if(mMaxMuonDca<1e4 && dca > mMaxMuonDca)            return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Int_t StMtdProdQAMaker::getMtdHitIndex(const StPicoTrack *track)
{
  Int_t index = -1;
  if(track->mtdPidTraitsIndex()>=0)
    {
      StPicoMtdPidTraits *mtdPid = mPicoDst->mtdPidTraits(track->mtdPidTraitsIndex());
      Int_t nMtdHits = mPicoDst->numberOfMtdHits();
      for(Int_t i=0; i<nMtdHits; i++)
	{
	  StPicoMtdHit *hit = mPicoDst->mtdHit(i);
	  if(!hit) continue;
	  if(mtdPid->backleg()==hit->backleg() &&
	     mtdPid->module()==hit->module() &&
	     mtdPid->cell()==hit->cell())
	    {
	      index = i;
	      break;
	    }
	}
    }
  return index;
}

//_____________________________________________________________________________
Int_t StMtdProdQAMaker::getMtdPidTraitsIndex(const StPicoMtdHit *hit)
{
  Int_t index = -1;
  Int_t nPidTraits = mPicoDst->numberOfMtdPidTraits();
  for(Int_t i=0; i<nPidTraits; i++)
    {
      StPicoMtdPidTraits *mtdPid = mPicoDst->mtdPidTraits(i);
      if(!mtdPid) continue;
      if(mtdPid->backleg()==hit->backleg() &&
	 mtdPid->module()==hit->module() &&
	 mtdPid->cell()==hit->cell())
	{
	  index = i;
	  break;
	}
    }
  return index;
}

//_____________________________________________________________________________
Int_t StMtdProdQAMaker::getMtdHitTHUB(const Int_t backleg)
{
  if(backleg>=1 && backleg<=15)        return 2;
  else if (backleg>=16 && backleg<=30) return 1;
  else return -1;
}


//_____________________________________________________________________________
Bool_t StMtdProdQAMaker::isValidTrack(StMuTrack *track) 
{
  if(!track) return kFALSE;
  StThreeVectorF mom = track->momentum();
  Float_t pt = mom.perp();
  Float_t eta = mom.pseudoRapidity();
  Float_t phi = mom.phi();
  if(phi<0) phi += 2*pi;

  if(pt < mMinTrkPt   || pt > mMaxTrkPt)             return kFALSE;
  if(eta < mMinTrkEta || eta > mMaxTrkEta)           return kFALSE;
  if(phi < mMinTrkPhi || phi > mMaxTrkPhi)           return kFALSE;
  if(track->nHitsFit(kTpcId)<mMinNHitsFit)           return kFALSE;
  if(track->nHitsDedx()<mMinNHitsDedx)               return kFALSE;
  if(mMaxDca<1e4 && track->dcaGlobal().mag()>mMaxDca)      return kFALSE;
  if(track->nHitsFit(kTpcId)/(1.0*track->nHitsPoss(kTpcId))<mMinNHitsFrac) return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t StMtdProdQAMaker::isValidTrack(StPicoTrack *track, const StThreeVectorF mom, StThreeVectorD vtxPos) 
{
  if(!track) return kFALSE;
  Float_t pt = mom.perp();
  Float_t eta = mom.pseudoRapidity();
  Float_t phi = mom.phi();
  if(phi<0) phi += 2*pi;

  if(pt < mMinTrkPt   || pt > mMaxTrkPt)       return kFALSE;
  if(eta < mMinTrkEta || eta > mMaxTrkEta)     return kFALSE;
  if(phi < mMinTrkPhi || phi > mMaxTrkPhi)     return kFALSE;
  if(track->nHitsFit()<mMinNHitsFit)           return kFALSE;
  if(track->nHitsDedx()<mMinNHitsDedx)         return kFALSE;
  if(track->nHitsFit()/(1.0*track->nHitsMax())<mMinNHitsFrac) return kFALSE;
  if(mMaxDca<1e4)
    {
      // double dca = track->helix().distance(vtxPos);
      //double dca = (track->dca()-vtxPos).mag();
      // double dca = (track->dcaPoint()-vtxPos).mag();
      double dca = track->gDCA(mPicoDst->event()->primaryVertex()).Mag();
      if(dca>mMaxDca) return kFALSE;
    }
  return kTRUE;
}

//_____________________________________________________________________________
float StMtdProdQAMaker::rotatePhi(float phi) const
{
  Double_t outPhi = phi;
  while(outPhi<0) outPhi += 2*pi;
  while(outPhi>2*pi) outPhi -= 2*pi;
  return outPhi;
}

// $Id: StMtdProdQAMaker.cxx,v 1.18 2019/01/02 21:10:49 marr Exp $
// $Log: StMtdProdQAMaker.cxx,v $
// Revision 1.18  2019/01/02 21:10:49  marr
// First commit of macros used to select cosmic events with at least two tracks
//
// Revision 1.17  2018/08/14 15:25:05  marr
// i) fill MtdVpdTacDiff for matched hits; ii) merge primary track loops
//
// Revision 1.16  2018/08/13 17:54:08  marr
// Modify to be compatible with latest picoDst code
//
// Revision 1.15  2016/03/30 18:17:10  marr
// local backup
//
//
