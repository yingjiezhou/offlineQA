#include "headers.h"
#include "StQAMaker.h"

ClassImp(StQAMaker)

//_____________________________________________________________________________
StQAMaker::StQAMaker(const char* name, StPicoDstMaker *picoMaker) : StMaker(name), mPrintConfig(1), mPrintMemory(0), mPrintCpu(0), mStreamName(""), fOutFile(0), mOutFileName(""), mEvtTree(0), mDefaultVtx(1), mSelectVtxRank(0), mMaxVtxR(2.), mMaxVtxZ(100.), mMaxVzDiff(3.), mMinTrkPt(0.2), mMaxTrkEta(1.), mMinNHitsFit(15), mMinNHitsFitRatio(0.52), mMinNHitsDedx(10), mMaxDca(3.), mMaxnSigmaE(2.0), mMaxBeta2TOF(0.03), mMinBemcPt(3.5), mMinAdc0(290), mMinPoverE(0.3), mMaxPoverE(1.5), mMaxZDist(3), mMaxPhiDist(0.02), mMinNEta(1), mMinNPhi(1)
{
  // default constructor
  mPicoDstMaker = picoMaker; //??
  mPicoDst = 0;
  mevent = 0;
  
  // run20 st_physics
  mTriggerIDs.clear();
  mTriggerIDs.push_back(350503);  // NPE_18
  mTriggerIDs.push_back(350513);  // NPE_18
  //mTriggerIDs.push_back(350504);  // NPE_25
  //mTriggerIDs.push_back(350514);  // NPE_25
  //mTriggerIDs.push_back(350501);  // NPE_25_nozdc
  //mTriggerIDs.push_back(350511);  // NPE_25_nozdc
}
//_____________________________________________________________________________
StQAMaker::~StQAMaker()
{
  // default destructor
}
//_____________________________________________________________________________
Int_t StQAMaker::Init()
{
  
  if(!mOutFileName.Length()){
    LOG_ERROR << "StQAMaker:: no output file specified for tree and histograms." << endm;
    return kStERR;
  }
  
  fOutFile = new TFile(mOutFileName.Data(),"recreate");
  LOG_INFO << "StQAMaker:: create the output file to store the tree and histograms: " << mOutFileName.Data() << endm;
  
  bookHistos();
  
  ifstream inData;
  //  inData.open("StRoot/StQAMaker/run20_9p2_stPhysics_runnumber_DD.dat");
  inData.open("StRoot/StQAMaker/0827_runnumber_26p5_2020_pico_fixed.dat");
  
  if(!inData.is_open()){
    LOG_ERROR << "Failed to runnumber list from loal file !" << endm;
    return kStErr;
  }
  
  cout << "Retrieving runnumber list from local file ...";
  mTotalRunId.clear();
  Int_t oldId;
  Int_t newId = 0;
  while(inData>>oldId){
    mTotalRunId[oldId] = newId;
    newId++;
  }
  cout<<" [OK]"<<endl;
  
  if(Debug()){
    for(map<Int_t,Int_t>::iterator iter=mTotalRunId.begin();iter!=mTotalRunId.end();iter++)
      cout<< setiosflags(ios::left) << "Run index: "<< setw(10) <<iter->second<< "runId: " << setw(10) <<iter->first<<endl;
    cout<<endl;
  }
  
  inData.close();
  
  return kStOK;
}
//_____________________________________________________________________________
Int_t StQAMaker::Finish()
{
  if(fOutFile){
    fOutFile->cd();
    fOutFile->Write();
    fOutFile->Close();
    LOG_INFO << "StQAMaker::Finish() -> write out tree in " << mOutFileName.Data() << endm;
  }
  
  if(mPrintConfig) printConfig();
  
  return kStOK;
}
//_____________________________________________________________________________
Int_t StQAMaker::Make()
{
  StTimer timer;
  if(mPrintMemory) StMemoryInfo::instance()->snapshot();
  if(mPrintCpu)    timer.start();
  
  if (!mPicoDstMaker)
  {
    LOG_WARN << " StQAMaker::Make() - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }
  StPicoDst const* mPicoDst = mPicoDstMaker->picoDst();
  if (!mPicoDst)
  {
    LOG_WARN << "StQAMaker::Make() - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  
  
  if(!processPicoDstEvent()) return kStOK;
  
  if(mPrintMemory){
    StMemoryInfo::instance()->snapshot();
    StMemoryInfo::instance()->print();
  }
  
  if(mPrintCpu){
    timer.stop();
    LOG_INFO << "CPU time for StQAMaker::Make(): "
    << timer.elapsedTime() << "sec " << endm;
  }
  
  return kStOK;
}
//_____________________________________________________________________________
Bool_t StQAMaker::processPicoDstEvent()
{
  hEvent_DefVtx->Fill(0.5);
  hEvent->Fill(0.5);
  
  mevent=mPicoDst->event();
  if(!mevent){
    LOG_WARN<<"No event level information !"<<endm;
    return kFALSE;
  }
  
  Bool_t NPE_18       = kFALSE;
  Bool_t NPE_25       = kFALSE;
  Bool_t NPE_25_nozdc = kFALSE;
  
  
  //
//  Bool_t validTrig = kFALSE;
  //  if(mTriggerIDs.size()==0){
  ////    for(Int_t i=0;i<64;i++){
  ////      Int_t trgId = mevent->triggerIdCollection().nominal().triggerId(i);
  ////      if(trgId>0){
  //        if(mevent->isTrigger(610001) || mevent->isTrigger(610011) || mevent->isTrigger(610021) || mevent->isTrigger(610031) || mevent->isTrigger(610041) || mevent->isTrigger(610051))
  //        validTrig = kTRUE;
  //        //        if(350503==trgId || 350513==trgId) NPE_18       = kTRUE;
  //        //        if(350504==trgId || 350514==trgId) NPE_25       = kTRUE;
  //        //        if(350501==trgId || 350511==trgId) NPE_25_nozdc = kTRUE;
  ////      }
  //    }
  //  }
  //  else{
  ////    for(Int_t i=0;i<mTriggerIDs.size();i++){
  ////      if(mevent->triggerIdCollection().nominal().isTrigger(mTriggerIDs[i])){
  //        if(mevent->isTrigger(610001) || mevent->isTrigger(610011) || mevent->isTrigger(610021) || mevent->isTrigger(610031) || mevent->isTrigger(610041) || mevent->isTrigger(610051))
  //        validTrig = kTRUE;
  //
  //        //        if(350503==mTriggerIDs[i] || 350513==mTriggerIDs[i]) NPE_18       = kTRUE;
  //        //        if(350504==mTriggerIDs[i] || 350514==mTriggerIDs[i]) NPE_25       = kTRUE;
  //        //        if(350501==mTriggerIDs[i] || 350511==mTriggerIDs[i]) NPE_25_nozdc = kTRUE;
  //        //        if(610001==mTriggerIDs[i] || 610001==mTriggerIDs[i] || 610001==mTriggerIDs[i]
  //        //           || 610001==mTriggerIDs[i] || 610001==mTriggerIDs[i] || 610001==mTriggerIDs[i])
  //        //if(i<2)      NPE_18       = kTRUE;
  //        //else if(i<4) NPE_25       = kTRUE;
  //        //else         NPE_25_nozdc = kTRUE;
  ////      }
  ////    }
  //  }
  
//  if(mevent->isTrigger(610001) || mevent->isTrigger(610011) || mevent->isTrigger(610021) || mevent->isTrigger(610031) || mevent->isTrigger(610041) || mevent->isTrigger(610051))
//    validTrig = kTRUE;
//  if(!validTrig){
//    if(Debug()) LOG_WARN<<"No valid interested triggers !"<<endm;
//    return kFALSE;
//  }
  //  Bool_t validTrig = kFALSE;
  //  if (mevent->isTrigger(780020)) validTrig = kTRUE;;
  //  if(!validTrig) return kFALSE;
  
  Int_t runIdx;
  Int_t runId = mevent->runId();
  map<Int_t,Int_t>::iterator iter = mTotalRunId.find(runId);
  if(iter != mTotalRunId.end())
    runIdx = iter->second;
  else{
    runIdx = -1;
    LOG_INFO<<"Can not find the runNumber in the runNumber list"<<endm;
  }
  
  //  if(NPE_18)       { hEvent_DefVtx->Fill(2.5); hEvent->Fill(2.5); hnHT2EvtsvsRun->Fill(runIdx); }
  //  if(NPE_25)       { hEvent_DefVtx->Fill(3.5); hEvent->Fill(3.5); }
  //  if(NPE_25_nozdc) { hEvent_DefVtx->Fill(4.5); hEvent->Fill(4.5); }
  
  hnEvtsvsRun->Fill(runIdx);
  
  Double_t vpdVz = mevent->vzVpd();
  TVector3 vtxPos = mevent->primaryVertex();
  Double_t vx      = vtxPos.x();
  Double_t vy      = vtxPos.y();
  Double_t vz      = vtxPos.z();
  Double_t vr      = sqrt(vx*vx + vy*vy);
  Double_t vzDiff  = vz - vpdVz;
  Double_t refMult = mevent->refMult();
  
  hVyvsVx_DefVtx->Fill(vx, vy);
  hVpdVzvsTpcVz_DefVtx->Fill(vz, vpdVz);
  hVzDiffvsTpcVz_DefVtx->Fill(vz, vzDiff);
  hVzDiffvsRefMult_DefVtx->Fill(refMult, vzDiff);
  
  if(TMath::Abs(vx)>=1.e-5 || TMath::Abs(vy)>=1.e-5 || TMath::Abs(vz)>=1.e-5){
    hEvent_DefVtx->Fill(7.5);
    if(vr<=mMaxVtxR){
      hEvent_DefVtx->Fill(8.5);
      if(TMath::Abs(vz)<=mMaxVtxZ){
        hEvent_DefVtx->Fill(9.5);
        hRefMult_VzVrCut_DefVtx->Fill(refMult);
        if(TMath::Abs(vzDiff)<=mMaxVzDiff){
          hEvent_DefVtx->Fill(10.5);
          hRefMult_EvtCut_DefVtx->Fill(refMult);
          //          if(NPE_18)       hEvent_DefVtx->Fill(12.5);
          //          if(NPE_25)       hEvent_DefVtx->Fill(13.5);
          //          if(NPE_25_nozdc) hEvent_DefVtx->Fill(14.5);
        }
      }
    }
  }
  
  hVyvsVx->Fill(vx, vy);
  hVpdVzvsTpcVz->Fill(vz, vpdVz);
  hTpcVzvsRefMult->Fill(refMult, vz);
  hRawVpdVzvsRefMult->Fill(refMult, vpdVz);
  if(TMath::Abs(vz)<=mMaxVtxZ) hVpdVzvsRefMult->Fill(refMult, vpdVz);
  hVzDiffvsTpcVz->Fill(vz, vzDiff);
  hVzDiffvsRefMult->Fill(refMult, vzDiff);
  
  //if(Debug()){
  //  LOG_INFO<<"vtxIdx: "<<vtxIdx<<" \tTPC Vx: "<<vx<<" \tTPC Vy: "<<vy<<" \tTPC Vz: "<<vz<<endm;
  //}
  
  Double_t bField  = mevent->bField()*kilogauss; // event->bField()).Mag()
  Double_t bbcRate = mevent->BBCx();
  Double_t zdcRate = mevent->ZDCx();
  
  hRefMultvsZdcX->Fill(zdcRate/1000., refMult);
  hRefMultvsBbcX->Fill(bbcRate/1000., refMult);
  
  hBFieldvsRun->Fill(runIdx, bField);
  hBbcXvsRun->Fill(runIdx, bbcRate/1000.);
  hZdcXvsRun->Fill(runIdx, zdcRate/1000.);
  hZdcXoverBbcXvsRun->Fill(runIdx, zdcRate/bbcRate);
  hTpcVxvsRun->Fill(runIdx, vx);
  hTpcVyvsRun->Fill(runIdx, vy);
  hTpcVzvsRun->Fill(runIdx, vz);
  hRawVpdVzvsRun->Fill(runIdx, vpdVz);
  if(TMath::Abs(vpdVz)<600){ // +/-570cm are the locations where VPDs sit
    hVpdVzvsRun->Fill(runIdx, vpdVz);
    hVzDiffvsRun->Fill(runIdx, vzDiff);
  }
  hRefMultvsRun->Fill(runIdx, refMult);
  hgRefMultvsRun->Fill(runIdx, mevent->grefMult());
  hTofMultvsRun->Fill(runIdx, mevent->btofTrayMultiplicity());
  hnTofMatchvsRun->Fill(runIdx, mevent->nBTOFMatch());
  hnBEMCMatchvsRun->Fill(runIdx, mevent->nBEMCMatch());
  
  //hVtxIdxvsRefMult->Fill(refMult, vtxIdx);
  //if(TMath::Abs(vzDiff)<3) hVtxIdxvsRefMult_VzDiffCut->Fill(refMult, vtxIdx);
  
  // begin add cut
  Double_t vtxRanking  = mevent->ranking();
  if(mSelectVtxRank && vtxRanking<=0) return kFALSE;
  //  hEvent->Fill(6.5);
  if(TMath::Abs(vx)<1.e-5 && TMath::Abs(vy)<1.e-5 && TMath::Abs(vz)<1.e-5) return kFALSE;
  hEvent->Fill(7.5);
  if(vr>mMaxVtxR) return kFALSE;
  hEvent->Fill(8.5);
  if(TMath::Abs(vz)>mMaxVtxZ) return kFALSE;
  hEvent->Fill(9.5);
  hRefMult_VzVrCut->Fill(refMult);
  if(TMath::Abs(vzDiff)>mMaxVzDiff)  return kFALSE;
  hEvent->Fill(10.5);
  hRefMult_EvtCut->Fill(refMult);
  
  //  if(NPE_18)       hEvent->Fill(12.5);
  //  if(NPE_25)       hEvent->Fill(13.5);
  //  if(NPE_25_nozdc) hEvent->Fill(14.5);
  
  Int_t nNodes = mPicoDst->numberOfTracks();
  if(Debug()){
    LOG_INFO<<"# of Tracks in PicoDst: "<<nNodes<<endm;
  }
  
  Int_t nTrks      = 0;
  Int_t nMthTrks   = 0;
  Int_t nTrigTrks  = 0;
  Int_t nBEMCeCans = 0;
  for(Int_t i=0;i<nNodes;i++){
    StPicoTrack* trk = mPicoDst->track(i);
    if(!trk) continue;
    
    if(!isValidTrack(trk)) continue;
    
    TVector3 pMom = trk->pMom();
    TVector3 gmom = trk->gMom();
    Double_t pt         = pMom.Perp();
    Double_t eta        = pMom.PseudoRapidity();
    Double_t phi        = pMom.Phi();
    Int_t    charge     = trk->charge();
    Double_t nHitsFit   = trk->nHitsFit();
    Double_t nHitsPoss  = trk->nHitsPoss(); // = nHitsMax, Possible number of hits (in TPC)
    Double_t nHitsDedx  = trk->nHitsDedx();
    Double_t dca        = trk->gDCA(vtxPos).Mag();
    Double_t dEdx       = trk->dEdx();
    Double_t nSigmaE    = trk->nSigmaElectron();
    Double_t nSigmaPi   = trk->nSigmaPion();
    Double_t nSigmaK    = trk->nSigmaKaon();
    Double_t nSigmaP    = trk->nSigmaProton();
    
    hNHitsFit->Fill(charge*nHitsFit);
    hNHitsPoss->Fill(charge*nHitsPoss);
    hNHitsDedx->Fill(charge*nHitsDedx);
    hEtavsPt->Fill(charge*pt, eta);
    if(bField>0) hFFPhivsPt->Fill(charge*pt, phi);
    else         hRFFPhivsPt->Fill(charge*pt, phi);
    if(charge>0) hPosTrkEtavsPhi->Fill(phi, eta);
    else         hNegTrkEtavsPhi->Fill(phi, eta);
    hDcavsPt->Fill(charge*pt, dca);
    hdEdxvsP->Fill(charge*pMom.Mag(), dEdx);
    hdEdxvsPhi->Fill(phi, dEdx);
    hdEdxvsEta->Fill(eta, dEdx);
    
    Double_t beta       = getTofBeta(trk);;
    bool tofmatch = (beta!=std::numeric_limits<float>::quiet_NaN()) && beta>0;
    if(tofmatch){
      hBetavsP->Fill(charge*pMom.Mag(), 1/beta);
      hMassvsP->Fill(charge*pMom.Mag(),pMom.Mag()*sqrt(1-beta*beta)*1.0/beta);
    }
    
    hPtvsRun->Fill(runIdx, pt);
    hEtavsRun->Fill(runIdx, eta);
    hPhivsRun->Fill(runIdx, phi);
    hDcavsRun->Fill(runIdx, dca);
    hNHitsFitvsRun->Fill(runIdx, nHitsFit);
    hNHitsPossvsRun->Fill(runIdx, nHitsPoss);
    hNHitsDedxvsRun->Fill(runIdx, nHitsDedx);
    hDedxvsRun->Fill(runIdx, dEdx);
    hNSigmaEvsRun->Fill(runIdx, nSigmaE);
    hNSigmaPivsRun->Fill(runIdx, nSigmaPi);
    hNSigmaKvsRun->Fill(runIdx, nSigmaK);
    hNSigmaPvsRun->Fill(runIdx, nSigmaP);
    if(beta>0) hBetavsRun->Fill(runIdx, 1/beta);
    
    getBemcInfo(trk, runIdx, nMthTrks, nTrigTrks, nBEMCeCans);
    
    getMtdInfo(trk, runIdx); // zyj
    
    nTrks++;
  }
  
  hNTrksvsRun->Fill(runIdx, nTrks);
  hNMthTrksvsRun->Fill(runIdx, nMthTrks);
  hNTrigTrksvsRun->Fill(runIdx, nTrigTrks);
  hNBemcEsvsRun->Fill(runIdx, nBEMCeCans);
  
  if(Debug()){
    LOG_INFO<<"# of good primary tracks: "<<nTrks<<endm;
    LOG_INFO<<"# of BEMC matched Tracks: "<<nMthTrks<<endm;
    LOG_INFO<<"# of BEMC triggers Tracks: "<<nTrigTrks<<endm;
    LOG_INFO<<"# of BEMC electron candidates: "<<nBEMCeCans<<endm;
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t StQAMaker::isValidTrack(StPicoTrack *trk) const
{
  Float_t pt  = trk->pMom().Perp();
  Float_t eta = trk->pMom().PseudoRapidity();
  Float_t dca = trk->gDCA(mevent->primaryVertex()).Mag();
  
  if(pt<mMinTrkPt)                            return kFALSE;
  if(TMath::Abs(eta)>mMaxTrkEta)              return kFALSE;
  if(trk->nHitsFit()<mMinNHitsFit) return kFALSE;
  if(trk->nHitsFit()*1./trk->nHitsPoss()<mMinNHitsFitRatio)  return kFALSE;
  if(trk->nHitsDedx()<mMinNHitsDedx)     return kFALSE;
  if(dca>mMaxDca)                             return kFALSE;
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t StQAMaker::getBemcInfo(StPicoTrack *trk, const Int_t runIdx, Int_t &nMthTrks, Int_t &nTrigTrks, Int_t &nBEMCeCans)
{
  Float_t maxtowerE = -999., energy = 0.;
  Float_t zdist = -999., phidist = -999., mindist = 999.;
  Int_t mod = -1, eta=-1, sub=-1;
  Int_t neta = -1, nphi=-1;
  UInt_t maxadc = 0;
  Int_t dsmadc   = -999;
  
  Bool_t bemcMatchFlag = kFALSE;
  Int_t bemcPidTraitsIndex = trk->bemcPidTraitsIndex();
  if(bemcPidTraitsIndex>=0) bemcMatchFlag = kTRUE;
  if(!bemcMatchFlag) return kFALSE;
  nMthTrks++;
  
  StPicoBEmcPidTraits *bemcPidTraits   =  mPicoDst->bemcPidTraits(bemcPidTraitsIndex);
  int TOWId = bemcPidTraits->btowId();
  energy     = bemcPidTraits->bemcE();
  maxtowerE = bemcPidTraits->bemcE0();
  float EMCpovere = (maxtowerE>0) ? trk->gMom().Mag()/maxtowerE : -999.; // p/E0(max)
  zdist    = bemcPidTraits->bemcZDist();
  phidist  = bemcPidTraits->bemcPhiDist();
  neta  = bemcPidTraits->bemcSmdNEta();
  nphi  = bemcPidTraits->bemcSmdNPhi();
  maxadc  = bemcPidTraits->bemcAdc0();
  
  int nEmcTrigger = mPicoDst->numberOfEmcTriggers();
  for(int nEmc=0;nEmc<nEmcTrigger;nEmc++){
    StPicoEmcTrigger *emcTrg = (StPicoEmcTrigger*)mPicoDst->emcTrigger(nEmc);
    int emcTrgID=emcTrg->id();
    if(emcTrgID==TOWId){ // trigger Id = bemc tow id
      dsmadc = emcTrg->adc();
      break;
    }
  }
  
  Int_t    charge     = trk->charge();
  Double_t trkp       = trk->pMom().Mag();
  Double_t trkpt      = trk->pMom().Perp();
  Double_t trketa     = trk->pMom().PseudoRapidity();
  Double_t trkphi     = trk->pMom().Phi();
  Double_t nSigmaE    = trk->nSigmaElectron();
  Double_t beta       = getTofBeta(trk);;
  
  hMthTrkPtvsRun->Fill(runIdx, trkpt);
  hMthTrkEtavsRun->Fill(runIdx, trketa);
  hMthTrkPhivsRun->Fill(runIdx, trkphi);
  hMthTrkNSigmaEvsRun->Fill(runIdx, nSigmaE);
  if(beta>0) hMthTrkBetavsRun->Fill(runIdx, 1/beta);
  hMthTrkAdc0vsRun->Fill(runIdx, maxadc);
  hMthTrkE0vsRun->Fill(runIdx, maxtowerE);
  hMthTrkEvsRun->Fill(runIdx, energy);
  hMthTrkZDistvsRun->Fill(runIdx, zdist);
  hMthTrkPhiDistvsRun->Fill(runIdx, phidist);
  hMthTrkNEtavsRun->Fill(runIdx, neta);
  hMthTrkNPhivsRun->Fill(runIdx, nphi);
  
  if(trkpt<mMinBemcPt)  return kFALSE;
  if(maxadc<mMinAdc0)   return kFALSE;
  nTrigTrks++;
  
  hTrigTrkPtvsRun->Fill(runIdx, trkpt);
  hTrigTrkEtavsRun->Fill(runIdx, trketa);
  hTrigTrkPhivsRun->Fill(runIdx, trkphi);
  hTrigTrkNSigmaEvsRun->Fill(runIdx, nSigmaE);
  hTrigTrkAdc0vsRun->Fill(runIdx, maxadc);
  hTrigTrkE0vsRun->Fill(runIdx, maxtowerE);
  hTrigTrkEvsRun->Fill(runIdx, energy);
  hTrigTrkZDistvsRun->Fill(runIdx, zdist);
  hTrigTrkPhiDistvsRun->Fill(runIdx, phidist);
  hTrigTrkNEtavsRun->Fill(runIdx, neta);
  hTrigTrkNPhivsRun->Fill(runIdx, nphi);
  
  if(beta>0 && TMath::Abs(1-1/beta)>mMaxBeta2TOF)      return kFALSE;
  if(TMath::Abs(nSigmaE)>mMaxnSigmaE)                  return kFALSE;
  if(trkp/energy<mMinPoverE || trkp/energy>mMaxPoverE) return kFALSE;
  if(TMath::Abs(zdist)>mMaxZDist)                      return kFALSE;
  if(TMath::Abs(phidist)>mMaxPhiDist)                  return kFALSE;
  if(neta<mMinNEta)                                    return kFALSE;
  if(nphi<mMinNPhi)                                    return kFALSE;
  nBEMCeCans++;
  
  hBEMCeEtavsPt->Fill(charge*trkpt, trketa);
  hBEMCePhivsPt->Fill(charge*trkpt, trkphi);
  hBEMCeEtavsPhi->Fill(trkphi, trketa);
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t StQAMaker::getMtdInfo(StPicoTrack *trk, const Int_t runIdx)
{
  
  int pmtdId = trk->mtdPidTraitsIndex();
  if(pmtdId>=0){
    StPicoMtdPidTraits * pmtdtrait = mPicoDst->mtdPidTraits(pmtdId);
    if(pmtdtrait){
      float pmtdDeltaY = pmtdtrait->deltaY();
      hMthTrkDeltaYvsRun->Fill(runIdx, pmtdDeltaY);
      float pmtdDeltaZ = pmtdtrait->deltaZ();
      hMthTrkDeltaZvsRun->Fill(runIdx, pmtdDeltaZ);
      float pmtdDeltaTOF = pmtdtrait->deltaTimeOfFlight();
      hMthTrkDeltaTOFvsRun->Fill(runIdx, pmtdDeltaTOF);
    }
  }
  /// MTD hits
  int nMtdHits = mPicoDst->numberOfMtdHits();
  int nMthHit = 0;
  for(int i=0; i<nMtdHits; i++)
  {
    StPicoMtdHit *hit = mPicoDst->mtdHit(i);
    if(!hit) continue;
    int backleg = hit->backleg();
    int module  = hit->module();
    int cell    = hit->cell();
    mhMtdHitMap->Fill(backleg, (module-1)*12+cell+1);
    
    bool isMth  = false;
    int index2traits = getMtdPidTraitsIndex(hit);
    if(index2traits>-1)
    {
      StPicoMtdPidTraits *mtdPid = mPicoDst->mtdPidTraits(index2traits);
      int index2track = mtdPid->trackIndex();
      StPicoTrack* gTrack = mPicoDst->track(index2track);
      if(gTrack)
      {
        if(isValidTrack(gTrack))isMth = true;
      }
    }
    if(isMth) {
      mhMtdMthHitMap->Fill(backleg, (module-1)*12+cell+1);
      nMthHit ++;
    }
  }
  
  mhNMtdHitsVsRun->Fill(runIdx, nMtdHits);
  mhNMtdMthHitsVsRun->Fill(runIdx, nMthHit);
  
  
  return kTRUE;
}
//_____________________________________________________________________________
void StQAMaker::bookHistos()
{
  // default vertex
  hEvent_DefVtx = new TH1D("hEvent_DefVtx","Event statistics",20,0,20);
  hEvent_DefVtx->GetXaxis()->SetBinLabel(1, "All events");
  //  hEvent_DefVtx->GetXaxis()->SetBinLabel(3, "NPE_18");
  //  hEvent_DefVtx->GetXaxis()->SetBinLabel(4, "NPE_25");
  //  hEvent_DefVtx->GetXaxis()->SetBinLabel(5, "NPE_25_nozdc");
  hEvent_DefVtx->GetXaxis()->SetBinLabel(7, "mRanking>0");
  hEvent_DefVtx->GetXaxis()->SetBinLabel(8, "None-Zero Vertex");
  hEvent_DefVtx->GetXaxis()->SetBinLabel(9,  Form("|V_{r}|<%1.2f cm",mMaxVtxR));
  hEvent_DefVtx->GetXaxis()->SetBinLabel(10, Form("|V_{z}|<%1.2f cm",mMaxVtxZ));
  hEvent_DefVtx->GetXaxis()->SetBinLabel(11, Form("|V_{z}Diff|<%1.2f cm",mMaxVzDiff));
  //  hEvent_DefVtx->GetXaxis()->SetBinLabel(13, "NPE_18");
  //  hEvent_DefVtx->GetXaxis()->SetBinLabel(14, "NPE_25");
  //  hEvent_DefVtx->GetXaxis()->SetBinLabel(15, "NPE_25_nozdc");
  hVyvsVx_DefVtx       = new TH2D("hVyvsVx_DefVtx","hVyvsVx_DefVtx; V_{x} (cm); V_{y} (cm)",1200,-6,6,1200,-6,6);
  hVpdVzvsTpcVz_DefVtx = new TH2D("hVpdVzvsTpcVz_DefVtx","hVpdVzvsTpcVz_DefVtx; TPC V_{z} (cm); VPD V_{z} (cm)",1600,-800,800,1600,-800,800);
  hVzDiffvsTpcVz_DefVtx   = new TH2D("hVzDiffvsTpcVz_DefVtx","hVzDiffvsTpcVz_DefVtx; TPC Vz (cm); Vz_{TPC} - Vz_{VPD} (cm)",1600,-800,800,32000,-1600,1600);
  hVzDiffvsRefMult_DefVtx = new TH2D("hVzDiffvsRefMult_DefVtx","hVzDiffvsRefMult_DefVtx; refMult; Vz_{TPC} - Vz_{VPD} (cm)",900,0,900,32000,-1600,1600);
  hRefMult_VzVrCut_DefVtx = new TH1D("hRefMult_VzVrCut_DefVtx","hRefMult_VzVrCut_DefVtx; refMult",900,0,900);
  hRefMult_EvtCut_DefVtx  = new TH1D("hRefMult_EvtCut_DefVtx","hRefMult_EvtCut_DefVtx; refMult",900,0,900);
  
  // selected vertex
  hEvent = new TH1D("hEvent","Event statistics",20,0,20);
  hEvent->GetXaxis()->SetBinLabel(1, "All events");
  //  hEvent->GetXaxis()->SetBinLabel(3, "NPE_18");
  //  hEvent->GetXaxis()->SetBinLabel(4, "NPE_25");
  //  hEvent->GetXaxis()->SetBinLabel(5, "NPE_25_nozdc");
  hEvent->GetXaxis()->SetBinLabel(7, "mRanking>0");
  hEvent->GetXaxis()->SetBinLabel(8, "None-Zero Vertex");
  hEvent->GetXaxis()->SetBinLabel(9,  Form("|V_{r}|<%1.2f cm",mMaxVtxR));
  hEvent->GetXaxis()->SetBinLabel(10, Form("|V_{z}|<%1.2f cm",mMaxVtxZ));
  hEvent->GetXaxis()->SetBinLabel(11, Form("|V_{z}Diff|<%1.2f cm",mMaxVzDiff));
  //  hEvent->GetXaxis()->SetBinLabel(13, "NPE_18");
  //  hEvent->GetXaxis()->SetBinLabel(14, "NPE_25");
  //  hEvent->GetXaxis()->SetBinLabel(15, "NPE_25_nozdc");
  hVyvsVx            = new TH2D("hVyvsVx","hVyvsVx; V_{x} (cm); V_{y} (cm)",1200,-6,6,1200,-6,6);
  hVpdVzvsTpcVz      = new TH2D("hVpdVzvsTpcVz","hVpdVzvsTpcVz; TPC V_{z} (cm); VPD V_{z} (cm)",1600,-800,800,1600,-800,800);
  hTpcVzvsRefMult    = new TH2D("hTpcVzvsRefMult","hTpcVzvsRefMult; refMult; TPC V_{z} (cm)",900,0,900,1600,-800,800);
  hRawVpdVzvsRefMult = new TH2D("hRawVpdVzvsRefMult","hRawVpdVzvsRefMult; refMult; VPD V_{z} (cm)",900,0,900,1600,-800,800);
  hVpdVzvsRefMult    = new TH2D("hVpdVzvsRefMult","hVpdVzvsRefMult; refMult; VPD V_{z} (cm)",900,0,900,1600,-800,800);
  hVzDiffvsTpcVz     = new TH2D("hVzDiffvsTpcVz","hVzDiffvsTpcVz; TPC Vz (cm); Vz_{TPC} - Vz_{VPD} (cm)",1600,-800,800,32000,-1600,1600);
  hVzDiffvsRefMult   = new TH2D("hVzDiffvsRefMult","hVzDiffvsRefMult; refMult; Vz_{TPC} - Vz_{VPD} (cm)",900,0,900,32000,-1600,1600);
  hRefMult_VzVrCut   = new TH1D("hRefMult_VzVrCut","hRefMult_VzVrCut; refMult",900,0,900);
  hRefMult_EvtCut    = new TH1D("hRefMult_EvtCut","hRefMult_EvtCut; refMult",900,0,900);
  hVtxIdxvsRefMult             = new TH2D("hVtxIdxvsRefMult","hVtxIdxvsRefMult; refMult; Vertex Index",900,0,900,20,0,20);
  hVtxIdxvsRefMult_VzDiffCut = new TH2D("hVtxIdxvsRefMult_VzDiffCut","hVtxIdxvsRefMult_VzDiffCut; refMult; Vertex Index",900,0,900,20,0,20);
  hRefMultvsZdcX     = new TH2D("hRefMultvsZdcX","hRefMultvsZdcX; zdcRate (kHz); refMult",600,0,60,900,0,900);
  hRefMultvsBbcX     = new TH2D("hRefMultvsBbcX","hRefMultvsBbcX; bbcRate (kHz); refMult",600,0,60,900,0,900);
  
  /***********   inclusive QA   ***********/
  hNHitsFit       = new TH1D("hNHitsFit","hNHitsFit;charge*nHitsFit;Counts",100,-50,50);
  hNHitsPoss      = new TH1D("hNHitsPoss","hNHitsPoss;charge*nHitsPoss;Counts",100,-50,50);
  hNHitsDedx      = new TH1D("hNHitsDedx","hNHitsDedx;charge*nHitsDedx;Counts",100,-50,50);
  hEtavsPt        = new TH2D("hEtavsPt","hEtavsPt; charge*p_{T} (GeV/c); #eta",300,-15,15,200,-1,1);
  hFFPhivsPt      = new TH2D("hFFPhivsPt","hFFPhivsPt; charge*p_{T} (GeV/c); #phi",750,-15,15,360,-TMath::Pi(),TMath::Pi());
  hRFFPhivsPt     = new TH2D("hRFFPhivsPt","hRFFPhivsPt; charge*p_{T} (GeV/c); #phi",750,-15,15,360,-TMath::Pi(),TMath::Pi());
  hPosTrkEtavsPhi = new TH2D("hPosTrkEtavsPhi","hPosTrkEtavsPhi;#phi;#eta",720,-TMath::Pi(),TMath::Pi(),200,-1,1);
  hNegTrkEtavsPhi = new TH2D("hNegTrkEtavsPhi","hNegTrkEtavsPhi;#phi;#eta",720,-TMath::Pi(),TMath::Pi(),200,-1,1);
  hDcavsPt        = new TH2D("hDcavsPt","hDcavsPt;charge*p_{T} (GeV/c);dca (cm)",300,-15,15,300,0,3);
  hdEdxvsP        = new TH2D("hdEdxvsP","hdEdxvsP; charge*p (GeV/c); dE/dx (KeV/cm)",300,-15,15,400,0,20);
  hdEdxvsPhi      = new TH2D("hdEdxvsPhi","hdEdxvsPhi;#phi;dEdx (KeV/cm)",360,-TMath::Pi(),TMath::Pi(),400,0,20);
  hdEdxvsEta      = new TH2D("hdEdxvsEta","hdEdxvsEta;#eta;dEdx (KeV/cm)",200,-1,1,400,0,20);
  hBetavsP        = new TH2D("hBetavsP","hBetavsP; charge*p (GeV/c); 1/#beta",300,-15,15,800,0,4);
  hMassvsP        = new TH2D("hMassvsP","hMassvsP; charge*p (GeV/c); mass(GeV/c^{2})",300,-15,15,100,0,2);
  
  hBEMCeEtavsPt   = new TH2D("hBEMCeEtavsPt","hBEMCeEtavsPt; charge*p_{T} (GeV/c); #eta",300,-15,15,200,-1,1);
  hBEMCePhivsPt   = new TH2D("hBEMCePhivsPt","hBEMCePhivsPt; charge*p_{T} (GeV/c); #phi",300,-15,15,360,-TMath::Pi(),TMath::Pi());
  hBEMCeEtavsPhi  = new TH2D("hBEMCeEtavsPhi","hBEMCeEtavsPhi;#phi;#eta",360,-TMath::Pi(),TMath::Pi(),200,-1,1);
  
  /***********   run by run QA   ***********/
  //event level QA
  hnEvtsvsRun = new TH1D("hnEvtsvsRun","hnEvtsvsRun;Run index;# of events",mTotalRuns, -0.5, mTotalRuns-0.5);
  hBFieldvsRun   = new TProfile("hBFieldvsRun","hBFieldvsRun;Run index;Magnetic field (kGauss)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hZdcXvsRun     = new TProfile("hZdcXvsRun","hZdcXvsRun;Run index;zdcRate (KHz)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hBbcXvsRun     = new TProfile("hBbcXvsRun","hBbcXvsRun;Run index;bbcRate (KHz)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hZdcXoverBbcXvsRun = new TProfile("hZdcXoverBbcXvsRun","hZdcXoverBbcXvsRun;Run index;zdcRate/bbcRate",mTotalRuns, -0.5, mTotalRuns-0.5);
  hTpcVxvsRun    = new TProfile("hTpcVxvsRun","hTpcVxvsRun;Run index;TPC V_{x} (cm)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hTpcVyvsRun    = new TProfile("hTpcVyvsRun","hTpcVyvsRun;Run index;TPC V_{y} (cm)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hTpcVzvsRun    = new TProfile("hTpcVzvsRun","hTpcVzvsRun;Run index;TPC V_{z} (cm)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hRawVpdVzvsRun = new TProfile("hRawVpdVzvsRun","hRawVpdVzvsRun;Run index;VPD V_{z} (cm)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hVpdVzvsRun    = new TProfile("hVpdVzvsRun","hVpdVzvsRun;Run index;VPD V_{z} (cm)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hVzDiffvsRun   = new TProfile("hVzDiffvsRun","hVzDiffvsRun;Run index; TPC V_{z} - VPD V_{z} (cm)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hRefMultvsRun  = new TProfile("hRefMultvsRun","hRefMultvsRun;Run index; refMult",mTotalRuns, -0.5, mTotalRuns-0.5);
  hgRefMultvsRun  = new TProfile("hgRefMultvsRun","hgRefMultvsRun;Run index; grefMult",mTotalRuns, -0.5, mTotalRuns-0.5);
  hTofMultvsRun  = new TProfile("hTofMultvsRun","hTofMultvsRun;Run index; tofMult",mTotalRuns, -0.5, mTotalRuns-0.5);
  hnTofMatchvsRun  = new TProfile("hnTofMatchvsRun","hnTofMatchvsRun;Run index; tofMatch",mTotalRuns, -0.5, mTotalRuns-0.5);
  hnBEMCMatchvsRun  = new TProfile("hnBEMCMatchvsRun","hnBEMCMatchvsRun;Run index; BEMCMatch",mTotalRuns, -0.5, mTotalRuns-0.5);
  
  // primary tracks
  hNTrksvsRun     = new TProfile("hNTrksvsRun","hNTrksvsRun;Run index;# of primary tracks",mTotalRuns, -0.5, mTotalRuns-0.5);
  hPtvsRun        = new TProfile("hPtvsRun","hPtvsRun;Run index;p_{T} (GeV/c)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hEtavsRun       = new TProfile("hEtavsRun","hEtavsRun;Run index;#eta",mTotalRuns, -0.5, mTotalRuns-0.5);
  hPhivsRun       = new TProfile("hPhivsRun","hPhivsRun;Run index;#phi (rad)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hDcavsRun       = new TProfile("hDcavsRun","hDcavsRun;Run index;dca (cm)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hNHitsFitvsRun  = new TProfile("hNHitsFitvsRun","hNHitsFitvsRun;Run index;nHitsFit",mTotalRuns, -0.5, mTotalRuns-0.5);
  hNHitsPossvsRun = new TProfile("hNHitsPossvsRun","hNHitsPossvsRun;Run index;nHitsPoss",mTotalRuns, -0.5, mTotalRuns-0.5);
  hNHitsDedxvsRun = new TProfile("hNHitsDedxvsRun","hNHitsDedxvsRun;Run index;nHitsDedx",mTotalRuns, -0.5, mTotalRuns-0.5);
  hDedxvsRun      = new TProfile("hDedxvsRun","hDedxvsRun;Run index;dE/dx (KeV/cm)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hNSigmaEvsRun   = new TProfile("hNSigmaEvsRun","hNSigmaEvsRun;Run index;n#sigma_{e}",mTotalRuns, -0.5, mTotalRuns-0.5);
  hNSigmaPivsRun  = new TProfile("hNSigmaPivsRun","hNSigmaPivsRun;Run index;n#sigma_{#pi}",mTotalRuns, -0.5, mTotalRuns-0.5);
  hNSigmaKvsRun   = new TProfile("hNSigmaKvsRun","hNSigmaKvsRun;Run index;n#sigma_{k}",mTotalRuns, -0.5, mTotalRuns-0.5);
  hNSigmaPvsRun   = new TProfile("hNSigmaPvsRun","hNSigmaPvsRun;Run index;n#sigma_{p}",mTotalRuns, -0.5, mTotalRuns-0.5);
  hBetavsRun      = new TProfile("hBetavsRun","hBetavsRun;Run index;1/#beta",mTotalRuns, -0.5, mTotalRuns-0.5);
  
  // BEMC match tracks
  hNMthTrksvsRun        = new TProfile("hNMthTrksvsRun","hNMthTrksvsRun;Run index;# of BEMC matched tracks",mTotalRuns, -0.5, mTotalRuns-0.5);
  hMthTrkPtvsRun        = new TProfile("hMthTrkPtvsRun","hMthTrkPtvsRun;Run index;p_{T} (GeV/c)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hMthTrkEtavsRun       = new TProfile("hMthTrkEtavsRun","hMthTrkEtavsRun;Run index;#eta",mTotalRuns, -0.5, mTotalRuns-0.5);
  hMthTrkPhivsRun       = new TProfile("hMthTrkPhivsRun","hMthTrkPhivsRun;Run index;#phi (rad)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hMthTrkNSigmaEvsRun   = new TProfile("hMthTrkNSigmaEvsRun","hMthTrkNSigmaEvsRun;Run index;n#sigma_{e}",mTotalRuns, -0.5, mTotalRuns-0.5);
  hMthTrkBetavsRun      = new TProfile("hMthTrkBetavsRun","hMthTrkBetavsRun;Run index;1/#beta",mTotalRuns, -0.5, mTotalRuns-0.5);
  hMthTrkAdc0vsRun      = new TProfile("hMthTrkAdc0vsRun","hMthTrkAdc0vsRun;Run index;adc0",mTotalRuns, -0.5, mTotalRuns-0.5);
  hMthTrkE0vsRun        = new TProfile("hMthTrkE0vsRun","hMthTrkE0vsRun;Run index;e0",mTotalRuns, -0.5, mTotalRuns-0.5);
  hMthTrkEvsRun         = new TProfile("hMthTrkEvsRun","hMthTrkEvsRun;Run index;e",mTotalRuns, -0.5, mTotalRuns-0.5);
  hMthTrkZDistvsRun     = new TProfile("hMthTrkZDistvsRun","hMthTrkZDistvsRun;Run index;#DeltaZ (cm)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hMthTrkPhiDistvsRun   = new TProfile("hMthTrkPhiDistvsRun","hMthTrkPhiDistvsRun;Run index;#Delta#phi (rad)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hMthTrkNEtavsRun      = new TProfile("hMthTrkNEtavsRun","hMthTrkNEtavsRun;Run index;n_{#eta}",mTotalRuns, -0.5, mTotalRuns-0.5);
  hMthTrkNPhivsRun      = new TProfile("hMthTrkNPhivsRun","hMthTrkNPhivsRun;Run index;n_{#phi}",mTotalRuns, -0.5, mTotalRuns-0.5);
  
  // BEMC trigger tracks
  hNTrigTrksvsRun        = new TProfile("hNTrigTrksvsRun","hNTrigTrksvsRun;Run index;# of BEMC triggered tracks",mTotalRuns, -0.5, mTotalRuns-0.5);
  hTrigTrkPtvsRun        = new TProfile("hTrigTrkPtvsRun","hTrigTrkPtvsRun;Run index;p_{T} (GeV/c)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hTrigTrkEtavsRun       = new TProfile("hTrigTrkEtavsRun","hTrigTrkEtavsRun;Run index;#eta",mTotalRuns, -0.5, mTotalRuns-0.5);
  hTrigTrkPhivsRun       = new TProfile("hTrigTrkPhivsRun","hTrigTrkPhivsRun;Run index;#phi (rad)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hTrigTrkNSigmaEvsRun   = new TProfile("hTrigTrkNSigmaEvsRun","hTrigTrkNSigmaEvsRun;Run index;n#sigma_{e}",mTotalRuns, -0.5, mTotalRuns-0.5);
  hTrigTrkAdc0vsRun      = new TProfile("hTrigTrkAdc0vsRun","hTrigTrkAdc0vsRun;Run index;adc0",mTotalRuns, -0.5, mTotalRuns-0.5);
  hTrigTrkE0vsRun        = new TProfile("hTrigTrkE0vsRun","hTrigTrkE0vsRun;Run index;e0",mTotalRuns, -0.5, mTotalRuns-0.5);
  hTrigTrkEvsRun         = new TProfile("hTrigTrkEvsRun","hTrigTrkEvsRun;Run index;e",mTotalRuns, -0.5, mTotalRuns-0.5);
  hTrigTrkZDistvsRun     = new TProfile("hTrigTrkZDistvsRun","hTrigTrkZDistvsRun;Run index;#DeltaZ (cm)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hTrigTrkPhiDistvsRun   = new TProfile("hTrigTrkPhiDistvsRun","hTrigTrkPhiDistvsRun;Run index;#Delta#phi (rad)",mTotalRuns, -0.5, mTotalRuns-0.5);
  hTrigTrkNEtavsRun      = new TProfile("hTrigTrkNEtavsRun","hTrigTrkNEtavsRun;Run index;n_{#eta}",mTotalRuns, -0.5, mTotalRuns-0.5);
  hTrigTrkNPhivsRun      = new TProfile("hTrigTrkNPhivsRun","hTrigTrkNPhivsRun;Run index;n_{#phi}",mTotalRuns, -0.5, mTotalRuns-0.5);
  
  // BEMC electron candidates
  hNBemcEsvsRun        = new TProfile("hNBemcEsvsRun","hNBemcEsvsRun;Run index;# BEMC electron candidates",mTotalRuns, -0.5, mTotalRuns-0.5);
  
  // MTD match tracks
  hMthTrkDeltaYvsRun = new TProfile("hMthTrkDeltaYvsRun", "hMthTrkDeltaYvsRun;Runindex;#DeltaY (cm)", mTotalRuns, -0.5, mTotalRuns-0.5);
  hMthTrkDeltaZvsRun = new TProfile("hMthTrkDeltaZvsRun", "hMthTrkDeltaZvsRun;Runindex;#DeltaZ (cm)", mTotalRuns, -0.5, mTotalRuns-0.5);
  hMthTrkDeltaTOFvsRun = new TProfile("hMthTrkDeltaTOFvsRun", "hMthTrkDeltaTOFvsRun;Runindex;#DeltaTOF", mTotalRuns, -0.5, mTotalRuns-0.5);
  
  mhNMtdHitsVsRun = new TProfile("mhNMtdHitsVsRun", "Number of MTD hits per event vs run; run index; N",mTotalRuns, -0.5, mTotalRuns-0.5);
  mhNMtdMthHitsVsRun = new TProfile("mhNMtdMthHitsVsRun", "Number of matched MTD hits per event vs run; run index; N",mTotalRuns, -0.5, mTotalRuns-0.5);
  mhMtdHitMap = new TH2F("mhMtdHitMap","Channel vs backleg of MTD hits;backleg;channel",30,0.5,30.5,60,0.5,60.5);
  mhMtdMthHitMap = new TH2F("mhMtdMthHitMap","Channel vs backleg of matched MTD hits;backleg;channel",30,0.5,30.5,60,0.5,60.5);
  
  
}

//_____________________________________________________________________________
void StQAMaker::printConfig()
{
  const char *decision[2] = {"no","yes"};
  printf("=== Configuration for StQAMaker ===\n");
  printf("Use default vertex: %s\n", decision[mDefaultVtx]);
  printf("Select positive vertex ranking: %s\n", decision[mSelectVtxRank]);
  printf("Maximum |Vr|: %1.2f\n", mMaxVtxR);
  printf("Maximum |Vz|: %1.2f\n", mMaxVtxZ);
  printf("Maximum |VzDiff|: %1.2f\n", mMaxVzDiff);
  printf("Minimum track pt: %1.2f\n", mMinTrkPt);
  printf("Maximum track |eta| : %1.2f\n", mMaxTrkEta);
  printf("Minimum number of fit hits: %d\n", mMinNHitsFit);
  printf("Minimum ratio of fit hits: %1.2f\n", mMinNHitsFitRatio);
  printf("Minimum number of dedx hits: %d\n", mMinNHitsDedx);
  printf("Maximum dca: %1.2f\n", mMaxDca);
  printf("Maximum |nSigmaE| for BEMCe: %1.2f\n", mMaxnSigmaE);
  printf("Maximum |1-1/beta| for BEMCe: %1.2f\n", mMaxBeta2TOF);
  printf("Minimum pt for BEMCe: %1.2f\n", mMinBemcPt);
  printf("Minimum adc0 for BEMCe: %d\n", mMinAdc0);
  printf("Minimum p/E for BEMCe: %1.2f\n", mMinPoverE);
  printf("Maximum p/E for BEMCe: %1.2f\n", mMaxPoverE);
  printf("Maximum |zDist| for BEMCe: %1.2f\n", mMaxZDist);
  printf("Maximum |phiDist| for BEMCe: %1.2f\n", mMaxPhiDist);
  printf("Minimum n_{#eta} for BEMCe: %d\n", mMinNEta);
  printf("Minimum n_{#phi} for BEMCe: %d\n", mMinNPhi);
  printf("=======================================\n");
}

float StQAMaker::getTofBeta(StPicoTrack const* const trk) const
{
  int index2tof = trk->bTofPidTraitsIndex();
  float beta = std::numeric_limits<float>::quiet_NaN();
  if (index2tof >= 0)
  {
    StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
    if (tofPid)
    {
      beta = tofPid->btofBeta();
      if (beta < 1e-4)
      {
        TVector3 const vtx3 = mevent->primaryVertex();
        StThreeVectorF vtx(vtx3.x(),vtx3.y(),vtx3.z());
        TVector3 const btofHitPos3 = tofPid->btofHitPos();
        // StThreeVectorF const btofHitPos = tofPid->btofHitPos();
        StThreeVectorF btofHitPos(btofHitPos3.x(),btofHitPos3.y(),btofHitPos3.z());
        StPicoPhysicalHelix helix = trk->helix(mevent->bField());
        float L = tofPathLength(&vtx, &btofHitPos, helix.curvature());
        float tof = tofPid->btof();
        if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
        else beta = std::numeric_limits<float>::quiet_NaN();
      }
    }
  }
  return beta;
}

Int_t StQAMaker::getMtdPidTraitsIndex(const StPicoMtdHit *hit)
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
