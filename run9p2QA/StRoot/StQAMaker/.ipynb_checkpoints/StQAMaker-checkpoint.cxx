#include "headers.h"
#include "StQAMaker.h"

ClassImp(StQAMaker)

//_____________________________________________________________________________
StQAMaker::StQAMaker(const Char_t *name) : StMaker(name), mPrintConfig(1), mPrintMemory(0), mPrintCpu(0), mStreamName(""), fOutFile(0), mOutFileName(""), mEvtTree(0), mDefaultVtx(1), mSelectVtxRank(0), mMaxVtxR(2.), mMaxVtxZ(100.), mMaxVzDiff(3.), mMinTrkPt(0.2), mMaxTrkEta(1.), mMinNHitsFit(15), mMinNHitsFitRatio(0.52), mMinNHitsDedx(10), mMaxDca(3.), mMaxnSigmaE(2.0), mMaxBeta2TOF(0.03), mMinBemcPt(3.5), mMinAdc0(290), mMinPoverE(0.3), mMaxPoverE(1.5), mMaxZDist(3), mMaxPhiDist(0.02), mMinNEta(1), mMinNPhi(1)
{
	// default constructor

	// run11 st_hlt 
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
	mEmcPosition = new StEmcPosition();
	for(Int_t i=0;i<4;i++){
		if(i==1) continue;
		mEmcGeom[i] = StEmcGeom::getEmcGeom(detname[i].Data());
	}

	if(!mOutFileName.Length()){
		LOG_ERROR << "StQAMaker:: no output file specified for tree and histograms." << endm;
		return kStERR;
	}
	fOutFile = new TFile(mOutFileName.Data(),"recreate");
	LOG_INFO << "StQAMaker:: create the output file to store the tree and histograms: " << mOutFileName.Data() << endm;

	bookHistos();

	ifstream inData;
	inData.open("StRoot/StQAMaker/run11_stHT_runnumber_Jamie.dat");
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

	mMuDstMaker = (StMuDstMaker *)GetMaker("MuDst");
	if(Debug()){
		LOG_INFO<<"MuDstMaker pointer: "<<mMuDstMaker<<endm;
	}

	if(mMuDstMaker){
		if(Debug()) LOG_INFO<<"Use MuDst file as input"<<endm;
		mMuDst = mMuDstMaker->muDst();
		if(!mMuDst){
			LOG_WARN<<"No MuDst !"<<endm;
			return kStOK;
		}
	}
	else{
		LOG_WARN<<"No StMuDstMaker !"<<endm;
		return kStOK;
	}

	if(!processMuDstEvent()) return kStOK;

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
Bool_t StQAMaker::processMuDstEvent()
{
	hEvent_DefVtx->Fill(0.5);
	hEvent->Fill(0.5);

	StMuEvent *mMuEvent = mMuDst->event();
	if(!mMuEvent){
		LOG_WARN<<"No event level information !"<<endm;
		return kFALSE;
	}

	Bool_t NPE_18       = kFALSE;
	Bool_t NPE_25       = kFALSE;
	Bool_t NPE_25_nozdc = kFALSE;

	Bool_t validTrig = kFALSE;
	if(mTriggerIDs.size()==0){
		for(Int_t i=0;i<64;i++){
			Int_t trgId = mMuEvent->triggerIdCollection().nominal().triggerId(i);
			if(trgId>0){
				validTrig = kTRUE;

				if(350503==trgId || 350513==trgId) NPE_18       = kTRUE;
				if(350504==trgId || 350514==trgId) NPE_25       = kTRUE;
				if(350501==trgId || 350511==trgId) NPE_25_nozdc = kTRUE;
			}
		}
	}
	else{
		for(Int_t i=0;i<mTriggerIDs.size();i++){
			if(mMuEvent->triggerIdCollection().nominal().isTrigger(mTriggerIDs[i])){
				validTrig = kTRUE;

				if(350503==mTriggerIDs[i] || 350513==mTriggerIDs[i]) NPE_18       = kTRUE;
				if(350504==mTriggerIDs[i] || 350514==mTriggerIDs[i]) NPE_25       = kTRUE;
				if(350501==mTriggerIDs[i] || 350511==mTriggerIDs[i]) NPE_25_nozdc = kTRUE;
				//if(i<2)      NPE_18       = kTRUE;
				//else if(i<4) NPE_25       = kTRUE;
				//else         NPE_25_nozdc = kTRUE;
			}
		}
	}

	if(!validTrig){
		if(Debug()) LOG_WARN<<"No valid interested triggers !"<<endm;
		return kFALSE;
	}

	Int_t runIdx;
	Int_t runId = mMuEvent->runId();
	map<Int_t,Int_t>::iterator iter = mTotalRunId.find(runId);
	if(iter != mTotalRunId.end())
		runIdx = iter->second;
	else{
		runIdx = -1;
		LOG_INFO<<"Can not find the runNumber in the runNumber list"<<endm;
	}

	if(NPE_18)       { hEvent_DefVtx->Fill(2.5); hEvent->Fill(2.5); hnHT2EvtsvsRun->Fill(runIdx); }
	if(NPE_25)       { hEvent_DefVtx->Fill(3.5); hEvent->Fill(3.5); }
	if(NPE_25_nozdc) { hEvent_DefVtx->Fill(4.5); hEvent->Fill(4.5); }

	Double_t vpdVz = -999; 
	StBTofHeader *mBTofHeader = mMuDst->btofHeader();
	if(mBTofHeader) vpdVz = mBTofHeader->vpdVz();

	StThreeVectorF vtxPos = mMuEvent->primaryVertexPosition();
	Double_t vx      = vtxPos.x();
	Double_t vy      = vtxPos.y();
	Double_t vz      = vtxPos.z();
	Double_t vr      = sqrt(vx*vx + vy*vy);
	Double_t vzDiff  = vz - vpdVz;
	Double_t refMult = mMuEvent->refMult();

	StMuPrimaryVertex *vtx = mMuDst->primaryVertex();
	if(vtx){
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
						if(NPE_18)       hEvent_DefVtx->Fill(12.5);
						if(NPE_25)       hEvent_DefVtx->Fill(13.5);
						if(NPE_25_nozdc) hEvent_DefVtx->Fill(14.5);
					}
				}
			}
		}
	}

	//select the right vertex using VPD
	Int_t vtxIdx = 0;
	for(UInt_t i=0;!mDefaultVtx&&i<mMuDst->numberOfPrimaryVertices();i++){
		vtx = mMuDst->primaryVertex(i);
		if(!vtx) continue;
		Float_t vz = vtx->position().z();
		if(fabs(vpdVz)<200 && fabs(vpdVz-vz)<3.){
			mMuDst->setVertexIndex(i); 
			vtxIdx = i;
			break;
		}
	}

	vtx = mMuDst->primaryVertex();
	if(!vtx) return kFALSE;

	vtxPos  = mMuEvent->primaryVertexPosition();
	vx      = vtxPos.x();
	vy      = vtxPos.y();
	vz      = vtxPos.z();
	vr      = sqrt(vx*vx + vy*vy);
	vzDiff  = vz - vpdVz;
	refMult = mMuEvent->refMult();

	hVyvsVx->Fill(vx, vy);
	hVpdVzvsTpcVz->Fill(vz, vpdVz);
	hTpcVzvsRefMult->Fill(refMult, vz);
	hRawVpdVzvsRefMult->Fill(refMult, vpdVz);
	if(TMath::Abs(vz)<=mMaxVtxZ) hVpdVzvsRefMult->Fill(refMult, vpdVz);
	hVzDiffvsTpcVz->Fill(vz, vzDiff);
	hVzDiffvsRefMult->Fill(refMult, vzDiff);

	if(Debug()){
		LOG_INFO<<"vtxIdx: "<<vtxIdx<<" \tTPC Vx: "<<vx<<" \tTPC Vy: "<<vy<<" \tTPC Vz: "<<vz<<endm;
	}

	Double_t bField  = mMuEvent->runInfo().magneticField();
	Double_t bbcRate = mMuEvent->runInfo().bbcCoincidenceRate();
	Double_t zdcRate = mMuEvent->runInfo().zdcCoincidenceRate();

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

	hVtxIdxvsRefMult->Fill(refMult, vtxIdx);
	if(TMath::Abs(vzDiff)<3) hVtxIdxvsRefMult_VzDiffCut->Fill(refMult, vtxIdx);

	Double_t vtxRanking  = mMuDst->primaryVertex()->ranking();
	if(mSelectVtxRank && vtxRanking<=0) return kFALSE;
	//hEvent->Fill(6.5);
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

	if(NPE_18)       hEvent->Fill(12.5);
	if(NPE_25)       hEvent->Fill(13.5);
	if(NPE_25_nozdc) hEvent->Fill(14.5);

	Int_t nNodes = mMuDst->numberOfPrimaryTracks();
	if(Debug()){
		LOG_INFO<<"# of primary Tracks in muDst: "<<nNodes<<endm;
	}

	Int_t nTrks      = 0;
	Int_t nMthTrks   = 0;
	Int_t nTrigTrks  = 0;
	Int_t nBEMCeCans = 0;
	for(Int_t i=0;i<nNodes;i++){
		StMuTrack* pMuTrack = mMuDst->primaryTracks(i);
		if(!pMuTrack) continue;

		if(!isValidTrack(pMuTrack)) continue;

		StThreeVectorF pMom = pMuTrack->p();
		Double_t pt         = pMom.perp();
		Double_t eta        = pMom.pseudoRapidity();
		Double_t phi        = pMom.phi();
		Int_t    charge     = pMuTrack->charge();
		Double_t nHitsFit   = pMuTrack->nHitsFit(kTpcId);
		Double_t nHitsPoss  = pMuTrack->nHitsPoss(kTpcId);
		Double_t nHitsDedx  = pMuTrack->nHitsDedx();
		Double_t dca        = pMuTrack->dcaGlobal().mag();
		Double_t dEdx       = pMuTrack->dEdx()*1.e6; 
		Double_t nSigmaE    = pMuTrack->nSigmaElectron();
		Double_t nSigmaPi   = pMuTrack->nSigmaPion();
		Double_t nSigmaK    = pMuTrack->nSigmaKaon();
		Double_t nSigmaP    = pMuTrack->nSigmaProton();

		hNHitsFit->Fill(charge*nHitsFit);
		hNHitsPoss->Fill(charge*nHitsPoss);
		hNHitsDedx->Fill(charge*nHitsDedx);
		hEtavsPt->Fill(charge*pt, eta);
		if(bField>0) hFFPhivsPt->Fill(charge*pt, phi);
		else         hRFFPhivsPt->Fill(charge*pt, phi);
		if(charge>0) hPosTrkEtavsPhi->Fill(phi, eta);
		else         hNegTrkEtavsPhi->Fill(phi, eta);
		hDcavsPt->Fill(charge*pt, dca);
		hdEdxvsP->Fill(charge*pMom.mag(), dEdx);
		hdEdxvsPhi->Fill(phi, dEdx);
		hdEdxvsEta->Fill(eta, dEdx);

		Double_t beta       = -999.;
		if( &(pMuTrack->btofPidTraits()) ){
			const StMuBTofPidTraits& btofPidTraits = pMuTrack->btofPidTraits();
			beta = btofPidTraits.beta();

			hBetavsP->Fill(charge*pMom.mag(), 1/beta);
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

		getBemcInfo(pMuTrack, runIdx, nMthTrks, nTrigTrks, nBEMCeCans);

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
Bool_t StQAMaker::isValidTrack(StMuTrack *pMuTrack) const
{
	Float_t pt  = pMuTrack->pt();
	Float_t eta = pMuTrack->eta(); 
	Float_t dca = pMuTrack->dcaGlobal().mag();

	if(pt<mMinTrkPt)                            return kFALSE;
	if(TMath::Abs(eta)>mMaxTrkEta)              return kFALSE;
	if(pMuTrack->nHitsFit(kTpcId)<mMinNHitsFit) return kFALSE;
	if(pMuTrack->nHitsFit(kTpcId)*1./pMuTrack->nHitsPoss(kTpcId)<mMinNHitsFitRatio)  return kFALSE;
	if(pMuTrack->nHitsDedx()<mMinNHitsDedx)     return kFALSE;
	if(dca>mMaxDca)                             return kFALSE;

	return kTRUE;
}

//_____________________________________________________________________________
Bool_t StQAMaker::getBemcInfo(StMuTrack *pMuTrack, const Int_t runIdx, Int_t &nMthTrks, Int_t &nTrigTrks, Int_t &nBEMCeCans)
{
	Float_t maxtowerE = -999., energy = 0.;
	Float_t zdist = -999., phidist = -999., mindist = 999.;
	Int_t mod = -1, eta=-1, sub=-1;
	Int_t neta = -1, nphi=-1;
	UInt_t maxadc = 0;

	mEmcCollection = mMuDst->emcCollection();
	if(!mEmcCollection) {
		LOG_WARN << " No Emc Collection for this event " << endm;
		return kFALSE;
	}

	StThreeVectorD position, momentum;
	StThreeVectorD positionBSMDE, momentumBSMDE;
	StThreeVectorD positionBSMDP, momentumBSMDP;

	Double_t bField = mMuDst->event()->runInfo().magneticField()/10.; //Tesla
	Bool_t ok       = kFALSE;
	Bool_t okBSMDE  = kFALSE;
	Bool_t okBSMDP  = kFALSE;
	if(mEmcPosition) {
		ok      = mEmcPosition->projTrack(&position, &momentum, pMuTrack, bField, mEmcGeom[0]->Radius());
		okBSMDE = mEmcPosition->projTrack(&positionBSMDE, &momentumBSMDE, pMuTrack, bField, mEmcGeom[2]->Radius());
		okBSMDP = mEmcPosition->projTrack(&positionBSMDP, &momentumBSMDP, pMuTrack, bField, mEmcGeom[3]->Radius());
	}
	if(!ok) {
		LOG_WARN << " Projection failed for this track ... " << endm;
		return kFALSE;
	}

	Bool_t bemcMatchFlag = kFALSE;
	if(ok && okBSMDE && okBSMDP){

		StSPtrVecEmcPoint& bEmcPoints = mEmcCollection->barrelPoints();
		mindist=1.e9;
		mEmcGeom[0]->getBin(positionBSMDP.phi(), positionBSMDE.pseudoRapidity(), mod, eta, sub); //project on SMD plan
		for(StSPtrVecEmcPointIterator it = bEmcPoints.begin(); it != bEmcPoints.end(); it++) {
			Bool_t associatedPoint = kFALSE;
			StPtrVecEmcCluster& bEmcClusters = (*it)->cluster(kBarrelEmcTowerId);
			if(bEmcClusters.size()==0 ) continue;
			if(bEmcClusters[0]==NULL) continue;
			for(StPtrVecEmcClusterIterator cIter = bEmcClusters.begin(); cIter != bEmcClusters.end(); cIter++){
				Bool_t associatedCluster = kFALSE;
				StPtrVecEmcRawHit& bEmcHits = (*cIter)->hit();
				for(StPtrVecEmcRawHitIterator hIter = bEmcHits.begin(); hIter != bEmcHits.end(); hIter++) {
					if(mod == (Int_t)(*hIter)->module() && eta == (Int_t)(*hIter)->eta() && sub == (Int_t)(*hIter)->sub()) {
						bemcMatchFlag = kTRUE;
						associatedPoint = kTRUE;
						associatedCluster = kTRUE;
						break;
					}
				}
				if(associatedCluster) {
					for(StPtrVecEmcRawHitIterator hitit=bEmcHits.begin(); hitit!=bEmcHits.end();hitit++) {
						if((*hitit)->energy()>maxtowerE) maxtowerE = (*hitit)->energy();
						if((*hitit)->adc()>maxadc) maxadc = (*hitit)->adc();
					}
				}
			}

			StPtrVecEmcCluster& smdeClusters = (*it)->cluster(kBarrelSmdEtaStripId);
			StPtrVecEmcCluster& smdpClusters = (*it)->cluster(kBarrelSmdPhiStripId);

			if(associatedPoint) {
				energy += (*it)->energy(); //use point's energy, not tower cluster's energy

				float deltaphi=(*it)->position().phi()-positionBSMDP.phi();
				if(deltaphi>=TMath::Pi()) deltaphi=deltaphi-TMath::TwoPi();
				if(deltaphi<-TMath::Pi()) deltaphi=deltaphi+TMath::TwoPi();

				float rsmdp=mEmcGeom[3]->Radius();
				float pointz=(*it)->position().z();
				float deltaz=pointz-positionBSMDE.z();
				if(sqrt(deltaphi*deltaphi*rsmdp*rsmdp+deltaz*deltaz)<mindist) {
					phidist=deltaphi;
					zdist  =deltaz;
					if(smdeClusters.size()>=1) neta=smdeClusters[0]->nHits();
					if(smdpClusters.size()>=1) nphi=smdpClusters[0]->nHits();
					mindist=sqrt(deltaphi*deltaphi*rsmdp*rsmdp+deltaz*deltaz);
				}
			}//associated
		}
	} // end if (ok && okBSMDE && okBSMDP)

	if(!bemcMatchFlag) return kFALSE;
	nMthTrks++;

	Int_t    charge     = pMuTrack->charge();
	Double_t trkp       = pMuTrack->p().mag();
	Double_t trkpt      = pMuTrack->pt();
	Double_t trketa     = pMuTrack->eta();
	Double_t trkphi     = pMuTrack->phi();
	Double_t nSigmaE    = pMuTrack->nSigmaElectron();
	Double_t beta       = -999.;
	if( &(pMuTrack->btofPidTraits()) ){
		const StMuBTofPidTraits& btofPidTraits = pMuTrack->btofPidTraits();
		beta = btofPidTraits.beta();
	}

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
void StQAMaker::bookHistos()
{
	// default vertex
	hEvent_DefVtx = new TH1D("hEvent_DefVtx","Event statistics",20,0,20);
	//hEvent_DefVtx->GetXaxis()->SetBinLabel(1, "All events");
	//hEvent_DefVtx->GetXaxis()->SetBinLabel(3, "NPE_18");
	//hEvent_DefVtx->GetXaxis()->SetBinLabel(4, "NPE_25");
	//hEvent_DefVtx->GetXaxis()->SetBinLabel(5, "NPE_25_nozdc");
	//hEvent_DefVtx->GetXaxis()->SetBinLabel(7, "mRanking>0");
	//hEvent_DefVtx->GetXaxis()->SetBinLabel(8, "None-Zero Vertex");
	//hEvent_DefVtx->GetXaxis()->SetBinLabel(9,  Form("|V_{r}|<%1.2f cm",mMaxVtxR));
	//hEvent_DefVtx->GetXaxis()->SetBinLabel(10, Form("|V_{z}|<%1.2f cm",mMaxVtxZ));
	//hEvent_DefVtx->GetXaxis()->SetBinLabel(11, Form("|V_{z}Diff|<%1.2f cm",mMaxVzDiff));
	//hEvent_DefVtx->GetXaxis()->SetBinLabel(13, "NPE_18");
	//hEvent_DefVtx->GetXaxis()->SetBinLabel(14, "NPE_25");
	//hEvent_DefVtx->GetXaxis()->SetBinLabel(15, "NPE_25_nozdc");
	hVyvsVx_DefVtx       = new TH2D("hVyvsVx_DefVtx","hVyvsVx_DefVtx; V_{x} (cm); V_{y} (cm)",600,-3,3,600,-3,3); 
	hVpdVzvsTpcVz_DefVtx = new TH2D("hVpdVzvsTpcVz_DefVtx","hVpdVzvsTpcVz_DefVtx; TPC V_{z} (cm); VPD V_{z} (cm)",400,-200,200,400,-200,200);
	hVzDiffvsTpcVz_DefVtx   = new TH2D("hVzDiffvsTpcVz_DefVtx","hVzDiffvsTpcVz_DefVtx; TPC Vz (cm); Vz_{TPC} - Vz_{VPD} (cm)",400,-200,200,400,-20,20);
	hVzDiffvsRefMult_DefVtx = new TH2D("hVzDiffvsRefMult_DefVtx","hVzDiffvsRefMult_DefVtx; refMult; Vz_{TPC} - Vz_{VPD} (cm)",800,0,800,400,-20,20);
	hRefMult_VzVrCut_DefVtx = new TH1D("hRefMult_VzVrCut_DefVtx","hRefMult_VzVrCut_DefVtx; refMult",800,0,800);
	hRefMult_EvtCut_DefVtx  = new TH1D("hRefMult_EvtCut_DefVtx","hRefMult_EvtCut_DefVtx; refMult",800,0,800);

	// selected vertex
	hEvent = new TH1D("hEvent","Event statistics",20,0,20);
	//hEvent->GetXaxis()->SetBinLabel(1, "All events");
	//hEvent->GetXaxis()->SetBinLabel(3, "NPE_18");
	//hEvent->GetXaxis()->SetBinLabel(4, "NPE_25");
	//hEvent->GetXaxis()->SetBinLabel(5, "NPE_25_nozdc");
	//hEvent->GetXaxis()->SetBinLabel(7, "mRanking>0");
	//hEvent->GetXaxis()->SetBinLabel(8, "None-Zero Vertex");
	//hEvent->GetXaxis()->SetBinLabel(9,  Form("|V_{r}|<%1.2f cm",mMaxVtxR));
	//hEvent->GetXaxis()->SetBinLabel(10, Form("|V_{z}|<%1.2f cm",mMaxVtxZ));
	//hEvent->GetXaxis()->SetBinLabel(11, Form("|V_{z}Diff|<%1.2f cm",mMaxVzDiff));
	//hEvent->GetXaxis()->SetBinLabel(13, "NPE_18");
	//hEvent->GetXaxis()->SetBinLabel(14, "NPE_25");
	//hEvent->GetXaxis()->SetBinLabel(15, "NPE_25_nozdc");
	hVyvsVx            = new TH2D("hVyvsVx","hVyvsVx; V_{x} (cm); V_{y} (cm)",600,-3,3,600,-3,3); 
	hVpdVzvsTpcVz      = new TH2D("hVpdVzvsTpcVz","hVpdVzvsTpcVz; TPC V_{z} (cm); VPD V_{z} (cm)",400,-200,200,400,-200,200);
	hTpcVzvsRefMult    = new TH2D("hTpcVzvsRefMult","hTpcVzvsRefMult; refMult; TPC V_{z} (cm)",800,0,800,600,-600,600);
	hRawVpdVzvsRefMult = new TH2D("hRawVpdVzvsRefMult","hRawVpdVzvsRefMult; refMult; VPD V_{z} (cm)",800,0,800,600,-600,600);
	hVpdVzvsRefMult    = new TH2D("hVpdVzvsRefMult","hVpdVzvsRefMult; refMult; VPD V_{z} (cm)",800,0,800,600,-600,600);
	hVzDiffvsTpcVz     = new TH2D("hVzDiffvsTpcVz","hVzDiffvsTpcVz; TPC Vz (cm); Vz_{TPC} - Vz_{VPD} (cm)",400,-200,200,400,-20,20);
	hVzDiffvsRefMult   = new TH2D("hVzDiffvsRefMult","hVzDiffvsRefMult; refMult; Vz_{TPC} - Vz_{VPD} (cm)",800,0,800,400,-20,20);
	hRefMult_VzVrCut   = new TH1D("hRefMult_VzVrCut","hRefMult_VzVrCut; refMult",800,0,800);
	hRefMult_EvtCut    = new TH1D("hRefMult_EvtCut","hRefMult_EvtCut; refMult",800,0,800);
	hVtxIdxvsRefMult             = new TH2D("hVtxIdxvsRefMult","hVtxIdxvsRefMult; refMult; Vertex Index",800,0,800,20,0,20);
	hVtxIdxvsRefMult_VzDiffCut = new TH2D("hVtxIdxvsRefMult_VzDiffCut","hVtxIdxvsRefMult_VzDiffCut; refMult; Vertex Index",800,0,800,20,0,20);
	hRefMultvsZdcX     = new TH2D("hRefMultvsZdcX","hRefMultvsZdcX; zdcRate (kHz); refMult",600,0,60,800,0,800);
	hRefMultvsBbcX     = new TH2D("hRefMultvsBbcX","hRefMultvsBbcX; bbcRate (kHz); refMult",600,0,60,800,0,800);

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
	hBEMCeEtavsPt   = new TH2D("hBEMCeEtavsPt","hBEMCeEtavsPt; charge*p_{T} (GeV/c); #eta",300,-15,15,200,-1,1);
	hBEMCePhivsPt   = new TH2D("hBEMCePhivsPt","hBEMCePhivsPt; charge*p_{T} (GeV/c); #phi",300,-15,15,360,-TMath::Pi(),TMath::Pi());
	hBEMCeEtavsPhi  = new TH2D("hBEMCeEtavsPhi","hBEMCeEtavsPhi;#phi;#eta",360,-TMath::Pi(),TMath::Pi(),200,-1,1);

	/***********   run by run QA   ***********/
	//event level QA
	hnHT2EvtsvsRun = new TH1D("hnHT2EvtsvsRun","hnHT2EvtsvsRun;Run index;# of HT2 events",mTotalRuns, -0.5, mTotalRuns-0.5);
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
