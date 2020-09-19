#include <map>
#include "StAnaTree.h"
#include "StEventHeader.h"
#include "StElectronTrack.h"
#include "StPartElectronTrack.h"
#include "StMuonTrack.h"
#include "StHadronTrack.h"
#include "StEEPair.h"
#include "StPhoEEPair.h"
#include "StEMuPair.h"
#include "StMuMuPair.h"
#include "StEmcTrigger.h"
#include "StMtdTrigger.h"
#include "StPicoAnaTreeMaker.h"
#include "StMessMgr.h"

TClonesArray** StAnaTree::anaTreeArrays       = 0;

//-----------------------------------------------------------------------
StAnaTree::StAnaTree() {
  /* no-op */
}
//-----------------------------------------------------------------------
void StAnaTree::unset() {
    anaTreeArrays   = 0;
}
//-----------------------------------------------------------------------
void StAnaTree::set(StPicoAnaTreeMaker* maker) {
  if (!maker) { LOG_WARN << "No PicoDstMaker!" << endm; return;}
  anaTreeArrays   = maker->mAnaTreeArrays;
}
//-----------------------------------------------------------------------
void StAnaTree::set(TClonesArray** theAnaTreeArrays)
{
  anaTreeArrays    = theAnaTreeArrays;
}
//-----------------------------------------------------------------------
void StAnaTree::Print(Option_t *option) const {
  LOG_INFO << endm;
  LOG_INFO << "=========== event header =============" << endm << endm;
  LOG_INFO << " run/event Id = " << event()->runId() << "/" << event()->eventId() << endm;
  LOG_INFO << "=====================================" << endm << endm;
}

//-----------------------------------------------------------------------
void StAnaTree::printTriggers()  {

  LOG_INFO << "=== triggers ===" << endm;
  LOG_INFO << endm;
  LOG_INFO << "+++++++++ trigger list ( " <<  " entries )" << endm << endm;
}
//-----------------------------------------------------------------------
void StAnaTree::printETracks()  {
  if (numberOfETracks() == 0) {
    LOG_INFO << "No electron found!" << endm;
    return;
  }
  LOG_INFO << endm;
  LOG_INFO << "+++++++++ electron list ( " << numberOfETracks() << " entries )" << endm << endm;
  for (UInt_t i_trk = 0; i_trk < numberOfETracks(); i_trk++) {      
    LOG_INFO << "+++ electron " << i_trk << endm;
    eTrack(i_trk)->Print();   
    LOG_INFO << endm;
  }
}
//-----------------------------------------------------------------------
void StAnaTree::printMuTracks()  {
  if (numberOfMuTracks() == 0) {
    LOG_INFO << "No muon found!" << endm;
    return;
  }
  LOG_INFO << endm;
  LOG_INFO << "+++++++++ muon list ( " << numberOfMuTracks() << " entries )" << endm << endm;
  for (UInt_t i_trk = 0; i_trk < numberOfMuTracks(); i_trk++) {      
    LOG_INFO << "+++ muon " << i_trk << endm;
    muTrack(i_trk)->Print();   
    LOG_INFO << endm;
  }
}
ClassImp(StAnaTree)
