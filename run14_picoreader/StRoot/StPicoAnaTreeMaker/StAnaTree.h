#ifndef StAnaTree_h
#define StAnaTree_h

class StPicoAnaTreeMaker;
class StEventHeader;
class StElectronTrack;
class StPartElectronTrack;
class StMuonTrack;
class StHadronTrack;
class StEEPair;
class StPhoEEPair;
class StEMuPair;
class StMuMuPair;
class StEmcTrigger;
class StMtdTrigger;

#include "TObject.h"
#include "TClonesArray.h"
#include "StAnaTreeArrays.h"

class StAnaTree : public TObject {
public:
  /// constructor
  StAnaTree();
  /// set the pointers to the TClonesArrays
  static void set(StPicoAnaTreeMaker* maker);
  /// set the pointers to the TClonesArrays
  static void set(TClonesArray**);
  /// resets the pointers to the TClonesArrays to 0
  static void unset();

protected:
  /// array of TClonesArrays
  static TClonesArray** anaTreeArrays;
  /// array of TClonesArrays for the stuff inherited from the StPicoV0

public:
  /// returns pointer to the n-th TClonesArray 
  static TClonesArray* anaTreeArray(int type) { return anaTreeArrays[type]; }

  /// returns pointer to current StPicoEvent (class holding the event wise information)
  static StEventHeader* event() { return (StEventHeader*)anaTreeArrays[anaTreeEvent]->UncheckedAt(0); }

  /// return pointer to i-th track 
  static StElectronTrack* eTrack(int i) { return (StElectronTrack*)anaTreeArrays[anaTreeETrack]->UncheckedAt(i); }
  static StPartElectronTrack* partETrack(int i) { return (StPartElectronTrack*)anaTreeArrays[anaTreePartETrack]->UncheckedAt(i); }
  static StMuonTrack* muTrack(int i) { return (StMuonTrack*)anaTreeArrays[anaTreeMuTrack]->UncheckedAt(i); }
  static StHadronTrack* hTrack(int i) { return (StHadronTrack*)anaTreeArrays[anaTreeHTrack]->UncheckedAt(i); }

  /// return pointer to i-th pair 
  static StEEPair* eePair(int i) { return (StEEPair*)anaTreeArrays[anaTreeEEPair]->UncheckedAt(i); }
  static StPhoEEPair* phoEEPair(int i) { return (StPhoEEPair*)anaTreeArrays[anaTreePhoEEPair]->UncheckedAt(i); }
  static StEMuPair* eMuPair(int i) { return (StEMuPair*)anaTreeArrays[anaTreeEMuPair]->UncheckedAt(i); }
  static StMuMuPair* MuMuPair(int i) { return (StMuMuPair*)anaTreeArrays[anaTreeMuMuPair]->UncheckedAt(i); }
  static StEmcTrigger* emcTrigger(int i) { return (StEmcTrigger*)anaTreeArrays[anaTreeEmcTrigger]->UncheckedAt(i); }
  static StMtdTrigger* mtdTrigger(int i) { return (StMtdTrigger*)anaTreeArrays[anaTreeMtdTrigger]->UncheckedAt(i); }

  static unsigned int numberOfETracks() { return anaTreeArrays[anaTreeETrack]->GetEntries(); }
  static unsigned int numberOfMuTracks() { return anaTreeArrays[anaTreeMuTrack]->GetEntries(); }
  static unsigned int numberOfHTracks() { return anaTreeArrays[anaTreeHTrack]->GetEntries(); }
  static unsigned int numberOfEEPairs() {return anaTreeArrays[anaTreeEEPair]->GetEntries(); }
  static unsigned int numberOfPhoEEPairs() {return anaTreeArrays[anaTreePhoEEPair]->GetEntries(); }
  static unsigned int numberOfMuMuPairs() {return anaTreeArrays[anaTreeMuMuPair]->GetEntries(); }
  static unsigned int numberOfEMuPairs() {return anaTreeArrays[anaTreeEMuPair]->GetEntries(); }
  static unsigned int numberOfEmcTriggers() {return anaTreeArrays[anaTreeEmcTrigger]->GetEntries(); }
 
  virtual void Print(Option_t *option = "") const; ///< Print basic event info
  static void printTriggers();
  static void printETracks();
  static void printMuTracks();
       
  friend class StPicoAnaTreeMaker;

  ClassDef(StAnaTree,1)
};

#endif
