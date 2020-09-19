#ifndef StAnaTreeArrays_hh
#define StAnaTreeArrays_hh

enum anaTreeTypes {anaTreeEvent=0, anaTreeETrack, anaTreePartETrack, anaTreeMuTrack, anaTreeHTrack, 
    anaTreeEEPair, anaTreePhoEEPair, anaTreeEMuPair, anaTreeMuMuPair, anaTreeEmcTrigger, anaTreeMtdTrigger};
enum NAnaTreeArrays {
__NANATREEARRAYS__ = 11
};

class StAnaTreeArrays {
  public:
  StAnaTreeArrays();
///< names of the TBranches in the TTree/File 
  static const char*   anaTreeArrayNames[__NANATREEARRAYS__];

///< names of the classes, the TClonesArrays are arrays of this type
  static const char*   anaTreeArrayTypes[__NANATREEARRAYS__    ];

///< maximum sizes of the TClonesArrays
  static int           anaTreeArraySizes[__NANATREEARRAYS__    ];

///< number of entries in current event, currently not used
  static int           anaTreeArrayCounters[__NANATREEARRAYS__    ];

};

#endif
