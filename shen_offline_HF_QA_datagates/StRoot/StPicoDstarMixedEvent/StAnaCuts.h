#ifndef StAnaCuts_H
#define StAnaCuts_H

/* **************************************************
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include "Rtypes.h"
#include <string>
#include <array>

namespace anaCuts
{
   std::array<unsigned int, 6> const triggers = {
   610001,610011,610021,610031,610041,610051
   };    
   //cut before QA
   float const qavz = 200.0;// < cm.
   float const qaVerror = 1.0e-5; //
   float const qaVr = 2.0; //cm
   float const qavzVpdVz = 3; //cm

   // QA tracks cuts
   float const qaGPt = 0.15;
   int const qaNHitsFit = 15;
   int const qaNHitsDedx = 15;
   float const qaDca = 3;// < cm
   float const qaEta = 1.5; 
   float const qaTofPion=4;
   float const qaTpcPion=4;

   //cut 
   float const vz = 200.0;// < cm.
   float const Verror = 1.0e-5; //
   float const Vr = 2.0; //cm
   // float const vzVpdVz = 3; //cm
   float const vzVpdVz = 1e3; //cm

   // QA tracks cuts
   // float const GPt = 0.20;
   float const GPt = 0.15;
   int const NHitsFit = 15;
   int const NHitsDedx = 15;
   int const NHitsFit2Poss = 0.52;
   // float const Dca = 1;// < cm
   float const Dca = 1.5;// < cm
   float const Eta = 1.5; 
   float const TofPion=4;
   float const TpcPion=4;
}
#endif
