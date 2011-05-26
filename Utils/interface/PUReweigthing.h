//--------------------------------------------------------------------------------------------------
// $Id: 
//
// PUReweigthing
//
// Authors: M. Zanetti
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_PUREWEIGTHING_H
#define MITPHYSICS_UTILS_PUREWEIGTHING_H

#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include <TH1F.h>

namespace mithep {

  class PUReweigthing {

  public:
    static double reweightOOT(int inTimePUMultiplicity, int outOfTimePUMultiplicity, 
			      const char* fileName, const char* histoName);


    ClassDef(PUReweigthing, 0) // PUReweigthing
      };


}

#endif
