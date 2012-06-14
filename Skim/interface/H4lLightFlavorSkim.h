//--------------------------------------------------------------------------------------------------
// $Id: H4lLightFlavorSkim.h,v 1.1 2012/06/03 20:25:56 paus Exp $
//
// H4lLightFlavorSkim
//
// This module selects events to control the light flavor background in the Higgs to ZZ to 4 lepton
// analysis.
//
// Authors: D.Ralph, C.Paus
//--------------------------------------------------------------------------------------------------
#ifndef H4LLIGHTFLAVORSKIM_H
#define H4LLIGHTFLAVORSKIM_H

#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitPhysics/Skim/interface/BaseH4lSkim.h"

namespace mithep
{
  class H4lLightFlavorSkim : public BaseH4lSkim
  {
  public:
    H4lLightFlavorSkim(const char *name  = "H4lLightFlavorSkim",
		       const char *title = "Light Flavor Control Skim");
    ~H4lLightFlavorSkim();
    
  protected:
    void Begin();
    void BeginRun();
    void EndRun();
    void SlaveBegin();
    void SlaveTerminate();
    void Terminate();
    void Process();

    ClassDef(H4lLightFlavorSkim,1)
  };
}
#endif
