//--------------------------------------------------------------------------------------------------
// $Id: H4lZPlusFakeSkim.h,v 1.1 2012/06/03 20:25:56 paus Exp $
//
// H4lZPlusFakeSkim
//
// This module selects events with Z plus at least one fake electron or muon. It is used to 
// compute fake rates, among other things.
//
// Authors: D.Ralph, C.Paus
//--------------------------------------------------------------------------------------------------
#ifndef H4LZPLUSFAKESKIM_H
#define H4LZPLUSFAKESKIM_H

#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/TriggerMask.h"
#include "MitPhysics/Skim/interface/BaseH4lSkim.h"
#include "TLorentzVector.h"

namespace mithep
{  
  class H4lZPlusFakeSkim : public BaseH4lSkim
  {    
  public:
    H4lZPlusFakeSkim(const char *name="H4lZPlusFakeSkim",
		     const char *title="Z Plus Fake Skim");
    ~H4lZPlusFakeSkim();	    

  protected:
    void Begin();
    void BeginRun();
    void EndRun();
    void SlaveBegin();
    void SlaveTerminate();
    void Terminate();
    void Process();

    ClassDef(H4lZPlusFakeSkim,1)
  };
}

#endif
