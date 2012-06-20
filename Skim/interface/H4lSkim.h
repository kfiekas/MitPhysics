//--------------------------------------------------------------------------------------------------
// $Id: H4lSkim.h,v 1.2 2012/06/03 20:25:56 paus Exp $
//
// H4lSkim
//
// This module selects Higgs to ZZ to 4 lepton events for skimming purposes.
//
// Authors: D.Ralph                   (but C.Paus stole it, changed the name and wrote the doc :-) )
//--------------------------------------------------------------------------------------------------
#ifndef H4LSKIM_H
#define H4LSKIM_H

#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityFwd.h"
#include "MitAna/DataTree/interface/BaseVertex.h"
#include "MitAna/DataTree/interface/TriggerMask.h"
#include "MitPhysics/Skim/interface/BaseH4lSkim.h"

namespace mithep
{  
  class H4lSkim : public BaseH4lSkim
  {    
  public:
    H4lSkim(const char *name="H4lSkim", const char *title="Four Lepton Skim");
    ~H4lSkim();	    
    
  protected:
    void Begin();
    void BeginRun();
    void EndRun();
    void SlaveBegin();
    void SlaveTerminate();
    void Terminate();
    void Process();
    
    ClassDef(H4lSkim,1)
  };
}
#endif
