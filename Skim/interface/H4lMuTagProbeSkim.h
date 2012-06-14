#ifndef H4LMUTAGPROBESKIM_H
#define H4LMUTAGPROBESKIM_H

#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/TriggerMask.h"
#include "MitPhysics/Skim/interface/BaseH4lSkim.h"
#include "TLorentzVector.h"

namespace mithep
{  
  class H4lMuTagProbeSkim : public BaseH4lSkim
  {    
  public:
    H4lMuTagProbeSkim(const char *name="H4lMuTagProbeSkim", const char *title="Muon Tag and Probe Skim");
    ~H4lMuTagProbeSkim();	    

  protected:
    void Begin();
    void BeginRun();
    void EndRun();
    void SlaveBegin();
    void SlaveTerminate();
    void Terminate();
    void Process();

    ClassDef(H4lMuTagProbeSkim,1)
  };
}

#endif
