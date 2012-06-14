#ifndef H4LELETAGPROBESKIM_H
#define H4LELETAGPROBESKIM_H

#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/TriggerMask.h"
#include "MitPhysics/Skim/interface/BaseH4lSkim.h"
#include "TLorentzVector.h"

namespace mithep
{  
  class H4lEleTagProbeSkim : public BaseH4lSkim
  {    
  public:
    H4lEleTagProbeSkim(const char *name="H4lEleTagProbeSkim", const char *title="Electron Tag and Probe Skim");
    ~H4lEleTagProbeSkim();	    

  protected:
    void Begin();
    void BeginRun();
    void EndRun();
    void SlaveBegin();
    void SlaveTerminate();
    void Terminate();
    void Process();

    ClassDef(H4lEleTagProbeSkim,1)
  };
}

#endif
