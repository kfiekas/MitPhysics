//--------------------------------------------------------------------------------------------------
// $Id: MCParticlesValMod.h,v 1.2 2009/03/23 22:17:06 loizides Exp $
//
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_VALIDATION_MCPARTICLESVALMOD_H
#define MITPHYSICS_VALIDATION_MCPARTICLESVALMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MCParticleFwd.h"

class TH1D;

namespace mithep 
{
  class MCParticlesValMod : public BaseMod 
  {
    public:
      MCParticlesValMod(const char *name="MCParticlesValMod", 
                        const char *title="MCParticles analysis module");

      const char              *GetPartName()              const { return fPartName; }
      void                     SetPartName(const char *n)       { fPartName=n;      }

    protected:
      void                     Process();
      void                     SlaveBegin();
      void                     SlaveTerminate();

      TString                  fPartName;   //branch name of MCParticle collection
      const MCParticleCol     *fParticles;  //!pointer to generated particle branch
      TH1D                    *fHs[100];    //!histograms

      ClassDef(MCParticlesValMod, 1) // MCParticles analysis module
  };
}
#endif
