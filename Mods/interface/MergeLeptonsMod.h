//--------------------------------------------------------------------------------------------------
// $Id: MergeLeptonsMod.h,v 1.3 2008/12/09 10:18:33 loizides Exp $
//
// MergeLeptonsMod
//
// This module merges muon and electron collections. (Note if need be this can easily be
// generalized).
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_MERGELEPTONSMOD_H
#define MITPHYSICS_MODS_MERGELEPTONSMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"

namespace mithep 
{
  class MergeLeptonsMod : public BaseMod 
  {
    public:
      MergeLeptonsMod(const char *name="MergeLeptonsMod", 
                      const char *title="Merging leptons module");
      ~MergeLeptonsMod() {}

      const char              *GetElectronsName()           const { return fElName;     }
      const char              *GetMergedName()              const { return fMergedName; }
      const char              *GetMuonsName()               const { return fMuName;     }
      void                     SetElectronsName(const char *n)    { fElName=n;          }
      void                     SetMergedName(const char *n)       { fMergedName=n;      }
      void                     SetMuonsName(const char *n)        { fMuName=n;          }

    protected:
      void                     Process();

      TString                  fElName;        //name of electrons collection
      TString                  fMuName;        //name of muons collection
      TString                  fMergedName;    //name of merged collection
      const ElectronCol       *fElIn;          //!pointer to electron collection (in) 
      const MuonCol           *fMuIn;          //!pointer to muon collection (in) 
      ParticleOArr            *fColOut;        //!pointer to merged collection (out)

    ClassDef(MergeLeptonsMod, 1) // Merging leptons module
  };
}
#endif
