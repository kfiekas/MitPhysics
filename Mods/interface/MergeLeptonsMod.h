//--------------------------------------------------------------------------------------------------
// $Id: MergeLeptonsMod.h,v 1.2 2008/12/10 17:28:22 loizides Exp $
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

      const char              *GetElectronsName()           const { return fElName;         }
      const char              *GetMergedName()              const { return fMergedName;     }
      const char              *GetMuonsName()               const { return fMuName;         }
      const char              *GetOutputName()              const { return GetMergedName(); }
      void                     SetElectronsName(const char *n)    { fElName=n;              }
      void                     SetMergedName(const char *n)       { fMergedName=n;          }
      void                     SetMuonsName(const char *n)        { fMuName=n;              }
      void                     SetOutputName(const char *n)       { SetMergedName(n);       }

    protected:
      void                     Process();

      TString                  fElName;        //name of electrons collection (input)
      TString                  fMuName;        //name of muons collection (input)
      TString                  fMergedName;    //name of merged collection (output)
      const ElectronCol       *fElIn;          //!pointer to electron collection
      const MuonCol           *fMuIn;          //!pointer to muon collection 
      ParticleOArr            *fColOut;        //!pointer to merged collection

    ClassDef(MergeLeptonsMod, 1) // Merging leptons module
  };
}
#endif
