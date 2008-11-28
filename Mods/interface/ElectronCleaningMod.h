//--------------------------------------------------------------------------------------------------
// $Id: ElectronCleaningMod.h,v 1.2 2008/11/27 16:30:26 loizides Exp $
//
// ElectronCleaningMod
//
// This Module performs cleaning of electrons, ie. it removes duplicate objects and good muons 
// from the good electrons.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_ELECTRONCLEANINGMOD_H
#define MITPHYSICS_MODS_ELECTRONCLEANINGMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"

namespace mithep 
{
  class ElectronCleaningMod : public BaseMod
  {
    public:
      ElectronCleaningMod(const char *name="ElectronCleaningMod", 
                     const char *title="Electron cleaning module");
      ~ElectronCleaningMod() {}

      void               SetGoodElectronsName(const char *name)  { fGoodElectronsName  = name; }
      void               SetCleanMuonsName(const char *name)     { fCleanMuonsName     = name; }
      void               SetCleanElectronsName(const char *name) { fCleanElectronsName = name; }
    protected:
      TString            fGoodElectronsName;  //name of good electrons (input)
      TString            fCleanMuonsName;     //name of clean muons (input)
      TString            fCleanElectronsName; //name of clean electrons (output)

      void               Process();
    
      ClassDef(ElectronCleaningMod,1) // Electron cleaning module
  };
}
#endif
