//--------------------------------------------------------------------------------------------------
// $Id: JetCleaningMod.h,v 1.2 2008/11/27 16:30:26 loizides Exp $
//
// JetCleaningMod
//
// This Module performs cleaning of jets, ie it removes jets which point 
// in the same direction as a clean isolated electrons.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_JETCLEANINGMOD_H
#define MITPHYSICS_MODS_JETCLEANINGMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"

namespace mithep 
{
  class JetCleaningMod : public BaseMod
  {
    public:
      JetCleaningMod(const char *name="JetCleaningMod", 
                     const char *title="Jet cleaning module");
      ~JetCleaningMod() {}

      void               SetCleanElectronsName(const char *name) { fCleanElectronsName = name; }
      void               SetGoodJetsName(const char *name)       { fGoodJetsName       = name; }  
      void               SetCleanJetsName(const char *name)      { fCleanJetsName      = name; }
 
    protected:
      TString            fCleanElectronsName; //name of clean electrons (input)
      TString            fGoodJetsName;       //name of good jets (input)
      TString            fCleanJetsName;      //name of clean jets (output)

      void               Process();
   
      ClassDef(JetCleaningMod,1) // Jet cleaning module
  };
}
#endif
