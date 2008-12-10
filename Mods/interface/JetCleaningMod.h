//--------------------------------------------------------------------------------------------------
// $Id: JetCleaningMod.h,v 1.5 2008/12/04 13:53:33 loizides Exp $
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

      void               SetCleanElectronsName(const char *name)  { fCleanElectronsName = name; }
      void               SetCleanPhotonsName(const char *name)    { fCleanPhotonsName = name; }
      void               SetGoodJetsName(const char *name)        { fGoodJetsName       = name; }  
      void               SetCleanJetsName(const char *name)       { fCleanJetsName      = name; }
      void               SetMinDeltaRToElectron(const Double_t x) { fMinDeltaRToElectron   = x; }
      void               SetMinDeltaRToPhoton(const Double_t x)   { fMinDeltaRToPhoton   = x;   }

    protected:
      void               Process();

      TString            fCleanElectronsName;   //name of clean electrons (input)
      TString            fCleanPhotonsName;     //name of clean photons   (input)
      TString            fGoodJetsName;         //name of good jets       (input)
      TString            fCleanJetsName;        //name of clean jets      (output)
      Double_t           fMinDeltaRToElectron;  //delta R for separating electrons from jets
      Double_t           fMinDeltaRToPhoton;    //delta R for separating photons from jets
   
    ClassDef(JetCleaningMod,1) // Jet cleaning module
  };
}
#endif
