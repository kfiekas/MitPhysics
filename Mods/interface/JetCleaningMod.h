//--------------------------------------------------------------------------------------------------
// $Id: JetCleaningMod.h,v 1.9 2009/04/02 12:40:16 ceballos Exp $
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

      const char        *GetCleanElectronsName()  const { return fCleanElectronsName;  }
      const char        *GetCleanMuonsName()      const { return fCleanMuonsName;  }
      const char        *GetCleanJetsName()       const { return fCleanJetsName;       }
      const char        *GetCleanName()           const { return GetCleanJetsName();   }
      const char        *GetCleanPhotonsName()    const { return fCleanPhotonsName;    }
      const char        *GetGoodJetsName()        const { return fGoodJetsName;        }  
      Double_t           GetMinDeltaRToElectron() const { return fMinDeltaRToElectron; }
      Double_t           GetMinDeltaRToMuon()     const { return fMinDeltaRToMuon; }
      Double_t           GetMinDeltaRToPhoton()   const { return fMinDeltaRToPhoton;   }
      const char        *GetOutputName()          const { return GetCleanJetsName();   }
      void               SetCleanElectronsName(const char *name) { fCleanElectronsName  = name; }
      void               SetCleanJetsName(const char *name)      { fCleanJetsName       = name; }
      void               SetCleanMuonsName(const char *name)     { fCleanMuonsName      = name; }
      void               SetCleanName(const char *name)          { SetCleanJetsName(name);      }
      void               SetCleanPhotonsName(const char *name)   { fCleanPhotonsName    = name; }
      void               SetGoodJetsName(const char *name)       { fGoodJetsName        = name; }  
      void               SetMinDeltaRToElectron(Double_t dr)     { fMinDeltaRToElectron = dr;   }
      void               SetMinDeltaRToMuon(Double_t dr)         { fMinDeltaRToMuon     = dr;   }
      void               SetMinDeltaRToPhoton(Double_t dr)       { fMinDeltaRToPhoton   = dr;   }
      void               SetOutputName(const char *name)         { SetCleanJetsName(name);      }

    protected:
      void               Process();

      TString            fCleanElectronsName;   //name of clean electrons (input)
      TString            fCleanMuonsName;       //name of clean muons (input)
      TString            fCleanPhotonsName;     //name of clean photons   (input)
      TString            fGoodJetsName;         //name of good jets       (input)
      TString            fCleanJetsName;        //name of clean jets      (output)
      Double_t           fMinDeltaRToElectron;  //delta R for separating electrons from jets
      Double_t           fMinDeltaRToMuon;      //delta R for separating muons from jets
      Double_t           fMinDeltaRToPhoton;    //delta R for separating photons from jets
   
    ClassDef(JetCleaningMod, 1) // Jet cleaning module
  };
}
#endif
