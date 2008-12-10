//--------------------------------------------------------------------------------------------------
// $Id: PhotonCleaningMod.h,v 1.3 2008/12/10 11:44:33 loizides Exp $
//
// PhotonCleaningMod
//
// This Module performs cleaning of jets, ie it removes jets which point 
// in the same direction as a clean isolated electrons.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PHOTONCLEANINGMOD_H
#define MITPHYSICS_MODS_PHOTONCLEANINGMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"

namespace mithep 
{
  class PhotonCleaningMod : public BaseMod
  {
    public:
      PhotonCleaningMod(const char *name="PhotonCleaningMod", 
                     const char *title="Photon cleaning module");
      ~PhotonCleaningMod() {}

      const char      *GetCleanElectronsName()   const { return fCleanElectronsName;   }
      const char      *GetCleanName()            const { return GetCleanPhotonsName(); }
      const char      *GetCleanPhotonsName()     const { return fCleanPhotonsName;     }
      const char      *GetGoodPhotonsName()      const { return fGoodPhotonsName;      } 
      Double_t         GetMinDeltaRToElectron()  const { return fMinDeltaRToElectron;  }
      const char      *GetOutputName()           const { return GetCleanPhotonsName(); }
      void             SetCleanElectronsName(const char *name)  { fCleanElectronsName    = name; }
      void             SetCleanName(const char *name)           { SetCleanPhotonsName(name);     }
      void             SetCleanPhotonsName(const char *name)    { fCleanPhotonsName      = name; }
      void             SetGoodPhotonsName(const char *name)     { fGoodPhotonsName       = name; } 
      void             SetMinDeltaRToElectron(Double_t dr)      { fMinDeltaRToElectron   = dr;   }
      void             SetOutputName(const char *name)          { SetCleanPhotonsName(name);     }

    protected:
      void             Process();

      TString          fCleanElectronsName;   //name of clean electrons (input)
      TString          fGoodPhotonsName;      //name of good jets (input)
      TString          fCleanPhotonsName;     //name of clean jets (output)
      Double_t         fMinDeltaRToElectron;  //delta R threshold for separating electrons/photons
   
    ClassDef(PhotonCleaningMod, 1) // Photon cleaning module
  };
}
#endif
