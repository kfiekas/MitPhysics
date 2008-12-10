//--------------------------------------------------------------------------------------------------
// $Id: PhotonCleaningMod.h,v 1.2 2008/12/04 13:53:33 loizides Exp $
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

      void             SetCleanElectronsName(const char *name)    { fCleanElectronsName    = name; }
      void             SetGoodPhotonsName(const char *name)       { fGoodPhotonsName       = name; } 
      void             SetCleanPhotonsName(const char *name)      { fCleanPhotonsName      = name; }
      void             SetMinDeltaRToElectron(const Double_t x)   { fMinDeltaRToElectron   = x;    }

    protected:
      void             Process();

      TString          fCleanElectronsName;   //name of clean electrons (input)
      TString          fGoodPhotonsName;      //name of good jets (input)
      TString          fCleanPhotonsName;     //name of clean jets (output)
      Double_t         fMinDeltaRToElectron;  //delta R threshold for separating electrons/photons
   
    ClassDef(PhotonCleaningMod,1) // Photon cleaning module
  };
}
#endif
