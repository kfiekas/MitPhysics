//--------------------------------------------------------------------------------------------------
// $Id: PhotonIDMod.h,v 1.5 2008/11/28 09:13:50 loizides Exp $
//
// PhotonIDMod
//
// This module applies photon identification criteria and exports a pointer to a collection
// of "good photons" according to the specified identification scheme.
//
// Authors: S.Xie, C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PHOTONIDMOD_H
#define MITPHYSICS_MODS_PHOTONIDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"

namespace mithep 
{
  class PhotonIDMod : public BaseMod
  {
    public:
      PhotonIDMod(const char *name="PhotonIDMod", 
                    const char *title="Photon identification module");
      ~PhotonIDMod() {}

      void          SetPhotonBranchName(const char *n) { fPhotonBranchName       = n;    }   
      void          SetGoodPhotonsName(const char *n)  { fGoodPhotonsName        = n;    }   
      void          SetPhotonIDType(const char *type)  { fPhotonIDType           = type; }
      void          SetPhotonIsoType(const char *type) { fPhotonIsoType          = type; }
      void          SetPhotonPtMin(Double_t pt)        { fPhotonPtMin            = pt;   }


      enum EPhIdType {
        kIdUndef = 0,       //not defined
        kTight,             //"Tight"
        kLoose,             //"Loose"
        kLooseEM,        //"LooseEM"
        kCustomId           //"Custom"
      };
      enum EPhIsoType {
        kIsoUndef = 0,      //not defined        
        kNoIso,             //"NoIso"
        kCustomIso          //"Custom"
      };

    protected:
      TString       fPhotonBranchName;     //branch name of electron collection
      TString       fGoodPhotonsName;      //name of exported "good electrons" collection
      TString       fPhotonIDType;	     //type of electron ID we impose
      TString       fPhotonIsoType;	     //type of electron Isolation we impose
      Double_t      fPhotonPtMin;	     //min pt cut
      PhotonCol    *fPhotons;		     //!electron branch
      EPhIdType     fPhIdType;               //!identification scheme
      EPhIsoType    fPhIsoType;              //!isolation scheme

      void          Process();
      void          SlaveBegin();
    
      ClassDef(PhotonIDMod,1) // Photon identification module
  };
}
#endif
