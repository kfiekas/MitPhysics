//--------------------------------------------------------------------------------------------------
// $Id: PhotonIDMod.h,v 1.7 2009/03/23 14:23:06 loizides Exp $
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

      Bool_t              GetApplyPixelSeed()         const { return fApplyPixelSeed;      }
      const char         *GetGoodName()               const { return GetGoodPhotonsName(); }   
      const char         *GetGoodPhotonsName()        const { return fGoodPhotonsName;     }   
      Double_t            GetHadOverEmMax()           const { return fHadOverEmMax;        }
      const char         *GetIDType()                 const { return fPhotonIDType;        }
      const char         *GetInputName()              const { return fPhotonBranchName;    }   
      const char         *GetIsoType()                const { return fPhotonIsoType;       }
      const char         *GetOutputName()             const { return GetGoodPhotonsName(); }   
      Double_t            GetPtMin()                  const { return fPhotonPtMin;         }
      void                SetApplyPixelSeed(Bool_t b)       { fApplyPixelSeed  = b;        }
      void                SetGoodName(const char *n)        { SetGoodPhotonsName(n);       }   
      void                SetGoodPhotonsName(const char *n) { fGoodPhotonsName = n;        }   
      void                SetHadOverEmMax(Double_t hoe)     { fHadOverEmMax    = hoe;      }
      void                SetIDType(const char *type)       { fPhotonIDType    = type;     }
      void                SetInputName(const char *n)       { fPhotonBranchName= n;        }   
      void                SetIsoType(const char *type)      { fPhotonIsoType   = type;     }
      void                SetOutputName(const char *n)      { SetGoodPhotonsName(n);       }    
      void                SetPtMin(Double_t pt)             { fPhotonPtMin     = pt;       }
      void                SetR9Min(Double_t x)              { fPhotonR9Min     = x;       }

      enum EPhIdType {
        kIdUndef = 0,       //not defined
        kTight,             //"Tight"
        kLoose,             //"Loose"
        kLooseEM,           //"LooseEM"
        kCustomId           //"Custom"
      };
      enum EPhIsoType {
        kIsoUndef = 0,      //not defined        
        kNoIso,             //"NoIso"
        kCombinedIso,       //"CombinedIso"
        kCustomIso          //"Custom"
      };

    protected:
      void                Process();
      void                SlaveBegin();

      TString             fPhotonBranchName;     //name of photon collection (input)
      TString             fGoodPhotonsName;      //name of exported "good photon" collection
      TString             fPhotonIDType;         //type of photon identification we impose
      TString             fPhotonIsoType;        //type of photon isolation we impose
      Double_t            fPhotonPtMin;          //min pt cut
      Double_t            fHadOverEmMax;         //maximum of hadronic/em energy
      Bool_t              fApplyPixelSeed;       //=true then apply pixel seed constraint
      Double_t            fPhotonR9Min;          //min R9 value
      const PhotonCol    *fPhotons;              //!photon branch
      EPhIdType           fPhIdType;             //!identification scheme
      EPhIsoType          fPhIsoType;            //!isolation scheme
    
    ClassDef(PhotonIDMod, 1) // Photon identification module
  };
}
#endif
