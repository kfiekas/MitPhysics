//--------------------------------------------------------------------------------------------------
// $Id: PhotonIDMod.h,v 1.12 2009/12/06 14:59:43 ceballos Exp $
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
#include "MitAna/DataTree/interface/PhotonFwd.h"

namespace mithep 
{
  class PhotonIDMod : public BaseMod
  {
    public:
      PhotonIDMod(const char *name="PhotonIDMod", 
                  const char *title="Photon identification module");

      Bool_t              GetApplySpikeRemoval()      const { return fApplySpikeRemoval;   }
      Bool_t              GetApplyPixelSeed()         const { return fApplyPixelSeed;      }
      const char         *GetGoodName()               const { return GetGoodPhotonsName(); }   
      const char         *GetGoodPhotonsName()        const { return fGoodPhotonsName;     }   
      Double_t            GetHadOverEmMax()           const { return fHadOverEmMax;        }
      const char         *GetIDType()                 const { return fPhotonIDType;        }
      const char         *GetInputName()              const { return fPhotonBranchName;    }   
      const char         *GetIsoType()                const { return fPhotonIsoType;       }
      const char         *GetOutputName()             const { return GetGoodPhotonsName(); }   
      Double_t            GetPtMin()                  const { return fPhotonPtMin;         }
      Bool_t              GetApplyFiduciality()       const { return fFiduciality;         }
      Double_t            GetEtaWidthEB()	      const { return fEtaWidthEB;	   }
      Double_t            GetEtaWidthEE()	      const { return fEtaWidthEE;	   }
      Double_t            GetAbsEtaMax()	      const { return fAbsEtaMax;	   }
      void                SetApplySpikeRemoval(Bool_t b)    { fApplySpikeRemoval  = b;     }
      void                SetApplyPixelSeed(Bool_t b)       { fApplyPixelSeed  = b;        }
      void                SetGoodName(const char *n)        { SetGoodPhotonsName(n);       }   
      void                SetGoodPhotonsName(const char *n) { fGoodPhotonsName = n;        }   
      void                SetHadOverEmMax(Double_t hoe)     { fHadOverEmMax    = hoe;      }
      void                SetIDType(const char *type)       { fPhotonIDType    = type;     }
      void                SetInputName(const char *n)       { fPhotonBranchName= n;        }   
      void                SetIsoType(const char *type)      { fPhotonIsoType   = type;     }
      void                SetOutputName(const char *n)      { SetGoodPhotonsName(n);       }    
      void                SetPtMin(Double_t pt)             { fPhotonPtMin     = pt;       }
      void                SetR9Min(Double_t x)              { fPhotonR9Min     = x;        }
      void                SetEtaWidthEB(Double_t x)	    { fEtaWidthEB      = x;	   }
      void                SetEtaWidthEE(Double_t x)         { fEtaWidthEE      = x;	   }
      void                SetAbsEtaMax(Double_t x)          { fAbsEtaMax       = x;	   }

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
      Bool_t              fApplySpikeRemoval;    //whether apply spike removal      
      Bool_t              fApplyPixelSeed;       //=true then apply pixel seed constraint
      Double_t            fPhotonR9Min;          //min R9 value
      EPhIdType           fPhIdType;             //!identification scheme
      EPhIsoType          fPhIsoType;            //!isolation scheme
      Bool_t              fFiduciality;          //=true then apply fiducual requirement
      Double_t            fEtaWidthEB;  	 //max Eta Width in ECAL Barrel
      Double_t            fEtaWidthEE;  	 //max Eta Width in ECAL End Cap
      Double_t            fAbsEtaMax;  	         //max Abs Eta
      const PhotonCol    *fPhotons;              //!photon branch
    
    ClassDef(PhotonIDMod, 1) // Photon identification module
  };
}
#endif
