//--------------------------------------------------------------------------------------------------
// $Id: ElectronIDMod.h,v 1.7 2008/12/10 11:44:32 loizides Exp $
//
// ElectronIDMod
//
// This module applies electron identification criteria and exports a pointer to a collection
// of "good electrons" according to the specified identification scheme.
//
// Authors: S.Xie, C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_ELECTRONIDMOD_H
#define MITPHYSICS_MODS_ELECTRONIDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"

namespace mithep 
{
  class ElectronIDMod : public BaseMod
  {
    public:
      ElectronIDMod(const char *name="ElectronIDMod", 
                    const char *title="Electron identification module");
      ~ElectronIDMod() {}

      Double_t            GetCaloIsoCut()             const { return fCaloIsolationCut;       }
      Double_t            GetEcalJurIsoCut()          const { return fEcalJuraIsoCut;         }
      const char         *GetGoodName()               const { return GetGoodElectronsName();  }   
      const char         *GetGoodElectronsName()      const { return fGoodElectronsName;      }   
      Double_t            GetHcalIsoCut()             const { return fHcalIsolationCut;       }
      Double_t            GetIDLikelihoodCut()        const { return fIDLikelihoodCut;        }
      const char         *GetIDType()                 const { return fElectronIDType;         }
      const char         *GetInputName()              const { return fElectronBranchName;     }   
      const char         *GetIsoType()                const { return fElectronIsoType;        }
      const char         *GetOutputName()             const { return GetGoodElectronsName();  }
      Double_t            GetPtMin()                  const { return fElectronPtMin;          }
      Double_t            GetTrackIsoCut()            const { return fTrackIsolationCut;      }
      void                SetCaloIsoCut(Double_t cut)           { fCaloIsolationCut   = cut;  }
      void                SetEcalJurIsoCut(Double_t cut)        { fEcalJuraIsoCut     = cut;  }
      void                SetGoodName(const char *n)            { SetGoodElectronsName(n);    }   
      void                SetGoodElectronsName(const char *n)   { fGoodElectronsName  = n;    }   
      void                SetHcalIsoCut(Double_t cut)           { fHcalIsolationCut   = cut;  }
      void                SetIDLikelihoodCut(Double_t cut)      { fIDLikelihoodCut    = cut;  }
      void                SetIDType(const char *type)           { fElectronIDType     = type; }
      void                SetInputName(const char *n)           { fElectronBranchName = n;    }   
      void                SetIsoType(const char *type)          { fElectronIsoType    = type; }
      void                SetOutputName(const char *n)          { SetGoodElectronsName(n);    }   
      void                SetPtMin(Double_t pt)                 { fElectronPtMin      = pt;   }
      void                SetTrackIsoCut(Double_t cut)          { fTrackIsolationCut  = cut;  }

      enum EElIdType {
        kIdUndef = 0,       //not defined
        kTight,             //"Tight"
        kLoose,             //"Loose"
        kLikelihood,        //"Likelihood"
        kCustomId           //"Custom"
      };
      enum EElIsoType {
        kIsoUndef = 0,      //not defined
        kTrackCalo,         //"TrackCalo"
        kTrackJura,         //"TrackJura"
        kTrackJuraSliding,  //"TrackJuraSliding"
        kNoIso,             //"NoIso"
        kCustomIso          //"Custom"
      };

    protected:
      void                Process();
      void                SlaveBegin();

      TString             fElectronBranchName;     //name of electron collection (input)
      TString             fGoodElectronsName;      //name of exported "good electrons" collection
      TString             fElectronIDType;         //type of electron ID we impose
      TString             fElectronIsoType;        //type of electron Isolation we impose
      Double_t            fElectronPtMin;          //min pt cut
      Double_t            fIDLikelihoodCut;        //cut value for ID likelihood
      Double_t            fTrackIsolationCut;      //cut value for track isolation
      Double_t            fCaloIsolationCut;       //cut value for calo isolation
      Double_t            fEcalJuraIsoCut;         //cut value for ecal jurassic isolation
      Double_t            fHcalIsolationCut;       //cut value for hcal isolation
      const ElectronCol  *fElectrons;              //!electron collection
      EElIdType           fElIdType;               //!identification scheme
      EElIsoType          fElIsoType;              //!isolation scheme
    
    ClassDef(ElectronIDMod, 1) // Electron identification module
  };
}
#endif
