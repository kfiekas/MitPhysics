//--------------------------------------------------------------------------------------------------
// $Id: ElectronIDMod.h,v 1.25 2009/09/29 19:41:13 loizides Exp $
//
// ElectronIDMod
//
// This module applies electron identification criteria and exports a pointer to a collection
// of "good electrons" according to the specified identification scheme.
//
// See http://indico.cern.ch/contributionDisplay.py?contribId=1&confId=42251
//
// Authors: S.Xie, C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_ELECTRONIDMOD_H
#define MITPHYSICS_MODS_ELECTRONIDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/DecayParticleFwd.h"

namespace mithep 
{
  class ElectronIDMod : public BaseMod
  {
    public:
      ElectronIDMod(const char *name="ElectronIDMod", 
                    const char *title="Electron identification module");

     enum EElIdType {
        kIdUndef = 0,       //not defined
        kTight,             //"Tight"
        kLoose,             //"Loose"
        kLikelihood,        //"Likelihood"
        kNoId,              //"NoId"
        kZeeId,             //"ZeeId"
        kCustomIdLoose,     //"CustomLoose"
        kCustomIdTight      //"CustomTight"
      };

      enum EElIsoType {
        kIsoUndef = 0,      //not defined
        kTrackCalo,         //"TrackCalo"
        kTrackJura,         //"TrackJura"
        kTrackJuraSliding,  //"TrackJuraSliding"
        kNoIso,             //"NoIso"
        kZeeIso,            //"ZeeIso"
        kCustomIso          //"Custom"
      };

      Bool_t              GetApplyConversionFilter()  const { return fApplyConvFilter;        }
      Bool_t              GetApplyD0Cut()             const { return fApplyD0Cut;             }
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
      Bool_t              GetReverseD0Cut()           const { return fReverseD0Cut;           }
      Bool_t              GetReverseIsoCut()          const { return fReverseIsoCut;          }
      Double_t            GetTrackIsoCut()            const { return fTrackIsolationCut;      }
      Bool_t              GetChargeFilter()           const { return fChargeFilter;           }
      Bool_t              PassChargeFilter(const Electron *el) const;
      Bool_t              PassConversionFilter(const Electron *el, const DecayParticleCol *conversions) const;
      Bool_t              PassD0Cut(const Electron *el, const VertexCol *vertices) const;
      Bool_t              PassIDCut(const Electron *el, EElIdType idType) const;
      Bool_t              PassIsolationCut(const Electron *el, EElIsoType isoType) const;
      void                Setup();
      void                SetApplyConversionFilter(Bool_t b)    { fApplyConvFilter    = b;    }
      void                SetApplyD0Cut(Bool_t b)               { fApplyD0Cut         = b;    }
      void                SetCaloIsoCut(Double_t cut)           { fCaloIsolationCut   = cut;  }
      void                SetD0Cut(Double_t cut)                { fD0Cut = cut;               }
      void                SetEcalJurIsoCut(Double_t cut)        { fEcalJuraIsoCut     = cut;  }
      void                SetGoodElectronsName(const char *n)   { fGoodElectronsName  = n;    }   
      void                SetGoodName(const char *n)            { SetGoodElectronsName(n);    }   
      void                SetHcalIsoCut(Double_t cut)           { fHcalIsolationCut   = cut;  }
      void                SetIDLikelihoodCut(Double_t cut)      { fIDLikelihoodCut    = cut;  }
      void                SetIDType(const char *type)           { fElectronIDType     = type; }
      void                SetInputName(const char *n)           { fElectronBranchName = n;    }   
      void                SetIsoType(const char *type)          { fElectronIsoType    = type; }
      void                SetOutputName(const char *n)          { SetGoodElectronsName(n);    }   
      void                SetPtMin(Double_t pt)                 { fElectronPtMin      = pt;   }
      void                SetReverseD0Cut(Bool_t b)             { fReverseD0Cut = b;          }
      void                SetReverseIsoCut(Bool_t b)            { fReverseIsoCut = b;         }
      void                SetTrackIsoCut(Double_t cut)          { fTrackIsolationCut  = cut;  }
      void                SetChargeFilter(Bool_t b)             { fChargeFilter = b;          }
      void                SetWrongHitsRequirement(Bool_t b)     { fWrongHitsRequirement = b;  }



    protected:
      Bool_t              PassCustomID(const Electron *el) const;
      void                Process();
      void                SetCustomIDCuts(EElIdType idt);
      void                SlaveBegin();


      TString                 fElectronBranchName;     //name of electron collection (input)
      TString                 fConversionBranchName;   //name of electron collection (input)
      TString                 fGoodElectronsName;      //name of exported "good electrons" col
      TString                 fVertexName;	       //name of vertex collection
      TString                 fElectronIDType;         //type of electron ID we impose
      TString                 fElectronIsoType;        //type of electron Isolation we impose
      Double_t                fElectronPtMin;	       //min pt cut
      Double_t                fIDLikelihoodCut;        //cut value for ID likelihood
      Double_t                fTrackIsolationCut;      //cut value for track isolation
      Double_t                fCaloIsolationCut;       //cut value for calo isolation
      Double_t                fEcalJuraIsoCut;         //cut value for ecal jurassic isolation
      Double_t                fHcalIsolationCut;       //cut value for hcal isolation
      Bool_t                  fApplyConvFilter;        //whether remove conversions
      Bool_t                  fWrongHitsRequirement;   //whether to use wrong hits req for conversion removal
      Bool_t                  fApplyD0Cut;             //whether apply d0 cut
      Bool_t                  fChargeFilter;           //whether apply GSF and CFT equal requirement
      Double_t                fD0Cut;                  //max d0
      Bool_t                  fReverseIsoCut;          //apply reversion iso cut (default=0)
      Bool_t                  fReverseD0Cut;           //apply reversion d0 cut (default=0)
      EElIdType               fElIdType;               //!identification scheme
      EElIsoType              fElIsoType;              //!isolation scheme
      const ElectronCol      *fElectrons;              //!electron collection
      const DecayParticleCol *fConversions;            //!conversion collection
      const VertexCol        *fVertices;               //!vertices branches
      Double_t                fCuts[6][8];             //!custom id cuts

    ClassDef(ElectronIDMod, 1) // Electron identification module
  };
}
#endif
