//--------------------------------------------------------------------------------------------------
// $Id: ElectronIDMod.h,v 1.19 2009/06/02 05:30:44 loizides Exp $
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
      Double_t            GetTrackIsoCut()            const { return fTrackIsolationCut;      }
      Bool_t              GetReverseIsoCut()          const { return fReverseIsoCut;          }
      Bool_t              GetReverseD0Cut()           const { return fReverseD0Cut;           }
      void                SetApplyConversionFilter(Bool_t b)    { fApplyConvFilter    = b;    }
      void                SetApplyD0Cut(Bool_t b)               { fApplyD0Cut         = b;    }
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
      void                SetD0Cut(Double_t cut)                { fD0Cut = cut;               }
      void                SetReverseIsoCut(Bool_t b)            { fReverseIsoCut = b;         }
      void                SetReverseD0Cut(Bool_t b)             { fReverseD0Cut = b;          }

      enum EElIdType {
        kIdUndef = 0,       //not defined
        kTight,             //"Tight"
        kLoose,             //"Loose"
        kLikelihood,        //"Likelihood"
        kNoId,              //"NoId"
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
      void                    Process();
      void                    SlaveBegin();

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
      Bool_t                  fApplyD0Cut;             //whether apply d0 cut
      Double_t                fD0Cut;                  //max d0
      Bool_t                  fReverseIsoCut;          //apply reversion iso cut (default=0)
      Bool_t                  fReverseD0Cut;           //apply reversion d0 cut (default=0)
      EElIdType               fElIdType;               //!identification scheme
      EElIsoType              fElIsoType;              //!isolation scheme
      const ElectronCol      *fElectrons;              //!electron collection
      const DecayParticleCol *fConversions;            //!conversion collection
      const VertexCol        *fVertices;               //!vertices branches

    ClassDef(ElectronIDMod, 1) // Electron identification module
  };
}
#endif
