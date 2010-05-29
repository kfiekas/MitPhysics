//--------------------------------------------------------------------------------------------------
// $Id: ElectronIDMod.h,v 1.35 2010/05/29 11:56:17 bendavid Exp $
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
#include "MitPhysics/Utils/interface/ElectronTools.h"

namespace mithep 
{
  class ElectronIDMod : public BaseMod
  {
    public:
      ElectronIDMod(const char *name="ElectronIDMod", 
                    const char *title="Electron identification module");

      Bool_t              GetApplyConversionFilterType1() const { return fApplyConvFilterType1;	  }
      Bool_t              GetApplyConversionFilterType2() const { return fApplyConvFilterType2;	  }
      Bool_t              GetApplySpikeRemoval()      	  const { return fApplySpikeRemoval;	  }
      Bool_t              GetApplyD0Cut()             	  const { return fApplyD0Cut;		  }
      Double_t            GetCaloIsoCut()             	  const { return fCaloIsolationCut;	  }
      Double_t            GetEcalJurIsoCut()          	  const { return fEcalJuraIsoCut;	  }
      Double_t            GetCombIsoCut()             	  const { return fCombIsolationCut;	  }
      const char         *GetGoodName()                   const { return GetGoodElectronsName();  }
      const char         *GetGoodElectronsName()          const { return fGoodElectronsName;	  }
      Double_t            GetHcalIsoCut()             	  const { return fHcalIsolationCut;	  }
      Double_t            GetIDLikelihoodCut()        	  const { return fIDLikelihoodCut;	  }
      const char         *GetIDType()                 	  const { return fElectronIDType;	  }
      const char         *GetInputName()                  const { return fElectronBranchName;	  }
      const char         *GetIsoType()                	  const { return fElectronIsoType;	  }
      const char         *GetOutputName()             	  const { return GetGoodElectronsName();  }
      Double_t            GetPtMin()                  	  const { return fElectronPtMin;	  }
      Bool_t              GetReverseD0Cut()           	  const { return fReverseD0Cut; 	  }
      Bool_t              GetReverseIsoCut()          	  const { return fReverseIsoCut;	  }
      Bool_t              GetApplyTriggerMatching()       const { return fApplyTriggerMatching;   }
      Double_t            GetTrackIsoCut()            	  const { return fTrackIsolationCut;	  }
      Bool_t              GetChargeFilter()           	  const { return fChargeFilter; 	  }
      Bool_t              PassIDCut(const Electron *el, ElectronTools::EElIdType idType) const;
      Bool_t              PassIsolationCut(const Electron *el, 
                                           ElectronTools::EElIsoType isoType) const;
      Bool_t              GetCombinedIdCut()               const { return fCombinedIdCut;      }
      void                SetApplyConversionFilterType1(Bool_t b){ fApplyConvFilterType1 = b;  }
      void                SetApplyConversionFilterType2(Bool_t b){ fApplyConvFilterType2 = b;  }
      void                SetNExpectedHitsInnerCut(Double_t cut) {fNExpectedHitsInnerCut = cut;}
      void                SetApplySpikeRemoval(Bool_t b)         { fApplySpikeRemoval  = b;    }
      void                SetApplyD0Cut(Bool_t b)                { fApplyD0Cut         = b;    }
      void                SetCaloIsoCut(Double_t cut)            { fCaloIsolationCut   = cut;  }
      void                SetCombIsoCut(Double_t cut)            { fCombIsolationCut  = cut;   }
      void                SetD0Cut(Double_t cut)                 { fD0Cut = cut;	       }
      void                SetEcalJurIsoCut(Double_t cut)         { fEcalJuraIsoCut     = cut;  }
      void                SetGoodElectronsName(const char *n)    { fGoodElectronsName  = n;    }  
      void                SetGoodName(const char *n)             { SetGoodElectronsName(n);    }  
      void                SetHcalIsoCut(Double_t cut)            { fHcalIsolationCut   = cut;  }
      void                SetIDLikelihoodCut(Double_t cut)       { fIDLikelihoodCut    = cut;  }
      void                SetIDType(const char *type)            { fElectronIDType     = type; }
      void                SetInputName(const char *n)            { fElectronBranchName = n;    }  
      void                SetIsoType(const char *type)           { fElectronIsoType    = type; }
      void                SetTriggerObjectsName(const char *n)   { fTrigObjectsName = n;       }
      void                SetOutputName(const char *n)           { SetGoodElectronsName(n);    }  
      void                SetPtMin(Double_t pt)                  { fElectronPtMin      = pt;   }
      void                SetReverseD0Cut(Bool_t b)              { fReverseD0Cut = b;	       }
      void                SetReverseIsoCut(Bool_t b)             { fReverseIsoCut = b;         }
      void                SetApplyTriggerMatching(Bool_t b)     { fApplyTriggerMatching = b;  }
      void                SetTrackIsoCut(Double_t cut)           { fTrackIsolationCut  = cut;  }
      void                SetChargeFilter(Bool_t b)              { fChargeFilter = b;	       }
      void                SetWrongHitsRequirement(Bool_t b)      { fWrongHitsRequirement = b;  }
      void                SetCombinedIdCut(Bool_t b)             { fCombinedIdCut  = b;        }
      void                Setup();

    protected:
      void                Process();
      void                SlaveBegin();


      TString                   fElectronBranchName;     //name of electron collection (input)
      TString                   fConversionBranchName;   //name of electron collection (input)
      TString                   fGoodElectronsName;      //name of exported "good electrons" col
      TString                   fVertexName;	         //name of vertex collection
      TString                   fElectronIDType;         //type of electron ID we impose
      TString                   fElectronIsoType;        //type of electron Isolation we impose
      TString                   fTrigObjectsName;        //name of trigger object collection
      Double_t                  fElectronPtMin;	         //min pt cut
      Double_t                  fIDLikelihoodCut;        //cut value for ID likelihood
      Double_t                  fTrackIsolationCut;      //cut value for track isolation
      Double_t                  fCaloIsolationCut;       //cut value for calo isolation
      Double_t                  fEcalJuraIsoCut;         //cut value for ecal jurassic isolation
      Double_t                  fHcalIsolationCut;       //cut value for hcal isolation
      Double_t                  fCombIsolationCut;       //cut value for combined isolation
      Bool_t                    fApplyConvFilterType1;   //whether remove conversions using fit method
      Bool_t                    fApplyConvFilterType2;   //whether remove conversions using DCotTheta method
      Bool_t                    fWrongHitsRequirement;   //whether to use wrong hits req 
                                                         //for conversion removal
      Double_t                  fNExpectedHitsInnerCut;  //cut value for NExpectedHitsInner maximum
      Bool_t                    fCombinedIdCut;          //whether to use full combined id
      Bool_t                    fApplySpikeRemoval;      //whether apply d0 cut
      Bool_t                    fApplyD0Cut;             //whether apply d0 cut
      Bool_t                    fChargeFilter;           //whether apply GSF and CFT equal requirement
      Double_t                  fD0Cut;                  //max d0
      Bool_t                    fReverseIsoCut;          //apply reversion iso cut (default=0)
      Bool_t                    fReverseD0Cut;           //apply reversion d0 cut (default=0)
      Bool_t                    fApplyTriggerMatching;   //match to hlt electron (default=0)
      ElectronTools::EElIdType  fElIdType;              //!identification scheme
      ElectronTools::EElIsoType fElIsoType;              //!isolation scheme
      const ElectronCol        *fElectrons;              //!electron collection
      const DecayParticleCol   *fConversions;            //!conversion collection
      const VertexCol          *fVertices;               //!vertices branches
 
    ClassDef(ElectronIDMod, 1) // Electron identification module
  };
}
#endif
