//--------------------------------------------------------------------------------------------------
// $Id: ElectronIDMod.h,v 1.44 2011/03/11 15:13:13 ceballos Exp $
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
#include "MitAna/DataTree/interface/TrackFwd.h"
#include "MitAna/DataTree/interface/DecayParticleFwd.h"
#include "MitAna/DataTree/interface/PFCandidateFwd.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodSwitches.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodMeasurements.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include <TFile.h>
#include <TDirectory.h>

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
      Bool_t              GetApplyDZCut()             	  const { return fApplyDZCut;		  }
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
      Double_t            GetEtaMax()                  	  const { return fElectronEtaMax;	  }
      Bool_t              GetReverseD0Cut()           	  const { return fReverseD0Cut; 	  }
      Bool_t              GetReverseIsoCut()          	  const { return fReverseIsoCut;	  }
      Bool_t              GetApplyTriggerMatching()       const { return fApplyTriggerMatching;   }
      Double_t            GetTrackIsoCut()            	  const { return fTrackIsolationCut;	  }
      Bool_t              GetChargeFilter()           	  const { return fChargeFilter; 	  }
      Bool_t              Likelihood(const Electron *ele) const;
      Bool_t              PassIDCut(const Electron *el, ElectronTools::EElIdType idType) const;
      Bool_t              PassIsolationCut(const Electron *el, ElectronTools::EElIsoType isoType,
                                           const TrackCol *tracks, const Vertex *vertex, 
					   const Double_t rho) const;
      Bool_t              GetCombinedIdCut()               const { return fCombinedIdCut;      }
      void                SetApplyConversionFilterType1(Bool_t b){ fApplyConvFilterType1 = b;  }
      void                SetApplyConversionFilterType2(Bool_t b){ fApplyConvFilterType2 = b;  }
      void                SetNExpectedHitsInnerCut(Double_t cut) {fNExpectedHitsInnerCut = cut;}
      void                SetApplySpikeRemoval(Bool_t b)         { fApplySpikeRemoval  = b;    }
      void                SetApplyD0Cut(Bool_t b)                { fApplyD0Cut         = b;    }
      void                SetApplyDZCut(Bool_t b)                { fApplyDZCut         = b;    }
      void                SetCaloIsoCut(Double_t cut)            { fCaloIsolationCut   = cut;  }
      void                SetCombIsoCut(Double_t cut)            { fCombIsolationCut  = cut;   }
      void                SetD0Cut(Double_t cut)                 { fD0Cut = cut;	       }
      void                SetDZCut(Double_t cut)                 { fDZCut = cut;	       }
      void                SetWhichVertex(Int_t d) 		 { fWhichVertex = d; 	       }
      void                SetEcalJurIsoCut(Double_t cut)         { fEcalJuraIsoCut     = cut;  }
      void                SetGoodElectronsName(const char *n)    { fGoodElectronsName  = n;    }  
      void                SetOldMuonsName(const char *n)         { fNonIsolatedMuonsName  = n;    }  
      void                SetOldElectronsName(const char *n)     { fNonIsolatedElectronsName  = n;}  
      void                SetGoodName(const char *n)             { SetGoodElectronsName(n);    }  
      void                SetHcalIsoCut(Double_t cut)            { fHcalIsolationCut   = cut;  }
      void                SetIDLikelihoodCut(Double_t cut)       { fIDLikelihoodCut    = cut;  }
      void                SetIDType(const char *type)            { fElectronIDType     = type; }
      void                SetInputName(const char *n)            { fElectronBranchName = n;    }  
      void                SetIsoType(const char *type)           { fElectronIsoType    = type; }
      void                SetTriggerObjectsName(const char *n)   { fTrigObjectsName = n;       }
      void                SetOutputName(const char *n)           { SetGoodElectronsName(n);    }  
      void                SetPtMin(Double_t pt)                  { fElectronPtMin      = pt;   }
      void                SetEtMin(Double_t et)                  { fElectronEtMin      = et;   }      
      void                SetEtaMax(Double_t eta)                { fElectronEtaMax     = eta;  }
      void                SetReverseD0Cut(Bool_t b)              { fReverseD0Cut = b;	       }
      void                SetReverseIsoCut(Bool_t b)             { fReverseIsoCut = b;         }
      void                SetApplyTriggerMatching(Bool_t b)      { fApplyTriggerMatching = b;  }
      void                SetTrackIsoCut(Double_t cut)           { fTrackIsolationCut  = cut;  }
      void                SetChargeFilter(Bool_t b)              { fChargeFilter = b;	       }
      void                SetNWrongHitsMax(UInt_t n)             { fNWrongHitsMax = n;         }
      void                SetCombinedIdCut(Bool_t b)             { fCombinedIdCut  = b;        }
      void                SetApplyEcalFiducial(Bool_t b)         { fApplyEcalFiducial = b;     }
      void                SetApplyEcalSeeded(Bool_t b)           { fApplyEcalSeeded = b;       }      
      void                SetApplyCombinedIso(Bool_t b)          { fApplyCombinedIso = b;      }      
      void                SetElectronsFromBranch(Bool_t b)       { fElectronsFromBranch = b;   }      
      void                SetLH(ElectronLikelihood *l)	         { fLH = l;  	               }
      void                Setup();

    protected:
      void                Process();
      void                SlaveBegin();
      void                Terminate();

      TString                   fElectronBranchName;     //name of electron collection (input)
      TString                   fConversionBranchName;   //name of electron collection (input)
      TString                   fGoodElectronsName;      //name of exported "good electrons" col
      TString                   fNonIsolatedMuonsName;    //name of imported "old muon" collection
      TString                   fNonIsolatedElectronsName;//name of imported "old electron" collection
      TString                   fVertexName;	         //name of vertex collection
      TString                   fBeamSpotName;           //name of beamspot collection
      TString                   fTrackName;              //name of track collection
      TString                   fPFCandidatesName;       //name of pfcandidates collection
      TString                   fElectronIDType;         //type of electron ID we impose
      TString                   fElectronIsoType;        //type of electron Isolation we impose
      TString                   fTrigObjectsName;        //name of trigger object collection
      Double_t                  fElectronPtMin;	         //min pt cut
      Double_t                  fElectronEtMin;          //min pt cut            
      Double_t                  fElectronEtaMax;         //max eta cut
      Double_t                  fIDLikelihoodCut;        //cut value for ID likelihood
      Double_t                  fTrackIsolationCut;      //cut value for track isolation
      Double_t                  fCaloIsolationCut;       //cut value for calo isolation
      Double_t                  fEcalJuraIsoCut;         //cut value for ecal jurassic isolation
      Double_t                  fHcalIsolationCut;       //cut value for hcal isolation
      Double_t                  fCombIsolationCut;       //cut value for combined isolation
      Bool_t                    fApplyConvFilterType1;   //whether remove conversions using fit method
      Bool_t                    fApplyConvFilterType2;   //whether remove conversions using DCotTheta method
      UInt_t                    fNWrongHitsMax;          //whether to use wrong hits req 
                                                         //for conversion removal
      Double_t                  fNExpectedHitsInnerCut;  //cut value for NExpectedHitsInner maximum
      Bool_t                    fCombinedIdCut;          //whether to use full combined id
      Bool_t                    fApplySpikeRemoval;      //whether spike removal
      Bool_t                    fApplyD0Cut;             //whether apply d0 cut
      Bool_t                    fApplyDZCut;             //whether apply dz cut
      Bool_t                    fChargeFilter;           //whether apply GSF and CFT equal requirement
      Double_t                  fD0Cut;                  //max d0
      Double_t                  fDZCut;                  //max dz
      Int_t                     fWhichVertex;            //vertex to use (-2: beamspot, -1: closest in Z)
      Bool_t                    fReverseIsoCut;          //apply reversion iso cut (default=0)
      Bool_t                    fReverseD0Cut;           //apply reversion d0 cut (default=0)
      Bool_t                    fApplyTriggerMatching;   //match to hlt electron (default=0)
      Bool_t                    fApplyEcalSeeded;        //require ecal seeded flag
      Bool_t                    fApplyCombinedIso;       //apply combined isolation
      Bool_t                    fApplyEcalFiducial;      //apply ecal fiducial cuts on supercluster eta
      Bool_t                    fElectronsFromBranch;    //where to get input electrons
      ElectronTools::EElIdType  fElIdType;              //!identification scheme
      ElectronTools::EElIsoType fElIsoType;              //!isolation scheme
      const ElectronCol        *fElectrons;              //!electron collection
      const DecayParticleCol   *fConversions;            //!conversion collection
      const VertexCol          *fVertices;               //!vertices branches
      const BeamSpotCol        *fBeamSpot;               //!beamspot branch
      const TrackCol           *fTracks;                 //!Track branch     
      const PFCandidateCol     *fPFCandidates;           //!pfcandidate branch

      MuonCol  	               *fNonIsolatedMuons;	 //!pointer to old muon collection 
      ElectronCol	       *fNonIsolatedElectrons;	 //!pointer to old electron collection
      ElectronLikelihood       *fLH;                    //LH
      TString                   fPileupEnergyDensityName;
      const PileupEnergyDensityCol *fPileupEnergyDensity;

    ClassDef(ElectronIDMod, 1) // Electron identification module
  };
}
#endif
