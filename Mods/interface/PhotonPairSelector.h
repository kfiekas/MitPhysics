//--------------------------------------------------------------------------------------------------
// $Id: PhotonPairSelector.h,v 1.2 2011/07/15 17:24:37 fabstoec Exp $
//
// PhotonPairSelector
//
// Authors: F. Stoeckli
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PHOTONPAIRSELECTOR_H
#define MITPHYSICS_MODS_PHOTONPAIRSELECTOR_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PhotonFwd.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"

#include "MitPhysics/Utils/interface/PhotonTools.h"

class TNtuple;
class TRandom3;

namespace mithep 
{
  class PhotonPairSelectorPhoton : public TObject
  {
    public:  
      void SetVars(const Photon *p, const MCParticle *m = 0);

    private:  
      Float_t e;
      Float_t pt;
      Float_t eta;
      Float_t phi;
      Float_t r9;
      Float_t e5x5;
      Float_t sce;
      Float_t scrawe;
      Float_t scpse;
      Float_t sceta;
      Float_t scphi;
      Bool_t isbarrel;
      Bool_t isr9reco;
      Bool_t isr9cat;
      UChar_t phcat;
      Bool_t ispromptgen;
      Float_t gene;
      Float_t genpt;
      Float_t geneta;
      Float_t genphi;
      
      ClassDef(PhotonPairSelectorPhoton, 1)

      
  };
  
  class PhotonPairSelectorDiphotonEvent : public TObject
  {
    public:
      
      PhotonPairSelectorDiphotonEvent() : photons(PhotonPairSelectorPhoton::Class(),2)
      {
        new (photons[0]) PhotonPairSelectorPhoton;
        new (photons[1]) PhotonPairSelectorPhoton;
      }
      ~PhotonPairSelectorDiphotonEvent() { photons.Clear(); }
    
      Float_t rho;
      Float_t genHiggspt;
      Float_t genHiggsZ;
      Float_t gencostheta;
      Float_t vtxZ;
      Int_t   numPU;
      Int_t   numPUminus;
      Int_t   numPUplus;
      Float_t mass;
      Float_t ptgg;
      Float_t costheta;
      UInt_t  evt;
      UInt_t  run;
      UInt_t  lumi;
      UChar_t evtcat;
      
      TClonesArray photons;
      
      ClassDef(PhotonPairSelectorDiphotonEvent, 1)
    
  };
  
  class PhotonPairSelector : public BaseMod
  {
  public:
    PhotonPairSelector(const char *name ="PhotonPairSelector", 
		       const char *title="Selecting PhotonPairs");
    
    ~PhotonPairSelector();

    enum PhotonSelection {
      kNoPhSelection = 0,
      kCiCPhSelection,
      kMITPhSelection
    };
    enum VertexSelection {
      kStdVtxSelection = 0,
      kCiCVtxSelection,
      kMITVtxSelection
    };

    // setting all the input Names
    void                SetInputPhotonsName(const char *n){ fPhotonBranchName= n;        }
    void                SetPhotonsFromBranch(bool b)      { fPhotonsFromBranch = b;      }
    void                SetTrackName(const char *n)       { fTrackBranchName = n;        }
    void                SetElectronName(const char *n)    { fElectronName = n;           }
    void                SetConversionName(const char *n)  { fConversionName = n;         }
    void                SetPUDensityName(const char *n)   { fPileUpDenName = n;          }
    void                SetPVName(const char *n)          { fPVName = n;                 }
    void                SetPVFromBranch(bool b)           { fPVFromBranch = b;           }
    void                SetMCParticle(const char *n)      { fMCParticleName = n;         }
    void                SetPUInfoName(const char *n)      { fPileUpName = n;             }
    void                SetBeamspotName(const char *n)    { fBeamspotName = n;           }
    void                SetPFCandName(const char *n)      { fPFCandName = n;             }

    // set the type of selection
    void                SetPhotonSelType(const char *type){ fPhotonSelType    = type;    }
    void                SetVertexSelType(const char *type){ fVertexSelType    = type;    }

    // get/set the Names for the output Photon Collection
    const char         *GetOutputName()             const { return fGoodPhotonsName;     }   
    void                SetOutputName(const char *n)      { fGoodPhotonsName=n;          }    

    // set basic Cut variables (FOR PRE-SELECTION)
    void                SetPtMin(Double_t pt)             { fPhotonPtMin     = pt;       }
    void                SetAbsEtaMax(Double_t eta)        { fPhotonEtaMax    = eta;	 }
    
    void                SetLeadingPtMin(Double_t pt)      { fLeadingPtMin = pt;          }
    void                SetTrainingPtMin(Double_t pt)     { fTrailingPtMin = pt;         }

    // is Data Or Not?
    void                SetIsData (Bool_t b) { fIsData = b;};
    
    // methods to set the MC smearing/energy correction values
    void                AddEnCorrPerRun( UInt_t sRun, UInt_t eRun,
					 Double_t corr_EB_hR9,
					 Double_t corr_EB_lR9,
					 Double_t corr_EE_hR9,
					 Double_t corr_EE_lR9) {

      fDataEnCorr_EB_hR9.push_back(corr_EB_hR9);
      fDataEnCorr_EB_lR9.push_back(corr_EB_lR9);
      fDataEnCorr_EE_hR9.push_back(corr_EE_hR9);
      fDataEnCorr_EE_lR9.push_back(corr_EE_lR9);      
      fRunStart.push_back         (sRun);
      fRunEnd.push_back           (eRun);
    };

    void                SetMCSmearFactors(Double_t _EB_hR9, 
					  Double_t _EB_lR9, 
					  Double_t _EE_hR9, 
					  Double_t _EE_lR9) {
      fMCSmear_EB_hR9 = _EB_hR9;
      fMCSmear_EB_lR9 = _EB_lR9;
      fMCSmear_EE_hR9 = _EE_hR9;
      fMCSmear_EE_lR9 = _EE_lR9;
    };

    void                ApplyEleVeto(bool a)            { fApplyEleVeto  = a; }
    void                DoDataEneCorr(bool a)           { fDoDataEneCorr = a; }
    void                DoMCSmear(bool a)               { fDoMCSmear     = a; }

    void                SetTupleName(const char* c)     { fTupleName     = c; }

  protected:
    void                Process();
    void                SlaveBegin();

    // private auxiliary methods...
    void                FindHiggsPtAndZ(Float_t& pt, Float_t& z);
    Int_t               FindRunRangeIdx(UInt_t run);
    Double_t            GetDataEnCorr(Int_t runRange, PhotonTools::CiCBaseLineCats cat);
    Double_t            GetMCSmearFac(PhotonTools::CiCBaseLineCats cat);
    Float_t             GetEventCat(PhotonTools::CiCBaseLineCats cat1, PhotonTools::CiCBaseLineCats cat2);
    const MCParticle   *MatchMC(const Photon *ph) const;

    // Names for the input Collections
    TString             fPhotonBranchName;
    TString             fElectronName;
    TString             fConversionName;
    TString             fTrackBranchName;
    TString             fPileUpDenName;    
    TString             fPVName;
    TString             fBeamspotName;
    TString             fPFCandName;
    TString             fMCParticleName;
    TString             fPileUpName;
    
    TString             fGoodPhotonsName;      //name of exported "good photon" collection
    
    // Selection Types
    TString             fPhotonSelType;
    TString             fVertexSelType;
    PhotonSelection     fPhSelType;
    VertexSelection     fVtxSelType;    

    // Basic Pre-Selection kinematics
    Double_t            fPhotonPtMin;          // min pt cut fro PRE-SELECTION!
    Double_t            fPhotonEtaMax;         // max eta cut for PRE-SELECTION!
    
    Double_t            fLeadingPtMin;
    Double_t            fTrailingPtMin;

    // is it Data or MC?
    Bool_t              fIsData;
    
    // in case there's some PV pre-selection
    Bool_t              fPhotonsFromBranch;
    Bool_t              fPVFromBranch;

    const PhotonCol              *fPhotons;
    const ElectronCol            *fElectrons;
    const DecayParticleCol       *fConversions;
    const TrackCol               *fTracks;
    const PileupEnergyDensityCol *fPileUpDen;
    const VertexCol              *fPV;
    const BeamSpotCol            *fBeamspot;
    const PFCandidateCol         *fPFCands;
    const MCParticleCol          *fMCParticles;
    const PileupInfoCol          *fPileUp;    

    // Vectroes to hols smeraring/correction factors
    std::vector<Double_t> fDataEnCorr_EB_hR9;
    std::vector<Double_t> fDataEnCorr_EB_lR9;
    std::vector<Double_t> fDataEnCorr_EE_hR9;
    std::vector<Double_t> fDataEnCorr_EE_lR9;
    
    std::vector<UInt_t> fRunStart;
    std::vector<UInt_t> fRunEnd;
    
    Double_t fMCSmear_EB_hR9;
    Double_t fMCSmear_EB_lR9;
    Double_t fMCSmear_EE_hR9;
    Double_t fMCSmear_EE_lR9;

    // pointer to RNG ionstance for smearing
    TRandom3* rng;
        
    // --------------------------------
    // some streagin flags, not adjustable yet (FIX-ME)
    bool fDoDataEneCorr;
    bool fDoMCSmear;
    bool fDoVtxSelection;
    bool fApplyEleVeto;

    // --------------------------------
    // validation Tuple
    TString fTupleName;
    PhotonPairSelectorDiphotonEvent* fDiphotonEvent;
    TTree* hCiCTuple;

    ClassDef(PhotonPairSelector, 1) // Photon identification module
  };
}
#endif
