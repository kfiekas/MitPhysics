//--------------------------------------------------------------------------------------------------
// $Id: PhotonPairSelector.h,v 1.5 2011/07/06 13:59:40 fabstoec Exp $
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
    
    
  protected:
    void                Process();
    void                SlaveBegin();

    // private auxiliary methods...
    void                FindHiggsPtAndZ(Float_t& pt, Float_t& z);
    Int_t               FindRunRangeIdx(UInt_t run);
    Double_t            GetDataEnCorr(Int_t runRange, PhotonTools::CiCBaseLineCats cat);
    Double_t            GetMCSmearFac(PhotonTools::CiCBaseLineCats cat);
    Float_t             GetEventCat(PhotonTools::CiCBaseLineCats cat1, PhotonTools::CiCBaseLineCats cat2);

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

    // --------------------------------
    // validation Tuple
    TNtuple* hCiCTuple;

    ClassDef(PhotonPairSelector, 1) // Photon identification module
  };
}
#endif
