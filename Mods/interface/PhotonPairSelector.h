//--------------------------------------------------------------------------------------------------
// M.Yang 2011/10/12
// $Id: PhotonPairSelector.h,v 1.30 2012/08/02 12:59:31 fabstoec Exp $
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
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"

#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/EGEnergyCorrector.h"

#include "MitPhysics/Utils/interface/MVATools.h"
#include "MitPhysics/Utils/interface/VertexTools.h"

#include "MitPhysics/Utils/interface/MVAMet.h"
#include "MitPhysics/Utils/interface/RhoUtilities.h"


class TNtuple;
class TRandom3;
class TH1D;

namespace mithep 
{  
  class PhotonPairSelector : public BaseMod
  {
  public:
    PhotonPairSelector(const char *name  = "PhotonPairSelector", 
		       const char *title = "Selecting PhotonPairs");
    
    ~PhotonPairSelector();

    enum PhotonSelection {
      kNoPhSelection = 0,
      kCiCPhSelection,
      kCiCPFPhSelection,      
      kMVAPhSelection, //MVA
      kMITPhSelection,
      kMITPFPhSelection
    };
    enum VertexSelection {
      kStdVtxSelection = 0,
      kCiCVtxSelection,
      kMITVtxSelection,
      kCiCMVAVtxSelection,
      kMetSigVtxSelection
    };

    // outsourced to MVATools (IdMVAType)  (fab)
/*     enum IdMVA { */
/*       k2011IdMVA = 0, */
/*       k2012IdMVA_globe */
/*     }; */
   

    // outsourced to PhotonTools (fab)
/*     enum ShowerShape { */
/*       k2011ShowerShape  = 0, */
/*       k2012ShowerShape */
/*     }; */

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
    void                SetIdMVAType(const char *type)    { fIdMVAType        = type;    } 
    void                SetShowerShapeType(const char *type) { fShowerShapeType        = type;    } 

    // get/set the Names for the output Photon Collection
    const char         *GetOutputName()             const { return fGoodPhotonsName;     }   
    void                SetOutputName(const char *n)      { fGoodPhotonsName=n;          }    

    const char         *GetOutputVtxName()          const { return fChosenVtxName;       }
    void                SetOutputVtxName(const char* n)   { fChosenVtxName = n;          }


    // set basic Cut variables (FOR PRE-SELECTION)
    void                SetPtMin(Double_t pt)             { fPhotonPtMin     = pt;       }
    void                SetAbsEtaMax(Double_t eta)        { fPhotonEtaMax    = eta;	 }
    
    void                SetLeadingPtMin(Double_t pt)      { fLeadingPtMin = pt;          }
    void                SetTrailingPtMin(Double_t pt)     { fTrailingPtMin = pt;         }

    // is Data Or Not?
    void                SetIsData (Bool_t b) { fIsData = b;}

    void                Set2012HCP (Bool_t b) { f2012HCP = b;};

    
    // methods to set the MC smearing/energy correction values
    void                AddEnCorrPerRun( UInt_t sRun, UInt_t eRun,
					 Double_t corr_EBlowEta_hR9central,
                                         Double_t corr_EBlowEta_hR9gap,
					 Double_t corr_EBlowEta_lR9,
                                         Double_t corr_EBhighEta_hR9,
                                         Double_t corr_EBhighEta_lR9,                                         
					 Double_t corr_EElowEta_hR9,
					 Double_t corr_EElowEta_lR9,
                                         Double_t corr_EEhighEta_hR9,
                                         Double_t corr_EEhighEta_lR9) {

      fDataEnCorr_EBlowEta_hR9central.push_back(corr_EBlowEta_hR9central);
      fDataEnCorr_EBlowEta_hR9gap.push_back(corr_EBlowEta_hR9gap);
      fDataEnCorr_EBlowEta_lR9.push_back(corr_EBlowEta_lR9);
      fDataEnCorr_EBhighEta_hR9.push_back(corr_EBhighEta_hR9);
      fDataEnCorr_EBhighEta_lR9.push_back(corr_EBhighEta_lR9);      
      fDataEnCorr_EElowEta_hR9.push_back(corr_EElowEta_hR9);
      fDataEnCorr_EElowEta_lR9.push_back(corr_EElowEta_lR9);      
      fDataEnCorr_EEhighEta_hR9.push_back(corr_EEhighEta_hR9);
      fDataEnCorr_EEhighEta_lR9.push_back(corr_EEhighEta_lR9);       
      fRunStart.push_back         (sRun);
      fRunEnd.push_back           (eRun);
    };

    void                SetMCSmearFactors(Double_t _EBlowEta_hR9central,
                                          Double_t _EBlowEta_hR9gap, 
					  Double_t _EBlowEta_lR9,
                                          Double_t _EBhighEta_hR9, 
                                          Double_t _EBhighEta_lR9,
					  Double_t _EElowEta_hR9, 
					  Double_t _EElowEta_lR9,
                                          Double_t _EEhighEta_hR9,
                                          Double_t _EEhighEta_lR9) {
      fMCSmear_EBlowEta_hR9central = _EBlowEta_hR9central;
      fMCSmear_EBlowEta_hR9gap = _EBlowEta_hR9gap;
      fMCSmear_EBlowEta_lR9 = _EBlowEta_lR9;
      fMCSmear_EBhighEta_hR9 = _EBhighEta_hR9;
      fMCSmear_EBhighEta_lR9 = _EBhighEta_lR9;      
      fMCSmear_EElowEta_hR9 = _EElowEta_hR9;
      fMCSmear_EElowEta_lR9 = _EElowEta_lR9;
      fMCSmear_EEhighEta_hR9 = _EEhighEta_hR9;
      fMCSmear_EEhighEta_lR9 = _EEhighEta_lR9;      
    };


    void                AddEnCorrPerRun2012HCP( UInt_t sRun, UInt_t eRun,
						Double_t corr_EBlowEta_hR9central,
						Double_t corr_EBlowEta_hR9gap,
						Double_t corr_EBlowEta_lR9central,
						Double_t corr_EBlowEta_lR9gap,
						Double_t corr_EBhighEta_hR9,
						Double_t corr_EBhighEta_lR9,                                         
						Double_t corr_EElowEta_hR9,
						Double_t corr_EElowEta_lR9,
						Double_t corr_EEhighEta_hR9,
						Double_t corr_EEhighEta_lR9) {
      fDataEnCorr_EBlowEta_hR9central.push_back(corr_EBlowEta_hR9central);
      fDataEnCorr_EBlowEta_hR9gap.push_back(corr_EBlowEta_hR9gap);
      fDataEnCorr_EBlowEta_lR9central.push_back(corr_EBlowEta_lR9central);
      fDataEnCorr_EBlowEta_lR9gap.push_back(corr_EBlowEta_lR9gap);
      fDataEnCorr_EBhighEta_hR9.push_back(corr_EBhighEta_hR9);
      fDataEnCorr_EBhighEta_lR9.push_back(corr_EBhighEta_lR9);      
      fDataEnCorr_EElowEta_hR9.push_back(corr_EElowEta_hR9);
      fDataEnCorr_EElowEta_lR9.push_back(corr_EElowEta_lR9);      
      fDataEnCorr_EEhighEta_hR9.push_back(corr_EEhighEta_hR9);
      fDataEnCorr_EEhighEta_lR9.push_back(corr_EEhighEta_lR9);       
      fRunStart.push_back         (sRun);
      fRunEnd.push_back           (eRun);
    };
    
    void SetMCSmearFactors2012HCP(Double_t _EBlowEta_hR9central,
				  Double_t _EBlowEta_hR9gap, 
				  Double_t _EBlowEta_lR9central,
				  Double_t _EBlowEta_lR9gap,
				  Double_t _EBhighEta_hR9, 
				  Double_t _EBhighEta_lR9,
				  Double_t _EElowEta_hR9, 
				  Double_t _EElowEta_lR9,
				  Double_t _EEhighEta_hR9,
				  Double_t _EEhighEta_lR9) {
      fMCSmear_EBlowEta_hR9central = _EBlowEta_hR9central;
      fMCSmear_EBlowEta_hR9gap = _EBlowEta_hR9gap;
      fMCSmear_EBlowEta_lR9central = _EBlowEta_lR9central;
      fMCSmear_EBlowEta_lR9gap = _EBlowEta_lR9gap;
      fMCSmear_EBhighEta_hR9 = _EBhighEta_hR9;
      fMCSmear_EBhighEta_lR9 = _EBhighEta_lR9;      
      fMCSmear_EElowEta_hR9 = _EElowEta_hR9;
      fMCSmear_EElowEta_lR9 = _EElowEta_lR9;
      fMCSmear_EEhighEta_hR9 = _EEhighEta_hR9;
      fMCSmear_EEhighEta_lR9 = _EEhighEta_lR9;      
    };
        
    void                SetApplyEleVeto(bool a)            { fApplyEleVeto  = a; }
    void                SetInvertElectronVeto(Bool_t b)   { fInvertElectronVeto = b;     }          
    void                DoDataEneCorr(bool a)           { fDoDataEneCorr = a; }
    void                DoMCSmear(bool a)               { fDoMCSmear     = a; }

    void                SetGoodElectronsFromBranch(Bool_t b) { fGoodElectronsFromBranch = b; }
    void                SetGoodElectronName(TString name) { fGoodElectronName = name; }
    void                SetUseSingleLegConversions(Bool_t b) { fUseSingleLegConversions = b; }
    void                SetDoRegression(Bool_t b)         { fDoRegression = b; }
    void                SetEtaCorrections(const TH1D *h)  { fEtaCorrections = h; }
    void                SetBdtCutBarrel(Float_t x)        { fbdtCutBarrel = x; }
    void                SetBdtCutEndcap(Float_t x)        { fbdtCutEndcap = x; }
 
    void                SetDoShowerShapeScaling(Bool_t b)  { fDoShowerShapeScaling = b; }
    
    void                SetJetsName(const char *n)       { fJetsName = n;              }    
        
    void                SetRhoType(RhoUtilities::RhoType type) { fRhoType = type ; }

    void                SetApplyLeptonTag(bool a)            {  fApplyLeptonTag = a; }

    void                SetLeptonTagElectronsName(TString name) { fLeptonTagElectronsName = name; }
    void                SetLeptonTagMuonsName    (TString name) { fLeptonTagMuonsName     = name; }

  protected:
    void                Process();
    void                SlaveBegin();

    // private auxiliary methods...
    void                FindHiggsPtAndZ(Float_t& pt, Float_t& z, Float_t& mass);
    Int_t               FindRunRangeIdx(UInt_t run);
    Double_t            GetDataEnCorr(Int_t runRange, PhotonTools::eScaleCats cat);
    Double_t            GetMCSmearFac(PhotonTools::eScaleCats cat);
    Double_t            GetDataEnCorrHCP(Int_t runRange, PhotonTools::eScaleCats cat);
    Double_t            GetMCSmearFacHCP(PhotonTools::eScaleCats cat);
    Float_t             GetEventCat(PhotonTools::CiCBaseLineCats cat1, PhotonTools::CiCBaseLineCats cat2);

    // Names for the input Collections
    TString             fPhotonBranchName;
    TString             fElectronName;
    TString             fGoodElectronName;
    TString             fConversionName;
    TString             fPFConversionName;    
    TString             fTrackBranchName;
    TString             fPileUpDenName;    
    TString             fPVName;
    TString             fBeamspotName;
    TString             fPFCandName;
    TString             fMCParticleName;
    TString             fPileUpName;
    TString             fJetsName;
    TString             fPFMetName;
    
    TString             fGoodPhotonsName;      //name of exported "good photon" collection
    TString             fChosenVtxName;        //name of exported "chosen Vtx"  collection

    TString             fLeptonTagElectronsName;
    TString             fLeptonTagMuonsName;
    
    // Selection Types
    TString             fPhotonSelType;
    TString             fVertexSelType;
    PhotonSelection     fPhSelType;
    VertexSelection     fVtxSelType;    

    // Id Type
    TString             fIdMVAType;
    MVATools::IdMVAType fIdType;
    //IdMVA               fIdType;

    // showershape
    TString                                fShowerShapeType;
    PhotonTools::ShowerShapeScales         fSSType;

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
    Bool_t              fGoodElectronsFromBranch;
    Bool_t              fUseSingleLegConversions;

    Bool_t              f2012HCP;

    const PhotonCol              *fPhotons;
    const ElectronCol            *fElectrons;
    const ElectronCol            *fGoodElectrons;    
    const DecayParticleCol       *fConversions;
    const DecayParticleCol       *fPFConversions;
    const TrackCol               *fTracks;
    const PileupEnergyDensityCol *fPileUpDen;
    const VertexCol              *fPV;
    const BeamSpotCol            *fBeamspot;
    const PFCandidateCol         *fPFCands;
    const MCParticleCol          *fMCParticles;
    const PileupInfoCol          *fPileUp;    
    const JetCol                 *fJets;
    const PFMetCol               *fPFMet;
    const ElectronCol            *fLeptonTagElectrons;
    const MuonCol                *fLeptonTagMuons;
    
    // Vectroes to hols smeraring/correction factors
    std::vector<Double_t> fDataEnCorr_EBlowEta_hR9central;
    std::vector<Double_t> fDataEnCorr_EBlowEta_hR9gap;
    std::vector<Double_t> fDataEnCorr_EBlowEta_lR9;
    std::vector<Double_t> fDataEnCorr_EBlowEta_lR9central;
    std::vector<Double_t> fDataEnCorr_EBlowEta_lR9gap;
    std::vector<Double_t> fDataEnCorr_EBhighEta_hR9;
    std::vector<Double_t> fDataEnCorr_EBhighEta_lR9;    
    std::vector<Double_t> fDataEnCorr_EElowEta_hR9;
    std::vector<Double_t> fDataEnCorr_EElowEta_lR9;
    std::vector<Double_t> fDataEnCorr_EEhighEta_hR9;
    std::vector<Double_t> fDataEnCorr_EEhighEta_lR9;
    
    std::vector<UInt_t>   fRunStart;
    std::vector<UInt_t>   fRunEnd;
    
    Double_t              fMCSmear_EBlowEta_hR9central;
    Double_t              fMCSmear_EBlowEta_hR9gap;
    Double_t              fMCSmear_EBlowEta_lR9;
    Double_t              fMCSmear_EBlowEta_lR9central;
    Double_t              fMCSmear_EBlowEta_lR9gap;
    Double_t              fMCSmear_EBhighEta_hR9;
    Double_t              fMCSmear_EBhighEta_lR9;    
    Double_t              fMCSmear_EElowEta_hR9;
    Double_t              fMCSmear_EElowEta_lR9;
    Double_t              fMCSmear_EEhighEta_hR9;
    Double_t              fMCSmear_EEhighEta_lR9;    

    // pointer to RNG ionstance for smearing
    TRandom3*             fRng;
    EGEnergyCorrector     fEgCor;
    Bool_t                fDoRegression;
    TString               fPhFixString;
    TString               fPhFixFile;
    TString               fRegWeights; 
  		         
    const TH1D           *fEtaCorrections;
    
    // --------------------------------
    // some streagin flags, not adjustable yet (FIX-ME)
    bool                  fDoDataEneCorr;
    bool                  fDoMCSmear;
    bool                  fDoVtxSelection;
    bool                  fApplyEleVeto;
    Bool_t                fInvertElectronVeto; //=true (invert ele veto, for cic sel only atm)
    
    // --------------------------------
    bool                  fApplyLeptonTag;   

    //MVA
    int                   fVariableType_2011;
    TString               fEndcapWeights_2011;
    TString               fBarrelWeights_2011;
    int                   fVariableType_2012_globe;
    TString               fEndcapWeights_2012_globe;
    TString               fBarrelWeights_2012_globe;  
    MVATools              fTool;
    Float_t               fbdtCutBarrel;
    Float_t               fbdtCutEndcap;

    VertexTools           fVtxTools;
    
    Bool_t                fDoShowerShapeScaling; 
       	                  
    Bool_t                fDoMCErrScaling;
    Double_t              fMCErrScaleEB;
    Double_t              fMCErrScaleEE;    
    UInt_t                fRegressionVersion;
    	                  
    Bool_t                fRelativePtCuts;
    
    MVAMet                fMVAMet;
    
    RhoUtilities::RhoType fRhoType;

    ClassDef(PhotonPairSelector, 1) // Photon identification module
  };
}
#endif
