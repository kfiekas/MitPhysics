//--------------------------------------------------------------------------------------------------
// $Id: LeptonPairPhotonTreeWriter.h,v 1.0 2012/06/23 21:25:01 auhess ksingh
//
// LeptonPairPhotonTreeWriter
//
// Authors: A. Hess & K. Singh
//--------------------------------------------------------------------------------------------------


#ifndef MITPHYSICS_MODS_LEPTONPAIRPHOTONTREEWRITER_H
#define MITPHYSICS_MODS_LEPTONPAIRPHOTONTREEWRITER_H

#include "MitPhysics/Utils/interface/MVATools.h"
#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PhotonFwd.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitPhysics/Utils/interface/PhotonFix.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/MVAMet.h"
#include "MitPhysics/Utils/interface/ElectronEnergyRegression.h"

#include "MitPhysics/Utils/interface/rochcor.h"



//
// correction header is passed in as a define so as not to make 
// others checkout UserCode. Set the env var using this syntax : 
//     export PHOSPHOR_CORRECTIONS_HEADER='\"MyPhosphorDir/PhosphorCorrectorFunctor.hh\"'
// proper quoting is important ...
#ifdef PHOSPHOR_CORRECTIONS_HEADER
#include PHOSPHOR_CORRECTIONS_HEADER
using namespace zgamma;
#endif


class TNtuple;
class TRandom3;


namespace mithep 
{
    // This class holds all of the ZgllTuple variables. Be sure to include it in MitPhysicsModLinks.h
    class LeptonPairPhotonEvent
  {
    public:
    Float_t electronZmass;
    Float_t mllg;// 3-body mass
    Float_t mllgCorr;// 3-body mass
    Bool_t muonZgVeto;
    Float_t muonZmass;
    Int_t year;   
    UInt_t run;
    UInt_t lumi;
    UInt_t event;

    Float_t ele1MVA;// Electron 1 MVA Value
    Float_t ele1charge;// Electron 1 Charge
    Float_t ele1energy;// Electron 1 Energy
    Float_t ele1px;// Electron 1 Px
    Float_t ele1py;// Electron 1 Py
    Float_t ele1pz;// Electron 1 Pz
    Float_t ele1pt;// Electron 1 Pt
    Float_t ele1eta;// Electron 1 Eta
    Float_t ele1SCeta;// Electron 1 SC Eta
    Float_t ele1mass;// Electron 1 Mass
    Float_t ele1phi;// Electron 1 Phi
    Float_t ele1dEtaIn;
    Float_t ele1dPhiIn;
    Float_t ele1sigmaIEtaIEta;
    Float_t ele1HadOverEm;
    Float_t ele1D0;
    Float_t ele1DZ;
    Float_t ele1OneOverEMinusOneOverP;
    Float_t ele1PFIsoOverPt;
    Bool_t  ele1Conversion;
    Float_t ele1missinghits;
    Float_t ele1RegressionEnergyV0;
    Float_t ele1RegressionEnergyV1;
    Float_t ele1RegressionEnergyV2;
    Float_t ele1RegressionEnergyErrorV0;
    Float_t ele1RegressionEnergyErrorV1;
    Float_t ele1RegressionEnergyErrorV2;

    // scale/res corrected
    Float_t ele1energyCorr;// Electron 1 Energy
    Float_t ele1pxCorr;// Electron 1 Px
    Float_t ele1pyCorr;// Electron 1 Py
    Float_t ele1pzCorr;// Electron 1 Pz
    Float_t ele1ptCorr;// Electron 1 Pt


    Float_t ele2MVA;// Electron 2 MVA Value
    Float_t ele2charge;// Electron 2 Charge
    Float_t ele2energy;// Electron 2 Energy
    Float_t ele2px;// Electron 2 Px
    Float_t ele2py;// Electron 2 Py
    Float_t ele2pz;// Electron 2 Pz
    Float_t ele2pt;// Electron 2 Pt
    Float_t ele2eta;// Electron 2 Eta
    Float_t ele2SCeta;// Electron 2 SC Eta
    Float_t ele2mass;// Electron 2 Mass
    Float_t ele2phi;// Electron 2 Phi
    Float_t ele2dEtaIn;
    Float_t ele2dPhiIn;
    Float_t ele2sigmaIEtaIEta;
    Float_t ele2HadOverEm;
    Float_t ele2D0;
    Float_t ele2DZ;
    Float_t ele2OneOverEMinusOneOverP;
    Float_t ele2PFIsoOverPt;
    Bool_t  ele2Conversion;
    Float_t ele2missinghits;
    Float_t ele2RegressionEnergyV0;
    Float_t ele2RegressionEnergyV1;
    Float_t ele2RegressionEnergyV2;
    Float_t ele2RegressionEnergyErrorV0;
    Float_t ele2RegressionEnergyErrorV1;
    Float_t ele2RegressionEnergyErrorV2;

    // scale/res corrected
    Float_t ele2energyCorr;// Electron 1 Energy
    Float_t ele2pxCorr;// Electron 1 Px
    Float_t ele2pyCorr;// Electron 1 Py
    Float_t ele2pzCorr;// Electron 1 Pz
    Float_t ele2ptCorr;// Electron 1 Pt

    Float_t chargediso_ele1;
    Float_t gammaiso_ele1;
    Float_t neutraliso_ele1;
    Float_t rho_ele1;
    Float_t effectivearea_ele1;
    Float_t chargediso_ele2;
    Float_t gammaiso_ele2;
    Float_t neutraliso_ele2;
    Float_t rho_ele2;
    Float_t effectivearea_ele2;

    Float_t photonidmva;// Photon MVA Value
    Float_t photonenergy;// Photon Energy
    Float_t photonpx;// Photon Px
    Float_t photonpy;// Photon Py
    Float_t photonpz;// Photon Pz
    Float_t photonpt;// Photon Pt
    Float_t photoneta;// Photon Eta
    Float_t photonSCeta;// Photon SC Eta
    Float_t photonmass;// Photon Mass??? photon->Mass() 
    Float_t photonphi;// Photon Phi
    Float_t photonr9;// Photon R9
    Float_t photonenergyerror;// Photon Energy Error

    Float_t photonenergyCorr;// Photon Energy
    Float_t photonpxCorr;// Photon Px
    Float_t photonpyCorr;// Photon Py
    Float_t photonpzCorr;// Photon Pz
    Float_t photonptCorr;// Photon Pt

    Float_t m1E;// Muon 1 Energy
    Float_t m1Pt;// Muon 1 Pt
    Float_t m1Mass;// Muon 1 Mass
    Float_t m1Px;// Muon 1 Px
    Float_t m1Py;// Muon 1 Py
    Float_t m1Pz;// Muon 1 Pz
    Float_t m1Eta;// Muon 1 Eta
    Float_t m1Phi;// Muon 1 Phi
    Float_t m1Charge;// Muon 1 Charge
    Float_t m1PtErr; // Muon 1 Pt Error
    Float_t m1ECorr;// Muon 1 Energy
    Float_t m1PtCorr;// Muon 1 Pt
    Float_t m1PxCorr;// Muon 1 Px
    Float_t m1PyCorr;// Muon 1 Py
    Float_t m1PzCorr;// Muon 1 Pz

    Float_t m2E;// Muon 2 Energy
    Float_t m2Pt;// Muon 2 Pt
    Float_t m2Mass;// Muon 2 Mass
    Float_t m2Px;// Muon 2 Px
    Float_t m2Py;// Muon 2 Py
    Float_t m2Pz;// Muon 2 Pz
    Float_t m2Eta;// Muon 2 Eta
    Float_t m2Phi;// Muon 2 Phi
    Float_t m2Charge;// Muon 2 Charge
    Float_t m2PtErr; // Muon 1 Pt Error
    Float_t m2ECorr;// Muon 1 Energy
    Float_t m2PtCorr;// Muon 1 Pt
    Float_t m2PxCorr;// Muon 1 Px
    Float_t m2PyCorr;// Muon 1 Py
    Float_t m2PzCorr;// Muon 1 Pz

    Float_t NPu; //Number of Pileup events
    Float_t NPuPlus;//Number of Pileup events in next signal readout
    Float_t NPuMinus;//Number of Pileup events in previous signal readout
    Float_t photonmatchmc;   

    Float_t costheta_lm_electrons;
    Float_t costheta_lp_electrons;
    Float_t phi_electrons;
    Float_t cosTheta_electrons;
    Float_t cosThetaG_electrons;
    Float_t costheta_lm_muons;
    Float_t costheta_lp_muons;
    Float_t phi_muons;
    Float_t cosTheta_muons;
    Float_t cosThetaG_muons;
  };

  class LeptonPairPhotonTreeWriter : public BaseMod
  {
  public:
    LeptonPairPhotonTreeWriter(const char *name ="LeptonPairPhotonTreeWriter", 
			       const char *title="Selecting PhotonPairs");
    
    ~LeptonPairPhotonTreeWriter();

    // setting all the input Names
    void                SetInputPhotonsName(const char *n){ fPhotonBranchName= n;        }
    void                SetPhotonsFromBranch(bool b)      { fPhotonsFromBranch = b;      }

    void                SetGoodElectronName(TString name)    { fGoodElectronName = name; }
    void                SetGoodElectronsFromBranch(Bool_t b) { fGoodElectronsFromBranch = b; }

    void                SetPVName(const char *n)          { fPVName = n;                 }
    void                SetPVFromBranch(bool b)           { fPVFromBranch = b;           }
    void                SetPFCandName(const char *n)          { fPFCandName = n;                 }
    void                SetPFCandsFromBranch(bool b)           { fPFCandsFromBranch = b;           }
    void                SetTrackName(const char *n)          { fTrackName = n;                 }
    void                SetTracksFromBranch(bool b)           { fTracksFromBranch = b;           }
    void                SetPileUpDenName(const char *n)          { fPileUpDenName = n;                 }
    void                SetPileUpDenFromBranch(bool b)           { fPileUpDenFromBranch = b;           }
    void 		SetGoodMuonName(TString name) { fGoodMuonName = name; }
    void 		SetGoodMuonsFromBranch(Bool_t b) { fGoodMuonsFromBranch = b; }
    void                SetPUInfoName(const char *n)      { fPileUpName = n;             }
    void                SetApplyElectronVeto(Bool_t b)   { fApplyElectronVeto = b;     }

    // is Data Or Not?
    void                SetIsData (Bool_t b)                 { fIsData = b; }
    
    void                SetTupleName(const char* c)          { fTupleName = c; }
    void		SetYear(Float_t y)		     {YEAR = y;} 
    void		SetVerbose(bool b)		     {verbose = b;} 
    void		SetDoElectronChannel(bool b)	     {_do_ElectronChannel = b;}
    void		SetDoMuonChannel(bool b)	     {_do_MuonChannel = b;} 

    void                SetPhosphorDataFile(TString s)       {phosphorDataFile = s;}

    // wrappers to keep things readable
    void fillEle1Variables(const mithep::Electron * ele1);  
    void fillEle2Variables(const mithep::Electron * ele2);  
    void resetTreeVariables();
    void regressEle1(const mithep::Electron              *ele1,
		     const mithep::PileupEnergyDensityCol *fPileUpDen,
		     const mithep::VertexCol              *fPV   );    
    void regressEle2(const mithep::Electron              *ele2,
		     const mithep::PileupEnergyDensityCol *fPileUpDen,
		     const mithep::VertexCol              *fPV   );    

    struct PhotonPtComparison { 
      bool operator()( const Photon* const  &l1,  const Photon* const &l2 ) { 
	if( l1->Pt() > l2->Pt() ) return true;
	else return false;
      }
    };

  protected:
    void                Process();
    void                SlaveBegin();
    // Names for the input Collections
    TString             fPhotonBranchName;
    TString             fGoodElectronName;
    TString		fGoodMuonName;
    TString             fPVName;
    TString             fPFCandName;
    TString             fTrackName;
    TString             fPileUpDenName; 
    TString             fPileUpName;
    TString 		fMCParticleName;
    TString	    	fConversionName;
    TString             fBeamSpotName;

    // Is it Data or MC?
    Bool_t              fIsData;
    Float_t 		YEAR; 
    
    // Determines whether the input comes from a branch or not
    Bool_t              fPhotonsFromBranch;
    Bool_t              fPVFromBranch;
    Bool_t              fGoodElectronsFromBranch;
    Bool_t		fGoodMuonsFromBranch;
    Bool_t		fPFCandsFromBranch;
    Bool_t		fTracksFromBranch;
    Bool_t		fPileUpDenFromBranch;
    Bool_t		fApplyElectronVeto;
    //Collection Names
    const PhotonCol               *fPhotons;
    const ElectronCol             *fGoodElectrons;    
    const VertexCol               *fPV;
    const MuonCol		  *fGoodMuons;
    const PFCandidateCol	  *fPFCands;
    const TrackCol		  *fTracks;
    const PileupEnergyDensityCol  *fPileUpDen;
    const PileupInfoCol           *fPileUp;
    const MCParticleCol		  *fMCParticles;
    const DecayParticleCol	  *fConversions;    
    const BeamSpotCol             *fBeamSpot;

// --------------------------------
    TString                        fTupleName;
    LeptonPairPhotonEvent*         fLeptonPairPhotonEvent;
    TTree*                         ZgllTuple;
    //Photon MVA initialization variables
    int                   fVariableType_2011;
    TString               fEndcapWeights_2011;
    TString               fBarrelWeights_2011;
    int                   fVariableType_2012_globe;
    TString               fEndcapWeights_2012_globe;
    TString               fBarrelWeights_2012_globe;
    MVATools              fTool;
    ElectronEnergyRegression *eleRegressionEvaluator_V0;
    ElectronEnergyRegression *eleRegressionEvaluator_V1;
    ElectronEnergyRegression *eleRegressionEvaluator_V2;

    bool verbose;
    bool _do_ElectronChannel, _do_MuonChannel;

    TRandom3 rand;
    rochcor * rmcor;

    TString phosphorDataFile;
#ifdef PHOSPHOR_CORRECTIONS_HEADER
    PhosphorCorrectionFunctor* phosphor;
#endif

    ClassDef(LeptonPairPhotonTreeWriter, 1) // Lepton Pair + Photon identification module
  };


}
#endif

