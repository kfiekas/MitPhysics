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

class TNtuple;
class TRandom3;

namespace mithep 
{
    // This class holds all of the ZgllTuple variables. Be sure to include it in MitPhysicsModLinks.h
    class LeptonPairPhotonEvent
  {
    public:
    Float_t Meeg;// 3-body mass with two electrons
    Float_t Mmmg;// 3-body mass with two muons
    Float_t Mee;// 2-body mass with two electrons
    Float_t Mmm;// 2-body mass with two muons
    Float_t ele1MVA;// Electron 1 MVA Value
    Float_t ele1charge;// Electron 1 Charge
    Float_t ele1energy;// Electron 1 Energy
    Float_t ele1px;// Electron 1 Px
    Float_t ele1py;// Electron 1 Py
    Float_t ele1pz;// Electron 1 Pz
    Float_t ele1pt;// Electron 1 Pt
    Float_t ele1eta;// Electron 1 Eta
    Float_t ele1mass;// Electron 1 Mass
    Float_t ele1phi;// Electron 1 Phi
   
    Float_t ele2MVA;// Electron 2 MVA Value
    Float_t ele2charge;// Electron 2 Charge
    Float_t ele2energy;// Electron 2 Energy
    Float_t ele2px;// Electron 2 Px
    Float_t ele2py;// Electron 2 Py
    Float_t ele2pz;// Electron 2 Pz
    Float_t ele2pt;// Electron 2 Pt
    Float_t ele2eta;// Electron 2 Eta
    Float_t ele2mass;// Electron 2 Mass
    Float_t ele2phi;// Electron 2 Phi

    Float_t photonisgen;// Photon MVA Value
    Float_t photonidmva;// Photon MVA Value
    Float_t photonenergy;// Photon Energy
    Float_t photonpx;// Photon Px
    Float_t photonpy;// Photon Py
    Float_t photonpz;// Photon Pz
    Float_t photonpt;// Photon Pt
    Float_t photoneta;// Photon Eta
    Float_t photonmass;// Photon Mass??? photon->Mass() 
    Float_t photonphi;// Photon Phi
    Float_t photonr9;// Photon R9

    Float_t m1E;// Muon 1 Energy
    Float_t m1Pt;// Muon 1 Pt
    Float_t m1Mass;// Muon 1 Mass
    Float_t m1Px;// Muon 1 Px
    Float_t m1Py;// Muon 1 Py
    Float_t m1Pz;// Muon 1 Pz
    Float_t m1Eta;// Muon 1 Eta
    Float_t m1Phi;// Muon 1 Phi
    Float_t m1Charge;// Muon 1 Charge
    Float_t m2E;// Muon 2 Energy
    Float_t m2Pt;// Muon 2 Pt
    Float_t m2Mass;// Muon 2 Mass
    Float_t m2Px;// Muon 2 Px
    Float_t m2Py;// Muon 2 Py
    Float_t m2Pz;// Muon 2 Pz
    Float_t m2Eta;// Muon 2 Eta
    Float_t m2Phi;// Muon 2 Phi
    Float_t m2Charge;// Muon 2 Charge

    Int_t NPu; //Number of Pileup events
    Int_t NPuPlus;//Number of Pileup events in next signal readout
    Int_t NPuMinus;//Number of Pileup events in previous signal readout
   

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

    // is Data Or Not?
    void                SetIsData (Bool_t b)                 { fIsData    = b; }
    
    void                SetDoMCCheck (Bool_t b)              { fDoMCCheck = b; }

    void                SetTupleName(const char* c)          { fTupleName = c; }
    
    // The following flag configures LeptonPairPhotonTreeWriter to select either dielectron or dimuon events
    void		SetAreElectronsOrMuons (Int_t i)    { fElectronMuonFlag = i;} 

  protected:
    void                Process();
    void                SlaveBegin();

    // Names for the input Collections
    TString             fPhotonBranchName;
    TString             fGoodElectronName;
    TString		fGoodMuonName;
    TString             fPVName;
    TString             fVertexName;
    TString             fPFCandName;
    TString             fTrackName;
    TString             fPileUpDenName; 
    TString             fPileUpName;

    // Is it Data or MC?
    Bool_t              fIsData;
    Bool_t              fDoMCCheck;
    
    // Dielectron or dimuon flag
    Int_t		fElectronMuonFlag;    
    
    // Determines whether the input comes from a branch or not
    Bool_t              fPhotonsFromBranch;
    Bool_t              fPVFromBranch;
    Bool_t              fGoodElectronsFromBranch;
    Bool_t		fGoodMuonsFromBranch;
    Bool_t		fPFCandsFromBranch;
    Bool_t		fTracksFromBranch;
    Bool_t		fPileUpDenFromBranch;

    //Collection Names
    const PhotonCol               *fPhotons;
    const ElectronCol             *fGoodElectrons;    
    const VertexCol               *fPV;
    const MuonCol		  *fGoodMuons;
    const VertexCol               *fVertices;
    const PFCandidateCol	  *fPFCands;
    const TrackCol		  *fTracks;
    const PileupEnergyDensityCol  *fPileUpDen;
    const PileupInfoCol           *fPileUp;

    const MCParticleCol*           fMCParticles;

    // --------------------------------
    TString                        fTupleName;
    LeptonPairPhotonEvent*         fLeptonPairPhotonEvent;
    TTree*                         ZgllTuple;
    
    //Photon MVA initialization variables
/*     int                   fVariableType_2011; */
/*     TString               fEndcapWeights_2011; */
/*     TString               fBarrelWeights_2011; */
/*     int                   fVariableType_2012_globe; */
/*     TString               fEndcapWeights_2012_globe; */
/*     TString               fBarrelWeights_2012_globe; */
    //MVATools              fTool;
    ClassDef(LeptonPairPhotonTreeWriter, 1) // Lepton Pair + Photon identification module
  };

}
#endif
