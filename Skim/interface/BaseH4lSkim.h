//--------------------------------------------------------------------------------------------------
// $Id: BaseH4lSkim.h,v 1.3 2012/06/21 21:37:46 dkralph Exp $
//
// BaseH4lSkim
//
// This module defines a common set of basic skimming actions for the H4l skims
//
// Authors: D.Ralph, C.Paus
//--------------------------------------------------------------------------------------------------
#ifndef BASEH4LSKIM_H
#define BASEH4LSKIM_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"

using namespace std;

namespace mithep
{  
  class BaseH4lSkim : public BaseMod
  {    
  public:
    BaseH4lSkim(const char *name  = "BaseH4lSkim",
		const char *title = "Skim H4l Base Class");
    ~BaseH4lSkim();	    

    unsigned makePFnoPUArray();
    bool    SetBestPv         ();
    Float_t ComputePfMuonIso  (const Muon *muon, const Double_t dRMax);
    Float_t ComputePfElecIso  (const Electron *electron, const Double_t dRMax);
    bool    PassMuonPreselNoIp(const mithep::Muon *mu);
    bool    PassElecPreselNoIp(const mithep::Electron *ele);
    bool    muon2012CutBasedIDTight(const mithep::Muon *mu);
    bool    muon2012CutBasedIDTightForTagProbe(const mithep::Muon *mu);
    bool    electron2012CutBasedIDMedium(const mithep::Electron *ele);
    bool    electron2012CutBasedIDMediumForTagProbe(const mithep::Electron *ele);
  
  protected:
    void    Begin();
    void    BeginRun();
    void    EndRun();
    void    SlaveBegin();
    void    SlaveTerminate();
    void    Terminate();
    void    Process();

    UInt_t                        fTotal,fSelected;

    TString                       fMuonName;             // muon collection name
    TString                       fElectronName;         // electron collection name
    TString                       fPrimVtxName;          // primary vertex collection name
    TString                       fBeamSpotName;         // pointer to beam spot branch
    TString                       fPuEnergyDensityName;  // Fastjet correction info name
    TString                       fPfCandidateName;      // particle flow candidates collection name
    TString                       fTracksName;

    const MuonCol                *fMuons;           // muon collection handle
    const ElectronCol            *fElectrons;       // electron collection handle
    const VertexCol              *fPrimVerts;       // primary vertex collection handle
    const BeamSpotCol            *fBeamSpot;        // pointer to beam spot branch
    const PileupEnergyDensityCol *fPuEnergyDensity; // Fastjet correction info handle
    const PFCandidateCol         *fPfCandidates;    // particle flow candidates collection handle
    const TrackCol               *fTracks;          // tracks collection handle
    const Vertex                 *fBestPv;          // best primary vertex in the event

    vector<bool> pfNoPileUpflag;               //flag which says whether each pf candidate is in the pf no pileup collection


    ClassDef(BaseH4lSkim,1)
  };
}
#endif
