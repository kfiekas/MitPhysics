//$Id:MakeNtuple.h,v 1.1 2011/07/12 Mingming Yang 

#ifndef MITHGG_MODS_MAKENTUPLE_H
#define MITHGG_MODS_MAKENTUPLE_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/MVATools.h"

#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitPhysics/Utils/interface/VertexTools.h"

class TH1D;
class MCEventInfo;
class IClassifierReader;
class ClusteringModuleCl3d;

class TNtuple;

namespace mithep 
{
  class MakeNtuple : public BaseMod
  {
  public:
    MakeNtuple(const char *name = "MakeNtuple", const char *title="Make Ntuple");

  protected:
    void Process();
    void SlaveBegin();
    
    void SetTrigObjsName(const char *name) { fTrigObjsName = name; }
    void SetPhotonName(TString name)       { fPhotonName=name; }
    void SetPhotonsFromBranch(Bool_t b)    { fPhotonsFromBranch = b; }
    void SetIsData(Bool_t b)               { fIsData = b; }
    void SetIsSignal(Bool_t b)             { fIsSignal = b; }
    void SetOverlapCut(double c)           { fOverlapCut = c; }
    int MatchRecPhotonsToGenPhotonsReal(const Photon *photonRec);
    void FillPhotonTree(const Photon *p,const Photon *p_accompany,const Vertex *SelVtx,const PFCandidateCol *fPFCands,Float_t GenPhotonID,Float_t GenPhotonMotherID,Float_t VtxProb, Float_t EventNum,Float_t _decayZ);
    void FindHiggsPtAndZ(Float_t& pt, Float_t& decayZ, Float_t& mass);

    Bool_t                      fPhotonsFromBranch;
    Bool_t                      fIsData;                   // looking at Data (or MC)
    Bool_t                      fIsSignal;  
    Bool_t                      fApplyElectronVeto; 
    Bool_t                      finvertElectronVeto; 
    TString                     fTrigObjsName;             // name of trigger objects
    TString                     fPhotonName;   
    TString                     fElectronName;       
    TString                     fTrackBranchName;
    TString                     fPVName;
    TString                     fPileUpDenName;
    TString                     fMcEventInfoName;          
    TString                     fMcParticleName;           
    TString                     fPileupInfoName;     
    TString                     fBeamSpotName;
    TString                     fConversionName;
    TString                     fPFCandName;
    const PhotonCol              *fPhotons;                  //! photon from data stream
    const ElectronCol            *fElectrons; 
    const TrackCol               *fTracks;
    const VertexCol              *fPV;
    const PileupEnergyDensityCol *fPileUpDen;
    const MCEventInfo          *fMcEventInfo;              //! MC event information branch
    const PileupInfoCol        *fPileupInfos;              //  
    const MCParticleCol        *fMcParticles;                
    const BeamSpotCol          *fBeamSpot;
    const DecayParticleCol     *fConversions;
    const PFCandidateCol         *fPFCands;
    double                      fPhotonPtMin;                 
    double                      fOverlapCut;               // cut to allow for rejection of overlap
    Float_t                     fprocessid;
    Float_t                     fNVertexesGenPile;
    Float_t                     fMatchedGenPhotonID;
    Float_t                     fMatchedGenPhotonMotherID;
    TNtuple                    *hPhotonNtuple;
    VertexTools                 fVtxTools;
    ClassDef(MakeNtuple, 1) // Make Ntuple module
      };
}
#endif
