//--------------------------------------------------------------------------------------------------
// $Id: $
//
// H4lSkim
//
// This module selects Higgs to ZZ to 4 lepton events for skimming purposes.
//
// Authors: D.Ralph                   (but C.Paus stole it, changed the name and wrote the doc :-) )
//--------------------------------------------------------------------------------------------------
#ifndef H4LSKIM_H
#define H4LSKIM_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityFwd.h"
#include "MitAna/DataTree/interface/BaseVertex.h"
#include "MitAna/DataTree/interface/TriggerMask.h"

#include <TLorentzVector.h>

namespace mithep
{  
  class H4lSkim : public BaseMod
  {    
  public:
    H4lSkim(const char *name="H4lSkim", const char *title="four lepton skim");
    ~H4lSkim();	    
    
    ULong64_t MatchHLT(const Double_t pt, const Double_t eta, const Double_t phi);
    ULong64_t MatchHLT(const Double_t eta, const Double_t phi);
    
    void AddTrigger(const char* name, const ULong64_t id,
		    const char* objName1="", const ULong64_t objId1=0, const Double_t minPt1=0,
		    const char* objName2="", const ULong64_t objId2=0, const Double_t minPt2=0) {
      fTriggerNamesv.push_back(name);
      fTriggerIdsv.push_back(id);
      
      fTriggerObjNames1v.push_back(objName1);
      fTriggerObjIds1v.push_back(objId1);
      fTriggerObjMinPt1v.push_back(minPt1);
        
      fTriggerObjNames2v.push_back(objName2);
      fTriggerObjIds2v.push_back(objId2);
      fTriggerObjMinPt2v.push_back(minPt2);
    };
      
  protected:
    void Begin();
    void BeginRun();
    void EndRun();
    void SlaveBegin();
    void SlaveTerminate();
    void Terminate();
    void Process();
    
    Bool_t                        fSkipIfHLTFail;   // flag: skip event proc. if HLT not accept
    UInt_t                        flpt_gt_5;
    UInt_t                        flpt_gt_10;
    UInt_t                        fselected;
    
    
    // Compute PF isolation
    Bool_t                        isMuFO(const Muon *mu);
    Bool_t                        isLooseEleFO(const Electron *ele);
    Float_t                       computePFEleIso(const Electron *electron, const Double_t dRMax);
    
    TString                       fMuonName;             // muon collection name
    TString                       fElectronName;         // electron collection name
    TString                       fPrimVtxName;          // primary vertex collection name
    TString                       fBeamSpotName;         // pointer to beam spot branch
    TString                       fPUEnergyDensityName;  // Fastjet correction info name
    TString                       fPFCandidateName;      // particle flow candidates collection name
    TString                       fTrigMaskName;         // trigger mask name
    
    
    const MuonCol                *fMuons;           // muon collection handle
    const ElectronCol            *fElectrons;       // electron collection handle
    const VertexCol              *fPrimVerts;       // primary vertex collection handle
    const BeamSpotCol            *fBeamSpot;        // pointer to beam spot branch
    const PileupEnergyDensityCol *fPUEnergyDensity; // Fastjet correction info handle
    const PFCandidateCol         *fPFCandidates;    // particle flow candidates collection handle 
    const TriggerMask            *fTrigMask;        // trigger mask handle
    
    BaseVertex                    fVertex;          // best primary vertex in the event
    
    vector<TString>               fTriggerNamesv;       // names of triggers we're interested in 
    vector<ULong64_t>             fTriggerIdsv;         // corresponding ETriggerBit value
    
    vector<TString>               fTriggerObjNames1v;
    vector<ULong64_t>             fTriggerObjIds1v;
    vector<Double_t>              fTriggerObjMinPt1v;
    
    vector<TString>               fTriggerObjNames2v;
    vector<ULong64_t>             fTriggerObjIds2v;
    vector<Double_t>              fTriggerObjMinPt2v;
    
    ClassDef(H4lSkim,1)
  };
}

#endif
