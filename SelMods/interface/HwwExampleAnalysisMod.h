//--------------------------------------------------------------------------------------------------
// $Id $
//
// HwwExampleAnalysisMod
//
// A Module for Selecting H->WW events
// and produces some distributions.
//
//
// Authors: Si Xie
//--------------------------------------------------------------------------------------------------

#ifndef ANA_SELMODS_HWWEXAMPLEANALYSISMOD_H
#define ANA_SELMODS_HWWEXAMPLEANALYSISMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitPhysics/Utils/interface/MuonTools.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class HwwExampleAnalysisMod : public BaseMod
  {
    public:
    HwwExampleAnalysisMod(const char *name="HwwExampleAnalysisMod", 
		 const char *title="Example analysis module with all branches");
      ~HwwExampleAnalysisMod() {}

      const char *GetMuonBranchName()                const { return fMuonBranchName;               }
      const char *GetMetName()                       const { return fMetName;                      }
      const char *GetCleanJetsName()                 const { return fCleanJetsName;                }
      const char *GetCleanJetsNoPtCutName()          const { return fCleanJetsNoPtCutName;         }

      void  SetMuonBranchName(const char *name)            { fMuonBranchName               = name; }
      void  SetMetName(const char *name)                   { fMetName                      = name; }
      void  SetCleanJetsName(const char *name)             { fCleanJetsName                = name; }
      void  SetCleanJetsNoPtCutName(const char *name)      { fCleanJetsNoPtCutName         = name; }

    protected:
      TString                  fMuonBranchName;	         //name of muon branch
      TString                  fMetName;                 //name of met collection
      TString                  fCleanJetsName;           //name of clean central jets collection
      TString                  fCleanJetsNoPtCutName;    //name of clean all jets collection
      TString                  fVertexName;              //name of vertex collection
      TString                  fPFCandidatesName;        //name of PFCandidate collection
      const MuonCol	      *fMuons;	                 //!Muon branch
      const MetCol            *fMet;                     //!Missing Et branch
      const VertexCol         *fVertices;                //!Vertex branch
      const PFCandidateCol    *fPFCandidates;            //!PFCandidates branch
      Int_t                    fNEventsSelected;         //selected events

      TH1D                    *fHWWSelection;
      TH1D                    *fHWWToEESelection;
      TH1D                    *fHWWToMuMuSelection;
      TH1D                    *fHWWToEMuSelection;
      TH1D                    *fHWWToMuESelection;

      TH1D                    *fLeptonEta;
      TH1D                    *fLeptonPtMax;
      TH1D                    *fLeptonPtMin;
      TH1D                    *fMetPtHist;
      TH1D                    *fMetPhiHist;
      TH1D                    *fUncorrMetPtHist;
      TH1D                    *fUncorrMetPhiHist;
      TH1D                    *fDeltaPhiLeptons;
      TH1D                    *fDeltaEtaLeptons;
      TH1D                    *fDileptonMass;

      TH1D                    *fLeptonEta_NMinusOne;   
      TH1D                    *fLeptonPtMax_NMinusOne;    
      TH1D                    *fLeptonPtMin_NMinusOne;    
      TH1D                    *fMetPtHist_NMinusOne;	      
      TH1D                    *fMetPhiHist_NMinusOne;        
      TH1D                    *fMETdeltaPhilEtHist_NMinusOne;     
      TH1D                    *fNCentralJets_NMinusOne;   
      TH1D                    *fNSoftMuonsHist_NMinusOne;
      TH1D                    *fDeltaPhiLeptons_NMinusOne;
      TH1D                    *fDeltaEtaLeptons_NMinusOne;
      TH1D                    *fDileptonMass_NMinusOne;   
      TH1D                    *fMinDeltaPhiLeptonMet_NMinusOne; 

      TH1D                    *fMinDeltaPhiLeptonMet_afterCuts;                                     
      TH1D                    *fMtLepton1_afterCuts;
      TH1D                    *fMtLepton2_afterCuts;
      TH1D                    *fMtHiggs_afterCuts;
      TH1D                    *fLeptonPtPlusMet_afterCuts;      

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      ClassDef(HwwExampleAnalysisMod,1) // TAM example analysis module
  };
}
#endif
