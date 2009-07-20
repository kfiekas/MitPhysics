//--------------------------------------------------------------------------------------------------
// $Id: FakeLeptonExampleAnaMod.h,v 1.2 2009/07/13 11:27:13 loizides Exp $
//
// FakeLeptonExampleAnaMod
//
// TODO
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITHIGGS_HWWMODS_FAKELEPTONEXAMPLEANAMOD_H
#define MITHIGGS_HWWMODS_FAKELEPTONEXAMPLEANAMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class FakeLeptonExampleAnaMod : public BaseMod
  {
    public:
      FakeLeptonExampleAnaMod(const char *name="FakeLeptonExampleAnaMod", 
                              const char *title="Example analysis module with all branches");

      Bool_t      GetUseMCFake()                         { return fUseMCFake;                    }
      Bool_t      GetPerformFakeMuonMetCorrection()      { return fPerformFakeMuonMetCorrection; }
      const char *GetSampleName()                  const { return fSampleName;                   }
      const char *GetFakeEventHeaderName()         const { return fFakeEventHeaderName;          }
      const char *GetElectronFakeableObjectsName() const { return fElectronFakeableObjectsName;  }
      const char *GetMuonFakeableObjectsName()     const { return fMuonFakeableObjectsName;      }
      const char *GetMCPartBranchName()            const { return fMCPartBranchName;             }
      const char *GetGenJetBranchName()            const { return fGenJetBranchName;             }
      const char *GetTrackBranchName()             const { return fTrackBranchName;              }
      const char *GetMuonBranchName()              const { return fMuonBranchName;               }
      const char *GetMetName()                     const { return fMetName;                      }
      const char *GetTriggerObjectsName()          const { return fTriggerObjectsName;           }

      void   SetUseMCFake(Bool_t b)                        { fUseMCFake                    = b;  }
      void   SetPerformFakeMuonMetCorrection(Bool_t b)     { fPerformFakeMuonMetCorrection = b;  }
      void   SetSampleName(const char *s)                  { fSampleName                   = s;  }
      void   SetFakeEventHeaderName (const char *s)        { fFakeEventHeaderName          = s;  }
      void   SetElectronFakeableObjectsName(const char *s) { fElectronFakeableObjectsName  = s;  }
      void   SetMuonFakeableObjectsName(const char *s)     { fMuonFakeableObjectsName      = s;  }
      void   SetCleanJetsName (TString s)                  { fCleanJetsName                = s;  }
      void   SetMCPartBranchName (TString s)               { fMCPartBranchName             = s;  }
      void   SetGenJetBranchName (TString s)               { fGenJetBranchName             = s;  }
      void   SetTrackBranchName (TString s)                { fTrackBranchName              = s;  }
      void   SetMuonBranchName (TString s)                 { fMuonBranchName               = s;  }
      void   SetMetName (TString s)                        { fMetName                      = s;  }
      void   SetTriggerObjectsName(const char *name)       { fTriggerObjectsName         = name; }

    protected:
      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      Bool_t                fUseMCFake;                    //whether to use MC simulation fakes
      Bool_t                fPerformFakeMuonMetCorrection; //whether to perform fake muon 
                                                           //met correction
      TString               fSampleName;                   //name of Sample
      TString               fFakeEventHeaderName;          //name of fake event header       (input)
      TString               fElectronFakeableObjectsName;  //name of electron denominators   (input)
      TString               fMuonFakeableObjectsName;      //name of muon denominators       (input)
      TString               fMCPartBranchName;             //name of particle collection
      TString               fGenJetBranchName;             //name of genjet collection
      TString               fTrackBranchName;              //name of track collection
      TString               fMuonBranchName;	           //name of muon collection
      TString               fMetName;                      //name of met collection
      TString               fCleanJetsName;                //name of clean central jets collection
      TString               fTriggerObjectsName;           //name of trigger objects
      const MCParticleCol  *fParticles;                    //!GenParticle branch
      const GenJetCol      *fGenJets;                      //!GenJet branch
      const TrackCol	   *fTracks;                       //!Track branch      
      const MuonCol	   *fMuons;	                   //!Muon branch
      const MetCol         *fMet;                          //!Missing Et
  

      TH1D                    *fDileptonCharge;
      TH1D                    *fLeptonPtMax;
      TH1D                    *fLeptonPtMin;
      TH1D                    *fMetPtHist;
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
      TH1D                    *fNDirtyMuonsHist_NMinusOne;
      TH1D                    *fNCleanExtraTracksHist_NMinusOne;
      TH1D                    *fDeltaPhiLeptons_NMinusOne;
      TH1D                    *fDeltaEtaLeptons_NMinusOne;
      TH1D                    *fDileptonMass_NMinusOne;   
      TH1D                    *fMinDeltaPhiLeptonMet_NMinusOne; 

      TH1D                    *fMinDeltaPhiLeptonMet_afterCuts;                                     
      TH1D                    *fMtLepton1_afterCuts;
      TH1D                    *fMtLepton2_afterCuts;
      TH1D                    *fMtHiggs_afterCuts;

      ClassDef(FakeLeptonExampleAnaMod,1) // Fake lepton analysis example
  };
}
#endif
