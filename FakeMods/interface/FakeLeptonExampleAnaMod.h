//--------------------------------------------------------------------------------------------------
// $Id: FakeLeptonExampleAnaMod.h,v 1.3 2009/07/20 19:05:04 loizides Exp $
//
// FakeLeptonExampleAnaMod
//
// This is an example analysis module which makes use of the FakeEventHeader objects in order to 
// loop over all possible fake lepton combinations. This module is meant to illustrate how
// the FakeEventHeader objects are to be used for analysis involving fake leptons.
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
      Bool_t                fPerformFakeMuonMetCorrection; //whether to do fake muon met correction
      TString               fSampleName;                   //name of Sample
      TString               fFakeEventHeaderName;          //name of fake event header       (input)
      TString               fElectronFakeableObjectsName;  //name of electron denominators   (input)
      TString               fMuonFakeableObjectsName;      //name of muon denominators       (input)
      TString               fMCPartBranchName;             //name of particle collection     (input)
      TString               fGenJetBranchName;             //name of genjet collection       (input)
      TString               fTrackBranchName;              //name of track collection        (input)
      TString               fMuonBranchName;	           //name of muon collection         (input)
      TString               fMetName;                      //name of met collection          (input)
      TString               fCleanJetsName;                //name of clean central jets collection
      TString               fTriggerObjectsName;           //name of trigger objects
      const MCParticleCol  *fParticles;                    //!GenParticle branch
      const GenJetCol      *fGenJets;                      //!GenJet branch
      const TrackCol	   *fTracks;                       //!Track branch      
      const MuonCol	   *fMuons;	                   //!Muon branch
      const MetCol         *fMet;                          //!Missing Et
  

      TH1D                    *fDileptonCharge;                   //Dilepton Charge Histogram
      TH1D                    *fLeptonPtMax;                      //Lepton1 Pt Histogram
      TH1D                    *fLeptonPtMin;                      //Lepton2 Pt Histogram
      TH1D                    *fMetPtHist;                        //Met Histogram
      TH1D                    *fDeltaPhiLeptons;                  //DeltaPhi Histogram
      TH1D                    *fDeltaEtaLeptons;                  //DeltaEta Histogram
      TH1D                    *fDileptonMass;                     //Dilepton Mass Histogram

      TH1D                    *fLeptonEta_NMinusOne;              //Lepton Eta Histogram
      TH1D                    *fLeptonPtMax_NMinusOne;            //Lepton1 Pt Histogram
      TH1D                    *fLeptonPtMin_NMinusOne;            //Lepton2 Pt Histogram
      TH1D                    *fMetPtHist_NMinusOne;	          //Met Histogram
      TH1D                    *fMetPhiHist_NMinusOne;             //Met Phi Histogram
      TH1D                    *fMETdeltaPhilEtHist_NMinusOne;     //METdeltaPhilEtHist Histogram   
      TH1D                    *fNCentralJets_NMinusOne;           //# of Central Jets Histogram
      TH1D                    *fNDirtyMuonsHist_NMinusOne;        //# of Dirty Muons Histogram
      TH1D                    *fNCleanExtraTracksHist_NMinusOne;  //# of Extra Tracks Histogram
      TH1D                    *fDeltaPhiLeptons_NMinusOne;        //DeltaPhi Histogram
      TH1D                    *fDeltaEtaLeptons_NMinusOne;        //DeltaEta Histogram
      TH1D                    *fDileptonMass_NMinusOne;           //Dilepton Mass Histogram
      TH1D                    *fMinDeltaPhiLeptonMet_NMinusOne;   //MinDeltaPhiLeptonMet Histogram

      TH1D                    *fMinDeltaPhiLeptonMet_afterCuts;   //MinDeltaPhiLeptonMet Histogram
      TH1D                    *fMtLepton1_afterCuts;              //M_t1 Histogram Histogram
      TH1D                    *fMtLepton2_afterCuts;              //M_t2 Histogram Histogram
      TH1D                    *fMtHiggs_afterCuts;                //M_t Higgs Histogram

      ClassDef(FakeLeptonExampleAnaMod,1) // Fake lepton analysis example
  };
}
#endif
