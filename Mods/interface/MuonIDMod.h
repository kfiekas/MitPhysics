//--------------------------------------------------------------------------------------------------
// $Id: MuonIDMod.h,v 1.4 2008/11/11 21:22:58 ceballos Exp $
//
// MuonIDMod
//
// This Module applies Muon ID criteria and exports a pointer to a collection
// of Good Muons according to some specified ID scheme
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITANA_TREEMOD_MUONIDMOD_H
#define MITANA_TREEMOD_MUONIDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"
#include "MitPhysics/Utils/interface/MuonTools.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class MuonIDMod : public BaseMod
  {
    public:
      MuonIDMod(const char *name="MuonIDMod", 
                     const char *title="Example analysis module with all branches");
      ~MuonIDMod() {}
      void     SetPrintDebug(bool b)                { fPrintDebug	         = b;	}
      void     SetTrackIsolationCut(Double_t cut)   { fTrackIsolationCut         = cut; }
      void     SetCaloIsolationCut(Double_t cut)    { fCaloIsolationCut          = cut; }
      void     SetCombIsolationCut(Double_t cut)    { fCombIsolationCut          = cut; }
      void     SetMuonPtMin(double p)               { fMuonPtMin	         = p;	}
      void     SetTMOneStationLooseCut(bool b)      { fTMOneStationLooseCut      = b;	}
      void     SetTMOneStationTightCut(bool b)      { fTMOneStationTightCut      = b;	}
      void     SetTM2DCompatibilityLooseCut(bool b) { fTM2DCompatibilityLooseCut = b;   }
      void     SetTM2DCompatibilityTightCut(bool b) { fTM2DCompatibilityTightCut = b;   }
      void     SetMuonSlidingIso(bool b)            { fMuonSlidingIso	         = b;   }
      void     SetCleanMuonsName(TString s)         { fCleanMuonsName            = s;   }   
    protected:
      bool      fPrintDebug;                // flag for printing debug output
      TString   fMuonName;                  // name of muon collection
      TString   fCleanMuonsName ;           // name of good muons collection  
      MuonCol   *fMuons;                    // !Muon branch

      double    fTrackIsolationCut;         // Cut value for track isolation
      double    fCaloIsolationCut;          // Cut value for calo isolation
      double    fCombIsolationCut;          // Cut value for combined isolation
      bool      fTMOneStationLooseCut;	   // apply TMOneStationLooseCut?
      bool      fTMOneStationTightCut;	   // apply TMOneStationTightCut?
      bool      fTM2DCompatibilityLooseCut; // apply TM2DCompatibilityLooseCut?
      bool      fTM2DCompatibilityTightCut; // apply TM2DCompatibilityTightCut?
      bool      fMuonSlidingIso;            // apply sliding iso?
      double    fMuonPtMin;                 // min Pt requirement

      int       fNEventsProcessed;          // Number of events processed

      MuonTools myMuonTools;

      void      Begin();
      void      Process();
      void      SlaveBegin();
      void      SlaveTerminate();
      void      Terminate();
     
    
      ClassDef(MuonIDMod,1) // TAM example analysis module
  };
}
#endif
