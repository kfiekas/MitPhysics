//--------------------------------------------------------------------------------------------------
// $Id: MuonIDMod.h,v 1.2 2008/10/23 12:23:31 ceballos Exp $
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
      void     SetPrintDebug(bool b)              { fPrintDebug         = b;     }
      void     SetMuonIDType(TString type)        { fMuonIDType         = type;  }
      void     SetTrackIsolationCut(Double_t cut) { fTrackIsolationCut  = cut;   }
      void     SetCaloIsolationCut(Double_t cut)  { fCaloIsolationCut   = cut;   }
      void     SetMuonPtMin(double p)             { fMuonPtMin          = p;     }
      void     SetCleanMuonsName(TString s)       { fCleanMuonsName     = s;     }     
    protected:
      bool     fPrintDebug;               //flag for printing debug output
      TString  fMuonName;                 //name of muon collection
      TString  fCleanMuonsName ;           //name of good muons collection  
      TString  fMuonIDType;               //Type of Muon ID we impose
      TString  fMuonIsoType;              //Type of Muon Isolation we impose
      MuonCol  *fMuons;                   //!Muon branch

      double   fTrackIsolationCut;        //Cut value for track isolation
      double   fCaloIsolationCut;         //Cut value for calo isolation
      double   fMuonPtMin;                //min Pt requirement

      int      fNEventsProcessed;         // Number of events processed

      void     Begin();
      void     Process();
      void     SlaveBegin();
      void     SlaveTerminate();
      void     Terminate();
     
    
      ClassDef(MuonIDMod,1) // TAM example analysis module
  };
}
#endif
