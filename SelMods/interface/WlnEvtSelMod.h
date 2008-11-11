//------------------------------------------------------------------------------
// $Id: WlnEvtSelMod.h,v 1.1 2008/10/23 12:21:48 ceballos Exp $
//
// WlnEvtSelMod
//
// A Module for Selecting ttbar events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_WlnEvtSelMod_H
#define HWWMODS_WlnEvtSelMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"
#include "TRandom.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class WlnEvtSelMod : public BaseMod
  {
    public:
    WlnEvtSelMod(const char *name="WlnEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~WlnEvtSelMod() {}
      void         SetPrintDebug(bool b)         { fPrintDebug = b;   }	
      void         SetCleanJetsName (TString s)  { fCleanJetsName = s;}	
      void         SetTrigObjsName(const char *n){ fObjsName = n;     }
      void         SetMCLeptonsName(TString s)   { fMCLeptonsName = s;}     

    protected:
      bool         fPrintDebug;     // Debug info
      TString      fPlotType;	    // Type of histograms to make
      TString      fMetName;	    // name of met collection
      TString      fTrackName;      // name of track collection
      TString      fCaloTowerName;  // name of calotower collection
      TString      fCleanJetsName;  // name of clean central jets collection
      TString      fMuonName;	    // name of muon collection
      TString      fElectronName;   // name of electron collection
      TString      fMCLeptonsName;  // new lepton coll
      MetCol       *fMet;	    // Missing Et
      TrackCol     *fTracks;        // Track branch     
      CaloTowerCol *fCaloTowers;    // CaloTowers branch     
      MuonCol      *fMuons;     	    
      ElectronCol  *fElectrons; 		
      JetCol       *fJets;      	    

      // Trigger info
      TString      fObjsName;   //name of trigger objects

      TH1D         *hDWlnSel[300];
      TH2D         *hDWlnSel2D[10];
    
      int          fNEventsProcessed;

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      TRandom RND;

      void         getTrackIsolation(double theEta, double thePhi,
                                     vector<ChargedParticle*> leptons,
                                     mithep::TrackCol *iTracks,
				     double &sumPt, int &nTracks);
      void         getCaloIsolation(double theEta, double thePhi,
                                    vector<ChargedParticle*> leptons,
                                    mithep::CaloTowerCol *iTowers,
				    double lIso[4]);

      ClassDef(WlnEvtSelMod,1) // TAM example analysis module
  };
}
#endif
