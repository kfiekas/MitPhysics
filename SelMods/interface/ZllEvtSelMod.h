//------------------------------------------------------------------------------
// $Id: ZllEvtSelMod.h,v 1.1 2008/10/23 12:21:48 ceballos Exp $
//
// ZllEvtSelMod
//
// A Module for Selecting ttbar events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_ZLLEVTSELMOD_H
#define HWWMODS_ZLLEVTSELMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"
#include "TRandom.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class ZllEvtSelMod : public BaseMod
  {
    public:
    ZllEvtSelMod(const char *name="ZllEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~ZllEvtSelMod() {}
      void               SetPrintDebug(bool b)       { fPrintDebug = b;   }   
      void               SetCleanJetsName (TString s){ fCleanJetsName = s;}   

    protected:
      bool         fPrintDebug;     // Debug info
      TString      fPlotType;	    // Type of histograms to make
      TString      fMetName;	    // name of met collection
      TString      fTrackName;      // name of track collection
      TString      fCaloTowerName;  // name of calotower collection
      TString      fCleanJetsName;  // name of clean central jets collection
      TString      fMCLeptonsName;  // new lepton coll
      TString      fMCQuarksName;   // new Quarks coll
      MetCol       *fMet;	    // Missing Et
      TrackCol     *fTracks;        // Track branch     
      CaloTowerCol *fCaloTowers;    // CaloTowers branch     
      MuonCol      *fMuons;     	    
      ElectronCol  *fElectrons; 		
      JetCol       *fJets;      	    

      TH1D         *hDZllMET[300];
      TH1D         *hDZllIso[300];
      TH2D         *hDZllIsoEtaIsoTrack0;  
      TH2D         *hDZllIsoEtaIsoNTrack0; 
      TH2D         *hDZllIsoEtaIsoCaloAll0;
      TH2D         *hDZllIsoEtaIsoEem0;
      TH2D         *hDZllIsoEtaIsoHcal0;  
      TH2D         *hDZllIsoEtaIsoHo0;
      TH2D         *hDZllIsoEtaIsoTrack;  
      TH2D         *hDZllIsoEtaIsoNTrack; 
      TH2D         *hDZllIsoEtaIsoCaloAll;
      TH2D         *hDZllIsoEtaIsoEem;    
      TH2D         *hDZllIsoEtaIsoHcal;   
      TH2D         *hDZllIsoEtaIsoHo;
    
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

      ClassDef(ZllEvtSelMod,1) // TAM example analysis module
  };
}
#endif
