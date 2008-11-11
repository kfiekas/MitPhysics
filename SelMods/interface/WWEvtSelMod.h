//------------------------------------------------------------------------------
// $Id: WWEvtSelMod.h,v 1.2 2008/10/23 12:21:48 ceballos Exp $
//
// WWEvtSelMod
//
// A Module for Selecting ttbar events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_WWEVTSELMOD_H
#define HWWMODS_WWEVTSELMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class WWEvtSelMod : public BaseMod
  {
    public:
    WWEvtSelMod(const char *name="WWEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~WWEvtSelMod() {}
      void      SetPrintDebug(bool b)	       { fPrintDebug = b;      }
      void      SetCleanJetsName (TString s)   { fCleanJetsName = s;   }
      void      SetMCLeptonsName(TString s)    { fMCLeptonsName = s;   }
      void      SetMCNeutrinosName(TString s)  { fMCNeutrinosName = s; }

    protected:
      bool      fPrintDebug;	 // Debug info
      TString   fPlotType;	 // Type of histograms to make
      TString   fMetName;	 // name of met collection
      TString   fMuonName;	 // name of muon collection
      TString   fTrackName;	 // name of track collection
      TString   fVertexName;	 // name of vertex collection
      TString   fCleanJetsName;  // name of clean central jets collection
      TString   fMCLeptonsName;  // new lepton coll
      TString   fMCNeutrinosName;// new neutrinos coll
      MetCol    *fMet;  	 // Missing Et
      MuonCol   *fMuons;	 // Muon branch
      TrackCol  *fTracks;	 // Track branch     
      VertexCol *fVertices;	 // Vertices branches

      TH1D      *hDwwSel[300];
      TH1D      *hDwwPresel[50];
      TH2D      *hDwwSelAlphaEP;
      TH1D      *hDwwVert[5];
      TH2D      *hDwwSelD0Phi;

      int       fNEventsProcessed;

      void      Begin();
      void      Process();
      void      SlaveBegin();
      void      SlaveTerminate();
      void      Terminate();	  

      double    DecayXY(const mithep::Track *lTrack, mithep::Vertex *iVertex);
      double    DecayXY(const mithep::Track *lTrack, mithep::VertexCol *iVertices);

      ClassDef(WWEvtSelMod,1) // TAM example analysis module
  };
}
#endif
