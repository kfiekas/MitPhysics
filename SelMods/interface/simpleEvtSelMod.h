//------------------------------------------------------------------------------
// $Id: simpleEvtSelMod.h,v 1.2 2008/10/23 12:21:48 ceballos Exp $
//
// simpleEvtSelMod
//
// A Module for Selecting ttbar events
// and produces some distributions
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_SIMPLEEVTSELMOD_H
#define HWWMODS_SIMPLEEVTSELMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class simpleEvtSelMod : public BaseMod
  {
    public:
    simpleEvtSelMod(const char *name="simpleEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~simpleEvtSelMod() {}
      void      SetPrintDebug(bool b)	       { fPrintDebug = b;      }   
      void      SetCleanJetsName (TString s)   { fCleanJetsName = s;   }   

    protected:
      bool      fPrintDebug;	   // Debug info
      TString   fPlotType;	   // Type of histograms to make
      TString   fMetName;	   // name of met collection
      TString   fMuonName;	   // name of muon collection
      TString   fTrackName;	   // name of track collection
      TString   fVertexName;	   // name of vertex collection
      TString   fCleanJetsName;    // name of clean central jets collection
      MetCol    *fMet;  	   // Missing Et
      MuonCol	*fMuons;	   // Muon branch
      TrackCol  *fTracks;	   // Track branch     
      VertexCol *fVertices;	   // Vertices branches

      TH1D      *hDSelVert[4];
      TH2D      *hDSelD0Phi;

      int       fNEventsProcessed;

      void      Begin();
      void      Process();
      void      SlaveBegin();
      void      SlaveTerminate();
      void      Terminate();	  

      double    DecayXY(const mithep::Track *lTrack, mithep::Vertex *iVertex);
      double    DecayXY(const mithep::Track *lTrack, mithep::VertexCol *iVertices);

      ClassDef(simpleEvtSelMod,1) // TAM example analysis module
  };
}
#endif
