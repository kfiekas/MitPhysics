//--------------------------------------------------------------------------------------------------
// $Id: MuonIDMod.h,v 1.7 2008/11/28 09:13:50 loizides Exp $
//
// MuonIDMod
//
// This module applies muon identification criteria and exports a pointer to a collection
// of "good muons" according to the specified ID scheme.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_MUONIDMOD_H
#define MITPHYSICS_MODS_MUONIDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"
#include "MitPhysics/Utils/interface/MuonTools.h"

namespace mithep 
{
  class MuonIDMod : public BaseMod
  {
    public:
      MuonIDMod(const char *name="MuonIDMod", 
                const char *title="Muon identification module");
      ~MuonIDMod() {}

      const char   *GetMuonBranchName(const char *name) const { return fMuonBranchName;    }   
      const char   *GetCleanMuonsName(const char *name) const { return fCleanMuonsName;    }   
      const char   *GetMuonIDType(const char *type)     const { return fMuonIDType;	   }
      const char   *GetMuonIsoType(const char *type)    const { return fMuonIsoType;	   }
      const char   *GetMuonClassType(const char *type)  const { return fMuonClassType;     }
      Double_t      GetTrackIsolationCut(Double_t cut)  const { return fTrackIsolationCut; }
      Double_t      GetCaloIsolationCut(Double_t cut)   const { return fCaloIsolationCut;  }
      Double_t      GetCombIsolationCut(Double_t cut)   const { return fCombIsolationCut;  }
      Double_t      GetMuonPtMin(Double_t pt)           const { return fMuonPtMin;         }
      void          SetMuonBranchName(const char *name)  { fMuonBranchName    = name; }   
      void          SetCleanMuonsName(const char *name)  { fCleanMuonsName    = name; }   
      void          SetMuonIDType(const char *type)      { fMuonIDType	      = type; }
      void          SetMuonIsoType(const char *type)     { fMuonIsoType	      = type; }
      void          SetMuonClassType(const char *type)   { fMuonClassType     = type; }
      void          SetTrackIsolationCut(Double_t cut)   { fTrackIsolationCut = cut;  }
      void          SetCaloIsolationCut(Double_t cut)    { fCaloIsolationCut  = cut;  }
      void          SetCombIsolationCut(Double_t cut)    { fCombIsolationCut  = cut;  }
      void          SetMuonPtMin(Double_t pt)            { fMuonPtMin	      = pt;   }

      enum EMuIdType {
        kIdUndef = 0,       //not defined
        kTight,             //"Tight"
        kLoose,             //"Loose"
        kCustomId           //"Custom"
      };
      enum EMuIsoType {
        kIsoUndef = 0,      //not defined
        kTrackCalo,         //"TrackCalo"
        kTrackCaloCombined, //"TrackCaloCombined"
        kTrackCaloSliding,  //"TrackCaloSliding"
        kCustomIso,         //"Custom"
        kNoIso              //"NoIso"
      };
      enum EMuClassType {
        kClassUndef = 0,    //not defined
        kAll,               //"All"
        kGlobal,            //"Global"
        kSta,               //"Standalone"
        kTrackerOnly        //"TrackerOnly"
      };

    protected:
      TString       fMuonBranchName;            //name of muon collection (in branch)
      TString       fCleanMuonsName;            //name of exported "good muon" collection
      TString       fMuonIDType;                //type of muon id scheme we impose
      TString       fMuonIsoType;               //type of muon isolations scheme we impose
      TString       fMuonClassType;             //type of muon class we impose
      Double_t      fTrackIsolationCut;         //cut value for track isolation
      Double_t      fCaloIsolationCut;          //cut value for calo isolation
      Double_t      fCombIsolationCut;          //cut value for combined isolation
      Double_t      fMuonPtMin;                 //min muon pt
      MuonCol      *fMuons;                     //!muon branch
      MuonTools    *fMuonTools;                 //!muon tool
      EMuIdType     fMuIDType;                  //!muon id type (imposed)
      EMuIsoType    fMuIsoType;                 //!muon iso type (imposed)
      EMuClassType  fMuClassType;               //!muon class type (imposed)

      void          Process();
      void          SlaveBegin();
    
      ClassDef(MuonIDMod,1) // Muon identification module
  };
}
#endif
