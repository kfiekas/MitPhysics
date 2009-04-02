//--------------------------------------------------------------------------------------------------
// $Id: MuonIDMod.h,v 1.12 2008/12/11 10:55:43 loizides Exp $
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

      Double_t           GetCaloIsoCut()                const { return fCaloIsolationCut;   }
      const char        *GetClassType()                 const { return fMuonClassType;      }
      const char        *GetCleanName()                 const { return GetCleanMuonsName(); }   
      const char        *GetCleanMuonsName()            const { return fCleanMuonsName;     }   
      Double_t           GetCombIsoCut()                const { return fCombIsolationCut;   }
      const char        *GetIDType()                    const { return fMuonIDType;         }
      const char        *GetInputName()                 const { return fMuonBranchName;     }   
      const char        *GetIsoType()                   const { return fMuonIsoType;        }
      const char        *GetOutputName()                const { return GetCleanMuonsName(); }   
      Double_t           GetPtMin(Double_t pt)          const { return fMuonPtMin;          }
      Double_t           GetTrackIsoCut()               const { return fTrackIsolationCut;  }
      void               SetCaloIsoCut(Double_t cut)          { fCaloIsolationCut  = cut;   }
      void               SetCombIsoCut(Double_t cut)          { fCombIsolationCut  = cut;   }
      void               SetClassType(const char *type)       { fMuonClassType     = type;  }
      void               SetCleanName(const char *name)       { SetCleanMuonsName(name);    }   
      void               SetCleanMuonsName(const char *name)  { fCleanMuonsName    = name;  }   
      void               SetIDType(const char *type)          { fMuonIDType        = type;  }
      void               SetInputName(const char *name)       { fMuonBranchName    = name;  }   
      void               SetIsoType(const char *type)         { fMuonIsoType       = type;  }
      void               SetOutputName(const char *name)      { SetCleanMuonsName(name);    }   
      void               SetPtMin(Double_t pt)                { fMuonPtMin         = pt;    }
      void               SetTrackIsoCut(Double_t cut)         { fTrackIsolationCut = cut;   }
      void               SetD0Cut(Double_t cut)               { fD0Cut = cut;               }

      enum EMuIdType {
        kIdUndef = 0,       //not defined
        kTight,             //"Tight"
        kLoose,             //"Loose"
        kNoId,              //"NoId"
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
      void               Process();
      void               SlaveBegin();

      TString            fMuonBranchName;      //name of muon collection (input)
      TString            fCleanMuonsName;      //name of exported "good muon" collection
      TString            fVertexName;	       //name of vertex collection
      TString            fMuonIDType;          //type of muon id scheme we impose
      TString            fMuonIsoType;         //type of muon isolations scheme we impose
      TString            fMuonClassType;       //type of muon class we impose
      Double_t           fTrackIsolationCut;   //cut value for track isolation
      Double_t           fCaloIsolationCut;    //cut value for calo isolation
      Double_t           fCombIsolationCut;    //cut value for combined isolation
      Double_t           fMuonPtMin;           //min muon pt
      const MuonCol     *fMuons;               //!muon collection
      const VertexCol   *fVertices;            // Vertices branches
      MuonTools         *fMuonTools;           //!muon tool
      EMuIdType          fMuIDType;            //!muon id type (imposed)
      EMuIsoType         fMuIsoType;           //!muon iso type (imposed)
      EMuClassType       fMuClassType;         //!muon class type (imposed)
      Double_t           fD0Cut;               //max d0
    
    ClassDef(MuonIDMod, 1) // Muon identification module
  };
}
#endif
