//--------------------------------------------------------------------------------------------------
// $Id: MuonIDMod.h,v 1.17 2009/04/30 06:34:02 ceballos Exp $
//
// MuonIDMod
//
// This module applies muon identification criteria and exports a pointer to a collection
// of "good muons" according to the specified ID scheme.
//
// See http://indico.cern.ch/contributionDisplay.py?contribId=1&confId=45945
// See http://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=42229
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
      Bool_t             GetReverseIsoCut()             const { return fReverseIsoCut;      }
      void               SetCaloIsoCut(Double_t cut)          { fCaloIsolationCut  = cut;   }
      void               SetClassType(const char *type)       { fMuonClassType     = type;  }
      void               SetCleanMuonsName(const char *name)  { fCleanMuonsName    = name;  }   
      void               SetCleanName(const char *name)       { SetCleanMuonsName(name);    }   
      void               SetCombIsoCut(Double_t cut)          { fCombIsolationCut  = cut;   }
      void               SetD0Cut(Double_t cut)               { fD0Cut = cut;               }
      void               SetIDType(const char *type)          { fMuonIDType        = type;  }
      void               SetInputName(const char *name)       { fMuonBranchName    = name;  }   
      void               SetIsoType(const char *type)         { fMuonIsoType       = type;  }
      void               SetOutputName(const char *name)      { SetCleanMuonsName(name);    }   
      void               SetPtMin(Double_t pt)                { fMuonPtMin         = pt;    }
      void               SetReverseIsoCut(Bool_t b)           { fReverseIsoCut = b;         }
      void               SetTrackIsoCut(Double_t cut)         { fTrackIsolationCut = cut;   }

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
      Double_t           fD0Cut;               //max d0
      EMuIdType          fMuIDType;            //!muon id type (imposed)
      EMuIsoType         fMuIsoType;           //!muon iso type (imposed)
      EMuClassType       fMuClassType;         //!muon class type (imposed)
      const MuonCol     *fMuons;               //!muon collection
      const VertexCol   *fVertices;            //!vertices branch
      MuonTools         *fMuonTools;           //!muon tool
      Bool_t             fReverseIsoCut;       //!apply reversion iso cut
    
    ClassDef(MuonIDMod, 1) // Muon identification module
  };
}
#endif
