//--------------------------------------------------------------------------------------------------
// $Id: ElectronIDMod.h,v 1.5 2008/11/28 09:13:50 loizides Exp $
//
// ElectronIDMod
//
// This module applies electron identification criteria and exports a pointer to a collection
// of "good electrons" according to the specified identification scheme.
//
// Authors: S.Xie, C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_ELECTRONIDMOD_H
#define MITPHYSICS_MODS_ELECTRONIDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"

namespace mithep 
{
  class ElectronIDMod : public BaseMod
  {
    public:
      ElectronIDMod(const char *name="ElectronIDMod", 
                    const char *title="Electron identification module");
      ~ElectronIDMod() {}

      void                SetElectronBranchName(const char *n) { fElectronBranchName = n;    }   
      void                SetGoodElectronsName(const char *n)  { fGoodElectronsName  = n;    }   
      void                SetElectronIDType(const char *type)  { fElectronIDType     = type; }
      void                SetElectronIsoType(const char *type) { fElectronIsoType    = type; }
      void                SetElectronPtMin(Double_t pt)        { fElectronPtMin      = pt;   }
      void                SetIDLikelihoodCut(Double_t cut)     { fIDLikelihoodCut    = cut;  }
      void                SetTrackIsolationCut(Double_t cut)   { fTrackIsolationCut  = cut;  }
      void                SetCaloIsolationCut(Double_t cut)    { fCaloIsolationCut   = cut;  }
      void                SetEcalJurIsoCut(Double_t cut)       { fEcalJuraIsoCut     = cut;  }
      void                SetHcalIsolationCut(Double_t cut)    { fHcalIsolationCut   = cut;  }

      enum EElIdType {
        kIdUndef = 0,       //not defined
        kTight,             //"Tight"
        kLoose,             //"Loose"
        kLikelihood,        //"Likelihood"
        kCustomId           //"Custom"
      };
      enum EElIsoType {
        kIsoUndef = 0,      //not defined
        kTrackCalo,         //"TrackCalo"
        kTrackJura,         //"TrackJura"
        kTrackJuraSliding,  //"TrackJuraSliding"
        kNoIso,             //"NoIso"
        kCustomIso          //"Custom"
      };

    protected:
      TString             fElectronBranchName;     //branch name of electron collection
      TString             fGoodElectronsName;      //name of exported "good electrons" collection
      TString             fElectronIDType;         //type of electron ID we impose
      TString             fElectronIsoType;        //type of electron Isolation we impose
      Double_t            fElectronPtMin;          //min pt cut
      Double_t            fIDLikelihoodCut;        //cut value for ID likelihood
      Double_t            fTrackIsolationCut;      //cut value for track isolation
      Double_t            fCaloIsolationCut;       //cut value for calo isolation
      Double_t            fEcalJuraIsoCut;         //cut value for ecal jurassic isolation
      Double_t            fHcalIsolationCut;       //cut value for hcal isolation
      const ElectronCol  *fElectrons;              //!electron branch
      EElIdType           fElIdType;               //!identification scheme
      EElIsoType          fElIsoType;              //!isolation scheme

      void                Process();
      void                SlaveBegin();
    
      ClassDef(ElectronIDMod,1) // Electron identification module
  };
}
#endif
