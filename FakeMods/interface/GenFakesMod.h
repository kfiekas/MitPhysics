//--------------------------------------------------------------------------------------------------
// $Id: GenFakesMod.h,v 1.2 2009/07/02 12:17:32 phedex Exp $
//
// Genfakesmod
//
// This Module generates a collection of FakeEventHeaders containing information
// about possible fakes and their weight. The collection generated takes into account all possible
// faking combinatorics from the given set of fakable objects.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_FAKEMODS_GENFAKESMOD_H
#define MITPHYSICS_FAKEMODS_GENFAKESMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 

namespace mithep 
{
  class FakeRate;

  class GenFakesMod : public BaseMod
  {
    public:
      GenFakesMod(const char *name="GenFakesMod", 
                  const char *title="Fake Object Generation Module");

      const char   *GetCleanElectronsName()          const { return fCleanElectronsName;          }
      const char   *GetCleanMuonsName()              const { return fCleanMuonsName;              }
      const char   *GetCleanPhotonsName()            const { return fCleanPhotonsName;            }
      const char   *GetCleanJetsName()               const { return fCleanJetsName;               }
      const char   *GetElFakeableObjsName()          const { return fElFakeableObjsName;          }
      const char   *GetMuFakeableObjsName()          const { return fMuFakeableObjsName;          }
      const char   *GetMCLeptonsName()               const { return fMCLeptonsName;               }
      const char   *GetMCTausName()                  const { return fMCTausName;                  }
      const char   *GetFakeEventHeadersName()        const { return fFakeEventHeadersName;        }
      const char   *GetOutputName()                  const { return GetFakeEventHeadersName();    }
      const Bool_t  GetUse2DFakeRate()               const { return fUse2DFakeRate;               }
      const Bool_t  GetUseFitFunction()              const { return fUseFitFunction;              }

      void SetCleanElectronsName(const char *name)          { fCleanElectronsName          = name; }
      void SetCleanMuonssName(const char *name)             { fCleanMuonsName              = name; }
      void SetCleanPhotonsName(const char *name)            { fCleanPhotonsName            = name; }
      void SetCleanJetsName(const char *name)               { fCleanJetsName               = name; }
      void SetElFakeableObjsName(const char *name)          { fElFakeableObjsName          = name; }
      void SetMuFakeableObjsName(const char *name)          { fMuFakeableObjsName          = name; }
      void SetMCLeptonsName(const char *name)               { fMCLeptonsName               = name; }
      void SetMCTausName(const char *name)                  { fMCTausName                  = name; }
      void SetFakeEventHeadersName(const char *name)        { fFakeEventHeadersName        = name; }
      void SetOutputName(const char *name)                  { SetFakeEventHeadersName(name);       }
      void SetElectronFRFilename(const char *name)          { fElectronFRFilename          = name; }
      void SetMuonFRFilename(const char *name)              { fMuonFRFilename              = name; }
      void SetUse2DFakeRate(Bool_t b)                       { fUse2DFakeRate               = b;    }
      void SetUseFitFunction(Bool_t b)                      { fUseFitFunction              = b;    }
      void SetElectronFRFunctionName(const char *name)      { fElectronFRFunctionName      = name; }
      void SetMuonFRFunctionName(const char *name)          { fMuonFRFunctionName          = name; }
      void SetElectronFRHistName(const char *name)          { fElectronFRHistName          = name; }
      void SetMuonFRHistName(const char *name)              { fMuonFRHistName              = name; }

      void LoadFakeRate();
      
    protected:
      void               Process();

      FakeRate          *fFakeRate;                     //holds the fake probabilities
      TString            fElectronFRFilename;           //file containing electron fake rate
      TString            fMuonFRFilename;               //file containing muon fake rate
      Bool_t             fUse2DFakeRate;                //whether to use fit function or not
      Bool_t             fUseFitFunction;               //whether to use fit function or not
      TString            fElectronFRFunctionName;       //fit function containing electron fake rate
      TString            fMuonFRFunctionName;           //fit function containing muon fake rate
      TString            fElectronFRHistName;           //hist containing electron fake rate
      TString            fMuonFRHistName;               //hist containing muon fake rate
      TString            fCleanElectronsName;           //name of clean electrons           (input)
      TString            fCleanMuonsName;               //name of clean muons               (input)
      TString            fCleanPhotonsName;             //name of clean photons             (input)
      TString            fCleanJetsName;                //name of clean jets                (input)
      TString            fMCLeptonsName;                //name of MC leptons
      TString            fMCTausName;                   //name of MC taus
      TString            fElFakeableObjsName;           //name of electron fakeable objects (input)
      TString            fMuFakeableObjsName;           //name of muon fakeable objects     (input)
      TString            fFakeEventHeadersName;         //name of collection of FakeEventHeaders    (output)
   
    ClassDef(GenFakesMod, 1) // Jet cleaning module
  };
}
#endif
