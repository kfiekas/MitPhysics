//--------------------------------------------------------------------------------------------------
// $Id: JetIDMod.h,v 1.6 2008/12/10 11:44:33 loizides Exp $
//
// JetIDMod
//
// This module applies jet identification criteria and exports a pointer to a collection
// of "good jet" according to the specified identification scheme.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_JETIDMOD_H
#define MITPHYSICS_MODS_JETIDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"

namespace mithep 
{
  class JetIDMod : public BaseMod
  {
    public:
      JetIDMod(const char *name="JetIDMod", 
               const char *title="Jed identification module");
      ~JetIDMod() {}

      Double_t          GetEtCut()                     const { return fJetEtCut;          }     
      const char       *GetInputName()                 const { return fJetBranchName;     }   
      const char       *GetGoodName()                  const { return GetGoodJetsName();  }     
      const char       *GetGoodJetsName()              const { return fGoodJetsName;      }     
      const char       *GetOutputName()                const { return GetGoodJetsName();  }     
      Bool_t            GetUseCorrection()             const { return fUseJetCorrection;  }     
      void              SetEtCut(Double_t cut)               { fJetEtCut = cut;           }     
      void              SetGoodJetsName(const char *name)    { fGoodJetsName = name;      }     
      void              SetGoodName(const char *name)        { SetGoodJetsName(name);     }     
      void              SetInputName(const char *name)       { fJetBranchName = name;     }  
      void              SetOutputName(const char *name)      { SetGoodJetsName(name);     }     
      void              SetUseCorrection(Bool_t b)           { fUseJetCorrection = b;     }     

    protected:
      void              Process();
      void              SlaveBegin();

      TString           fJetBranchName;         //name of jet collection (input)
      TString           fGoodJetsName;          //name of good jets collection (output)
      Bool_t            fUseJetCorrection;      //=true then use corrected energy
      Double_t          fJetEtCut;              //jet et cut
      const JetCol     *fJets;                  //!jet collection

      ClassDef(JetIDMod, 1) // Jet identification module
  };
}
#endif
