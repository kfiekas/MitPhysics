//--------------------------------------------------------------------------------------------------
// $Id: JetIDMod.h,v 1.12 2009/05/09 18:35:57 ceballos Exp $
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

namespace mithep 
{
  class JetIDMod : public BaseMod
  {
    public:
      JetIDMod(const char *name="JetIDMod", 
               const char *title="Jet identification module");

      const char       *GetInputName()                 const { return fJetsName;     }   
      const char       *GetGoodName()                  const { return GetGoodJetsName();  }     
      const char       *GetGoodJetsName()              const { return fGoodJetsName;      }     
      const char       *GetOutputName()                const { return GetGoodJetsName();  }     
      Double_t          GetPtCut()                     const { return fJetPtCut;          }     
      Bool_t            GetUseCorrection()             const { return fUseJetCorrection;  }     
      Double_t          GetEtaMaxCut()                 const { return fJetEtaMaxCut;      }     
      void              SetGoodJetsName(const char *name)    { fGoodJetsName = name;      }     
      void              SetGoodName(const char *name)        { SetGoodJetsName(name);     }     
      void              SetInputName(const char *name)       { fJetsName = name;          }  
      void              SetOutputName(const char *name)      { SetGoodJetsName(name);     }     
      void              SetPtCut(Double_t cut)               { fJetPtCut = cut;           }     
      void              SetUseCorrection(Bool_t b)           { fUseJetCorrection = b;     }     
      void              SetEtaMaxCut(Double_t cut)           { fJetEtaMaxCut = cut;       }     

    protected:
      void              Process();

      TString           fJetsName;              //name of jet collection (input)
      TString           fGoodJetsName;          //name of good jets collection (output)
      Bool_t            fUseJetCorrection;      //=true then use corrected energy
      Double_t          fJetPtCut;              //jet pt cut
      Double_t          fJetEtaMaxCut;          //jet eta max cut

      ClassDef(JetIDMod, 1) // Jet identification module
  };
}
#endif
