//--------------------------------------------------------------------------------------------------
// $Id: JetIDMod.h,v 1.3 2008/11/27 16:30:26 loizides Exp $
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

      void              SetJetBranchName(const char *name) { fJetBranchName = name; }  
      void              SetGoodJetsName(const char *name)  { fGoodJetsName = name;  }     
      void              SetUseJetCorrection(bool b)        { fUseJetCorrection = b; }     
      void              SetJetEtCut(Double_t cut)          { fJetEtCut = cut;       }     

    protected:
      TString           fJetBranchName;         //name of muon collection
      TString           fGoodJetsName ;         //name of good jets collection  
      Bool_t            fUseJetCorrection;      //=true then use corrected energy
      Double_t          fJetEtCut;              //value of Jet Et cut
      JetCol           *fJets;                  //!jet branch

      void              Process();
      void              SlaveBegin();

      ClassDef(JetIDMod,1) // Jet identification module
  };
}
#endif
