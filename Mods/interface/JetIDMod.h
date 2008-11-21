//--------------------------------------------------------------------------------------------------
// $Id: JetIDMod.h,v 1.1 2008/10/15 06:04:59 loizides Exp $
//
// JetIDMod
//
// This Module applies electron ID criteria and exports a pointer to a collection
// of Good Electrons according to some specified ID scheme
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITANA_TREEMOD_JETIDMOD_H
#define MITANA_TREEMOD_JETIDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class JetIDMod : public BaseMod
  {
    public:
      JetIDMod(const char *name="JetIDMod", 
                     const char *title="Example analysis module with all branches");
      ~JetIDMod() {}
      void                     SetPrintDebug(bool b)             { fPrintDebug = b;          }
      void                     SetJetIDType(TString type)        { fJetIDType     = type;    }
      void                     SetJetName(TString s)             { fJetName = s;             }  
      void                     SetGoodJetsName(TString s)        { fGoodJetsName = s;        }     
      void                     SetUseJetCorrection(bool b)       { fUseJetCorrection = b;    }     
      void                     SetJetEtCut(Double_t x)           { fJetEtCut = x;            }     

    protected:
      bool                     fPrintDebug;               //flag for printing debug output
      TString                  fJetName;                  //name of muon collection
      TString                  fGoodJetsName ;            //name of good jets collection  
      TString                  fJetIDType;                //Type of ID we impose
      JetCol		      *fJets;                     //!Jet branch
      int                      fNEventsProcessed;         // Number of events processed
      bool                     fUseJetCorrection;         //bool to use JetCorrected pt
      Double_t                 fJetEtCut;                 //value of Jet Et cut

      void                     Begin();
      void                     Process();
      void                     SlaveBegin();
      void                     SlaveTerminate();
      void                     Terminate();

      ClassDef(JetIDMod,1) // TAM example analysis module
  };
}
#endif
