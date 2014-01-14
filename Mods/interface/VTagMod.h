//--------------------------------------------------------------------------------------------------
// $Id: $
//
// VTagMod
//
// This module analyses input jets and tries to identify them as jets corresponding to Vector
// bosons. See documentation here: 
//
// Authors: C.Paus
//--------------------------------------------------------------------------------------------------
#ifndef MITPHYSICS_MODS_VTAGMOD_H
#define MITPHYSICS_MODS_VTAGMOD_H

#include "fastjet/tools/Pruner.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ActiveAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"

namespace mithep
{
  class VTagMod : public BaseMod
  {
    public:
      VTagMod(const char *name  = "VTagMod",
              const char *title = "Vector Boson Tagging module");

      const char                   *GetGoodVTagsName()              const { return fGoodVTagsName;     }
      const char                   *GetOutputName()                 const { return GetGoodVTagsName(); }

      void                          SetGoodVTagsName(const char *n)       { fGoodVTagsName = n;        }
			            
    protected:		            
      void                          Process();
      void                          SlaveBegin();
      void                          Terminate();

    private:
      //float                         GetTau(fastjet::PseudoJet &iJet,int iN, float iKappa);
			            
      TString                       fJetsBranchName;        //(i) name of
      TString                       fPFCandidatesName;      //(i) name of PF candidates colleciont
      TString                       fGoodVTagsName;         //(o) name of VTags collection

      const PFCandidateCol         *fPFCandidates;          //particle flow candidates collection handle

      // Objects from fastjet we want to use
      double                        fConeSize;
      fastjet::Pruner              *fPruner;
      fastjet::JetDefinition       *fCAJetDef;
      fastjet::ActiveAreaSpec      *fActiveArea;
      fastjet::AreaDefinition      *fAreaDefinition;
      //fastjet::ClusterSequenceArea *fClustering;

    ClassDef(VTagMod, 1) // Vector boson tagging module
  };
}
#endif
