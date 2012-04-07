//--------------------------------------------------------------------------------------------------
// $Id: JetIDMod.h,v 1.16 2012/04/05 12:25:09 pharris Exp $
//
// MVAMetMod
//
// Example on how to call regressed MET
//
// Authors: P.Harris
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_MVAMETMOD_H
#define MITPHYSICS_MODS_MVAMETOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitPhysics/Utils/interface/MVAMet.h"

namespace mithep 
{
  class MVAMetMod : public BaseMod
  {
    public:
      MVAMetMod(const char *name="MVAMetMod", 
		const char *title="MVAMet example");

    protected:
      void                Process();
      TString             fMVAMetName;
      TString             fJetsName  ;
      TString             fPFCandName;
      TString             fVertexName;
      TString             fPFMetName ;
      const PFJetCol       *fJets;
      const PFCandidateCol *fCands;
      const VertexCol      *fVertices;
      const PFMetCol       *fPFMet;
      const MuonCol        *fMuons;

      MVAMet         *fMVAMet;
      ClassDef(MVAMetMod, 1) // Jet identification module
  };
}
#endif
