//--------------------------------------------------------------------------------------------------
// $Id: MVAMetMod.h,v 1.2 2012/04/07 11:52:02 ceballos Exp $
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
#include "MitPhysics/Utils/interface/MVAMet.h"

namespace mithep 
{
  class MVAMetMod : public BaseMod
  {
    public:
      MVAMetMod(const char *name="MVAMetMod", 
		const char *title="MVAMet example");
      void   SetJetsName(TString s) { fJetsName = s;}   

    protected:
      void                Process();
      void                SlaveBegin();
      void                SlaveTerminate();

      TString             fMVAMetName;
      TString             fJetsName  ;
      TString             fPFCandName;
      TString             fVertexName;
      TString             fPFMetName ;
      const JetCol         *fJets;
      const PFCandidateCol *fCands;
      const VertexCol      *fVertices;
      const PFMetCol       *fPFMet;

      MVAMet         *fMVAMet;
      ClassDef(MVAMetMod, 1) // Jet identification module
  };
}
#endif
