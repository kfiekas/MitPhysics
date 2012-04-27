//--------------------------------------------------------------------------------------------------
// $Id: SeparatePileUpMod.h,v 1.6 2009/06/15 15:00:21 loizides Exp $
//
// SeparatePileUpMod
//
// This module applies PFNoPU selection
//
// Authors: G.Ceballos
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_SEPARATEPILEUPMOD_H
#define MITPHYSICS_MODS_SEPARATEPILEUPMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/PFCandidateFwd.h"

namespace mithep 
{
  class SeparatePileUpMod : public BaseMod
  {
    public:
      SeparatePileUpMod(const char *name="SeparatePileUpMod", 
               const char *title="PFNoPU identification module");

      void                SetPFCandidatesName(const char *n)    { fPFCandidatesName  = n;  }  
      void                SetPFPileUpName(const char *n)	{ fPFPileUpName  = n;      }  
      void                SetPFNoPileUpName(const char *n)	{ fPFNoPileUpName  = n;    }  
      void                SetVertexName(const char *n)          { fVertexName = n;         }
      void                SetCheckClosestZVertex(Bool_t b)      { fCheckClosestZVertex = b;}

    protected:
      void                Process();
      void                SlaveBegin();

      TString               fPFCandidatesName;    //name of PF collection (input)
      TString               fPFPileUpName;        //name of exported PFPileUp collection (output)
      TString               fPFNoPileUpName;      //name of exported PFNoPileUp collection (output)
      TString               fVertexName;	  //name of vertex collection
      const PFCandidateCol *fPFCandidates;	  //!pfcandidate branch
      const VertexCol      *fVertices;  	  //!vertices branches
      Bool_t                fCheckClosestZVertex; //boolean to use the closest vertex approach

    ClassDef(SeparatePileUpMod, 1) // Tau identification module
  };
}
#endif
