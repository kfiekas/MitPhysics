//--------------------------------------------------------------------------------------------------
// $Id $
//
// VertexTools
//
// Helper Class for photon Identification decisions.
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_VERTEXTOOLS_H
#define MITPHYSICS_UTILS_VERTEXTOOLS_H

#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/Electron.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/BaseVertex.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/VertexMVA.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"

#include "TMVA/Reader.h"

namespace mithep {
  typedef std::vector<double> VertexZarray; 

  class VertexTools {
  public:
   
    static double NewMass(const Photon* ph1, const Photon* ph2, const BaseVertex* vert);

    static VertexZarray ExtractZarray(const VertexCol* vcol, float zmin=-30, float zmax = 30,
				      const BaseVertex  *fBeamSpot = NULL);
    static VertexZarray ExtractZarray(float zmin=-30, float zmax=30, float step=0.05);
    
    static const Vertex* BestVtx(const PFCandidateCol *fPFJets, const VertexCol *c,
				 const BaseVertex  *fBeamSpot, FourVector diboso);
    static double BestVtx(const PFCandidateCol *fPFJets, VertexZarray zcol, 
			  const BaseVertex  *fBeamSpot, FourVector diboso);
    static double Prob(const PFCandidateCol *fPFJets, double zpos,  
		       const BaseVertex  *fBeamSpot, FourVector diboso);
    
    static double VertexWidth(const Vertex*,  const BaseVertex* );
    
    static VertexTools* instance(const char* str){
      if(meobject == NULL){
	meobject = new VertexTools(str);
      }
      return meobject;
    }

    Float_t tmvar1, tmvar2, tmvar3, tmvar4, tmvar5, tmvar6;
    TMVA::Reader* reader;
    
  private:
    static VertexTools *meobject;

    VertexTools();
    VertexTools(const char* str);
    TString relname;

    ClassDef(VertexTools, 0) // Muon tools
      };
}

#endif
