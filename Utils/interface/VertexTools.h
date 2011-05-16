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


namespace mithep {
  class VertexTools {
    public:
      VertexTools();
        
      static const Vertex*      BestVtx(const VertexCol *c, const VertexMVA *mva);
      
    ClassDef(VertexTools, 0) // Muon tools
  };
}

#endif
