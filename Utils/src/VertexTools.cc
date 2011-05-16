// $Id: VertexTools.cc,v 1.3 2011/04/21 15:07:40 bendavid Exp $

#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include <TFile.h>

ClassImp(mithep::VertexTools)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
VertexTools::VertexTools()  
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
const Vertex *VertexTools::BestVtx(const VertexCol *c, const VertexMVA *mva) {

  if (!c || !c->GetEntries()) return 0;
  
  return c->At(0);
  
}