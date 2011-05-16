// $Id: VertexMVA.cc,v 1.3 2011/04/21 15:07:40 bendavid Exp $

#include "MitPhysics/Utils/interface/VertexMVA.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include <TFile.h>

ClassImp(mithep::VertexMVA)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
VertexMVA::VertexMVA()  
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
Double_t VertexMVA::GetProb(const Vertex *v) const  
{
  return 0.0;
}
