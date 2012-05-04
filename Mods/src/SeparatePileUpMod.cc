// $Id: SeparatePileUpMod.cc,v 1.2 2012/04/27 22:41:41 ceballos Exp $

#include "MitPhysics/Mods/interface/SeparatePileUpMod.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::SeparatePileUpMod)

//------------------------------------------------------------------------------
SeparatePileUpMod::SeparatePileUpMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fPFPileUpName("PFPileUp"),
  fPFNoPileUpName("PFNoPileUp"),
  fAllVertexName(Names::gkPVBrn),
  fVertexName(ModNames::gkGoodVertexesName),
  fPFCandidates(0),
  fAllVertices(0),
  fVertices(0),
  fCheckClosestZVertex(kTRUE),
  fUseAllVertices(kTRUE)
{
  // Constructor.
}

//------------------------------------------------------------------------------
void SeparatePileUpMod::Process()
{
  // Process entries of the tree. 

  LoadBranch(fPFCandidatesName);

  PFCandidateOArr *pfPileUp = new PFCandidateOArr;
  pfPileUp->SetName(fPFPileUpName);
  
  PFCandidateOArr *pfNoPileUp = new PFCandidateOArr;
  pfNoPileUp->SetName(fPFNoPileUpName);

  LoadBranch(fAllVertexName);
  
  if(fUseAllVertices == kTRUE) fVertices = fAllVertices;
  else                         fVertices = GetObjThisEvt<VertexOArr>(fVertexName);

  for(UInt_t i = 0; i < fPFCandidates->GetEntries(); i++) {
    const PFCandidate *pf = fPFCandidates->At(i);
    assert(pf);

    if(pf->PFType() == PFCandidate::eHadron) {
      if(pf->HasTrackerTrk() && 
         fVertices->At(0)->HasTrack(pf->TrackerTrk()) &&
         fVertices->At(0)->TrackWeight(pf->TrackerTrk()) > 0)
      {
        pfNoPileUp->Add(pf);
      }
      else {
        Bool_t vertexFound = kFALSE;
        const Vertex *closestVtx = 0;
        Double_t dzmin = 10000;

	for(UInt_t j = 0; j < fAllVertices->GetEntries(); j++) {
	  const Vertex *vtx = fAllVertices->At(j);
	  assert(vtx);

	  if(pf->HasTrackerTrk() && 
	     vtx->HasTrack(pf->TrackerTrk()) &&
	     vtx->TrackWeight(pf->TrackerTrk()) > 0) {
	    vertexFound = kTRUE;
	    closestVtx = vtx;
	    break;
	  }
	  Double_t dz = fabs(pf->SourceVertex().Z() - vtx->Z());
	  if(dz < dzmin) {
	    closestVtx = vtx;
	    dzmin = dz;
	  }
	}

	if(fCheckClosestZVertex) {
	  // Fallback: if track is not associated with any vertex,
	  // associate it with the vertex closest in z
	  if(vertexFound || closestVtx != fVertices->At(0))
	    pfPileUp->Add(pf);
	  else
	    pfNoPileUp->Add(pf);
	}
	else {
	  if(vertexFound && closestVtx != fVertices->At(0))
	    pfPileUp->Add(pf);
	  else
	    pfNoPileUp->Add(pf); // Ridiculous but that's how it is
	}
      }
    }
    else {
      pfNoPileUp->Add(pf);
    }
  } // Loop over PF candidates

  // add to event for other modules to use
  AddObjThisEvt(pfPileUp);  
  AddObjThisEvt(pfNoPileUp);  
}

//------------------------------------------------------------------------------
void SeparatePileUpMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the Tau collection branch.

  ReqBranch(fAllVertexName,    fAllVertices);
  ReqBranch(fPFCandidatesName, fPFCandidates);
}
