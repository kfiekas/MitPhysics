// $Id: MetTools.cc,v 1.14 2011/01/21 11:25:29 ceballos Exp $

#include "MitPhysics/Utils/interface/MetTools.h"
#include <TFile.h>

ClassImp(mithep::MetTools)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
MetTools::MetTools(const MuonCol *fMuons, const PFCandidateCol *fPFCandidates, const Vertex *fVertex, 
		   float deltaZCut, float ptCut, float etaCut) {

  float trackNumeratorX  =0, trackNumeratorY  =0;
  float neutralNumeratorX=0, neutralNumeratorY=0;

  // muons Pt
  for (UInt_t m = 0; m < fMuons->GetEntries(); ++m) {
    trackNumeratorX -= fMuons->At(m)->Px();
    trackNumeratorY -= fMuons->At(m)->Py();
  }

  // PF candidates pT
  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {

    // charged
    if (fPFCandidates->At(i)->TrackerTrk()){

      bool isMuonTrack = false;
      for (UInt_t m = 0; m < fMuons->GetEntries(); ++m) {
	if (fMuons->At(m)->TrackerTrk() == fPFCandidates->At(i)->TrackerTrk()) {
	  isMuonTrack = true;
	  break;
	}
      }      
      if (isMuonTrack) continue;
      
      if (fabs(fPFCandidates->At(i)->TrackerTrk()->DzCorrected(*fVertex)) < deltaZCut) {
	trackNumeratorX -= fPFCandidates->At(i)->Px();
	trackNumeratorY -= fPFCandidates->At(i)->Py();
      }
    }

    // neutral    
    if (fPFCandidates->At(i)->PFType()== PFCandidate::eNeutralHadron || fPFCandidates->At(i)->PFType()== PFCandidate::eGamma) {
      if (fPFCandidates->At(i)->Pt() > ptCut and fabs(fPFCandidates->At(i)->Eta()) < etaCut ) {
	neutralNumeratorX -= fPFCandidates->At(i)->Px();
	neutralNumeratorY -= fPFCandidates->At(i)->Py();
      }
    }
  }

  mithep::Met fCorrectedMet(trackNumeratorX+neutralNumeratorX, trackNumeratorY+neutralNumeratorY);

}


MetTools::MetTools(const ElectronCol *fElectrons, const PFCandidateCol *fPFCandidates, const Vertex *fVertex, 
		   float deltaZCut, float ptCut, float etaCut) {

  float trackNumeratorX  =0, trackNumeratorY  =0;
  float neutralNumeratorX=0, neutralNumeratorY=0;

  // electrons Pt
  for (UInt_t m = 0; m < fElectrons->GetEntries(); ++m) {
    trackNumeratorX -= fElectrons->At(m)->Px();
    trackNumeratorY -= fElectrons->At(m)->Py();
  }

  // PF candidates pT
  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {

    // charged
    if (fPFCandidates->At(i)->TrackerTrk()){
      bool isElectronTrack = false;
      for (UInt_t m = 0; m < fElectrons->GetEntries(); ++m) {
	if (fElectrons->At(m)->TrackerTrk() == fPFCandidates->At(i)->TrackerTrk() or
	    fElectrons->At(m)->GsfTrk() == fPFCandidates->At(i)->GsfTrk() ) {
	  isElectronTrack = true;
	  break;
	}
      }      
      if (isElectronTrack) continue;

      if (fabs(fPFCandidates->At(i)->TrackerTrk()->DzCorrected(*fVertex)) < deltaZCut) {
	trackNumeratorX -= fPFCandidates->At(i)->Px();
	trackNumeratorY -= fPFCandidates->At(i)->Py();
      }
    }

    // neutral    
    if (fPFCandidates->At(i)->PFType()== PFCandidate::eNeutralHadron || fPFCandidates->At(i)->PFType()== PFCandidate::eGamma) {
      if (fPFCandidates->At(i)->Pt() > ptCut and fabs(fPFCandidates->At(i)->Eta()) < etaCut ) {
	neutralNumeratorX -= fPFCandidates->At(i)->Px();
	neutralNumeratorY -= fPFCandidates->At(i)->Py();
      }
    }
  }
  mithep::Met fCorrectedMet(trackNumeratorX+neutralNumeratorX, trackNumeratorY+neutralNumeratorY);

}



Met MetTools::GetMimumMet(const Met *UncorrectedMet) {

  return (fCorrectedMet.Pt() < UncorrectedMet->Pt()) ?  fCorrectedMet : *UncorrectedMet; 
}


