//--------------------------------------------------------------------------------------------------
// $Id: 
//
// MetTools
//
// Authors: M. Zanetti
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_METTOOLS_H
#define MITPHYSICS_UTILS_METTOOLS_H

#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

namespace mithep {

  class MetTools {

  public:

    MetTools(const MuonCol *fMuons, const PFCandidateCol *fPFCandidates, 
             const Vertex *fVertex, float deltaZCut = 0.1, float ptCut=4.0, float etaCut = 3.0);
    MetTools(const ElectronCol *fElectrons, const PFCandidateCol *fPFCandidates, 
             const Vertex *fVertex, float deltaZCut = 0.1, float ptCut=4.0, float etaCut = 3.0);

    MetTools(const MuonCol *fMuons, const PFCandidateCol *fPFCandidates, const PFJetCol *fPFJets,
             const Vertex *fVertex, float deltaZCut = 0.1, float ptCut=4.0, float etaCut = 3.0);
    MetTools(const ElectronCol *fElectrons, const PFCandidateCol *fPFCandidates, const PFJetCol *fPFJets, 
             const Vertex *fVertex, float deltaZCut = 0.1, float ptCut=4.0, float etaCut = 3.0);

    MetTools(const MuonCol *fMuons, const ElectronCol *fElectrons, const PFCandidateCol *fPFCandidates, 
             const Vertex *fVertex, float deltaZCut = 0.1, float ptCut=4.0, float etaCut = 3.0);

    MetTools(const MuonCol *fMuons, const ElectronCol *fElectrons, const PFCandidateCol *fPFCandidates, const PFJetCol *fPFJets, 
             const Vertex *fVertex, float deltaZCut = 0.1, float ptCut=4.0, float etaCut = 3.0);


    ~MetTools() {}
    
    Met GetCorrectedMet() { return fCorrectedMet; }
    Met GetMinimumMet(const Met *UncorrectedMet);
    Met GetCorrectedTrackMet() { return fCorrectedTrackMet; }
    Met GetMinimumTrackMet(const Met *UncorrectedMet);

    template<class V>
    double GetProjectedMet(const V *fV, const Met *UncorrectedMet);
    template<class V>
    double GetProjectedMet(const V *fV);
    template<class V>
    double GetProjectedTrackMet(const V *fV);

  private:
    Met fCorrectedMet;
    Met fCorrectedTrackMet;
    
    ClassDef(MetTools, 0) // Met tools
      };

  template<class V>
  double MetTools::GetProjectedMet(const V *fV, const Met *UncorrectedMet) {
    double projectedMet = UncorrectedMet->Pt();
    double minDPhi = 999;
    int index = -1;
    for (UInt_t m = 0; m < fV->GetEntries(); ++m) {
      if (MathUtils::DeltaPhi(UncorrectedMet->Phi(), fV->At(m)->Phi()) < minDPhi) {
	minDPhi = MathUtils::DeltaPhi(UncorrectedMet->Phi(), fV->At(m)->Phi());
	index = m;
      }
    }
    if (minDPhi < TMath::Pi()/2.) return projectedMet = projectedMet * sin(minDPhi);
    return projectedMet;
  }

  template<class V>
  double MetTools::GetProjectedMet(const V *fV) {
    double projectedMet = fCorrectedMet.Pt();
    double minDPhi = 999;
    int index = -1;
    for (UInt_t m = 0; m < fV->GetEntries(); ++m) {
      if (MathUtils::DeltaPhi(fCorrectedMet.Phi(), fV->At(m)->Phi()) < minDPhi) {
	minDPhi = MathUtils::DeltaPhi(fCorrectedMet.Phi(), fV->At(m)->Phi());
	index = m;
      }
    }
    if (minDPhi < TMath::Pi()/2.) return projectedMet = projectedMet * sin(minDPhi);
    return projectedMet;
  }

  template<class V>
  double MetTools::GetProjectedTrackMet(const V *fV) {
    double projectedMet = fCorrectedTrackMet.Pt();
    double minDPhi = 999;
    int index = -1;
    for (UInt_t m = 0; m < fV->GetEntries(); ++m) {
      if (MathUtils::DeltaPhi(fCorrectedTrackMet.Phi(), fV->At(m)->Phi()) < minDPhi) {
	minDPhi = MathUtils::DeltaPhi(fCorrectedTrackMet.Phi(), fV->At(m)->Phi());
	index = m;
      }
    }
    if (minDPhi < TMath::Pi()/2.) return projectedMet = projectedMet * sin(minDPhi);
    return projectedMet;
  }

}

#endif
