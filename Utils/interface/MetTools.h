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
#include "MitCommon/MathTools/interface/MathUtils.h"

namespace mithep {

  class MetTools {

  public:

    MetTools(const MuonCol *fMuons, const PFCandidateCol *fPFCandidates, const Vertex *fVertex, float deltaZCut = 0.1, float ptCut=4.0, float etaCut = 3.0);
    MetTools(const ElectronCol *fElectrons, const PFCandidateCol *fPFCandidates, const Vertex *fVertex, float deltaZCut = 0.1, float ptCut=4.0, float etaCut = 3.0);

    ~MetTools() {}
    
    Met GetCorrectedMet() { return fCorrectedMet; }
    Met GetMimumMet(const Met *UncorrectedMet);

    template<class V>
    double GetProjectedMet(const V *fV);
    template<class V>
    double GetProjectedMet(const V *fV, const Met *UncorrectedMet);

  private:
    Met fCorrectedMet;
    
    ClassDef(MetTools, 0) // Met tools
      };

  template<class V>
  double MetTools::GetProjectedMet(const V *fV, const Met *UncorrectedMet) {
    double projectedMet = 0;
    double minDPhi = 999;
    int index = -1;
    for (UInt_t m = 0; m < fV->GetEntries(); ++m) {
      if (MathUtils::DeltaPhi(UncorrectedMet->Phi(), fV->At(m)->Phi()) < minDPhi) {
	minDPhi = MathUtils::DeltaPhi(UncorrectedMet->Phi(), fV->At(m)->Phi());
	index = m;
      }
    }
    if (index==-1) return projectedMet;
    if (minDPhi < TMath::Pi()/2.) projectedMet = projectedMet * sin(minDPhi);
  }


  template<class V>
  double MetTools::GetProjectedMet(const V *fV) {
    double projectedMet = 0;
    double minDPhi = 999;
    int index = -1;
    for (UInt_t m = 0; m < fV->GetEntries(); ++m) {
      if (MathUtils::DeltaPhi(fCorrectedMet.Phi(), fV->At(m)->Phi()) < minDPhi) {
	minDPhi = MathUtils::DeltaPhi(fCorrectedMet.Phi(), fV->At(m)->Phi());
	index = m;
      }
    }
    if (index==-1) return projectedMet;
    if (minDPhi < TMath::Pi()/2.) projectedMet = projectedMet * sin(minDPhi);
  }

}

#endif
