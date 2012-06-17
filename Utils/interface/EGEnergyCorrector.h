//--------------------------------------------------------------------------------------------------
// $Id $
//
// EGEnergyCorrector
//
// Helper Class for photon Identification decisions.
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_EGEnergyCorrector_H
#define MITPHYSICS_UTILS_EGEnergyCorrector_H

#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/Electron.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/BaseVertex.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/SuperCluster.h"
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/PhotonFix.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"


class TRandom3;
class GBRForest;
// namespace TMVA {
//   class Reader;
// }


namespace mithep {
  class EGEnergyCorrector {
    public:
      EGEnergyCorrector();
      ~EGEnergyCorrector(); 

      void Initialize(TString phfixstring, TString phfixfile, TString regweights);
      Bool_t IsInitialized() const { return fIsInitialized; }
      
      void CorrectEnergyWithError(Photon *p, const VertexCol *vtxs = 0, Double_t rho = 0., UInt_t version=1, Bool_t applyRescale = kFALSE);
      std::pair<double,double> CorrectedEnergyWithError(const Photon *p);
      std::pair<double,double> CorrectedEnergyWithErrorV2(const Photon *p, const VertexCol *vtxs);
      std::pair<double,double> CorrectedEnergyWithErrorV3(const Photon *p, const VertexCol *vtxs, Double_t rho, Bool_t applyRescale = kFALSE);

      
    protected:
      PhotonFix fPhFix;
      GBRForest *fReadereb;
      GBRForest *fReaderebvariance;
      GBRForest *fReaderee;
      GBRForest *fReadereevariance;      

      TString fMethodname;
      
      Bool_t fIsInitialized;
      Bool_t fIsMC;

      Float_t *fVals;
      
    ClassDef(EGEnergyCorrector, 0) // Muon tools
      };
}

#endif
