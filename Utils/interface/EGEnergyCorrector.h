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
namespace TMVA {
  class Reader;
}

namespace mithep {
  class EGEnergyCorrector {
    public:
      EGEnergyCorrector();
      ~EGEnergyCorrector(); 

      void Initialize(Bool_t ismc, TString phfixstring, TString phfixfile, TString ebweights, TString ebvarweights, TString eeweights, TString eevarweights);
      Bool_t IsInitialized() const { return fIsInitialized; }
      
      void CorrectEnergyWithError(Photon *p);
      std::pair<double,double> CorrectedEnergyWithError(const Photon *p);
      
    protected:
      PhotonFix fPhFix;
      TMVA::Reader *fReadereb;
      TMVA::Reader *fReaderebvariance;
      TMVA::Reader *fReaderee;
      TMVA::Reader *fReadereevariance;      

      TString fMethodname;
      
      Bool_t fIsInitialized;
      Bool_t fIsMC;
      
      Float_t r9;
      Float_t e5x5norm;
      Float_t scrawe;
      Float_t scpsenorm;
      Float_t sceta;
      Float_t scphi;
      Float_t scetawidth;
      Float_t scphiwidth;
      Float_t hovere;
      Float_t sigietaieta;
      Float_t etac;
      Float_t etas;
      Float_t etam;
      Float_t phic;
      Float_t phis;
      Float_t phim;  
      Float_t xz;
      Float_t xc;
      Float_t xs;
      Float_t xm;
      Float_t yz;
      Float_t yc;
      Float_t ys;
      Float_t ym; 
      
      Float_t fSpectator;
      
      
    ClassDef(EGEnergyCorrector, 0) // Muon tools
      };
}

#endif
