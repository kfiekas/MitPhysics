//--------------------------------------------------------------------------------------------------
// $Id $
//
// ElectronIDMVA
//
// Helper Class for Electron Identification MVA
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_ElectronIDMVA_H
#define MITPHYSICS_UTILS_ElectronIDMVA_H

#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/TrackFwd.h"
#include "MitAna/DataTree/interface/Electron.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodSwitches.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodMeasurements.h"


class TRandom3;
namespace TMVA {
  class Reader;
}

namespace mithep {
  class ElectronIDMVA {
    public:
      ElectronIDMVA();
      ~ElectronIDMVA(); 

      void     Initialize(TString methodName,
                          TString Subdet0Pt10To20Weights , 
                          TString Subdet1Pt10To20Weights , 
                          TString Subdet2Pt10To20Weights,
                          TString Subdet0Pt20ToInfWeights, 
                          TString Subdet1Pt20ToInfWeights, 
                          TString Subdet2Pt20ToInfWeights,
                          ElectronLikelihood *LH );
      
      Bool_t   IsInitialized() const { return fIsInitialized; }
      Double_t MVAValue(const Electron *ele, const Vertex *vertex);
      
    protected:      
      TMVA::Reader            *fTMVAReader[6];
      TString                  fMethodname;
      
      ElectronLikelihood       *fLH;                     //Likelihood
      Bool_t                    fIsInitialized;
      
      Float_t                   fMVAVar_EleSigmaIEtaIEta; 
      Float_t                   fMVAVar_EleDEtaIn; 
      Float_t                   fMVAVar_EleDPhiIn; 
      Float_t                   fMVAVar_EleHoverE; 
      Float_t                   fMVAVar_EleD0; 
      Float_t                   fMVAVar_EleDZ; 
      Float_t                   fMVAVar_EleFBrem; 
      Float_t                   fMVAVar_EleEOverP; 
      Float_t                   fMVAVar_EleESeedClusterOverPout; 
      Float_t                   fMVAVar_EleSigmaIPhiIPhi; 
      Float_t                   fMVAVar_EleNBrem; 
      Float_t                   fMVAVar_EleOneOverEMinusOneOverP; 
      Float_t                   fMVAVar_EleESeedClusterOverPIn; 
      Float_t                   fMVAVar_EleIP3d; 
      Float_t                   fMVAVar_EleIP3dSig; 
      Float_t                   fMVAVar_EleStandardLikelihood; 
      
      
    ClassDef(ElectronIDMVA, 0) // Muon tools
      };
}

#endif
