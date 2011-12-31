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
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

class TRandom3;
namespace TMVA {
  class Reader;
}

namespace mithep {
  class ElectronIDMVA {
    public:
      ElectronIDMVA();
      ~ElectronIDMVA(); 

      enum MVAType {
        kBaseline = 0,      // SigmaIEtaIEta, DEtaIn, DPhiIn, FBrem, SigmaIPhiIPhi, NBrem, OneOverEMinusOneOverP
        kNoIPInfo,          // kBaseline + EOverP, ESeedClusterOverPout, ESeedClusterOverPIn
        kWithIPInfo,        // kBaseline + d0 , IP3d, IP3dSig
        kIDIsoCombined      // ID variables , PFIso03 , PFIso04
      };


      void     Initialize(TString methodName,
                          TString Subdet0Pt10To20Weights , 
                          TString Subdet1Pt10To20Weights , 
                          TString Subdet2Pt10To20Weights,
                          TString Subdet0Pt20ToInfWeights, 
                          TString Subdet1Pt20ToInfWeights, 
                          TString Subdet2Pt20ToInfWeights,
                          ElectronIDMVA::MVAType type );
      
      Bool_t   IsInitialized() const { return fIsInitialized; }
      Double_t MVAValue(const Electron *ele, const Vertex *vertex);
      Double_t MVAValue(const Electron *ele, const Vertex *vertex, 
                        const PFCandidateCol *PFCands, const PileupEnergyDensityCol *PileupEnergyDensity);
      Double_t MVAValue(Double_t ElePt , Double_t EleEta,
                        Double_t EleSigmaIEtaIEta,
                        Double_t EleDEtaIn,
                        Double_t EleDPhiIn,
                        Double_t EleHoverE,
                        Double_t EleD0,
                        Double_t EleDZ,
                        Double_t EleFBrem,
                        Double_t EleEOverP,
                        Double_t EleESeedClusterOverPout,
                        Double_t EleSigmaIPhiIPhi,
                        Double_t EleNBrem,
                        Double_t EleOneOverEMinusOneOverP,
                        Double_t EleESeedClusterOverPIn,
                        Double_t EleIP3d,
                        Double_t EleIP3dSig);
      Double_t MVAValue(Double_t ElePt , Double_t EleEta,
                        Double_t PileupEnergyDensity,
                        Double_t EleSigmaIEtaIEta,
                        Double_t EleDEtaIn,
                        Double_t EleDPhiIn,
                        Double_t EleHoverE,
                        Double_t EleD0,
                        Double_t EleDZ,
                        Double_t EleFBrem,
                        Double_t EleEOverP,
                        Double_t EleESeedClusterOverPout,
                        Double_t EleSigmaIPhiIPhi,
                        Double_t EleNBrem,
                        Double_t EleOneOverEMinusOneOverP,
                        Double_t EleESeedClusterOverPIn,
                        Double_t EleIP3d,
                        Double_t EleIP3dSig,
                        Double_t EleGsfTrackChi2OverNdof,
                        Double_t EledEtaCalo,
                        Double_t EledPhiCalo,
                        Double_t EleR9,
                        Double_t EleSCEtaWidth,
                        Double_t EleSCPhiWidth,
                        Double_t EleCovIEtaIPhi,
                        Double_t ElePreShowerOverRaw,
                        Double_t EleChargedIso03,
                        Double_t EleNeutralHadronIso03,
                        Double_t EleGammaIso03,
                        Double_t EleChargedIso04,
                        Double_t EleNeutralHadronIso04,
                        Double_t EleGammaIso04 );


    protected:      
      TMVA::Reader            *fTMVAReader[6];
      TString                  fMethodname;
      
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
      Float_t                   fMVAVar_EleGsfTrackChi2OverNdof;
      Float_t                   fMVAVar_EledEtaCalo;
      Float_t                   fMVAVar_EledPhiCalo;
      Float_t                   fMVAVar_EleR9;
      Float_t                   fMVAVar_EleSCEtaWidth;
      Float_t                   fMVAVar_EleSCPhiWidth;
      Float_t                   fMVAVar_EleCovIEtaIPhi;
      Float_t                   fMVAVar_ElePreShowerOverRaw;
      Float_t                   fMVAVar_EleChargedIso03OverPt;
      Float_t                   fMVAVar_EleNeutralHadronIso03OverPt;
      Float_t                   fMVAVar_EleGammaIso03OverPt;
      Float_t                   fMVAVar_EleChargedIso04OverPt;
      Float_t                   fMVAVar_EleNeutralHadronIso04OverPt;
      Float_t                   fMVAVar_EleGammaIso04OverPt;
        
    ClassDef(ElectronIDMVA, 0) // Muon tools
      };
}

#endif
