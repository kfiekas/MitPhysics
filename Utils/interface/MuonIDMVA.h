//--------------------------------------------------------------------------------------------------
// $Id $
//
// ElectronIDMVA
//
// Helper Class for Muon Identification MVA
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_MuonIDMVA_H
#define MITPHYSICS_UTILS_MuonIDMVA_H

#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/TrackFwd.h"
#include "MitAna/DataTree/interface/Muon.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

class TRandom3;
namespace TMVA {
  class Reader;
}

namespace mithep {
  class MuonIDMVA {
    public:
      MuonIDMVA();
      ~MuonIDMVA(); 

      enum MVAType {
        kV2,
        kV3,
        kV8
      };


      void     Initialize(TString methodName,
                          TString Subdet0Pt10To14p5Weights , 
                          TString Subdet1Pt10To14p5Weights , 
                          TString Subdet0Pt14p5To20Weights,
                          TString Subdet1Pt14p5To20Weights, 
                          TString Subdet0Pt20ToInfWeights, 
                          TString Subdet1Pt20ToInfWeights,
                          MuonIDMVA::MVAType type);
      
      Bool_t   IsInitialized() const { return fIsInitialized; }
//       Double_t MVAValue(const Muon *mu, const Vertex *vertex);
      Double_t MVAValue( Double_t MuPt , Double_t MuEta,
                         Double_t                   MuTkNchi2, 
                         Double_t                   MuGlobalNchi2, 
                         Double_t                   MuNValidHits, 
                         Double_t                   MuNTrackerHits, 
                         Double_t                   MuNPixelHits, 
                         Double_t                   MuNMatches, 
                         Double_t                   MuD0, 
                         Double_t                   MuIP3d, 
                         Double_t                   MuIP3dSig, 
                         Double_t                   MuTrkKink, 
                         Double_t                   MuSegmentCompatibility, 
                         Double_t                   MuCaloCompatibility, 
                         Double_t                   MuHadEnergyOverPt, 
                         Double_t                   MuHoEnergyOverPt, 
                         Double_t                   MuEmEnergyOverPt, 
                         Double_t                   MuHadS9EnergyOverPt, 
                         Double_t                   MuHoS9EnergyOverPt, 
                         Double_t                   MuEmS9EnergyOverPt, 
                         Double_t                   MuChargedIso03OverPt,
                         Double_t                   MuNeutralIso03OverPt,
                         Double_t                   MuChargedIso04OverPt,
                         Double_t                   MuNeutralIso04OverPt
        );


    protected:      
      TMVA::Reader            *fTMVAReader[6];
      TString                  fMethodname;
      
      Bool_t                    fIsInitialized;
      
      Float_t                   fMVAVar_MuTkNchi2; 
      Float_t                   fMVAVar_MuGlobalNchi2; 
      Float_t                   fMVAVar_MuNValidHits; 
      Float_t                   fMVAVar_MuNTrackerHits; 
      Float_t                   fMVAVar_MuNPixelHits; 
      Float_t                   fMVAVar_MuNMatches; 
      Float_t                   fMVAVar_MuD0; 
      Float_t                   fMVAVar_MuIP3d; 
      Float_t                   fMVAVar_MuIP3dSig; 
      Float_t                   fMVAVar_MuTrkKink; 
      Float_t                   fMVAVar_MuSegmentCompatibility; 
      Float_t                   fMVAVar_MuCaloCompatibility; 
      Float_t                   fMVAVar_MuHadEnergyOverPt; 
      Float_t                   fMVAVar_MuHoEnergyOverPt; 
      Float_t                   fMVAVar_MuEmEnergyOverPt; 
      Float_t                   fMVAVar_MuHadS9EnergyOverPt; 
      Float_t                   fMVAVar_MuHoS9EnergyOverPt; 
      Float_t                   fMVAVar_MuEmS9EnergyOverPt; 
      Float_t                   fMVAVar_MuChargedIso03OverPt;
      Float_t                   fMVAVar_MuNeutralIso03OverPt;
      Float_t                   fMVAVar_MuChargedIso04OverPt;
      Float_t                   fMVAVar_MuNeutralIso04OverPt;

      
    ClassDef(MuonIDMVA, 0) // Muon MVA
      };
}

#endif
