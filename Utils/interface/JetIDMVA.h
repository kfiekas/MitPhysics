//--------------------------------------------------------------------------------------------------
// $Id $
//
// JetIDMVA
//
// Helper Class for Jet Id MVA
//
// Authors: P. Harris
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_JetIDMVA_H
#define MITPHYSICS_UTILS_JetIDMVA_H

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "MitAna/DataTree/interface/PFJetFwd.h"
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/TrackFwd.h"
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

class TRandom3;
namespace TMVA {
  class Reader;
}

namespace mithep {
  class JetIDMVA {
    public:
      JetIDMVA();
      ~JetIDMVA(); 

      enum MVAType {
        kBaseline = 0
      };


      void     Initialize(TString           iMethodName="JetIDMVA",
                          TString           iWeights="$CMSSW_BASE/src/MitPhysics/data/mva_JetID.weights.xml",
                          JetIDMVA::MVAType iType=kBaseline );
      
      Bool_t   IsInitialized() const { return fIsInitialized; }
      Double_t MVAValue(    
			Float_t iNPV    ,
			Float_t iJPt1   ,
			Float_t iJEta1  ,
			Float_t iJPhi1  ,
			Float_t iJD01   ,
			Float_t iJDZ1   ,
			Float_t iJM1    ,
			Float_t iNPart1 ,
			Float_t iLPt1   ,
			Float_t iLEta1  ,
			Float_t iLPhi1  ,
			Float_t iSPt1   ,
			Float_t iSEta1  ,
			Float_t iSPhi1  ,
			Float_t iNEPt1  ,
			Float_t iNEEta1 ,
			Float_t iNEPhi1 ,
			Float_t iEMPt1  ,
			Float_t iEMEta1 ,
			Float_t iEMPhi1 ,
			Float_t iChPt1  ,
			Float_t iChPhi1 ,
			Float_t iLFr1   ,
			Float_t iDRlC1  ,
			Float_t iDRLS1  ,
			Float_t iDRM1   ,
			Float_t iDRMNE1 ,
			Float_t iDREM1  ,
			Float_t iDRCH1  
			);

      //UNcorrected Jets
      Bool_t   pass(const PFJet *iJet,const Vertex *iVertex,
		    FactorizedJetCorrector *iJetCorrector,
		    const PileupEnergyDensityCol *iPileupEnergyDensity);
      
      //Corrected Jets
      Bool_t   pass(const PFJet *iJet,const Vertex *iVertex);
		    			
      //Uncorrected Jets
      Double_t MVAValue(const PFJet *iJet,const Vertex *iVertex,
			FactorizedJetCorrector *iJetCorrector,
			const PileupEnergyDensityCol *iPileupEnergyDensity,
			Bool_t printDebug=false);

      //Corrected Jets
      Double_t MVAValue(const PFJet *iJet,const Vertex *iVertex,
			Bool_t printDebug=false);


      double  correctedPt(const PFJet *iJet, FactorizedJetCorrector *iJetCorrector,
			  const PileupEnergyDensityCol *iPUEnergyDensity);

      Float_t                  fJetPtMin;

    protected:      
      TMVA::Reader            *fReader;
      TString                  fMethodName;
      MVAType                  fType;
      Bool_t                   fIsInitialized;
      
      Float_t fNPV;
      Float_t fJPt1;
      Float_t fJEta1;
      Float_t fJPhi1;
      Float_t fJD01 ;
      Float_t fJDZ1 ;
      Float_t fJM1  ;
      Float_t fNPart1;
      Float_t fLPt1 ;
      Float_t fLEta1;
      Float_t fLPhi1;
      Float_t fSPt1 ;
      Float_t fSEta1;
      Float_t fSPhi1;
      Float_t fNEPt1;
      Float_t fNEEta1;
      Float_t fNEPhi1;
      Float_t fEMPt1;
      Float_t fEMEta1;
      Float_t fEMPhi1;
      Float_t fChPt1;
      Float_t fChPhi1;
      Float_t fLFr1 ;
      Float_t fDRLC1;
      Float_t fDRLS1;
      Float_t fDRM1 ;
      Float_t fDRNE1;
      Float_t fDREM1;
      Float_t fDRCH1;
        
      ClassDef(JetIDMVA,0)
	};
}


#endif
