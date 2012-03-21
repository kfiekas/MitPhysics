#include "MitPhysics/Utils/interface/JetIDMVA.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include <TFile.h>
#include <TRandom3.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"


ClassImp(mithep::JetIDMVA)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
JetIDMVA::JetIDMVA() :
  fMethodName("JetIDMA"),
  fIsInitialized(kFALSE),
  fJetPtMin(10), //We need to lower this
  fNPV    (0),
  fJPt1   (0),
  fJEta1  (0),
  fJPhi1  (0),
  fJD01   (0),
  fJDZ1   (0),
  fJM1    (0),
  fNPart1 (0),
  fLPt1   (0),
  fLEta1  (0),
  fLPhi1  (0),
  fSPt1   (0),
  fSEta1  (0),
  fSPhi1  (0),
  fNEPt1  (0),
  fNEEta1 (0),
  fNEPhi1 (0),
  fEMPt1  (0),
  fEMEta1 (0),
  fEMPhi1 (0),
  fChPt1  (0),
  fChPhi1 (0),
  fLFr1   (0),
  fDRLC1  (0),
  fDRLS1  (0),
  fDRM1   (0),
  fDRNE1  (0),
  fDREM1  (0),
  fDRCH1  (0)
{    
  fReader = 0;
}
//--------------------------------------------------------------------------------------------------
JetIDMVA::~JetIDMVA()
{
  fReader = 0;
}

//--------------------------------------------------------------------------------------------------
void JetIDMVA::Initialize( TString iMethodName,
			   TString iWeights, 
			    JetIDMVA::MVAType iType) 
{ 
  
  fIsInitialized = kTRUE;
  fMethodName    = iMethodName;
  fType          = iType;
  fReader        = 0;
  fReader        = new TMVA::Reader( "!Color:!Silent:Error" );  
  if (fType == kBaseline) {
    fReader->AddSpectator( "npv"     , &fNPV    );
    fReader->AddVariable( "jspt_1"   , &fJPt1   );
    fReader->AddVariable( "jseta_1"  , &fJEta1  );
    fReader->AddVariable( "jsphi_1"  , &fJPhi1  );
    fReader->AddVariable( "jd0_1"    , &fJD01   );
    fReader->AddVariable( "jdZ_1"    , &fJDZ1   );
    fReader->AddVariable( "jm_1"     , &fJM1    );
    fReader->AddVariable( "npart_1"  , &fNPart1 );
    fReader->AddVariable( "lpt_1"    , &fLPt1   );
    fReader->AddVariable( "leta_1"   , &fLEta1  );
    fReader->AddVariable( "lphi_1"   , &fLPhi1  );
    fReader->AddVariable( "spt_1"    , &fSPt1   );
    fReader->AddVariable( "seta_1"   , &fSEta1  );
    fReader->AddVariable( "sphi_1"   , &fSPhi1  );
    fReader->AddVariable( "lnept_1"  , &fNEPt1  );
    fReader->AddVariable( "lneeta_1" , &fNEEta1 );
    fReader->AddVariable( "lnephi_1" , &fNEPhi1 );
    fReader->AddVariable( "lempt_1"  , &fEMPt1  );
    fReader->AddVariable( "lemeta_1" , &fEMEta1 );
    fReader->AddVariable( "lemphi_1" , &fEMPhi1 );
    fReader->AddVariable( "lchpt_1"  , &fChPt1  );
    fReader->AddVariable( "lchphi_1" , &fChPhi1 );
    fReader->AddVariable( "lLfr_1"   , &fLFr1   );
    fReader->AddVariable( "drlc_1"   , &fDRLC1  );
    fReader->AddVariable( "drls_1"   , &fDRLS1  );
    fReader->AddVariable( "drm_1"    , &fDRM1   );
    fReader->AddVariable( "drmne_1"  , &fDRNE1  );
    fReader->AddVariable( "drem_1"   , &fDREM1  );
    fReader->AddVariable( "drch_1"   , &fDRCH1  );
  }
  fReader->BookMVA(fMethodName , iWeights );
  std::cout << "Jet ID MVA Initialization\n";
  std::cout << "MethodName : " << fMethodName << " , type == " << fType << std::endl;
}

//--------------------------------------------------------------------------------------------------
Double_t JetIDMVA::MVAValue(    
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
			    Float_t iDRLC1  ,
			    Float_t iDRLS1  ,
			    Float_t iDRM1   ,
			    Float_t iDRNE1 ,
			    Float_t iDREM1  ,
			    Float_t iDRCH1  
			    ){
  
  if (!fIsInitialized) { 
    std::cout << "Error: JetIDMVA not properly initialized.\n"; 
    return -9999;
  }
  
  fNPV    = iNPV;
  fJPt1   = iJPt1;
  fJEta1  = iJEta1;
  fJPhi1  = fJPhi1;
  fJD01   = iJD01;
  fJDZ1   = iJDZ1;
  fJM1    = iJM1 ;
  fNPart1 = iNPart1;
  fLPt1   = iLPt1;
  fLEta1  = iLEta1;
  fLPhi1  = iLPhi1;
  fSPt1   = iSPt1;
  fSEta1  = iSEta1;
  fSPhi1  = iSPhi1;
  fNEPt1  = iNEPt1;
  fNEEta1 = iNEEta1;
  fNEPhi1 = iNEPhi1;
  fEMPt1  = iEMPt1;
  fEMEta1 = iEMEta1;
  fEMPhi1 = iEMPhi1;
  fChPt1  = iChPt1;
  fChPhi1 = iChPhi1;
  fLFr1   = iLFr1;
  fDRLC1  = iDRLC1;
  fDRLS1  = iDRLS1;
  fDRM1   = iDRM1;
  fDRNE1  = iDRNE1;
  fDREM1  = iDREM1;
  fDRCH1  = iDRCH1;  

  Double_t lMVA = -9999;  
  lMVA = fReader->EvaluateMVA( fMethodName );

  return lMVA;
}
//--------------------------------------------------------------------------------------------------
Bool_t JetIDMVA::pass(const PFJet *iJet,const Vertex *iVertex,
		      FactorizedJetCorrector *iJetCorrector,
		      const PileupEnergyDensityCol *iPileupEnergyDensity) { 
  if(!JetTools::passPFLooseId(iJet))                 return false;
  if(correctedPt(iJet,iJetCorrector,iPileupEnergyDensity) < fJetPtMin && iJet->TrackCountingHighEffBJetTagsDisc() == -100) return false; //This line is a bug in the Met training
  if( fabs(JetTools::impactParameter(iJet,iVertex,true)) < 0.2) return true;
  double lMVA = MVAValue(iJet,iVertex,iJetCorrector,iPileupEnergyDensity);
  if(lMVA < -0.8)                            return false;
  if(lMVA < -0.5 && fabs(iJet->Eta()) > 3.0) return false;
  return true;
}
//--------------------------------------------------------------------------------------------------
Double_t JetIDMVA::MVAValue(const PFJet *iJet,const Vertex *iVertex, //Vertex here is the PV
			    FactorizedJetCorrector *iJetCorrector,
			    const PileupEnergyDensityCol *iPileupEnergyDensity,
			    Bool_t printDebug) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: JetIDMVA not properly initialized.\n"; 
    return -9999;
  }
  if(!JetTools::passPFLooseId(iJet)) return -2.;

  //set all input variables
  fJPt1      = correctedPt(iJet,iJetCorrector,iPileupEnergyDensity);
  fJEta1     = iJet->RawMom().Eta();
  fJPhi1     = iJet->RawMom().Phi();
  fJM1       = iJet->Mass();

  const mithep::PFCandidate *lLead     = JetTools::leadCand(iJet,-1); 
  const mithep::PFCandidate *lSecond   = JetTools::leadCand(iJet,-1,true); 
  const mithep::PFCandidate *lLeadNeut = JetTools::leadCand(iJet ,5); 
  const mithep::PFCandidate *lLeadEm   = JetTools::leadCand(iJet ,4); 
  const mithep::PFCandidate *lLeadCh   = JetTools::leadCand(iJet ,1); 

  fJD01         = JetTools::impactParameter(iJet,iVertex);  
  fJDZ1         = JetTools::impactParameter(iJet,iVertex,true);
  fNPart1       = iJet->NPFCands();
  fLPt1         = lLead    ->Pt(); 
  fLEta1        = lLead    ->Eta(); 
  fLPhi1        = lLead    ->Phi(); 
  fSPt1         = lSecond  ->Pt(); 
  fSEta1        = lSecond  ->Eta(); 
  fSPhi1        = lSecond  ->Phi(); 
  fNEPt1        = lLeadNeut->Pt(); 
  fNEEta1       = lLeadNeut->Eta(); 
  fNEPhi1       = lLeadNeut->Phi(); 
  fEMPt1        = lLeadEm  ->Pt(); 
  fEMEta1       = lLeadEm  ->Eta(); 
  fEMPhi1       = lLeadEm  ->Phi(); 
  fChPt1        = lLeadCh  ->Pt(); 
  //fChEta1       = lLeadCh  ->Eta(); 
  fChPhi1       = lLeadCh  ->Phi(); 
  fLFr1         = lLead->Pt()/iJet->Pt();

  fDRLC1        = MathUtils::DeltaR(iJet->Mom(),lLead  ->Mom());
  fDRLS1        = MathUtils::DeltaR(iJet->Mom(),lSecond->Mom());
  fDRM1         = JetTools::dRMean (iJet,-1);
  fDRNE1        = JetTools::dRMean (iJet, 5);
  fDREM1        = JetTools::dRMean (iJet, 4);
  fDRCH1        = JetTools::dRMean (iJet, 1);

  double lMVA = fReader->EvaluateMVA( fMethodName );
   
  if (printDebug == kTRUE) {
    std::cout << "Debug Jet MVA: "
	      << fNPV    << " "
	      << fJPt1   << " "
	      << fJEta1  << " "
	      << fJPhi1  << " "
	      << fJD01   << " "
	      << fJDZ1   << " "
	      << fJM1    << " "
	      << fNPart1 << " "
	      << fLPt1   << " "
	      << fLEta1  << " "
	      << fLPhi1  << " "
	      << fSPt1   << " "
	      << fSEta1  << " "
	      << fSPhi1  << " "
	      << fNEPt1  << " "
	      << fNEEta1 << " "
	      << fNEPhi1 << " "
	      << fEMPt1  << " "
	      << fEMEta1 << " "
	      << fEMPhi1 << " "
	      << fChPt1  << " "
	      << fChPhi1 << " "
	      << fLFr1   << " "
	      << fDRLC1  << " "
	      << fDRLS1  << " "
	      << fDRM1   << " "
	      << fDRNE1 << " "
	      << fDREM1  << " "
	      << fDRCH1  << " "
              << " === : === "
              << lMVA << " "    
              << std::endl;
  }

  return lMVA;
}
Double_t JetIDMVA::correctedPt(const PFJet *iJet, FactorizedJetCorrector *iJetCorrector,const PileupEnergyDensityCol *iPUEnergyDensity) { 
  const FourVectorM rawMom = iJet->RawMom();
  iJetCorrector->setJetEta(rawMom.Eta());
  iJetCorrector->setJetPt (rawMom.Pt());
  iJetCorrector->setJetPhi(rawMom.Phi());
  iJetCorrector->setJetE  (rawMom.E());
  iJetCorrector->setRho   (iPUEnergyDensity->At(0)->RhoHighEta());
  iJetCorrector->setJetA  (iJet->JetArea());
  iJetCorrector->setJetEMF(-99.0);     
  Double_t correction = iJetCorrector->getCorrection();
  return rawMom.Pt()*correction;
}
