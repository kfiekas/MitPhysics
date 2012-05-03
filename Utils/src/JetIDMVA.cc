#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
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
  fJetPtMin(0.)  , //We need to lower this
  fDZCut   (0.2),
  fLowPtMethodName ("JetIDMVALowPt" ),
  fHighPtMethodName("JetIDMVAHighPt"),
  fIsInitialized(kFALSE),
  fNVtx     (0),
  fJPt1     (0),
  fJEta1    (0),
  fJPhi1    (0),
  fJD01     (0),
  fJDZ1     (0),
  fBeta     (0),
  fBetaStar (0),
  fNCharged (0),
  fNNeutrals(0),
  fDRMean   (0),
  fFrac01   (0),
  fFrac02   (0),
  fFrac03   (0),
  fFrac04   (0),
  fFrac05   (0)
{    
  fReader = 0;
}
//--------------------------------------------------------------------------------------------------
JetIDMVA::~JetIDMVA() {

  fReader = 0;
}

//--------------------------------------------------------------------------------------------------
void JetIDMVA::Initialize( JetIDMVA::CutType iCutType,
			   TString iLowPtWeights,
			   TString iHighPtWeights, 
			   JetIDMVA::MVAType iType,
			   TString iCutFileName) 
{ 
  
  fIsInitialized = kTRUE;
  fType          = iType;
  fCutType       = iCutType;
  fReader        = 0;
  fReader        = new TMVA::Reader( "!Color:!Silent:Error" );  
  if (fType == kBaseline) {
    fReader->AddVariable( "nvtx"     , &fNVtx      ); 
    fReader->AddVariable( "jetPt"    , &fJPt1      );  
    fReader->AddVariable( "jetEta"   , &fJEta1     );
    fReader->AddVariable( "jetPhi"   , &fJPhi1     );             
    fReader->AddVariable( "dZ"       , &fJDZ1      );
    fReader->AddVariable( "d0"       , &fJD01      );
    fReader->AddVariable( "beta"     , &fBeta      );
    fReader->AddVariable( "betaStar" , &fBetaStar  );
    fReader->AddVariable( "nCharged" , &fNCharged  );
    fReader->AddVariable( "nNeutrals", &fNNeutrals );
    fReader->AddVariable( "dRMean"   , &fDRMean    );
    fReader->AddVariable( "frac01"   , &fFrac01    );
    fReader->AddVariable( "frac02"   , &fFrac02    );
    fReader->AddVariable( "frac03"   , &fFrac03    );
    fReader->AddVariable( "frac04"   , &fFrac04    );
    fReader->AddVariable( "frac05"   , &fFrac05    );
  }
  fReader->BookMVA(fLowPtMethodName  , iLowPtWeights );
  fReader->BookMVA(fHighPtMethodName , iHighPtWeights );
  std::cout << "Jet ID MVA Initialization\n";
  std::cout << "MethodName : " << fLowPtMethodName << " , type == " << fType << std::endl;

  //Load Cut Matrix
  edm::ParameterSet lConfig = edm::readPSetsFrom(iCutFileName.Data())->getParameter<edm::ParameterSet>("JetIdParams");
  std::string lCutType = "Tight";
  if(fCutType == kMedium) lCutType = "Medium";
  if(fCutType == kLoose ) lCutType = "Loose";
  if(fCutType == kMET   ) lCutType = "MET";
  std::vector<double> lPt010  = lConfig.getParameter<std::vector<double> >(("Pt010_" +lCutType).c_str());
  std::vector<double> lPt1020 = lConfig.getParameter<std::vector<double> >(("Pt1020_"+lCutType).c_str());
  std::vector<double> lPt2030 = lConfig.getParameter<std::vector<double> >(("Pt2030_"+lCutType).c_str());
  std::vector<double> lPt3050 = lConfig.getParameter<std::vector<double> >(("Pt3050_"+lCutType).c_str());
  for(int i0 = 0; i0 < 4; i0++) fMVACut[0][i0] = lPt010 [i0];
  for(int i0 = 0; i0 < 4; i0++) fMVACut[1][i0] = lPt1020[i0];
  for(int i0 = 0; i0 < 4; i0++) fMVACut[2][i0] = lPt2030[i0];
  for(int i0 = 0; i0 < 4; i0++) fMVACut[3][i0] = lPt3050[i0];
  //std::cout << " Working Points : << " << std::endl;
  //for(int i0 = 0; i0 < 4; i0++) for(int i1 = 0; i1 < 4; i1++) 
  //  std::cout << " ==> " << i0 << " -- " << i1 << " -- " << fMVACut[i0][i1] << std::endl;
}

//--------------------------------------------------------------------------------------------------
Double_t JetIDMVA::MVAValue(    
			    Float_t iNPV    ,
			    Float_t iJPt1   ,
			    Float_t iJEta1  ,
			    Float_t iJPhi1  ,
			    Float_t iJD01   ,
			    Float_t iJDZ1   ,
			    Float_t iBeta   ,
			    Float_t iBetaStar,
			    Float_t iNCharged,
			    Float_t iNNeutrals,
			    Float_t iDRMean  ,
			    Float_t iFrac01  ,
			    Float_t iFrac02  ,
			    Float_t iFrac03  ,
			    Float_t iFrac04  ,
			    Float_t iFrac05  
			    ){
  
  if(!fIsInitialized) { 
    std::cout << "Error: JetIDMVA not properly initialized.\n"; 
    return -9999;
  }
  
  fNVtx      = iNPV;
  fJPt1      = iJPt1;
  fJEta1     = iJEta1;
  fJPhi1     = fJPhi1;
  fJD01      = iJD01;
  fJDZ1      = iJDZ1;
  fBeta      = iBeta;
  fBetaStar  = iBetaStar;
  fNCharged  = iNCharged;
  fNNeutrals = iNNeutrals;
  fDRMean    = iDRMean;
  fFrac01    = iFrac01;
  fFrac02    = iFrac02;
  fFrac03    = iFrac03;
  fFrac04    = iFrac04;
  fFrac05    = iFrac05;

  Double_t lMVA = -9999;  
  if(iJPt1 < 10) lMVA = fReader->EvaluateMVA( fLowPtMethodName  );
  if(iJPt1 > 10) lMVA = fReader->EvaluateMVA( fHighPtMethodName );

  return lMVA;
}
//--------------------------------------------------------------------------------------------------
Bool_t JetIDMVA::pass(const PFJet *iJet,const Vertex *iVertex,const VertexCol *iVertices,
		      FactorizedJetCorrector *iJetCorrector,
		      const PileupEnergyDensityCol *iPileupEnergyDensity) { 
  
  if(!JetTools::passPFLooseId(iJet))                 return false;
  if(fabs(iJet->Eta()) > 4.99)                       return false;
  
  double lMVA = MVAValue   (iJet,iVertex,iVertices,iJetCorrector,iPileupEnergyDensity);
  double lPt  = correctedPt(iJet,                  iJetCorrector,iPileupEnergyDensity);
  if(lPt < fJetPtMin)                         return false; 
  if(lPt > 50)                                return true; //==> we can raise this 
  
  int lPtId = 0; 
  if(lPt > 10 && lPt < 20) lPtId = 1;
  if(lPt > 20 && lPt < 30) lPtId = 2;
  if(lPt > 30                   ) lPtId = 3;
  
  int lEtaId = 0;
  if(fabs(iJet->Eta()) > 2.5  && fabs(iJet->Eta()) < 2.75) lEtaId = 1; 
  if(fabs(iJet->Eta()) > 2.75 && fabs(iJet->Eta()) < 3.0 ) lEtaId = 2; 
  if(fabs(iJet->Eta()) > 3.0  && fabs(iJet->Eta()) < 5.0 ) lEtaId = 3; 
  
  double lMVACut = fMVACut[lPtId][lEtaId];
  if(lMVA < lMVACut) return false;
  return true;
   //if( fabs(JetTools::impactParameter(iJet,iVertex,true)) < 0.2) return true;
  //if(correctedPt(iJet,iJetCorrector,iPileupEnergyDensity) < fJetPtMin && iJet->TrackCountingHighEffBJetTagsDisc() == -100) return false; 
  //This line is a bug in the Met training
  //if(lMVA < -0.8)                            return false;
  //if(lMVA < -0.5 && fabs(iJet->Eta()) > 3.0) return false;
}
//--------------------------------------------------------------------------------------------------
Bool_t JetIDMVA::pass(const PFJet *iJet,const Vertex *iVertex,const VertexCol *iVertices) { 
  if(!JetTools::passPFLooseId(iJet))                 return false;
  if(iJet->Pt()        < fJetPtMin) return false; 
  if(iJet->Pt()        > 50)        return true; //==> we can raise this 
  if(fabs(iJet->Eta()) > 4.99)      return false;
  double lMVA = MVAValue(iJet,iVertex,iVertices);
  
  int lPtId = 0; 
  if(iJet->Pt() > 10 && iJet->Pt() < 20) lPtId = 1;
  if(iJet->Pt() > 20 && iJet->Pt() < 30) lPtId = 2;
  if(iJet->Pt() > 30                   ) lPtId = 3;
  
  int lEtaId = 0;
  if(fabs(iJet->Eta()) > 2.5  && fabs(iJet->Eta()) < 2.75) lEtaId = 1; 
  if(fabs(iJet->Eta()) > 2.75 && fabs(iJet->Eta()) < 3.0 ) lEtaId = 2; 
  if(fabs(iJet->Eta()) > 3.0  && fabs(iJet->Eta()) < 5.0 ) lEtaId = 3; 
  
  double lMVACut = fMVACut[lPtId][lEtaId];
  if(lMVA < lMVACut) return false;
  return true;
  //if(lMVA < -0.8)                            return false;
  //if(lMVA < -0.5 && fabs(iJet->Eta()) > 3.0) return false;
  //if(iJet->Pt() < fJetPtMin && iJet->TrackCountingHighEffBJetTagsDisc() == -100) return false; //This line is a bug in the Met training
  //if( fabs(JetTools::impactParameter(iJet,iVertex,true)) < 0.2) return true;
}
//--------------------------------------------------------------------------------------------------
Double_t JetIDMVA::MVAValue(const PFJet *iJet,const Vertex *iVertex,const VertexCol *iVertices, //Vertex here is the PV
			    FactorizedJetCorrector *iJetCorrector,
			    const PileupEnergyDensityCol *iPileupEnergyDensity,
			    Bool_t printDebug) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: JetIDMVA not properly initialized.\n"; 
    return -9999;
  }
  if(!JetTools::passPFLooseId(iJet)) return -2.;

  //set all input variables
  fNVtx      = iVertices->GetEntries();
  fJPt1      = correctedPt(iJet,iJetCorrector,iPileupEnergyDensity);
  fJEta1     = iJet->RawMom().Eta();
  fJPhi1     = iJet->RawMom().Phi();
  fJD01      = JetTools::impactParameter(iJet,iVertex);  
  fJDZ1      = JetTools::impactParameter(iJet,iVertex,true);
  fBeta      = JetTools::Beta(iJet,iVertex,fDZCut);
  fBetaStar  = JetTools::betaStar(iJet,iVertex,iVertices,fDZCut);
  fNCharged  = iJet->ChargedMultiplicity();
  fNNeutrals = iJet->NeutralMultiplicity();

  fDRMean    = JetTools::dRMean(iJet,-1);
  fFrac01    = JetTools::frac  (iJet,0.1,0. ,-1);
  fFrac02    = JetTools::frac  (iJet,0.2,0.1,-1);
  fFrac03    = JetTools::frac  (iJet,0.3,0.2,-1);
  fFrac04    = JetTools::frac  (iJet,0.4,0.3,-1);
  fFrac05    = JetTools::frac  (iJet,0.5,0.4,-1);

  double lMVA = 0;
  if(fJPt1 < 10) lMVA = fReader->EvaluateMVA( fLowPtMethodName  );
  if(fJPt1 > 10) lMVA = fReader->EvaluateMVA( fHighPtMethodName );
  if (printDebug == kTRUE) {
    std::cout << "Debug Jet MVA: "
	      << fNVtx      << " "
	      << fJPt1      << " "
	      << fJEta1     << " "
	      << fJPhi1     << " "
	      << fJD01      << " "
	      << fJDZ1      << " "
	      << fBeta      << " "
	      << fBetaStar  << " "
	      << fNCharged  << " "
	      << fNNeutrals << " "
	      << fDRMean    << " "
	      << fFrac01    << " "
	      << fFrac02    << " "
	      << fFrac03    << " "
	      << fFrac04    << " "
	      << fFrac05    
              << " === : === "
              << lMVA << " "    
              << std::endl;
  }

  return lMVA;
}
Double_t JetIDMVA::MVAValue(const PFJet *iJet,const Vertex *iVertex, const VertexCol *iVertices,//Vertex here is the PV
			    Bool_t printDebug) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: JetIDMVA not properly initialized.\n"; 
    return -9999;
  }
  if(!JetTools::passPFLooseId(iJet)) return -2.;

  //set all input variables
  fNVtx      = iVertices->GetEntries();
  fJPt1      = iJet->Pt();
  fJEta1     = iJet->RawMom().Eta();
  fJPhi1     = iJet->RawMom().Phi();
  fJD01      = JetTools::impactParameter(iJet,iVertex);  
  fJDZ1      = JetTools::impactParameter(iJet,iVertex,true);
  fBeta      = JetTools::Beta(iJet,iVertex,fDZCut);
  fBetaStar  = JetTools::betaStar(iJet,iVertex,iVertices,fDZCut);
  fNCharged  = iJet->ChargedMultiplicity();
  fNNeutrals = iJet->NeutralMultiplicity();

  fDRMean    = JetTools::dRMean(iJet,-1);
  fFrac01    = JetTools::frac  (iJet,0.1,0. ,-1);
  fFrac02    = JetTools::frac  (iJet,0.2,0.1,-1);
  fFrac03    = JetTools::frac  (iJet,0.3,0.2,-1);
  fFrac04    = JetTools::frac  (iJet,0.4,0.3,-1);
  fFrac05    = JetTools::frac  (iJet,0.5,0.4,-1);

  double lMVA = 0;
  if(fJPt1 < 10) lMVA = fReader->EvaluateMVA( fLowPtMethodName  );
  if(fJPt1 > 10) lMVA = fReader->EvaluateMVA( fHighPtMethodName );
  
  if (printDebug == kTRUE) {
    std::cout << "Debug Jet MVA: "
	      << fNVtx      << " "
	      << fJPt1      << " "
	      << fJEta1     << " "
	      << fJPhi1     << " "
	      << fJD01      << " "
	      << fJDZ1      << " "
	      << fBeta      << " "
	      << fBetaStar  << " "
	      << fNCharged  << " "
	      << fNNeutrals << " "
	      << fDRMean    << " "
	      << fFrac01    << " "
	      << fFrac02    << " "
	      << fFrac03    << " "
	      << fFrac04    << " "
	      << fFrac05    
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
