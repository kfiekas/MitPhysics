#include "MitPhysics/Utils/interface/ElectronIDMVA.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include <TFile.h>
#include <TRandom3.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"


ClassImp(mithep::ElectronIDMVA)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
ElectronIDMVA::ElectronIDMVA() :
fMethodname("BDTG method"),
fIsInitialized(kFALSE)
{
  // Constructor.
  for(UInt_t i=0; i<6; ++i) {
    fTMVAReader[i] = 0;
  }
}


//--------------------------------------------------------------------------------------------------
ElectronIDMVA::~ElectronIDMVA()
{
  for(UInt_t i=0; i<6; ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMVA::Initialize( TString methodName,
                                TString Subdet0Pt10To20Weights , 
                                TString Subdet1Pt10To20Weights , 
                                TString Subdet2Pt10To20Weights,
                                TString Subdet0Pt20ToInfWeights,
                                TString Subdet1Pt20ToInfWeights, 
                                TString Subdet2Pt20ToInfWeights,
                                ElectronIDMVA::MVAType type) {

  fIsInitialized = kTRUE;
  
  fMethodname = methodName;
    
  for(UInt_t i=0; i<6; ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];

    fTMVAReader[i] = new TMVA::Reader( "!Color:!Silent:Error" );  
    fTMVAReader[i]->SetVerbose(kTRUE);

    if (type == kBaseline) {
      fTMVAReader[i]->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta            );
      fTMVAReader[i]->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                   );
      fTMVAReader[i]->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                   );
      fTMVAReader[i]->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                    );
      fTMVAReader[i]->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi            );
      fTMVAReader[i]->AddVariable( "NBrem",                 &fMVAVar_EleNBrem                    );
      fTMVAReader[i]->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP    );      
    }
    
    if (type == kNoIPInfo) {
      fTMVAReader[i]->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta            );
      fTMVAReader[i]->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                   );
      fTMVAReader[i]->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                   );
      fTMVAReader[i]->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                    );
      fTMVAReader[i]->AddVariable( "EOverP",                &fMVAVar_EleEOverP                   );
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPout",  &fMVAVar_EleESeedClusterOverPout     );
      fTMVAReader[i]->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi            );
      fTMVAReader[i]->AddVariable( "NBrem",                 &fMVAVar_EleNBrem                    );
      fTMVAReader[i]->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP    );      
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPIn",   &fMVAVar_EleESeedClusterOverPIn      );
    }
    if (type == kWithIPInfo) {
      fTMVAReader[i]->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta            );
      fTMVAReader[i]->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                   );
      fTMVAReader[i]->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                   );
      fTMVAReader[i]->AddVariable( "D0",                    &fMVAVar_EleD0                       );
      fTMVAReader[i]->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                    );
      fTMVAReader[i]->AddVariable( "EOverP",                &fMVAVar_EleEOverP                   );
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPout",  &fMVAVar_EleESeedClusterOverPout     );
      fTMVAReader[i]->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi            );
      fTMVAReader[i]->AddVariable( "NBrem",                 &fMVAVar_EleNBrem                    );
      fTMVAReader[i]->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP    );      
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPIn",   &fMVAVar_EleESeedClusterOverPIn      );
      fTMVAReader[i]->AddVariable( "IP3d",                  &fMVAVar_EleIP3d                     );
      fTMVAReader[i]->AddVariable( "IP3dSig",               &fMVAVar_EleIP3dSig                  );
    }
    if (type == kIDIsoCombined) {
      fTMVAReader[i]->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta            );
      fTMVAReader[i]->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                   );
      fTMVAReader[i]->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                   );
      fTMVAReader[i]->AddVariable( "D0",                    &fMVAVar_EleD0                       );
      fTMVAReader[i]->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                    );
      fTMVAReader[i]->AddVariable( "EOverP",                &fMVAVar_EleEOverP                   );
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPout",  &fMVAVar_EleESeedClusterOverPout     );
      fTMVAReader[i]->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi            );
      fTMVAReader[i]->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP    );      
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPIn",   &fMVAVar_EleESeedClusterOverPIn      );
      fTMVAReader[i]->AddVariable( "IP3d",                  &fMVAVar_EleIP3d                     );
      fTMVAReader[i]->AddVariable( "IP3dSig",               &fMVAVar_EleIP3dSig                  );

      fTMVAReader[i]->AddVariable( "GsfTrackChi2OverNdof",  &fMVAVar_EleGsfTrackChi2OverNdof     );
      fTMVAReader[i]->AddVariable( "dEtaCalo",              &fMVAVar_EledEtaCalo                 );
      fTMVAReader[i]->AddVariable( "dPhiCalo",              &fMVAVar_EledPhiCalo                 );
      fTMVAReader[i]->AddVariable( "R9",                    &fMVAVar_EleR9                       );
      fTMVAReader[i]->AddVariable( "SCEtaWidth",            &fMVAVar_EleSCEtaWidth               );
      fTMVAReader[i]->AddVariable( "SCPhiWidth",            &fMVAVar_EleSCPhiWidth               );
      fTMVAReader[i]->AddVariable( "CovIEtaIPhi",           &fMVAVar_EleCovIEtaIPhi              );
      if (i == 2 || i == 5) {
        fTMVAReader[i]->AddVariable( "PreShowerOverRaw",      &fMVAVar_ElePreShowerOverRaw       );
      }
      fTMVAReader[i]->AddVariable( "ChargedIso03",          &fMVAVar_EleChargedIso03OverPt       );
      fTMVAReader[i]->AddVariable( "NeutralHadronIso03",    &fMVAVar_EleNeutralHadronIso03OverPt );
      fTMVAReader[i]->AddVariable( "GammaIso03",            &fMVAVar_EleGammaIso03OverPt         );
      fTMVAReader[i]->AddVariable( "ChargedIso04",          &fMVAVar_EleChargedIso04OverPt       );
      fTMVAReader[i]->AddVariable( "NeutralHadronIso04",    &fMVAVar_EleNeutralHadronIso04OverPt );
      fTMVAReader[i]->AddVariable( "GammaIso04",            &fMVAVar_EleGammaIso04OverPt         );

    }

    if (i==0) fTMVAReader[i]->BookMVA(fMethodname , Subdet0Pt10To20Weights );
    if (i==1) fTMVAReader[i]->BookMVA(fMethodname , Subdet1Pt10To20Weights );
    if (i==2) fTMVAReader[i]->BookMVA(fMethodname , Subdet2Pt10To20Weights );
    if (i==3) fTMVAReader[i]->BookMVA(fMethodname , Subdet0Pt20ToInfWeights );
    if (i==4) fTMVAReader[i]->BookMVA(fMethodname , Subdet1Pt20ToInfWeights );
    if (i==5) fTMVAReader[i]->BookMVA(fMethodname , Subdet2Pt20ToInfWeights );

  }

  std::cout << "Electron ID MVA Initialization\n";
  std::cout << "MethodName : " << fMethodname << " , type == " << type << std::endl;
  std::cout << "Load weights file : " << Subdet0Pt10To20Weights << std::endl;
  std::cout << "Load weights file : " << Subdet1Pt10To20Weights << std::endl;
  std::cout << "Load weights file : " << Subdet2Pt10To20Weights << std::endl;
  std::cout << "Load weights file : " << Subdet0Pt20ToInfWeights << std::endl;
  std::cout << "Load weights file : " << Subdet1Pt20ToInfWeights << std::endl;
  std::cout << "Load weights file : " << Subdet2Pt20ToInfWeights << std::endl;

}

//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVA::MVAValue(Double_t ElePt , Double_t EleEta,
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
                                 Double_t EleIP3dSig
  ) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
    return -9999;
  }

  Int_t subdet = 0;
  if (fabs(EleEta) < 1.0) subdet = 0;
  else if (fabs(EleEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ElePt > 20.0) ptBin = 1;
  
  //set all input variables
  fMVAVar_EleSigmaIEtaIEta = EleSigmaIEtaIEta;
  fMVAVar_EleDEtaIn = EleDEtaIn;
  fMVAVar_EleDPhiIn = EleDPhiIn;
  fMVAVar_EleHoverE = EleHoverE;
  fMVAVar_EleD0 = EleD0;
  fMVAVar_EleDZ = EleDZ;
  fMVAVar_EleFBrem = EleFBrem;
  fMVAVar_EleEOverP = EleEOverP;
  fMVAVar_EleESeedClusterOverPout = EleESeedClusterOverPout;
  fMVAVar_EleSigmaIPhiIPhi = EleSigmaIPhiIPhi;
  fMVAVar_EleNBrem = EleNBrem;
  fMVAVar_EleOneOverEMinusOneOverP = EleOneOverEMinusOneOverP;
  fMVAVar_EleESeedClusterOverPIn = EleESeedClusterOverPIn;
  fMVAVar_EleIP3d = EleIP3d;
  fMVAVar_EleIP3dSig = EleIP3dSig;

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;
  assert(MVABin >= 0 && MVABin <= 5);
  reader = fTMVAReader[MVABin];
                                                
  mva = reader->EvaluateMVA( fMethodname );

  return mva;
}

//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVA::MVAValue(Double_t ElePt , Double_t EleEta, Double_t PileupEnergyDensity,
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
                                 Double_t EleGammaIso04
  ) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
    return -9999;
  }

  Int_t subdet = 0;
  if (fabs(EleEta) < 1.0) subdet = 0;
  else if (fabs(EleEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ElePt > 20.0) ptBin = 1;
  
  //set all input variables
  fMVAVar_EleSigmaIEtaIEta = EleSigmaIEtaIEta;
  fMVAVar_EleDEtaIn = EleDEtaIn;
  fMVAVar_EleDPhiIn = EleDPhiIn;
  fMVAVar_EleHoverE = EleHoverE;
  fMVAVar_EleD0 = EleD0;
  fMVAVar_EleDZ = EleDZ;
  fMVAVar_EleFBrem = EleFBrem;
  fMVAVar_EleEOverP = EleEOverP;
  fMVAVar_EleESeedClusterOverPout = EleESeedClusterOverPout;
  fMVAVar_EleSigmaIPhiIPhi = EleSigmaIPhiIPhi;
  fMVAVar_EleNBrem = EleNBrem;
  fMVAVar_EleOneOverEMinusOneOverP = EleOneOverEMinusOneOverP;
  fMVAVar_EleESeedClusterOverPIn = EleESeedClusterOverPIn;
  fMVAVar_EleIP3d = EleIP3d;
  fMVAVar_EleIP3dSig = EleIP3dSig;
  fMVAVar_EleGsfTrackChi2OverNdof = EleGsfTrackChi2OverNdof;
  fMVAVar_EledEtaCalo = EledEtaCalo;
  fMVAVar_EledPhiCalo = EledPhiCalo;
  fMVAVar_EleR9 = EleR9;
  fMVAVar_EleSCEtaWidth = EleSCEtaWidth;
  fMVAVar_EleSCPhiWidth = EleSCPhiWidth;
  fMVAVar_EleCovIEtaIPhi = EleCovIEtaIPhi;
  fMVAVar_ElePreShowerOverRaw = ElePreShowerOverRaw;
  fMVAVar_EleChargedIso03OverPt 
    = (EleChargedIso03 
       - PileupEnergyDensity * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleChargedIso03, EleEta)) / ElePt;
  fMVAVar_EleNeutralHadronIso03OverPt 
    = (EleNeutralHadronIso03
       - PileupEnergyDensity * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso03, EleEta) 
       + PileupEnergyDensity * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso007,EleEta)) / ElePt;
  fMVAVar_EleGammaIso03OverPt 
    = (EleGammaIso03 
       - PileupEnergyDensity * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIso03, EleEta) 
       + PileupEnergyDensity * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIsoVetoEtaStrip03,EleEta))/ElePt;
  fMVAVar_EleChargedIso04OverPt 
    = (EleChargedIso04 
       - PileupEnergyDensity * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleChargedIso04, EleEta))/ElePt;
  fMVAVar_EleNeutralHadronIso04OverPt
    = (EleNeutralHadronIso04 
       - PileupEnergyDensity * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso04, EleEta) 
       + PileupEnergyDensity * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso007,EleEta))/ElePt;
  fMVAVar_EleGammaIso04OverPt 
    = (EleGammaIso04 
       - PileupEnergyDensity * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIso04, EleEta) 
       + PileupEnergyDensity * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIsoVetoEtaStrip04,EleEta))/ElePt;
  
  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;
  assert(MVABin >= 0 && MVABin <= 5);
  reader = fTMVAReader[MVABin];
                                                
  mva = reader->EvaluateMVA( fMethodname );

  Bool_t printdebug = kTRUE;
  if (printdebug == kTRUE) {
    std::cout << "Debug Electron MVA: "
	 << ElePt << " " << EleEta << " " << " --> MVABin " << MVABin << " : "     
	 << fMVAVar_EleSigmaIEtaIEta << " " 
	 << fMVAVar_EleDEtaIn << " " 
	 << fMVAVar_EleDPhiIn << " " 
	 << fMVAVar_EleHoverE << " " 
	 << fMVAVar_EleD0 << " " 
	 << fMVAVar_EleDZ << " " 
	 << fMVAVar_EleFBrem << " " 
	 << fMVAVar_EleEOverP << " " 
	 << fMVAVar_EleESeedClusterOverPout << " " 
	 << fMVAVar_EleSigmaIPhiIPhi << " " 
	 << fMVAVar_EleNBrem << " " 
	 << fMVAVar_EleOneOverEMinusOneOverP << " " 
	 << fMVAVar_EleESeedClusterOverPIn << " " 
	 << fMVAVar_EleIP3d << " " 
	 << fMVAVar_EleIP3dSig << " " 
	 << fMVAVar_EleGsfTrackChi2OverNdof << " "
	 << fMVAVar_EledEtaCalo << " "
	 << fMVAVar_EledPhiCalo << " "
	 << fMVAVar_EleR9 << " "
	 << fMVAVar_EleSCEtaWidth << " "
	 << fMVAVar_EleSCPhiWidth << " "
	 << fMVAVar_EleCovIEtaIPhi << " "
	 << fMVAVar_ElePreShowerOverRaw << " "
	 << fMVAVar_EleChargedIso03OverPt  << " "
	 << fMVAVar_EleNeutralHadronIso03OverPt  << " "
	 << fMVAVar_EleGammaIso03OverPt  << " "
	 << fMVAVar_EleChargedIso04OverPt  << " "
	 << fMVAVar_EleNeutralHadronIso04OverPt  << " "
	 << fMVAVar_EleGammaIso04OverPt  << " "
	 << " === : === "
	 << mva 
	 << std::endl;
  }

  return mva;
}


//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVA::MVAValue(const Electron *ele, const Vertex *vertex, 
                                 const PFCandidateCol *PFCands, 
                                 const PileupEnergyDensityCol *PileupEnergyDensity) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
    return -9999;
  }

  Int_t subdet = 0;
  if (ele->SCluster()->AbsEta() < 1.0) subdet = 0;
  else if (ele->SCluster()->AbsEta() < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ele->Pt() > 20.0) ptBin = 1;
  
  //set all input variables
  fMVAVar_EleSigmaIEtaIEta = ele->CoviEtaiEta() ; 
  fMVAVar_EleDEtaIn = ele->DeltaEtaSuperClusterTrackAtVtx(); 
  fMVAVar_EleDPhiIn = ele->DeltaPhiSuperClusterTrackAtVtx(); 
  fMVAVar_EleHoverE = ele->HadronicOverEm(); 
  fMVAVar_EleD0 = ele->BestTrk()->D0Corrected(*vertex); 
  fMVAVar_EleDZ = ele->BestTrk()->DzCorrected(*vertex); 
  fMVAVar_EleFBrem = ele->FBrem(); 
  fMVAVar_EleEOverP = ele->ESuperClusterOverP(); 
  fMVAVar_EleESeedClusterOverPout = ele->ESeedClusterOverPout(); 
  if (!TMath::IsNaN(ele->SCluster()->Seed()->CoviPhiiPhi())) fMVAVar_EleSigmaIPhiIPhi = TMath::Sqrt(ele->SCluster()->Seed()->CoviPhiiPhi()); 
  else fMVAVar_EleSigmaIPhiIPhi = ele->CoviEtaiEta();
  fMVAVar_EleNBrem = ele->NumberOfClusters() - 1; 
  fMVAVar_EleOneOverEMinusOneOverP = (1.0/(ele->SCluster()->Energy())) - 1.0 / ele->BestTrk()->P(); 
  fMVAVar_EleESeedClusterOverPIn = ele->ESeedClusterOverPIn(); 
  fMVAVar_EleIP3d = ele->Ip3dPV(); 
  fMVAVar_EleIP3dSig = ele->Ip3dPVSignificance(); 
  fMVAVar_EleGsfTrackChi2OverNdof = ele->BestTrk()->Chi2() / ele->BestTrk()->Ndof();
  fMVAVar_EledEtaCalo =  ele->DeltaEtaSeedClusterTrackAtCalo();
  fMVAVar_EledPhiCalo = ele->DeltaPhiSeedClusterTrackAtCalo();
  fMVAVar_EleR9 = ele->SCluster()->R9();
  fMVAVar_EleSCEtaWidth = ele->SCluster()->EtaWidth();
  fMVAVar_EleSCPhiWidth = ele->SCluster()->PhiWidth();
  fMVAVar_EleCovIEtaIPhi = ele->SCluster()->Seed()->CoviEtaiPhi();
  fMVAVar_ElePreShowerOverRaw = ele->SCluster()->PreshowerEnergy() / ele->SCluster()->RawEnergy();
  fMVAVar_EleChargedIso03OverPt 
    = (IsolationTools::PFElectronIsolation(ele, PFCands, vertex, 0.1, 99999, 0.3, 0.0) 
       - PileupEnergyDensity->At(0)->Rho() * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleChargedIso03, ele->SCluster()->Eta())) / ele->Pt();
  fMVAVar_EleNeutralHadronIso03OverPt 
    = (IsolationTools::PFElectronIsolation(ele, PFCands, vertex, 0.1, 0.5, 0.3, 0.0, PFCandidate::eNeutralHadron) 
       - PileupEnergyDensity->At(0)->Rho() * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso03, ele->SCluster()->Eta()) 
       + PileupEnergyDensity->At(0)->Rho() * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso007,ele->SCluster()->Eta())) / ele->Pt();
  fMVAVar_EleGammaIso03OverPt 
    = (IsolationTools::PFElectronIsolation(ele, PFCands, vertex, 0.1, 0.5, 0.3, 0.0, PFCandidate::eGamma) 
       - PileupEnergyDensity->At(0)->Rho() * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIso03, ele->SCluster()->Eta()) 
       + PileupEnergyDensity->At(0)->Rho() * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIsoVetoEtaStrip03,ele->SCluster()->Eta())) / ele->Pt();
  fMVAVar_EleChargedIso04OverPt 
    = (IsolationTools::PFElectronIsolation(ele, PFCands, vertex, 0.1, 99999, 0.4, 0.0) 
       - PileupEnergyDensity->At(0)->Rho() * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleChargedIso04, ele->SCluster()->Eta())) / ele->Pt();
  fMVAVar_EleNeutralHadronIso04OverPt 
    = (IsolationTools::PFElectronIsolation(ele, PFCands, vertex, 0.1, 0.5, 0.4, 0.0, PFCandidate::eNeutralHadron) 
       - PileupEnergyDensity->At(0)->Rho() * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso04, ele->SCluster()->Eta()) 
       + PileupEnergyDensity->At(0)->Rho() * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso007,ele->SCluster()->Eta())) / ele->Pt() ;
  fMVAVar_EleGammaIso04OverPt 
    = (IsolationTools::PFElectronIsolation(ele, PFCands, vertex, 0.1, 0.5, 0.4, 0.0, PFCandidate::eGamma) 
       - PileupEnergyDensity->At(0)->Rho() * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIso04, ele->SCluster()->Eta()) 
       + PileupEnergyDensity->At(0)->Rho() * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIsoVetoEtaStrip04,ele->SCluster()->Eta())) / ele->Pt();
  
  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;
  assert(MVABin >= 0 && MVABin <= 5);
  reader = fTMVAReader[MVABin];
                                                
  mva = reader->EvaluateMVA( fMethodname );

  Bool_t printdebug = kFALSE;
  if (printdebug == kTRUE) {
    std::cout << "Debug Electron MVA: "
              << ele->Pt() << " " << ele->Eta() << " " << ele->Phi() << " : "
              << ele->Pt() << " " << ele->SCluster()->AbsEta() << " --> MVABin " << MVABin << " : "     
              << fMVAVar_EleSigmaIEtaIEta << " " 
              << fMVAVar_EleDEtaIn << " " 
              << fMVAVar_EleDPhiIn << " " 
              << fMVAVar_EleHoverE << " " 
              << fMVAVar_EleD0 << " " 
              << fMVAVar_EleDZ << " " 
              << fMVAVar_EleFBrem << " " 
              << fMVAVar_EleEOverP << " " 
              << fMVAVar_EleESeedClusterOverPout << " " 
              << fMVAVar_EleSigmaIPhiIPhi << " " 
              << fMVAVar_EleNBrem << " " 
              << fMVAVar_EleOneOverEMinusOneOverP << " " 
              << fMVAVar_EleESeedClusterOverPIn << " " 
              << fMVAVar_EleIP3d << " " 
              << fMVAVar_EleIP3dSig << " " 
              << fMVAVar_EleGsfTrackChi2OverNdof << " "
              << fMVAVar_EledEtaCalo << " "
              << fMVAVar_EledPhiCalo << " "
              << fMVAVar_EleR9 << " "
              << fMVAVar_EleSCEtaWidth << " "
              << fMVAVar_EleSCPhiWidth << " "
              << fMVAVar_EleCovIEtaIPhi << " "
              << fMVAVar_ElePreShowerOverRaw << " "
              << fMVAVar_EleChargedIso03OverPt  << " "
              << fMVAVar_EleNeutralHadronIso03OverPt  << " "
              << fMVAVar_EleGammaIso03OverPt  << " "
              << fMVAVar_EleChargedIso04OverPt  << " "
              << fMVAVar_EleNeutralHadronIso04OverPt  << " "
              << fMVAVar_EleGammaIso04OverPt  << " "
              << " === : === "
              << mva << " "    
              << std::endl;
    
  }

  return mva;
}

//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVA::MVAValue(const Electron *ele, const Vertex *vertex) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
    return -9999;
  }

  Int_t subdet = 0;
  if (ele->SCluster()->AbsEta() < 1.0) subdet = 0;
  else if (ele->SCluster()->AbsEta() < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ele->Pt() > 20.0) ptBin = 1;
  
  //set all input variables
  fMVAVar_EleSigmaIEtaIEta = ele->CoviEtaiEta() ; 
  fMVAVar_EleDEtaIn = ele->DeltaEtaSuperClusterTrackAtVtx(); 
  fMVAVar_EleDPhiIn = ele->DeltaPhiSuperClusterTrackAtVtx(); 
  fMVAVar_EleHoverE = ele->HadronicOverEm(); 
  fMVAVar_EleD0 = ele->BestTrk()->D0Corrected(*vertex); 
  fMVAVar_EleDZ = ele->BestTrk()->DzCorrected(*vertex); 
  fMVAVar_EleFBrem = ele->FBrem(); 
  fMVAVar_EleEOverP = ele->ESuperClusterOverP(); 
  fMVAVar_EleESeedClusterOverPout = ele->ESeedClusterOverPout(); 
  if (!TMath::IsNaN(ele->SCluster()->Seed()->CoviPhiiPhi())) fMVAVar_EleSigmaIPhiIPhi = TMath::Sqrt(ele->SCluster()->Seed()->CoviPhiiPhi()); 
  else fMVAVar_EleSigmaIPhiIPhi = ele->CoviEtaiEta();
  fMVAVar_EleNBrem = ele->NumberOfClusters() - 1; 
  fMVAVar_EleOneOverEMinusOneOverP = (1.0/(ele->SCluster()->Energy())) - 1.0 / ele->BestTrk()->P(); 
  fMVAVar_EleESeedClusterOverPIn = ele->ESeedClusterOverPIn(); 
  fMVAVar_EleIP3d = ele->Ip3dPV(); 
  fMVAVar_EleIP3dSig = ele->Ip3dPVSignificance(); 

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;
  assert(MVABin >= 0 && MVABin <= 5);
  reader = fTMVAReader[MVABin];
                                                
  mva = reader->EvaluateMVA( fMethodname );

  return mva;
}
