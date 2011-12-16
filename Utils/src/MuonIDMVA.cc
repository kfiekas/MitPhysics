#include "MitPhysics/Utils/interface/MuonIDMVA.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include <TFile.h>
#include <TRandom3.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"


ClassImp(mithep::MuonIDMVA)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
MuonIDMVA::MuonIDMVA() :
fMethodname("BDTG method"),
fIsInitialized(kFALSE)
{
  // Constructor.
  for(UInt_t i=0; i<6; ++i) {
    fTMVAReader[i] = 0;
  }
}


//--------------------------------------------------------------------------------------------------
MuonIDMVA::~MuonIDMVA()
{
  for(UInt_t i=0; i<6; ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
}

//--------------------------------------------------------------------------------------------------
void MuonIDMVA::Initialize( TString methodName,
                            TString Subdet0Pt10To14p5Weights , 
                            TString Subdet1Pt10To14p5Weights , 
                            TString Subdet0Pt14p5To20Weights,
                            TString Subdet1Pt14p5To20Weights, 
                            TString Subdet0Pt20ToInfWeights, 
                            TString Subdet1Pt20ToInfWeights,
                            MuonIDMVA::MVAType type) {
  
  fIsInitialized = kTRUE;
  
  fMethodname = methodName;
    
  for(UInt_t i=0; i<6; ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];

    fTMVAReader[i] = new TMVA::Reader( "!Color:!Silent:Error" );  
    fTMVAReader[i]->SetVerbose(kTRUE);

    if (type == kV3) {
      fTMVAReader[i]->AddVariable( "TkNchi2",              &fMVAVar_MuTkNchi2               );
      fTMVAReader[i]->AddVariable( "GlobalNchi2",          &fMVAVar_MuGlobalNchi2           );
      fTMVAReader[i]->AddVariable( "NValidHits",           &fMVAVar_MuNValidHits            );
      fTMVAReader[i]->AddVariable( "NTrackerHits",         &fMVAVar_MuNTrackerHits          );
      fTMVAReader[i]->AddVariable( "NPixelHits",           &fMVAVar_MuNPixelHits            );
      fTMVAReader[i]->AddVariable( "NMatches",             &fMVAVar_MuNMatches              );
      fTMVAReader[i]->AddVariable( "D0",                   &fMVAVar_MuD0                    );      
      fTMVAReader[i]->AddVariable( "IP3d",                 &fMVAVar_MuIP3d                  );      
      fTMVAReader[i]->AddVariable( "IP3dSig",              &fMVAVar_MuIP3dSig               );      
      fTMVAReader[i]->AddVariable( "TrkKink",              &fMVAVar_MuTrkKink               );      
      fTMVAReader[i]->AddVariable( "SegmentCompatibility", &fMVAVar_MuSegmentCompatibility  );      
      fTMVAReader[i]->AddVariable( "CaloCompatibility",    &fMVAVar_MuCaloCompatibility     );      
      fTMVAReader[i]->AddVariable( "HadEnergyOverPt",      &fMVAVar_MuHadEnergyOverPt       );      
      if (i==0 || i==2 || i==4) {
        fTMVAReader[i]->AddVariable( "HoEnergyOverPt",     &fMVAVar_MuHoEnergyOverPt        );      
      }
      fTMVAReader[i]->AddVariable( "EmEnergyOverPt",       &fMVAVar_MuEmEnergyOverPt        );      
      fTMVAReader[i]->AddVariable( "HadS9EnergyOverPt",    &fMVAVar_MuHadS9EnergyOverPt     );      
      if (i==0 || i==2 || i==4) {
        fTMVAReader[i]->AddVariable( "HoS9EnergyOverPt",   &fMVAVar_MuHoS9EnergyOverPt      );      
      }
      fTMVAReader[i]->AddVariable( "EmS9EnergyOverPt",     &fMVAVar_MuEmS9EnergyOverPt      );      
    }
    
    if (i==0) fTMVAReader[i]->BookMVA(fMethodname , Subdet0Pt10To14p5Weights );
    if (i==1) fTMVAReader[i]->BookMVA(fMethodname , Subdet1Pt10To14p5Weights );
    if (i==2) fTMVAReader[i]->BookMVA(fMethodname , Subdet0Pt14p5To20Weights );
    if (i==3) fTMVAReader[i]->BookMVA(fMethodname , Subdet1Pt14p5To20Weights );
    if (i==4) fTMVAReader[i]->BookMVA(fMethodname , Subdet0Pt20ToInfWeights  );
    if (i==5) fTMVAReader[i]->BookMVA(fMethodname , Subdet1Pt20ToInfWeights  );

  }

  std::cout << "Muon ID MVA Initialization\n";
  std::cout << "MethodName : " << fMethodname << " , type == " << type << std::endl;
  std::cout << "Load weights file : " << Subdet0Pt10To14p5Weights << std::endl;
  std::cout << "Load weights file : " << Subdet1Pt10To14p5Weights << std::endl;
  std::cout << "Load weights file : " << Subdet0Pt14p5To20Weights << std::endl;
  std::cout << "Load weights file : " << Subdet1Pt14p5To20Weights << std::endl;
  std::cout << "Load weights file : " << Subdet0Pt20ToInfWeights << std::endl;
  std::cout << "Load weights file : " << Subdet1Pt20ToInfWeights << std::endl;

}

//--------------------------------------------------------------------------------------------------
Double_t MuonIDMVA::MVAValue(Double_t MuPt , Double_t MuEta,
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
                         Double_t                   MuEmS9EnergyOverPt 
  ) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: MuonIDMVA not properly initialized.\n"; 
    return -9999;
  }

  Int_t subdet = 0;
  if (fabs(MuEta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (MuPt > 14.5) ptBin = 1;
  if (MuPt > 20.0) ptBin = 2;

  
  //set all input variables
  fMVAVar_MuTkNchi2              = MuTkNchi2; 
  fMVAVar_MuGlobalNchi2          = MuGlobalNchi2; 
  fMVAVar_MuNValidHits           = MuNValidHits; 
  fMVAVar_MuNTrackerHits         = MuNTrackerHits; 
  fMVAVar_MuNPixelHits           = MuNPixelHits;  
  fMVAVar_MuNMatches             = MuNMatches; 
  fMVAVar_MuD0                   = MuD0; 
  fMVAVar_MuIP3d                 = MuIP3d; 
  fMVAVar_MuIP3dSig              = MuIP3dSig; 
  fMVAVar_MuTrkKink              = MuTrkKink; 
  fMVAVar_MuSegmentCompatibility = MuSegmentCompatibility; 
  fMVAVar_MuCaloCompatibility    = MuCaloCompatibility; 
  fMVAVar_MuHadEnergyOverPt      = MuHadEnergyOverPt; 
  fMVAVar_MuHoEnergyOverPt       = MuHoEnergyOverPt; 
  fMVAVar_MuEmEnergyOverPt       = MuEmEnergyOverPt; 
  fMVAVar_MuHadS9EnergyOverPt    = MuHadS9EnergyOverPt; 
  fMVAVar_MuHoS9EnergyOverPt     = MuHoS9EnergyOverPt; 
  fMVAVar_MuEmS9EnergyOverPt     = MuEmS9EnergyOverPt; 


  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;
  if (subdet == 0 && ptBin == 2) MVABin = 4;
  if (subdet == 1 && ptBin == 2) MVABin = 5;
  assert(MVABin >= 0 && MVABin <= 5);
  reader = fTMVAReader[MVABin];
                                                
  mva = reader->EvaluateMVA( fMethodname );

  return mva;
}



// //--------------------------------------------------------------------------------------------------
// Double_t MuonIDMVA::MVAValue(const Muon *ele, const Vertex *vertex) {
  
//   if (!fIsInitialized) { 
//     std::cout << "Error: MuonIDMVA not properly initialized.\n"; 
//     return -9999;
//   }

//   Int_t subdet = 0;
//   if (ele->SCluster()->AbsEta() < 1.0) subdet = 0;
//   else if (ele->SCluster()->AbsEta() < 1.479) subdet = 1;
//   else subdet = 2;
//   Int_t ptBin = 0;
//   if (ele->Pt() > 20.0) ptBin = 1;
  
//   //set all input variables
//   fMVAVar_EleSigmaIEtaIEta = ele->CoviEtaiEta() ; 
//   fMVAVar_EleDEtaIn = ele->DeltaEtaSuperClusterTrackAtVtx(); 
//   fMVAVar_EleDPhiIn = ele->DeltaPhiSuperClusterTrackAtVtx(); 
//   fMVAVar_EleHoverE = ele->HadronicOverEm(); 
//   fMVAVar_EleD0 = ele->BestTrk()->D0Corrected(*vertex); 
//   fMVAVar_EleDZ = ele->BestTrk()->DzCorrected(*vertex); 
//   fMVAVar_EleFBrem = ele->FBrem(); 
//   fMVAVar_EleEOverP = ele->ESuperClusterOverP(); 
//   fMVAVar_EleESeedClusterOverPout = ele->ESeedClusterOverPout(); 
//   if (!TMath::IsNaN(ele->SCluster()->Seed()->CoviPhiiPhi())) fMVAVar_EleSigmaIPhiIPhi = TMath::Sqrt(ele->SCluster()->Seed()->CoviPhiiPhi()); 
//   else fMVAVar_EleSigmaIPhiIPhi = ele->CoviEtaiEta();
//   fMVAVar_EleNBrem = ele->NumberOfClusters() - 1; 
//   fMVAVar_EleOneOverEMinusOneOverP = (1.0/(ele->SCluster()->Energy())) - 1.0 / ele->BestTrk()->P(); 
//   fMVAVar_EleESeedClusterOverPIn = ele->ESeedClusterOverPIn(); 
//   fMVAVar_EleIP3d = ele->Ip3dPV(); 
//   fMVAVar_EleIP3dSig = ele->Ip3dPVSignificance(); 

//   Double_t mva = -9999;  
//   TMVA::Reader *reader = 0;
//   Int_t MVABin = -1;
//   if (subdet == 0 && ptBin == 0) MVABin = 0;
//   if (subdet == 1 && ptBin == 0) MVABin = 1;
//   if (subdet == 2 && ptBin == 0) MVABin = 2;
//   if (subdet == 0 && ptBin == 1) MVABin = 3;
//   if (subdet == 1 && ptBin == 1) MVABin = 4;
//   if (subdet == 2 && ptBin == 1) MVABin = 5;
//   assert(MVABin >= 0 && MVABin <= 5);
//   reader = fTMVAReader[MVABin];
                                                
//   mva = reader->EvaluateMVA( fMethodname );

//   return mva;
// }
