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

    if (type == kV2) {
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
    }

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

    if (type == kV8) {
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
      fTMVAReader[i]->AddVariable( "EmEnergyOverPt",       &fMVAVar_MuEmEnergyOverPt        );      
      fTMVAReader[i]->AddVariable( "HadS9EnergyOverPt",    &fMVAVar_MuHadS9EnergyOverPt     );      
      fTMVAReader[i]->AddVariable( "EmS9EnergyOverPt",     &fMVAVar_MuEmS9EnergyOverPt      );      
      fTMVAReader[i]->AddVariable( "ChargedIso03OverPt",   &fMVAVar_MuChargedIso03OverPt    );
      fTMVAReader[i]->AddVariable( "NeutralIso03OverPt",   &fMVAVar_MuNeutralIso03OverPt    );      
      fTMVAReader[i]->AddVariable( "ChargedIso04OverPt",   &fMVAVar_MuChargedIso04OverPt    );
      fTMVAReader[i]->AddVariable( "NeutralIso04OverPt",   &fMVAVar_MuNeutralIso04OverPt    );      
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
                             Double_t                   MuEmS9EnergyOverPt,
                             Double_t                   MuChargedIso03OverPt,
                             Double_t                   MuNeutralIso03OverPt,
                             Double_t                   MuChargedIso04OverPt,
                             Double_t                   MuNeutralIso04OverPt                             
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
  fMVAVar_MuChargedIso03OverPt   = MuChargedIso03OverPt; 
  fMVAVar_MuNeutralIso03OverPt   = MuNeutralIso03OverPt; 
  fMVAVar_MuChargedIso04OverPt   = MuChargedIso04OverPt; 
  fMVAVar_MuNeutralIso04OverPt   = MuNeutralIso04OverPt; 

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

  Bool_t printdebug = kTRUE;
  if (printdebug == kTRUE) {
    std::cout << "Debug Muon MVA: "
	 << MuPt << " " << MuEta << " --> MVABin " << MVABin << " : "     
	 << fMVAVar_MuTkNchi2              << " " 
	 << fMVAVar_MuGlobalNchi2          << " " 
	 << fMVAVar_MuNValidHits           << " " 
	 << fMVAVar_MuNTrackerHits         << " " 
	 << fMVAVar_MuNPixelHits           << " "  
	 << fMVAVar_MuNMatches             << " " 
	 << fMVAVar_MuD0                   << " " 
	 << fMVAVar_MuIP3d                 << " " 
	 << fMVAVar_MuIP3dSig              << " " 
	 << fMVAVar_MuTrkKink              << " " 
	 << fMVAVar_MuSegmentCompatibility << " " 
	 << fMVAVar_MuCaloCompatibility    << " " 
	 << fMVAVar_MuHadEnergyOverPt      << " " 
	 << fMVAVar_MuHoEnergyOverPt       << " " 
	 << fMVAVar_MuEmEnergyOverPt       << " " 
	 << fMVAVar_MuHadS9EnergyOverPt    << " " 
	 << fMVAVar_MuHoS9EnergyOverPt     << " " 
	 << fMVAVar_MuEmS9EnergyOverPt     << " " 
	 << fMVAVar_MuChargedIso03OverPt   << " " 
	 << fMVAVar_MuNeutralIso03OverPt   << " " 
	 << fMVAVar_MuChargedIso04OverPt   << " " 
	 << fMVAVar_MuNeutralIso04OverPt   << " " 
	 << " === : === "
	 << mva 
	 << std::endl;
  }

  return mva;
}



//--------------------------------------------------------------------------------------------------
Double_t MuonIDMVA::MVAValue(const Muon *mu, const Vertex *vertex, MuonTools *fMuonTools,
                             const PFCandidateCol *PFCands, 
                             const PileupEnergyDensityCol *PileupEnergyDensity) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: MuonIDMVA not properly initialized.\n"; 
    return -9999;
  }

  const Track *muTrk=0;
  if(mu->HasTrackerTrk())         { muTrk = mu->TrackerTrk();    }
  else if(mu->HasStandaloneTrk()) { muTrk = mu->StandaloneTrk(); } 
  
  Double_t muNchi2 = 0.0; 
  if(mu->HasGlobalTrk())          { muNchi2 = mu->GlobalTrk()->RChi2();     }
  else if(mu->HasStandaloneTrk()) { muNchi2 = mu->StandaloneTrk()->RChi2(); }
  else if(mu->HasTrackerTrk())    { muNchi2 = mu->TrackerTrk()->RChi2();    }

  Double_t ChargedIso03 = IsolationTools::PFMuonIsolation(mu, PFCands, vertex, 0.1, 99999, 0.3, 0.0, 0.0);
  Double_t NeutralIso03_05Threshold = IsolationTools::PFMuonIsolation(mu, PFCands, vertex, 0.0, 0.5, 0.3, 0.0, 0.0);
  Double_t ChargedIso04 = IsolationTools::PFMuonIsolation(mu, PFCands, vertex, 0.1, 99999, 0.4, 0.0, 0.0);
  Double_t NeutralIso04_05Threshold = IsolationTools::PFMuonIsolation(mu, PFCands, vertex, 0.0, 0.5, 0.4, 0.0, 0.0);

  Double_t Rho = 0;
  if (!(TMath::IsNaN(PileupEnergyDensity->At(0)->Rho()) || isinf(PileupEnergyDensity->At(0)->Rho()))) Rho = PileupEnergyDensity->At(0)->Rho();

  Int_t subdet = 0;
  if (fabs(muTrk->Eta()) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (muTrk->Pt() > 14.5) ptBin = 1;
  if (muTrk->Pt() > 20.0) ptBin = 2;

  //set all input variables
  fMVAVar_MuTkNchi2              = muTrk->RChi2();
  fMVAVar_MuGlobalNchi2          = muNchi2;
  fMVAVar_MuNValidHits           = mu->NValidHits();
  fMVAVar_MuNTrackerHits         = muTrk->NHits();
  fMVAVar_MuNPixelHits           = muTrk->NPixelHits();
  fMVAVar_MuNMatches             = mu->NMatches();
  fMVAVar_MuD0                   = muTrk->D0Corrected(*vertex);
  fMVAVar_MuIP3d                 = mu->Ip3dPV();
  fMVAVar_MuIP3dSig              = mu->Ip3dPVSignificance();
  fMVAVar_MuTrkKink              = mu->TrkKink();
  fMVAVar_MuSegmentCompatibility = fMuonTools->GetSegmentCompatability(mu);
  fMVAVar_MuCaloCompatibility    = fMuonTools->GetCaloCompatability(mu, kTRUE, kTRUE);
  fMVAVar_MuHadEnergyOverPt      = (mu->HadEnergy() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuHadEnergy,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuHoEnergyOverPt       = (mu->HoEnergy() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuHoEnergy,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuEmEnergyOverPt       = (mu->EmEnergy() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuEmEnergy,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuHadS9EnergyOverPt    = (mu->HadS9Energy() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuHadS9Energy,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuHoS9EnergyOverPt     = (mu->HoS9Energy() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuHoS9Energy,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuEmS9EnergyOverPt     = (mu->EmS9Energy() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuEmS9Energy,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuChargedIso03OverPt   = (ChargedIso03 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuChargedIso03,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuNeutralIso03OverPt   = (NeutralIso03_05Threshold - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralIso03,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuChargedIso04OverPt   = (ChargedIso04 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuChargedIso04,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuNeutralIso04OverPt   = (NeutralIso04_05Threshold - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralIso04,muTrk->Eta()))/muTrk->Pt();

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

  Bool_t printdebug = kTRUE;
  if (printdebug == kTRUE) {
    std::cout << "Debug Muon MVA: "
              << mu->Pt() << " " << mu->Eta() << " " << mu->Phi() << " : "
              << muTrk->Pt() << " " << muTrk->Eta() << " --> MVABin " << MVABin << " : "     
              << fMVAVar_MuTkNchi2              << " " 
              << fMVAVar_MuGlobalNchi2          << " " 
              << fMVAVar_MuNValidHits           << " " 
              << fMVAVar_MuNTrackerHits         << " " 
              << fMVAVar_MuNPixelHits           << " "  
              << fMVAVar_MuNMatches             << " " 
              << fMVAVar_MuD0                   << " " 
              << fMVAVar_MuIP3d                 << " " 
              << fMVAVar_MuIP3dSig              << " " 
              << fMVAVar_MuTrkKink              << " " 
              << fMVAVar_MuSegmentCompatibility << " " 
              << fMVAVar_MuCaloCompatibility    << " " 
              << fMVAVar_MuHadEnergyOverPt      << " " 
              << fMVAVar_MuHoEnergyOverPt       << " " 
              << fMVAVar_MuEmEnergyOverPt       << " " 
              << fMVAVar_MuHadS9EnergyOverPt    << " " 
              << fMVAVar_MuHoS9EnergyOverPt     << " " 
              << fMVAVar_MuEmS9EnergyOverPt     << " " 
              << fMVAVar_MuChargedIso03OverPt   << " " 
              << fMVAVar_MuNeutralIso03OverPt   << " " 
              << fMVAVar_MuChargedIso04OverPt   << " " 
              << fMVAVar_MuNeutralIso04OverPt   << " " 
              << " === : === "
              << mva 
              << std::endl;
  }

  return mva;
}
