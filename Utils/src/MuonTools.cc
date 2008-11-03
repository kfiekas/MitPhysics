#include "MitPhysics/Utils/interface/MuonTools.h"
using namespace mithep;

MuonTools::MuonTools() :
  fpion_em_etaEmi(0),
  fpion_had_etaEmi(0),
  fpion_had_etaTmi(0),
  fpion_em_etaB(0),
  fpion_had_etaB(0),
  fpion_ho_etaB(0),
  fpion_had_etaTpl(0),
  fpion_em_etaEpl(0),
  fpion_had_etaEpl(0),
  fmuon_em_etaEmi(0),
  fmuon_had_etaEmi(0),
  fmuon_had_etaTmi(0),
  fmuon_em_etaB(0),
  fmuon_had_etaB(0),
  fmuon_ho_etaB(0),
  fmuon_had_etaTpl(0),
  fmuon_em_etaEpl(0),
  fmuon_had_etaEpl(0) {}
MuonTools::~MuonTools() {
  fPion_templates->Close(); 
  fMuon_templates->Close();
  delete fpion_em_etaEmi;
  delete fpion_had_etaEmi;
  delete fpion_had_etaTmi;
  delete fpion_em_etaB;
  delete fpion_had_etaB;
  delete fpion_ho_etaB;
  delete fpion_had_etaTpl;
  delete fpion_em_etaEpl;
  delete fpion_had_etaEpl;
  delete fmuon_em_etaEmi;
  delete fmuon_had_etaEmi;
  delete fmuon_had_etaTmi;
  delete fmuon_em_etaB;
  delete fmuon_had_etaB;
  delete fmuon_ho_etaB;
  delete fmuon_had_etaTpl;
  delete fmuon_em_etaEpl;
  delete fmuon_had_etaEpl;
}

double MuonTools::getCaloCompatability(mithep::Muon* iMuon,bool iEMSpecial, bool iCorrectedHCAL) {
  if(fpion_em_etaEmi == 0) {
    TFile* fPion_templates = new TFile("$CMSSW_BASE/src/MitPhysics/Init/PionCaloTemplate.root","READ");
    TFile* fMuon_templates = new TFile("$CMSSW_BASE/src/MitPhysics/Init/MuonCaloTemplate.root","READ");
    fpion_em_etaEmi  = (TH2D*) fPion_templates->Get("em_etaEmi");
    fpion_had_etaEmi = (TH2D*) fPion_templates->Get("had_etaEmi");
    fpion_had_etaTmi = (TH2D*) fPion_templates->Get("had_etaTmi");
    fpion_em_etaB    = (TH2D*) fPion_templates->Get("em_etaB")   ;
    fpion_had_etaB   = (TH2D*) fPion_templates->Get("had_etaB")  ;
    fpion_ho_etaB    = (TH2D*) fPion_templates->Get("ho_etaB")   ;
    fpion_had_etaTpl = (TH2D*) fPion_templates->Get("had_etaTpl");
    fpion_em_etaEpl  = (TH2D*) fPion_templates->Get("em_etaEpl") ;
    fpion_had_etaEpl = (TH2D*) fPion_templates->Get("had_etaEpl");
    fmuon_em_etaEmi  = (TH2D*) fMuon_templates->Get("em_etaEmi") ;
    fmuon_had_etaEmi = (TH2D*) fMuon_templates->Get("had_etaEmi");
    fmuon_had_etaTmi = (TH2D*) fMuon_templates->Get("had_etaTmi");
    fmuon_em_etaB    = (TH2D*) fMuon_templates->Get("em_etaB")   ;
    fmuon_had_etaB   = (TH2D*) fMuon_templates->Get("had_etaB")  ;
    fmuon_ho_etaB    = (TH2D*) fMuon_templates->Get("ho_etaB")   ;
    fmuon_had_etaTpl = (TH2D*) fMuon_templates->Get("had_etaTpl");
    fmuon_em_etaEpl  = (TH2D*) fMuon_templates->Get("em_etaEpl");
    fmuon_had_etaEpl = (TH2D*) fMuon_templates->Get("had_etaEpl");
  }
  double lEta = -1.; double lP = -1;
  double lEM  = -5.;      double lHad = 0;      double lHO = 0;
  lEta = iMuon->Eta();
  lP   = iMuon->P(); 
  if(lP >= 2000.) lP = 1999.9;
  if(!iEMSpecial || iMuon->EmEnergy() != 0.) lEM  = iMuon->EmEnergy();
  lHad = iMuon->HadEnergy();
  lHO  = iMuon->HoEnergy();
  if(lP < 0. )           return 0.5; 
  if(fabs(lEta) >  2.5 ) return 0.5; 
  TH2D* lTMuonHad = NULL;
  TH2D* lTPionHad = NULL;
  TH2D* lTMuonHo  = NULL;
  TH2D* lTPionHo  = NULL;
  TH2D* lTMuonEm  = NULL;
  TH2D* lTPionEm  = NULL;
  
  if(fabs(lEta) >=  1.27) {
    if(iCorrectedHCAL) lHad *= 1.8/2.2;
    if(lEta > 0) {
      lTPionHad = fpion_had_etaEpl;
      lTMuonHad = fmuon_had_etaEpl;
    } else {
      lTPionHad = fpion_had_etaEmi;
      lTMuonHad = fmuon_had_etaEmi;
    }
  }
  if(fabs(lEta) <  1.27  && fabs(lEta) >=  1.1 ) {
    if(iCorrectedHCAL)    lHad *= (1.8/(-2.2*fabs(lEta)+5.5));
    if(lEta > 0) {
      lTPionHad  = fpion_had_etaTpl;
      lTMuonHad  = fmuon_had_etaTpl;
    } else {
      lTPionHad  = fpion_had_etaTmi;
      lTMuonHad  = fmuon_had_etaTmi;
    }
  }
  if(fabs(lEta) <  1.1) {
    if(iCorrectedHCAL)    lHad *= sin(2*atan(exp(iMuon->Eta())));
    lTPionHad  = fpion_had_etaB;
    lTMuonHad  = fmuon_had_etaB;
  }
  if(lEta >  1.479  ) {
    lTPionEm  = fpion_em_etaEpl;
    lTMuonEm  = fmuon_em_etaEpl;
  }
  if(fabs(lEta) <=  1.479) {
    lTPionEm  = fpion_em_etaB;
    lTMuonEm  = fmuon_em_etaB;
  }
  if(lEta < -1.479 ) {
    lTPionEm  = fpion_em_etaEmi;
    lTMuonEm  = fmuon_em_etaEmi;
  }
  if(fabs(lEta) < 1.28) {
    lTPionHo  = fpion_ho_etaB;
    lTMuonHo  = fmuon_ho_etaB;
  }
  
  double lPBX = 1.;     double lPSX = 1.; 
  double lPBY = 1.;     double lPSY = 1.; 
  double lPBZ = 1.;     double lPSZ = 1.; 
  if(!overflow(lTPionEm, lP,lEM))  lPBX =  lTPionEm ->GetBinContent(lTPionEm ->GetXaxis()->FindBin(lP),lTPionEm ->GetYaxis()->FindBin(lEM) );
  if(!overflow(lTPionHad,lP,lHad)) lPBY =  lTPionHad->GetBinContent(lTPionHad->GetXaxis()->FindBin(lP),lTPionHad->GetYaxis()->FindBin(lHad));
  if(!overflow(lTPionHo, lP,lHO))  lPBZ =  lTPionHo ->GetBinContent(lTPionHo ->GetXaxis()->FindBin(lP),lTPionHo ->GetYaxis()->FindBin(lHO) );
  if(!overflow(lTMuonEm, lP,lEM )) lPSX =  lTMuonEm ->GetBinContent(lTMuonEm ->GetXaxis()->FindBin(lP),lTMuonEm ->GetYaxis()->FindBin(lEM) );
  if(!overflow(lTMuonHad,lP,lHad)) lPSY =  lTMuonHad->GetBinContent(lTMuonHad->GetXaxis()->FindBin(lP),lTMuonHad->GetYaxis()->FindBin(lHad));
  if(!overflow(lTMuonHo ,lP,lHO))  lPSZ =  lTMuonHo ->GetBinContent(lTMuonHo ->GetXaxis()->FindBin(lP),lTMuonHo ->GetYaxis()->FindBin(lHO) );
  
  if(lPSX == 0. || lPBX == 0. || (lEM <= 0. && !iEMSpecial)) {lPSX = 1.; lPBX = 1.;} 
  if(lPSY == 0. || lPBY == 0. || lHad == 0.) {lPSY = 1.; lPBY = 1.;}
  if(lPSZ == 0. || lPBZ == 0. || lHO  == 0.) {lPSZ = 1.; lPBZ = 1.;}
  if((lPSX*lPSY*lPSZ+lPBX*lPBY*lPBZ) > 0.) return lPSX*lPSY*lPSZ/(lPSX*lPSY*lPSZ+lPBX*lPBY*lPBZ);
  return 0.5;
}

bool MuonTools::isGood(mithep::Muon *iMuon,selection iSelection) {
  double lVal = 0;
  switch(iSelection) {
  case AllArbitrated:
    if(iMuon->StandaloneTrk() != 0 || iMuon->GlobalTrk()!= 0)  return true;
    if(iMuon->NSegments() > 0) return true;
    return false;
    break;
  case PromptTight: 
    return promptTight(-1,iMuon);
    break;
   case TMOneStationLoose:
     return TMOneStation(iMuon,99999,999999);
     break;
  case TMOneStationTight:
    return TMOneStation(iMuon);
    break;
  case TMLastStationLoose: 
    return TMLastStation(iMuon,999999,999999);
    break;
  case TMLastStationTight: 
    return TMLastStation(iMuon);
    break;
  case TM2DCompatibilityLoose: 
    lVal             = 1.2*getSegmentCompatability(iMuon); 
    if(lVal/1.2 == 0.5) return false;
    lVal += 0.8*getCaloCompatability(iMuon,true,true);
    if(lVal > 0.7) return true;
    return false;
    break;
  case TM2DCompatibilityTight: 
    lVal             = 1.2*getSegmentCompatability(iMuon); 
    if(lVal/1.2 == 0.5) return false;
    lVal += 0.8*getCaloCompatability(iMuon,true,true);
    if(lVal > 1.0) return true;
    return false;
    break;
  default:
    return false;
  }
  return false;
}
