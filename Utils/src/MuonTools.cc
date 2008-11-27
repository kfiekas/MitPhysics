// $Id:$

#include "MitPhysics/Utils/interface/MuonTools.h"
#include <TFile.h>

using namespace mithep;

//--------------------------------------------------------------------------------------------------
MuonTools::MuonTools(const char *mutemp, const char *pitemp) :
  fIsInit(kFALSE),
  fmuon_em_etaEmi(0),
  fmuon_had_etaEmi(0),
  fmuon_had_etaTmi(0),
  fmuon_em_etaB(0),
  fmuon_had_etaB(0),
  fmuon_ho_etaB(0),
  fmuon_had_etaTpl(0),
  fmuon_em_etaEpl(0),
  fmuon_had_etaEpl(0),
  fpion_em_etaEmi(0),
  fpion_had_etaEmi(0),
  fpion_had_etaTmi(0),
  fpion_em_etaB(0),
  fpion_had_etaB(0),
  fpion_ho_etaB(0),
  fpion_had_etaTpl(0),
  fpion_em_etaEpl(0),
  fpion_had_etaEpl(0)
{
  // Constructor.

  if (mutemp && pitemp)
    Init(mutemp, pitemp);
}

//--------------------------------------------------------------------------------------------------
MuonTools::~MuonTools() 
{
  // Destructor.

  DeleteHistos();
}

//--------------------------------------------------------------------------------------------------
void MuonTools::DeleteHistos()
{
  // Delete histograms.

  if (fIsInit) {
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
    fpion_em_etaEmi  = 0;
    fpion_had_etaEmi = 0;
    fpion_had_etaTmi = 0;
    fpion_em_etaB    = 0;
    fpion_had_etaB   = 0;
    fpion_ho_etaB    = 0;
    fpion_had_etaTpl = 0;
    fpion_em_etaEpl  = 0;
    fpion_had_etaEpl = 0;
    fmuon_em_etaEmi  = 0;
    fmuon_had_etaEmi = 0;
    fmuon_had_etaTmi = 0;
    fmuon_em_etaB    = 0;
    fmuon_had_etaB   = 0;
    fmuon_ho_etaB    = 0;
    fmuon_had_etaTpl = 0;
    fmuon_em_etaEpl  = 0;
    fmuon_had_etaEpl = 0;
    fIsInit = kFALSE;
  }
}

//--------------------------------------------------------------------------------------------------
Double_t MuonTools::GetCaloCompatability(const Muon *iMuon,
                                         Bool_t iEMSpecial, Bool_t iCorrectedHCAL) const
{
  // todo

  Double_t lEta = iMuon->Eta();
  Double_t aEta = TMath::Abs(lEta);
  if (aEta > 2.5) 
    return 0.5; 

  Double_t lP = iMuon->P();
  if (lP >= 2000.) 
    lP = 1999.9;
  if(lP < 0. )           
    return 0.5; 

  Double_t lEM  = -5.;      
  if (!iEMSpecial || iMuon->EmEnergy() != 0.) 
    lEM  = iMuon->EmEnergy();

  Double_t lHad = iMuon->HadEnergy();
  Double_t lHO = iMuon->HoEnergy();;

  TH2D *lTMuonHad = 0;
  TH2D *lTPionHad = 0;
  TH2D *lTMuonHo  = 0;
  TH2D *lTPionHo  = 0;
  TH2D *lTMuonEm  = 0;
  TH2D *lTPionEm  = 0;
    
  if (aEta >= 1.27) {
  if (iCorrectedHCAL) 
      lHad *= 1.8/2.2;
    if (lEta > 0) {
      lTPionHad = fpion_had_etaEpl;
      lTMuonHad = fmuon_had_etaEpl;
    } else {
      lTPionHad = fpion_had_etaEmi;
      lTMuonHad = fmuon_had_etaEmi;
    }
  }

  if (aEta < 1.27 && aEta >= 1.1) {
    if (iCorrectedHCAL)    
      lHad *= (1.8/(-2.2*aEta+5.5));
    if (lEta > 0) {
      lTPionHad  = fpion_had_etaTpl;
      lTMuonHad  = fmuon_had_etaTpl;
    } else {
      lTPionHad  = fpion_had_etaTmi;
      lTMuonHad  = fmuon_had_etaTmi;
    }
  }

  if (aEta < 1.1) {
    if(iCorrectedHCAL)    
      lHad *= TMath::Sin(2*TMath::ATan(TMath::Exp(lEta))); //todo ask!
    lTPionHad  = fpion_had_etaB;
    lTMuonHad  = fmuon_had_etaB;
  }
  if (lEta > 1.479) {
    lTPionEm = fpion_em_etaEpl;
    lTMuonEm = fmuon_em_etaEpl;
  }
  if (aEta <= 1.479) {
    lTPionEm  = fpion_em_etaB;
    lTMuonEm  = fmuon_em_etaB;
  }
  if (lEta < -1.479) {
    lTPionEm  = fpion_em_etaEmi;
    lTMuonEm  = fmuon_em_etaEmi;
  }
  if (aEta < 1.28) {
    lTPionHo  = fpion_ho_etaB;
    lTMuonHo  = fmuon_ho_etaB;
  }
  
  Double_t lPBX = 1.;     
  Double_t lPSX = 1.; 
  Double_t lPBY = 1.;     
  Double_t lPSY = 1.; 
  Double_t lPBZ = 1.;     
  Double_t lPSZ = 1.; 
  if (!Overflow(lTPionEm, lP,lEM))  
    lPBX = lTPionEm ->GetBinContent(lTPionEm ->GetXaxis()->FindBin(lP),
                                     lTPionEm ->GetYaxis()->FindBin(lEM));
  if (!Overflow(lTPionHad,lP,lHad)) 
    lPBY = lTPionHad->GetBinContent(lTPionHad->GetXaxis()->FindBin(lP),
                                     lTPionHad->GetYaxis()->FindBin(lHad));
  if (!Overflow(lTPionHo, lP,lHO))  
    lPBZ = lTPionHo ->GetBinContent(lTPionHo ->GetXaxis()->FindBin(lP),
                                     lTPionHo ->GetYaxis()->FindBin(lHO));
  if (!Overflow(lTMuonEm, lP,lEM )) 
    lPSX = lTMuonEm ->GetBinContent(lTMuonEm ->GetXaxis()->FindBin(lP),
                                     lTMuonEm ->GetYaxis()->FindBin(lEM));
  if (!Overflow(lTMuonHad,lP,lHad)) 
    lPSY = lTMuonHad->GetBinContent(lTMuonHad->GetXaxis()->FindBin(lP),
                                    lTMuonHad->GetYaxis()->FindBin(lHad));
  if (!Overflow(lTMuonHo ,lP,lHO))  
    lPSZ =  lTMuonHo ->GetBinContent(lTMuonHo ->GetXaxis()->FindBin(lP),
                                     lTMuonHo ->GetYaxis()->FindBin(lHO));
  
  if (lPSX == 0. || lPBX == 0. || (lEM <= 0. && !iEMSpecial)) {
    lPSX = 1.; 
    lPBX = 1.;
  } 
  if (lPSY == 0. || lPBY == 0. || lHad == 0.) {
    lPSY = 1.; 
    lPBY = 1.;
  }
  if (lPSZ == 0. || lPBZ == 0. || lHO  == 0.) {
    lPSZ = 1.; 
    lPBZ = 1.;
  }
  if ((lPSX*lPSY*lPSZ+lPBX*lPBY*lPBZ) > 0.) 
    return lPSX*lPSY*lPSZ/(lPSX*lPSY*lPSZ+lPBX*lPBY*lPBZ);

  return 0.5;
}

//--------------------------------------------------------------------------------------------------
Bool_t MuonTools::Init(const char *mutemp, const char *pitemp)
{
  // Read histograms from given files.

  if (fIsInit) {
    DeleteHistos();
  }

  TDirectory::TContext context(0);

  TFile *muon_templates = TFile::Open(mutemp);
  if (!muon_templates) {
    Fatal("Init", "Could not open file %s", mutemp);
    return kFALSE;
  }
  fmuon_em_etaEmi  = LoadHisto("em_etaEmi",  muon_templates);
  fmuon_had_etaEmi = LoadHisto("had_etaEmi", muon_templates);
  fmuon_had_etaTmi = LoadHisto("had_etaTmi", muon_templates);
  fmuon_em_etaB    = LoadHisto("em_etaB",    muon_templates);
  fmuon_had_etaB   = LoadHisto("had_etaB",   muon_templates);
  fmuon_ho_etaB    = LoadHisto("ho_etaB",    muon_templates);
  fmuon_had_etaTpl = LoadHisto("had_etaTpl", muon_templates);
  fmuon_em_etaEpl  = LoadHisto("em_etaEpl",  muon_templates);
  fmuon_had_etaEpl = LoadHisto("had_etaEpl", muon_templates);
  muon_templates->Close();
  delete muon_templates;

  TFile *pion_templates = TFile::Open(pitemp);
  if (!pion_templates) {
    Fatal("Init", "Could not open file %s", pitemp);
    return kFALSE;
  }

  fpion_em_etaEmi  = LoadHisto("em_etaEmi",  pion_templates);
  fpion_had_etaEmi = LoadHisto("had_etaEmi", pion_templates);
  fpion_had_etaTmi = LoadHisto("had_etaTmi", pion_templates);
  fpion_em_etaB    = LoadHisto("em_etaB",    pion_templates);
  fpion_had_etaB   = LoadHisto("had_etaB",   pion_templates);
  fpion_ho_etaB    = LoadHisto("ho_etaB",    pion_templates);
  fpion_had_etaTpl = LoadHisto("had_etaTpl", pion_templates);
  fpion_em_etaEpl  = LoadHisto("em_etaEpl",  pion_templates);
  fpion_had_etaEpl = LoadHisto("had_etaEpl", pion_templates);
  pion_templates->Close();
  delete pion_templates;

  fIsInit = kTRUE;
  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t MuonTools::IsGood(const mithep::Muon *iMuon, ESelType iSel) const
{
  // Return true if given muon qualifies given selection criterium.

  Double_t tm2dcut = 0.;

  switch(iSel) {
    case kAllArbitrated:
      if (iMuon->StandaloneTrk() != 0 || iMuon->GlobalTrk()!= 0)  
        return kTRUE;
      if (iMuon->NSegments() > 0) 
        return kTRUE;
      break;
    case kPromptTight: 
      return iMuon->PromptTight(Muon::kAny);
      break;
    case kTMOneStationLoose:
      return iMuon->TMOneStation(99999,999999);
      break;
    case kTMOneStationTight:
      return iMuon->TMOneStation();
      break;
    case kTMLastStationLoose: 
      return iMuon->TMLastStation(999999,999999);
      break;
    case kTMLastStationTight: 
      return iMuon->TMLastStation();
      break;
    case kTM2DCompatibilityLoose: 
      tm2dcut = 0.7;
      break;
    case kTM2DCompatibilityTight: 
      tm2dcut = 1.0;
      break;
    default:
      return kFALSE;
      break;
  }

  Double_t lVal = 1.2*GetSegmentCompatability(iMuon); 
  if (lVal/1.2 == 0.5) 
    return kFALSE;

  lVal += 0.8*GetCaloCompatability(iMuon,kTRUE,kTRUE);
  if (lVal > tm2dcut) 
    return kTRUE;
  return kFALSE;
}

//--------------------------------------------------------------------------------------------------
Double_t MuonTools::GetSegmentCompatability(const mithep::Muon *iMuon) const
{
  // todo phil

  Int_t lNStationsCrossed = 0;
  Int_t lNStationsSegment = 0;
  
  Int_t lStSegmentmatch[8];
  Int_t lStCrossed[8];
  Double_t lStBoundary[8];

  Double_t lWeight  = 0.;
  Bool_t lAdjust    = kTRUE;
  for (Int_t i0 = 0; i0 < 8; ++i0) {
    lStBoundary[i0] = 0.;
    if(iMuon->GetTrackDist(i0) < 999999. ) { 
      lNStationsCrossed++;
      lStCrossed[i0]  = 1;
      if (iMuon->GetTrackDist(i0) > -10. ) 
        lStBoundary[i0] = iMuon->GetTrackDist(i0); 
    } else
      lStCrossed[i0]  = 0;

    if(iMuon->GetDX(i0) < 999999.) { //Use iMuon->GetSegmentX--> CHECK
      lNStationsSegment++;
      lStSegmentmatch[i0] = 1;
    } else
      lStSegmentmatch[i0] = 0;

    if(iMuon->GetDY(i0) < 999999.) 
      lAdjust = kFALSE;
  }

  if (lNStationsCrossed == 0)
    return 0.5;

  Double_t lStWeight[8];
  Int_t lPCross = -1;
  const Double_t lAtWeight = 0.5; 
  for (Int_t i0 = 0; i0< 8; ++i0) { 
    lStWeight[i0] = 0;
    if (lStCrossed[i0] > 0) { 
      lPCross++;

      switch (lNStationsCrossed) { 
        case 1 : 
          lStWeight[i0] =  1.;
          break;
        case 2 :
          if (lPCross == 0 ) 
            lStWeight[i0] = 0.33;
          else 
            lStWeight[i0] = 0.67;
          break;
        case 3 : 
          if (lPCross == 0) 
            lStWeight[i0] = 0.23;
          else if (lPCross == 1) 
            lStWeight[i0] = 0.33;
          else 
            lStWeight[i0] = 0.44;
          break;
        case 4 : 
          if (lPCross == 0) 
            lStWeight[i0] = 0.10;
          else if (lPCross == 1) 
            lStWeight[i0] = 0.20;
          else if (lPCross == 2) 
            lStWeight[i0] = 0.30;
          else                    
            lStWeight[i0] = 0.40;
          break;
        default : 
          lStWeight[i0] = 1./lNStationsCrossed;
      }
    
      if (lStSegmentmatch[i0] <= 0 && lStBoundary[i0] != 0.)
        lStWeight[i0] *= lAtWeight*0.5*(TMath::Erf(lStBoundary[i0]/6.)+1.); 
      else if (lStSegmentmatch[i0] <= 0 && lStBoundary[i0] == 0) 
        lStWeight[i0] = 0.;
      
      if (lStSegmentmatch[i0] > 0) { 
        Double_t lP2X = TMath::Power(iMuon->GetPullX(i0),2.);
        Double_t lP2Y = TMath::Power(iMuon->GetPullY(i0),2.);
        Double_t lD2X = TMath::Power(iMuon->GetDX(i0),2.);
        Double_t lD2Y = TMath::Power(iMuon->GetDY(i0),2.);
        if (iMuon->GetDY(i0) < 999999 && iMuon->GetDX(i0) < 999999) 
          lStWeight[i0] *= SigWeight(TMath::Sqrt(lD2X+lD2Y),TMath::Sqrt(lP2X+lP2Y));
        else if (iMuon->GetDY(i0) >= 999999 && i0 < 4)
          lStWeight[i0] *= SigWeight(iMuon->GetDX(i0),iMuon->GetPullX(i0));
        else if(i0 < 4)
          lStWeight[i0] *= SigWeight(iMuon->GetDY(i0),iMuon->GetPullY(i0));
      }
    } 
    lWeight += lStWeight[i0];
  }

  return lWeight;
}

//--------------------------------------------------------------------------------------------------
TH2D *MuonTools::LoadHisto(const char *name, TFile *file) const
{
  // Load histogram with given name from given file and return it.

  TH2D *ret = dynamic_cast<TH2D*>(file->Get(name));
  if (!ret) {
    Fatal("LoadHisto", "Could not load histogram %s from file %s", name, file->GetName());
    return 0;
  }
  ret->SetDirectory(0);
  return ret;
}

