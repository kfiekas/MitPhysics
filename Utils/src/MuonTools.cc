// $Id: MuonTools.cc,v 1.21 2011/10/02 10:05:42 ceballos Exp $

#include "MitPhysics/Utils/interface/MuonTools.h"
#include <TFile.h>

ClassImp(mithep::MuonTools)

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
  // Get calo compatibility value for given muon based on calorimeter templates.
  // If iEMSpecial is true, then a use different arrangement of ECAL for compatibility.

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
      lHad *= TMath::Sin(2*TMath::ATan(TMath::Exp(lEta)));
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
      return iMuon->TMOneStation(999999,999999);
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

  Double_t lVal = GetSegmentCompatability(iMuon); 
  if (lVal == 0.5) // exclude this border case
    return kFALSE;

  lVal *= 1.2;
  lVal += 0.8*GetCaloCompatability(iMuon,kTRUE,kTRUE);
  if (lVal > tm2dcut) 
    return kTRUE;

  return kFALSE;
}

//--------------------------------------------------------------------------------------------------
Double_t MuonTools::GetSegmentCompatability(const mithep::Muon *iMuon) const
{
  // Get segment compatability for given muon based on likelihood of well defined 
  // track through chambers.

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

    if(iMuon->GetDX(i0) < 999999.) {
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
//--------------------------------------------------------------------------------------------------
Bool_t MuonTools::PassD0Cut(const Muon *mu, const VertexCol *vertices, Double_t fD0Cut, Int_t nVertex) 
{
  Bool_t d0cut = kFALSE;
  const Track *mt = mu->BestTrk();
  if (!mt) return kFALSE;

  Double_t d0_real = 1e30;
  if(nVertex >= 0) d0_real = TMath::Abs(mt->D0Corrected(*vertices->At(nVertex)));
  else            {
    Double_t distVtx = 999.0;
    Int_t closestVtx = 0;
    for(UInt_t nv=0; nv<vertices->GetEntries(); nv++){
      double dz = TMath::Abs(mt->DzCorrected(*vertices->At(nv)));
      if(dz < distVtx) {
	distVtx    = dz;
        closestVtx = nv;
      }
    }
    d0_real = TMath::Abs(mt->D0Corrected(*vertices->At(closestVtx)));
  }
  if(d0_real < fD0Cut) d0cut = kTRUE;
  
  return d0cut;
}

//--------------------------------------------------------------------------------------------------
Bool_t MuonTools::PassD0Cut(const Muon *mu, const BeamSpotCol *beamspots, Double_t fD0Cut) 
{
  Bool_t d0cut = kFALSE;
  const Track *mt = mu->BestTrk();
  if (!mt) return kFALSE;

  // d0 cut
  Double_t d0_real = 99999;
  for(UInt_t i0 = 0; i0 < beamspots->GetEntries(); i0++) {
    Double_t pD0 = mt->D0Corrected(*beamspots->At(i0));
    if(TMath::Abs(pD0) < TMath::Abs(d0_real)) d0_real = TMath::Abs(pD0);
  }
  if(d0_real < fD0Cut) d0cut = kTRUE;
  
  return d0cut;
}

//--------------------------------------------------------------------------------------------------
Bool_t MuonTools::PassDZCut(const Muon *mu, const VertexCol *vertices, Double_t fDZCut, Int_t nVertex) 
{
  Bool_t dzcut = kFALSE;
  const Track *mt = mu->BestTrk();
  if (!mt) return kFALSE;

  Double_t distVtx = 999.0;
  if(nVertex >= 0) distVtx = TMath::Abs(mt->DzCorrected(*vertices->At(nVertex)));
  else {
    for(UInt_t nv=0; nv<vertices->GetEntries(); nv++){
      double dz = TMath::Abs(mt->DzCorrected(*vertices->At(nv)));
      if(dz < distVtx) {
        distVtx = dz;
      }
    }
  }

  if(distVtx < fDZCut) dzcut = kTRUE;
  
  return dzcut;
}

//--------------------------------------------------------------------------------------------------
Bool_t MuonTools::PassSoftMuonCut(const Muon *mu, const VertexCol *vertices, const Double_t fDZCut,
                                  const Bool_t applyIso) 
{
  if(mu->Pt() <= 3.0) return kFALSE;

  if(!mu->IsTrackerMuon()) return kFALSE;
  
  if(!mu->Quality().Quality(MuonQuality::TMLastStationAngTight)) return kFALSE;

  if(mu->BestTrk()->NHits() <= 10) return kFALSE;

  if(!PassD0Cut(mu, vertices, 0.2, 0)) return kFALSE;

  if(!PassDZCut(mu, vertices, fDZCut, 0)) return kFALSE;

  if(applyIso == kTRUE){
    Double_t totalIso = 1.0 * mu->IsoR03SumPt() + 
                        1.0 * mu->IsoR03EmEt() + 
                        1.0 * mu->IsoR03HadEt();
    if (totalIso < (mu->Pt()*0.10) && mu->Pt() > 20.0) return kFALSE;
  }

  return kTRUE;
}
