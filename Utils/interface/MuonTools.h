#ifndef MITPHYSICS_UTIL_MUONTOOLS_H
#define MITPHYSICS_UTIL_MUONTOOLS_H
#include <vector>
#include <iostream>

#include "TH2D.h"
#include "TFile.h"
#include "MitAna/DataTree/interface/Collections.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

namespace mithep {
  class MuonTools {
  public:
    MuonTools();
    ~MuonTools();
    double getCaloCompatability(mithep::Muon* iMuon,bool iEMSpecial, bool iCorrectedHCAL);
    enum selection { 
      AllArbitrated,
      PromptTight,
      TMLastStationLoose,
      TMLastStationTight,
      TMOneStationLoose,
      TMOneStationTight,
      TM2DCompatibilityLoose,
      TM2DCompatibilityTight
    };
    static inline double sigWeight(double iVal0,double iVal1) {
      if(iVal1 < 1.) return 1;
      if(iVal0 < 3. && iVal1 > 3.) {
	double lVal = TMath::Max(iVal0,1.);
	return 1/TMath::Power(lVal,0.25);
      }
      double lVal = TMath::Max(iVal1,1.);
      return 1/TMath::Power(lVal,0.25);
    }
    static inline int lastHit(mithep::Muon *iMuon) {
      int lId = -1;
      for(int i0 = 0; i0 < 8; i0++) {
	if(iMuon->GetDX(i0) < 99999 || iMuon->GetDY(i0) < 99999) lId = i0;
      }
      return lId;
    }
    static inline int maxChamberId(mithep::Muon* iMuon,double iMaxD,double iMaxP) {
      int lId = -1;
      for(int i0 = 0; i0 < 8; i0++) {
	if(iMuon->GetTrackDist(i0)                            < iMaxD &&
	   iMuon->GetTrackDist(i0)/iMuon->GetTrackDistErr(i0) < iMaxP) {
	  if(iMuon->GetDX(i0) >= 999999) {lId = i0;} else {lId = -1;}
	}
      }
      return lId;
    }
    static inline int lastStation(mithep::Muon* iMuon,double iMaxD,double iMaxP) {
      int lId = -1;
      for(int i0 = 0; i0 < 8; i0++) {
	if((lId % 4) > (i0 % 4)) continue;
	if(iMuon->GetTrackDist(i0)                            < iMaxD &&
	   iMuon->GetTrackDist(i0)/iMuon->GetTrackDistErr(i0) < iMaxP) lId = i0;
      }
      return lId;
    }
    static inline int lastStation(mithep::Muon* iMuon,int iMax=8) {
      int lId = -1; if(iMax > 8) iMax = 8;
      for(int i0 = 0; i0 < iMax; i0++) {if(iMuon->StationBit(i0) && ((lId % 4) < (i0 % 4)))lId = i0;}
      return lId;
    }
    static inline bool overflow(TH2D *iHist,double lVal0,double lVal1) {
      if(iHist == 0) return true;
      if(iHist ->GetXaxis()->FindBin(lVal0) == 0               || 
	 iHist ->GetXaxis()->FindBin(lVal0) >  iHist->GetNbinsX() || 
	 iHist ->GetYaxis()->FindBin(lVal0) == 0               || 
	 iHist ->GetYaxis()->FindBin(lVal0) >  iHist->GetNbinsY()) return true;
      return false;
    }
    static inline MCParticle* etaPhiMatch(MCParticleCol* iMCParts,Muon *iMuon) {
      mithep::MCParticle *lMC = 0; double lDR = 1.; double lPt = -1.;
      for(unsigned int i0 = 0; i0 < iMCParts->GetEntries(); i0++) {
	mithep::MCParticle* pMC = iMCParts->At(i0);
	if(pMC->Status() != 1) continue;
	double pDR = mithep::MathUtils::DeltaR(pMC->Mom(),iMuon->Mom());
	if(pDR > lDR && pDR > 0.01) continue;
	if(fabs(pMC->Pt()) < lPt)   continue;
	lDR = pDR;
	lMC = pMC;
	lPt = pMC->Pt();
      }
      return lMC;
    }
    static inline Muon* match(MuonCol* iMuons,Track *iTrack) {
      mithep::Muon *lMuon = 0; double lDR = 1.; double lPt = -1.;
      for(unsigned int i0 = 0; i0 < iMuons->GetEntries(); i0++) {
	mithep::Muon* pMuon = iMuons->At(i0);
	double pDR = mithep::MathUtils::DeltaR(pMuon->Mom(),iTrack->Mom4(0.109));
	if(pDR > lDR && pDR > 0.01) continue;
	if(fabs(pMuon->Pt()) < lPt)   continue;
	lDR   = pDR;
	lMuon = pMuon;
	lPt   = pMuon->Pt();
      }
      return lMuon;
    }
    static inline int NSegments(mithep::Muon *iMuon) {
      int lNSegs = 0;
      for(int i0 = 0; i0 < 8; i0++) {
	if(iMuon->StationBit(i0)) lNSegs++;
      }
      return lNSegs;
    }
    static inline bool promptTight(int iId,Muon *iMuon) {
      const mithep::Track *lTrack = 0;
      if(iId == 0                              ) lTrack = iMuon->GlobalTrk();
      if(iId == 1 || (iId == -1 && lTrack == 0)) lTrack = iMuon->TrackerTrk();
      if(iId == 2 || (iId == -1 && lTrack == 0)) lTrack = iMuon->StandaloneTrk();
      if(lTrack == 0) return 0;
      if(lTrack->NHits() < 11)                                        return false;
      if(lTrack->Chi2()/lTrack->Ndof() > 10.)                         return false;
      if(lTrack->D0() > 0.2)                                          return false;
      if((lastHit(iMuon)% 4) == 0)                                    return false;
      return true;
    }
    static inline bool TMLastStation(Muon *iMuon,double iDYMin = 3., double iPYMin = 3.,double iDXMin = 3., double iPXMin = 3.,int iN = 2) {
      if(iMuon->NSegments()       < iN)               return false; //2Last 1One
      int lLast = lastStation(iMuon,-3.,-3.);                       //Last Required Station
      if(lLast < 0)                                   return false;
      if(iMuon->GetDX(lLast)  > 9999.)                return false;
      lLast = lastStation(iMuon);                                   //No Requirements Imply StationMask (Segment and Track Arbitrated) 
      if(lLast < 0)                                   return false;
      if(!(fabs(iMuon->GetDX(lLast))      < iDXMin  ||
	   fabs(iMuon->GetPullX(lLast))   < iPXMin))  return false;
      if(lLast == 3) lLast = lastStation(iMuon,3);
      if(lLast < 0)                                   return false;
      if(!(fabs(iMuon->GetDY(lLast))      < iDYMin ||               //Old Code it was 9999
	   fabs(iMuon->GetPullY(lLast))   < iPYMin))  return false;
      return true;
    }
    static inline bool TMOneStation(Muon *iMuon,double iDYMin = 3., double iPYMin = 3.,double iDXMin = 3., double iPXMin = 3.,int iN = 1) {
      if(iMuon->NSegments()       < iN)               return false; //2Last 1One
      bool pGoodX  = false; bool pBadY = false;
      for(int i0 = 0; i0 < 8; i0++) {
	if((fabs(iMuon->GetDX(i0))      < iDXMin  ||
	    fabs(iMuon->GetPullX(i0))   < iPXMin)) pGoodX = true; 
	if(pGoodX                                 &&
	   (fabs(iMuon->GetDY(i0))      < iDYMin  || 
	    fabs(iMuon->GetPullY(i0))   < iPYMin))  return true;
        if(fabs(iMuon->GetDY(i0)) < 999999) pBadY = true;
	if(i0 == 3 && pGoodX && !pBadY)             return true;
      }
      return false;
    }
    static inline double getCaloCompatabilitySlowly(mithep::Muon* iMuon,bool iEMSpecial, bool iCorrectedHCAL) {
      std::cout << "Warning Loading Compatability Slowly" << std::endl;
      TFile* lPion_templates = new TFile("MitCommon/MuonTools/data/PionCaloTemplate.root","READ");
      TFile* lMuon_templates = new TFile("MitCommon/MuonTools/data/MuonCaloTemplate.root","READ");
      TH2D* lpion_em_etaEmi  = (TH2D*) lPion_templates->Get("em_etaEmi");
      TH2D* lpion_had_etaEmi = (TH2D*) lPion_templates->Get("had_etaEmi");
      TH2D* lpion_had_etaTmi = (TH2D*) lPion_templates->Get("had_etaTmi");
      TH2D* lpion_em_etaB    = (TH2D*) lPion_templates->Get("em_etaB")   ;
      TH2D* lpion_had_etaB   = (TH2D*) lPion_templates->Get("had_etaB")  ;
      TH2D* lpion_ho_etaB    = (TH2D*) lPion_templates->Get("ho_etaB")   ;
      TH2D* lpion_had_etaTpl = (TH2D*) lPion_templates->Get("had_etaTpl");
      TH2D* lpion_em_etaEpl  = (TH2D*) lPion_templates->Get("em_etaEpl") ;
      TH2D* lpion_had_etaEpl = (TH2D*) lPion_templates->Get("had_etaEpl");
      TH2D* lmuon_em_etaEmi  = (TH2D*) lMuon_templates->Get("em_etaEmi") ;
      TH2D* lmuon_had_etaEmi = (TH2D*) lMuon_templates->Get("had_etaEmi");
      TH2D* lmuon_had_etaTmi = (TH2D*) lMuon_templates->Get("had_etaTmi");
      TH2D* lmuon_em_etaB    = (TH2D*) lMuon_templates->Get("em_etaB")   ;
      TH2D* lmuon_had_etaB   = (TH2D*) lMuon_templates->Get("had_etaB")  ;
      TH2D* lmuon_ho_etaB    = (TH2D*) lMuon_templates->Get("ho_etaB")   ;
      TH2D* lmuon_had_etaTpl = (TH2D*) lMuon_templates->Get("had_etaTpl");
      TH2D* lmuon_em_etaEpl  = (TH2D*) lMuon_templates->Get("em_etaEpl");
      TH2D* lmuon_had_etaEpl = (TH2D*) lMuon_templates->Get("had_etaEpl");

      double lEta = -1.; double lP = -1;
      double lEM  = -5.;      double lHad = 0;      double lHO = 0;
      lEta = iMuon->Eta();
      lP   = iMuon->P(); 
      if(lP >= 2000.) lP = 1999.9;
      if(!iEMSpecial || iMuon->EmEnergy() != 0.) lEM  = iMuon->EmEnergy();
      lHad = iMuon->HadEnergy();
      lHO  = iMuon->HoEnergy();
      if(lP < 0. )	     return 0.5; 
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
          lTPionHad = lpion_had_etaEpl;
          lTMuonHad = lmuon_had_etaEpl;
        } else {
          lTPionHad = lpion_had_etaEmi;
          lTMuonHad = lmuon_had_etaEmi;
        }
      }
      if(fabs(lEta) <  1.27  && fabs(lEta) >=  1.1 ) {
        if(iCorrectedHCAL)    lHad *= (1.8/(-2.2*fabs(lEta)+5.5));
        if(lEta > 0) {
          lTPionHad  = lpion_had_etaTpl;
          lTMuonHad  = lmuon_had_etaTpl;
        } else {
          lTPionHad  = lpion_had_etaTmi;
          lTMuonHad  = lmuon_had_etaTmi;
        }
      }
      if(fabs(lEta) <  1.1) {
        if(iCorrectedHCAL)    lHad *= sin(2*atan(exp(iMuon->Eta())));
        lTPionHad  = lpion_had_etaB;
        lTMuonHad  = lmuon_had_etaB;
      }
      if(lEta >  1.479  ) {
        lTPionEm  = lpion_em_etaEpl;
        lTMuonEm  = lmuon_em_etaEpl;
      }
      if(fabs(lEta) <=  1.479) {
        lTPionEm  = lpion_em_etaB;
        lTMuonEm  = lmuon_em_etaB;
      }
      if(lEta < -1.479 ) {
        lTPionEm  = lpion_em_etaEmi;
        lTMuonEm  = lmuon_em_etaEmi;
      }
      if(fabs(lEta) < 1.28) {
        lTPionHo  = lpion_ho_etaB;
        lTMuonHo  = lmuon_ho_etaB;
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
      lPion_templates->Close();
      lMuon_templates->Close();
      if((lPSX*lPSY*lPSZ+lPBX*lPBY*lPBZ) > 0.) return lPSX*lPSY*lPSZ/(lPSX*lPSY*lPSZ+lPBX*lPBY*lPBZ);
      return 0.5;
    }
    static inline float getSegmentCompatability(const mithep::Muon* iMuon) {
      int lNStationsCrossed = 0;
      int lNStationsSegment = 0;
      
      std::vector<int>   lStTrack(8);
      std::vector<int>   lStSegmentmatch(8);
      std::vector<int>   lStCrossed(8);
      std::vector<float> lStBoundary(8);
      std::vector<float> lStWeight(8);
      float lWeight    = 0.;
      bool  lAdjust    = true;
      for(int i0 = 0; i0 < 8; i0++) {
	if(iMuon->GetTrackDist(i0) < 999999. ) { 
	  lNStationsCrossed++;
	  lStCrossed[i0]  = 1;
	  lStBoundary[i0] = 0.;
	  if(iMuon->GetTrackDist(i0) > -10. ) lStBoundary[i0] = iMuon->GetTrackDist(i0); 
	}
	if(iMuon->GetDX(i0) < 999999.) { //Use iMuon->GetSegmentX--> CHECK
	  lNStationsSegment++;
	  lStSegmentmatch[i0] = 1;
	}
	if(iMuon->GetDY(i0) < 999999.) lAdjust = false;
      }
      int lPCross = -1;
      const float lAtWeight = 0.5; 
      for(int i0 = 0; i0< 8; i0++)     { 
	lStWeight[i0] = 0;
	if(lStCrossed[i0] > 0 )        { lPCross++;
	switch ( lNStationsCrossed ) { 
	case 1 : 
	  lStWeight[i0] =  1.;
	  break;
	case 2 :
	  if     ( lPCross  == 0 ) lStWeight[i0] =  0.33;
	  else   lStWeight[i0] =  0.67;
	  break;
	case 3 : 
	  if     ( lPCross == 0 ) lStWeight[i0] =  0.23;
	  else if( lPCross == 1 ) lStWeight[i0] =  0.33;
	  else                    lStWeight[i0] =  0.44;
	  break;
	case 4 : 
	  if     ( lPCross == 0 ) lStWeight[i0] =  0.10;
	  else if( lPCross == 1 ) lStWeight[i0] =  0.20;
	  else if( lPCross == 2 ) lStWeight[i0] =  0.30;
	  else                    lStWeight[i0] =  0.40;
	  break;
          
	default : 
	  lStWeight[i0] = 1./lNStationsCrossed;
	}
	
	if(lStSegmentmatch[i0] <= 0 && lStBoundary[i0] != 0. ) {
	  lStWeight[i0] *= lAtWeight*0.5*(TMath::Erf(lStBoundary[i0]/6.)+1.); 
	} else if(lStSegmentmatch[i0] <= 0 && lStBoundary[i0] == 0) {lStWeight[i0] = 0.;}
	
	if( lStSegmentmatch[i0] > 0) { 
	  double lP2X = TMath::Power(iMuon->GetPullX(i0),2.);
	  double lP2Y = TMath::Power(iMuon->GetPullY(i0),2.);
	  double lD2X = TMath::Power(iMuon->GetDX(i0),2.);
	  double lD2Y = TMath::Power(iMuon->GetDY(i0),2.);
	  if( iMuon->GetDY(i0) < 999999 && iMuon->GetDX(i0) < 999999) { 
	    lStWeight[i0] *= sigWeight(sqrt(lD2X+lD2Y),sqrt(lP2X+lP2Y));
	  } else if (iMuon->GetDY(i0) >= 999999 && i0 < 4) { 
	    lStWeight[i0] *= sigWeight(iMuon->GetDX(i0),iMuon->GetPullX(i0));
	  } else if(i0 < 4) {
	    lStWeight[i0] *= sigWeight(iMuon->GetDY(i0),iMuon->GetPullY(i0));
	  }
	}
	} 
	lWeight += lStWeight[i0];
      }
      if(lNStationsCrossed == 0) lWeight = 0.5;
      return lWeight;
    }
    static inline bool isGoodMuon(mithep::Muon *iMuon,selection iSelection) {
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
      case TMOneStationTight:
        return TMOneStation(iMuon);
      case TMLastStationLoose: 
	return TMLastStation(iMuon,999999,999999);
	break;
      case TMLastStationTight: 
	return TMLastStation(iMuon);
	break;
      case TM2DCompatibilityLoose: 
	lVal             = 1.2*getSegmentCompatability(iMuon); 
	if(lVal/1.2 == 0.5) return false;
	//if(fClass ) lVal += 0.8*getCaloCompatability(iMuon,true,true);
	lVal += 0.8*getCaloCompatabilitySlowly(iMuon,true,true);
	if(lVal > 0.7) return true;
	return false;
	break;
      case TM2DCompatibilityTight: 
	lVal             = 1.2*getSegmentCompatability(iMuon); 
	if(lVal/1.2 == 0.5) return false;
	//if(fClass ) lVal += 0.8*getCaloCompatability(iMuon,true,true);
	lVal += 0.8*getCaloCompatabilitySlowly(iMuon,true,true);
	if(lVal > 1.0) return true;
	return false;
	break;
      default:
	return false;
      }
      return false;
    }
    bool isGood(mithep::Muon *iMuon,selection iSelection);
  private:
    TFile *fMuon_templates;
    TFile *fPion_templates;
    TH2D* fpion_em_etaEmi;
    TH2D* fpion_had_etaEmi;
    TH2D* fpion_had_etaTmi;
    TH2D* fpion_em_etaB;
    TH2D* fpion_had_etaB;
    TH2D* fpion_ho_etaB;
    TH2D* fpion_had_etaTpl;
    TH2D* fpion_em_etaEpl;
    TH2D* fpion_had_etaEpl;
    TH2D* fmuon_em_etaEmi;
    TH2D* fmuon_had_etaEmi;
    TH2D* fmuon_had_etaTmi;
    TH2D* fmuon_em_etaB;
    TH2D* fmuon_had_etaB;
    TH2D* fmuon_ho_etaB;
    TH2D* fmuon_had_etaTpl;
    TH2D* fmuon_em_etaEpl;
    TH2D* fmuon_had_etaEpl;
  };
}

#endif
