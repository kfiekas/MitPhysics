//--------------------------------------------------------------------------------------------------
// $Id: IsolationTools.h,v 1.1 2008/10/15 06:02:05 loizides Exp $
//
// MuonTools
//
// Isolation functions to compute various kinds of isolation.
// http://cmslxr.fnal.gov/lxr/source/RecoMuon/MuonIdentification/
// http://cmslxr.fnal.gov/lxr/source/DataFormats/MuonReco/
// Authors: P.Harris, C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTIL_MUONTOOLS_H
#define MITPHYSICS_UTIL_MUONTOOLS_H

#include "MitAna/DataTree/interface/Muon.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "TH2D.h"

namespace mithep {
  class MuonTools {
    public:
      MuonTools(const char *mutemp="$CMSSW_BASE/src/MitPhysics/data/MuonCaloTemplate.root", 
                const char *pitemp="$CMSSW_BASE/src/MitPhysics/data/PionCaloTemplate.root");
      ~MuonTools();

      enum ESelType { 
        kAllArbitrated,
        kPromptTight,
        kTMLastStationLoose,
        kTMLastStationTight,
        kTMOneStationLoose,
        kTMOneStationTight,
        kTM2DCompatibilityLoose,
        kTM2DCompatibilityTight
      };

      Bool_t      Init(const char *mutemp, const char *pitemp);
      Bool_t      IsGood(const mithep::Muon *iMuon, ESelType iSel) const;

    protected:
      void        DeleteHistos();
      Double_t    GetCaloCompatability(const mithep::Muon *iMuon, 
                                       Bool_t iEMSpecial, Bool_t iCorrectedHCAL) const;
      Double_t    GetSegmentCompatability(const mithep::Muon *iMuon)             const;
      Bool_t      Overflow(const TH2D *iHist, Double_t lVal0, Double_t lVal1)    const;
      Double_t    SigWeight(Double_t iVal0, Double_t iVal1)                      const;

    private:
      Bool_t      fIsInit;              //!
      TH2D       *fmuon_em_etaEmi;      //!
      TH2D       *fmuon_had_etaEmi;     //!
      TH2D       *fmuon_had_etaTmi;     //!
      TH2D       *fmuon_em_etaB;        //!
      TH2D       *fmuon_had_etaB;       //!
      TH2D       *fmuon_ho_etaB;        //!
      TH2D       *fmuon_had_etaTpl;     //!
      TH2D       *fmuon_em_etaEpl;      //!
      TH2D       *fmuon_had_etaEpl;     //!
      TH2D       *fpion_em_etaEmi;      //!
      TH2D       *fpion_had_etaEmi;     //!
      TH2D       *fpion_had_etaTmi;     //!
      TH2D       *fpion_em_etaB;        //!
      TH2D       *fpion_had_etaB;       //!
      TH2D       *fpion_ho_etaB;        //!
      TH2D       *fpion_had_etaTpl;     //!
      TH2D       *fpion_em_etaEpl;      //!
      TH2D       *fpion_had_etaEpl;     //!

      TH2D       *LoadHisto(const char *fname, TFile *file)                      const;
  };
}

//--------------------------------------------------------------------------------------------------
inline Double_t mithep::MuonTools::SigWeight(Double_t iVal0, Double_t iVal1) const
{
  //

  if (iVal1 < 1.) 
    return 1.;
  if (iVal0 < 3. && iVal1 > 3.) {
    Double_t lVal = TMath::Max(iVal0,1.);
    return 1./TMath::Power(lVal,0.25);
  }

  Double_t lVal = TMath::Max(iVal1,1.);
  return 1./TMath::Power(lVal,0.25);
}

//--------------------------------------------------------------------------------------------------
inline Bool_t mithep::MuonTools::Overflow(const TH2D *iHist, Double_t lVal0, Double_t lVal1) const
{
  // Check if values are in overflow bins of given histogram.

  if(iHist == 0)
    return kTRUE;

  if (iHist ->GetXaxis()->FindBin(lVal0) == 0                  || 
      iHist ->GetXaxis()->FindBin(lVal0) >  iHist->GetNbinsX() || 
      iHist ->GetYaxis()->FindBin(lVal0) == 0                  || 
      iHist ->GetYaxis()->FindBin(lVal0) >  iHist->GetNbinsY()) {
    return kTRUE;
  }
  return kFALSE;
}
#endif

#if 0
    static inline MCParticle* etaPhiMatch(MCParticleCol* iMCParts,Muon *iMuon) {
      mithep::MCParticle *lMC = 0; Double_t lDR = 1.; Double_t lPt = -1.;
      for(unsigned Int_t i0 = 0; i0 < iMCParts->GetEntries(); i0++) {
	mithep::MCParticle* pMC = iMCParts->At(i0);
	if(pMC->Status() != 1) continue;
	Double_t pDR = mithep::MathUtils::DeltaR(pMC->Mom(),iMuon->Mom());
	if(pDR > lDR && pDR > 0.01) continue;
	if(fabs(pMC->Pt()) < lPt)   continue;
	lDR = pDR;
	lMC = pMC;
	lPt = pMC->Pt();
      }
      return lMC;
    }
    static inline Muon* match(MuonCol* iMuons,Track *iTrack) {
      mithep::Muon *lMuon = 0; Double_t lDR = 1.; Double_t lPt = -1.;
      for(unsigned Int_t i0 = 0; i0 < iMuons->GetEntries(); i0++) {
	mithep::Muon* pMuon = iMuons->At(i0);
	Double_t pDR = mithep::MathUtils::DeltaR(pMuon->Mom(),iTrack->Mom4(0.109));
	if(pDR > lDR && pDR > 0.01) continue;
	if(fabs(pMuon->Pt()) < lPt)   continue;
	lDR   = pDR;
	lMuon = pMuon;
	lPt   = pMuon->Pt();
      }
      return lMuon;
    }
#endif


