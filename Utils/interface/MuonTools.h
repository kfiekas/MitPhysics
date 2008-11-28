//--------------------------------------------------------------------------------------------------
// $Id: MuonTools.h,v 1.3 2008/11/27 16:28:58 loizides Exp $
//
// MuonTools
//
// This class allows you to classify a given muon according to defined criteria. For this purpose
// is loads histograms from two ROOT files, specified in the constructor. The main function then
// is "IsGood(*muon, selection)" which returns true if the given muon fulfills the selection
// criteria. 
// 
// Logically, the code has been put together by Phil who took most of the ideas from
//  http://cmslxr.fnal.gov/lxr/source/RecoMuon/MuonIdentification/
//  http://cmslxr.fnal.gov/lxr/source/DataFormats/MuonReco/
//
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
        kAllArbitrated,          // todo by phil
        kPromptTight,            //
        kTMLastStationLoose,     //
        kTMLastStationTight,     //
        kTMOneStationLoose,      //
        kTMOneStationTight,      //
        kTM2DCompatibilityLoose, //
        kTM2DCompatibilityTight  //
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
      Bool_t      fIsInit;              //!=true if histograms are loaded
      TH2D       *fmuon_em_etaEmi;      //!todo by phil
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
  // todo

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
