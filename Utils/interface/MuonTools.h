//--------------------------------------------------------------------------------------------------
// $Id: MuonTools.h,v 1.7 2009/04/07 15:37:09 loizides Exp $
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

#ifndef MITPHYSICS_UTILS_MUONTOOLS_H
#define MITPHYSICS_UTILS_MUONTOOLS_H

#include "MitAna/DataTree/interface/Muon.h"
#include "MitAna/DataTree/interface/Collections.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "TH2D.h"

namespace mithep {
  class MuonTools {
    public:
      MuonTools(const char *mutemp="$CMSSW_BASE/src/MitPhysics/data/MuonCaloTemplate.root", 
                const char *pitemp="$CMSSW_BASE/src/MitPhysics/data/PionCaloTemplate.root");
      ~MuonTools();

      enum ESelType { 
        kAllArbitrated,          //All arbitration (DT/CSC/RPC Hits) put on at least one 
                                 //  segments given a global Muon
        kPromptTight,            //Standard global muon identification
        kTMLastStationLoose,     //Loose matching requirements on lastmost muon station of reco 
        kTMLastStationTight,     //Tight matching requirements on lastmost muon station of reco 
        kTMOneStationLoose,      //Loose matching requirements on at least one muon station of reco 
        kTMOneStationTight,      //Tight matching requirements on at least one muon station of reco 
        kTM2DCompatibilityLoose, //Loose requirement on sum of compatabiliity variables 
                                 //  ===> 1.2 Segment compatability + 0.8 calo compatability > 0.8
        kTM2DCompatibilityTight  //Tight requirement on sum of compatabiliity variables 
                                 //  ===> 1.2 Segment compatability + 0.8 calo compatability > 1.2
      };

      Bool_t      Init(const char *mutemp, const char *pitemp);
      Bool_t      IsGood(const mithep::Muon *iMuon, ESelType iSel) const;
      Double_t    GetCaloCompatability(const mithep::Muon *iMuon,
                                       Bool_t iEMSpecial, Bool_t iCorrectedHCAL) const; 
      Double_t    GetSegmentCompatability(const mithep::Muon *iMuon)             const;

    protected:
      void        DeleteHistos();
      Bool_t      Overflow(const TH2D *iHist, Double_t lVal0, Double_t lVal1)    const; 
      Double_t    SigWeight(Double_t iVal0, Double_t iVal1)                      const; 
   
    private:
      Bool_t      fIsInit;              //!=true if histograms are loaded
      TH2D       *fmuon_em_etaEmi;      //!Neg Endcap EM       Calo Deposit Template for Muons
      TH2D       *fmuon_had_etaEmi;     //!Neg Endcap Hadronic Calo Deposit Template for Muons
      TH2D       *fmuon_had_etaTmi;     //!Neg Transition Hadronic Calo Deposit Template for Muons
      TH2D       *fmuon_em_etaB;        //!Barrel EM       Calo Deposit Template for Muons
      TH2D       *fmuon_had_etaB;       //!Barrel Hadronic Calo Deposit Template for Muons
      TH2D       *fmuon_ho_etaB;        //!Barrel HO       Calo Deposit Template for Muons
      TH2D       *fmuon_had_etaTpl;     //!Plus Transition Hadronic Calo Deposit Template for Muons
      TH2D       *fmuon_em_etaEpl;      //!Plus Endcap EM       Calo Deposit Template for Muons
      TH2D       *fmuon_had_etaEpl;     //!Plus Endcap Hadronic Calo Deposit Template for Muons
      TH2D       *fpion_em_etaEmi;      //!Neg  Endcap EM       Calo Deposit Template for Pions
      TH2D       *fpion_had_etaEmi;     //!Neg  Endcap Hadronic Calo Deposit Template for Pions
      TH2D       *fpion_had_etaTmi;     //!Neg Transition Hadronic Calo Deposit Template for Pions
      TH2D       *fpion_em_etaB;        //!Barrel EM       Calo Deposit Template for Pions
      TH2D       *fpion_had_etaB;       //!Barrel Hadronic Calo Deposit Template for Pions
      TH2D       *fpion_ho_etaB;        //!Barrel HO       Calo Deposit Template for Pions
      TH2D       *fpion_had_etaTpl;     //!Plus Transition Hadronic Calo Deposit Template for Pions
      TH2D       *fpion_em_etaEpl;      //!Plus Endcap EM       Calo Deposit Template for Pions
      TH2D       *fpion_had_etaEpl;     //!Plus Endcap Hadronic Calo Deposit Template for Pions

      TH2D       *LoadHisto(const char *fname, TFile *file)                      const;
  };
}

//--------------------------------------------------------------------------------------------------
inline Double_t mithep::MuonTools::SigWeight(Double_t iVal0, Double_t iVal1) const
{
  // Returns weighted uncertainty given segment matching uncertainty (iVal0) and
  // segment matching pull (iVal1).

  if (iVal1 < 1.)  //if pull defined and within range
    return 1.;
  if (iVal0 < 3. && iVal1 > 3.) {       //if pull not well defined and uncertainty defined
    Double_t lVal = TMath::Max(iVal0,1.);
    return 1./TMath::Power(lVal,0.25);
  }

  Double_t lVal = TMath::Max(iVal1,1.); //if pull > 1 and pull < 3 return 1/pull^4
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
