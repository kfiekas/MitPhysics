//--------------------------------------------------------------------------------------------------
// $Id: FakeRate.h,v 1.1 2009/06/30 10:47:17 loizides Exp $
//
// FakeRate
//
// Class for storing the fake rates.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_FAKEMODS_FAKERATE_H
#define MITPHYSICS_FAKEMODS_FAKERATE_H
 
#include "MitAna/DataTree/interface/DataObject.h"

class TH1F;
class TH2F;
class TF1;
class TF2;

namespace mithep 
{
  class FakeRate : public DataObject
  {
    public:
      FakeRate() {}
      FakeRate(TString eleFRFile, TString muonFRFile, TString eleFRFunctionName, 
               TString muonFRFunctionName, TString eleFRHistName, TString muonFRHistName, 
               Bool_t use2DFakeRate, 
               Bool_t useFitFunction) : fElectronFRFilename(eleFRFile), 
                                        fMuonFRFilename(muonFRFile),
                                        fElectronFRFunctionName(eleFRFunctionName),
                                        fMuonFRFunctionName(muonFRFunctionName),
                                        fElectronFRHistName(eleFRHistName),
                                        fMuonFRHistName(muonFRHistName),
                                        fUse2DFakeRate (use2DFakeRate),
                                        fUseFitFunction(useFitFunction) {        
        fIsInit = Init();
      }
      
      Bool_t       Init();
      Double_t     ElectronFakeRate(Double_t et, Double_t eta, Double_t phi);
      Double_t     MuonFakeRate(Double_t pt, Double_t eta, Double_t phi);
      Double_t     ElectronFakeRateError(Double_t et, Double_t eta, Double_t phi);
      Double_t     MuonFakeRateError(Double_t pt, Double_t eta, Double_t phi);
      const char  *GetElectronFRFilename()                     { return fElectronFRFilename;     }
      const char  *GetMuonFRFilename()                         { return fMuonFRFilename;         }
      const char  *GetElectronFRFunctionName()                 { return fElectronFRFunctionName; }
      const char  *GetMuonFRFunctionName()                     { return fMuonFRFunctionName;     }
      const char  *GetElectronFRHistName()                     { return fElectronFRHistName;     }
      const char  *GetMuonFRHistName()                         { return fMuonFRHistName;         }
      Bool_t       GetUse2DFakeRate()                          { return fUse2DFakeRate;          }
      Bool_t       GetUseFitFunction()                         { return fUseFitFunction;         }

      void         SetElectronFRFilename(const char *name)     { fElectronFRFilename     = name; }
      void         SetMuonFRFilename(const char *name)         { fMuonFRFilename         = name; }
      void         SetElectronFRFunctionName(const char *name) { fElectronFRFunctionName = name; }
      void         SetMuonFRFunctionName(const char *name)     { fMuonFRFunctionName     = name; }
      void         SetElectronFRHistName(const char *name)     { fElectronFRHistName     = name; }
      void         SetMuonFRHistName(const char *name)         { fMuonFRHistName         = name; }
      void         SetUse2DFakeRate(Bool_t b)                  { fUse2DFakeRate          = b;    }
      void         SetUseFitFunction(Bool_t b)                 { fUseFitFunction         = b;    }

    protected:
      void         DeleteHistos();
      TString      fElectronFRFilename;       //!filename of file containing electron fake rate
      TString      fMuonFRFilename;           //!filename of file containing muon fake rate
      TString      fElectronFRFunctionName;   //!name of electron fake rate function
      TString      fMuonFRFunctionName;       //!name of muon fake rate function
      TString      fElectronFRHistName;       //!name of histogram containing electron fake rate
      TString      fMuonFRHistName;           //!name of histogram containing muon fake rate
      Bool_t       fUse2DFakeRate;            //!whether to use 2D pt-eta fake rate
      Bool_t       fUseFitFunction;           //!whether to use fit function or not

    private:
      Bool_t       fIsInit;                       //!=true if histograms are loaded
      TH2F        *fElectronFakeRateHist_PtEta;   //!2D Fake Rate for electrons
      TH2F        *fMuonFakeRateHist_PtEta;       //!2D Fake Rate for muons
      TH1F        *fElectronFakeRateHist_Pt;      //!2D Fake Rate for electrons
      TH1F        *fMuonFakeRateHist_Pt;          //!2D Fake Rate for electrons
      TF2         *fElectronFakeRateFit_PtEta;    //!2D Fake Rate Fit for electrons
      TF2         *fMuonFakeRateFit_PtEta;        //!2D Fake Rate Fit for muons
      TF1         *fElectronFakeRateFit_Pt;       //!2D Fake Rate Fit for electrons
      TF1         *fMuonFakeRateFit_Pt;           //!2D Fake Rate Fit for electrons

      TH2F        *fElectronFakeRateHist_PtEta_sysError; //!2D Fake Rate Error for electrons
      TH2F        *fMuonFakeRateHist_PtEta_sysError;     //!2D Fake Rate Error for muons


    ClassDef(FakeRate, 1) // Charged particle class
  };
}
#endif
