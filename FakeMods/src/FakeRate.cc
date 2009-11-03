// $Id: FakeRate.cc,v 1.5 2009/08/11 11:19:40 phedex Exp $

#include "MitPhysics/FakeMods/interface/FakeRate.h"
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include <TMath.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TF2.h>

using namespace mithep;
using namespace std;

ClassImp(mithep::FakeRate)

//--------------------------------------------------------------------------------------------------
Bool_t FakeRate::Init()
{
  // Load all fake rate histogram files.

  // get the root file which is storing the fake rates
  TFile *ElectronFRFile = new TFile(fElectronFRFilename);  
  if (!ElectronFRFile) {
    cerr << "The Electron FakeRate file : " << fElectronFRFilename << " could not be opened!\n";
    return kFALSE;
  }
  TFile *MuonFRFile = new TFile(fMuonFRFilename);  
  if (!MuonFRFile) {
    cerr << "The Muon FakeRate file : " << fMuonFRFilename << " could not be opened!\n";
    return kFALSE;
  }
  
  // load histogram/fit functions from the file
  if (fUse2DFakeRate) {

    fElectronFakeRateHist_PtEta = (TH2DAsymErr*)(ElectronFRFile->Get(fElectronFRHistName));
    fMuonFakeRateHist_PtEta = (TH2DAsymErr*)(MuonFRFile->Get(fMuonFRHistName));

    if (!fElectronFakeRateHist_PtEta) {
      cout << "Error: Histogram " << fElectronFRHistName << " cannot be loaded from file "
           << fElectronFRFilename << endl;
    }
    if (!fMuonFakeRateHist_PtEta) {
      cout << "Error: Histogram " << fMuonFRHistName << " cannot be loaded. from file"
           << fMuonFRFilename << endl;
    }

    fElectronFakeRateHist_PtEta->SetDirectory(0);
    fMuonFakeRateHist_PtEta->SetDirectory(0);

    if (fUseFitFunction) {
      //Currently unsupported
      cerr << "Error: Using 2D Fake Rates with Fit Function is not currently supported.\n";
      return kFALSE;

      //Once supported, below code will load the functions.
      TF2 *ElectronFakeRateFit_PtEta_temp = (TF2*)(ElectronFRFile->Get(fElectronFRFunctionName));
      TF2 *MuonFakeRateFit_PtEta_temp = (TF2*)(MuonFRFile->Get(fMuonFRFunctionName));
      if (!ElectronFakeRateFit_PtEta_temp) {
        cout << "Error: Function " << fElectronFRFunctionName << " cannot be loaded from file "
             << fElectronFRFilename << endl;
      }
      if (!MuonFakeRateFit_PtEta_temp) {
      cout << "Error: Function " << fMuonFRFunctionName << " cannot be loaded. from file"
           << fMuonFRFilename << endl;
      }
      fElectronFakeRateFit_PtEta = (TF2*)(ElectronFakeRateFit_PtEta_temp->Clone());
      fMuonFakeRateFit_PtEta = (TF2*)(MuonFakeRateFit_PtEta_temp->Clone());
    }

  } else {

    fElectronFakeRateHist_Pt = (TH1F*)(ElectronFRFile->Get(fElectronFRHistName)); 
    fMuonFakeRateHist_Pt = (TH1F*)(MuonFRFile->Get(fMuonFRHistName)); 
    if (!fElectronFakeRateHist_Pt) {
      cout << "Error: Histogram " << fElectronFRHistName << " cannot be loaded from file "
           << fElectronFRFilename << endl;
    }
    if (!fMuonFakeRateHist_Pt) {
      cout << "Error: Histogram " << fMuonFRHistName << " cannot be loaded. from file"
           << fMuonFRFilename << endl;
    }
    fElectronFakeRateHist_Pt->SetDirectory(0);
    fMuonFakeRateHist_Pt->SetDirectory(0);

    if (fUseFitFunction) {
      TF1 *ElectronFakeRateFit_Pt_temp = (TF1*)(ElectronFRFile->Get(fElectronFRFunctionName));
      TF1 *MuonFakeRateFit_Pt_temp = (TF1*)(MuonFRFile->Get(fMuonFRFunctionName));
      if (!ElectronFakeRateFit_Pt_temp) {
        cout << "Error: Function " << fElectronFRFunctionName << " cannot be loaded from file "
             << fElectronFRFilename << endl;
      }
      if (!MuonFakeRateFit_Pt_temp) {
      cout << "Error: Function " << fMuonFRFunctionName << " cannot be loaded. from file"
           << fMuonFRFilename << endl;
      }
      fElectronFakeRateFit_Pt = (TF1*)(ElectronFakeRateFit_Pt_temp->Clone());
      fMuonFakeRateFit_Pt = (TF1*)(MuonFakeRateFit_Pt_temp->Clone());
    }
  }

  ElectronFRFile->Close();
  MuonFRFile->Close();
  delete ElectronFRFile;
  delete MuonFRFile;

  return true;
}


//--------------------------------------------------------------------------------------------------
Double_t FakeRate::ElectronFakeRate(Double_t pt, Double_t eta, Double_t phi)
{
  // Calculate the electron fake rate given pt, eta, and phi.

  Double_t prob = 0.0;

  if (fIsInit) {
    if (fUse2DFakeRate) {
      if (fUseFitFunction) {
        cerr << "Error: Using 2D Fake Rates with Fit Function is not currently supported.\n";      
      } else {
        if (fElectronFakeRateHist_PtEta) {
          Int_t ptbin = fElectronFakeRateHist_PtEta->GetXaxis()->FindFixBin(pt);
          Int_t etabin = fElectronFakeRateHist_PtEta->GetYaxis()->FindFixBin(eta);    
          prob = fElectronFakeRateHist_PtEta->GetBinContent(ptbin,etabin);
        } else {
          Fatal("ElectronFakeRate","Error: fElectronFakeRateHist_PtEta was not loaded properly.");
        }
      }
    } else {
      if (fUseFitFunction) {
        if (fElectronFakeRateFit_Pt) {
          prob = fElectronFakeRateFit_Pt->Eval(pt);
        } else {
          Fatal("ElectronFakeRate","Error: fElectronFakeRateFit_Pt was not loaded properly.");
        }
      } else {
        if (fElectronFakeRateHist_Pt) {
          Int_t ptbin = fElectronFakeRateHist_Pt->GetXaxis()->FindFixBin(pt);
          prob = fElectronFakeRateHist_Pt->GetBinContent(ptbin);
        } else {
          Fatal("ElectronFakeRate","Error: fElectronFakeRateHist_Pt was not loaded properly.");
        }
      }
    }
  } else {
    Fatal("ElectronFakeRate","Error: FakeRate was not properly initialized.");
  }
  if (prob > 1) {
    cerr << "Fake Rate = " << prob << " is larger than 1.0" << endl;
  }
  return prob;
}

//--------------------------------------------------------------------------------------------------
Double_t FakeRate::ElectronFakeRateError(Double_t pt, Double_t eta, Double_t phi, 
                                         TH2DAsymErr::EErrType errorType)
{
  // Calculate the electron fake rate error given pt, eta, and phi

  Double_t error = 0.0;

  if (fIsInit) {
    if (fUse2DFakeRate) {
      if (fUseFitFunction) {
      } else {
        if (fElectronFakeRateHist_PtEta) { 
          return fElectronFakeRateHist_PtEta->GetError(pt, eta, errorType);        
        } else {
          Fatal("ElectronFakeRate","Error: fElectronFakeRateHist_PtEta_sysError was not loaded properly.");
        }
      }
    } 
  } else {
    Fatal("ElectronFakeRate","Error: FakeRate was not properly initialized.");
  }
  return error;
}

//--------------------------------------------------------------------------------------------------
Double_t FakeRate::ElectronFakeRateStatErrorLow(Double_t pt, Double_t eta, Double_t phi)
{
  // Calculate the electron fake rate lower statistical error given pt, eta, and phi
  return ElectronFakeRateError(pt, eta, phi, mithep::TH2DAsymErr::kStatErrLow);
}

//--------------------------------------------------------------------------------------------------
Double_t FakeRate::ElectronFakeRateStatErrorHigh(Double_t pt, Double_t eta, Double_t phi)
{
  // Calculate the electron fake rate upper statistical error given pt, eta, and phi
  return ElectronFakeRateError(pt, eta, phi, TH2DAsymErr::kStatErrHigh);
}

//--------------------------------------------------------------------------------------------------
Double_t FakeRate::ElectronFakeRateSysErrorLow(Double_t pt, Double_t eta, Double_t phi)
{
  // Calculate the electron fake rate  lower systematic error given pt, eta, and phi
  return ElectronFakeRateError(pt, eta, phi, TH2DAsymErr::kSysErrLow);
}

//--------------------------------------------------------------------------------------------------
Double_t FakeRate::ElectronFakeRateSysErrorHigh(Double_t pt, Double_t eta, Double_t phi)
{
  // Calculate the electron fake rate upper systematic error given pt, eta, and phi
  return ElectronFakeRateError(pt, eta, phi, TH2DAsymErr::kSysErrHigh);
}

//--------------------------------------------------------------------------------------------------
Double_t FakeRate::ElectronFakeRateErrorLow(Double_t pt, Double_t eta, Double_t phi)
{
  // Calculate the electron fake rate total lower error given pt, eta, and phi
  return TMath::Sqrt( ElectronFakeRateError(pt, eta, phi, TH2DAsymErr::kSysErrLow)*
                      ElectronFakeRateError(pt, eta, phi, TH2DAsymErr::kSysErrLow) +
                      ElectronFakeRateError(pt, eta, phi, TH2DAsymErr::kStatErrLow)*
                      ElectronFakeRateError(pt, eta, phi, TH2DAsymErr::kStatErrLow)
    );
}

//--------------------------------------------------------------------------------------------------
Double_t FakeRate::ElectronFakeRateErrorHigh(Double_t pt, Double_t eta, Double_t phi)
{
  // Calculate the electron fake rate total upper error given pt, eta, and phi
  return TMath::Sqrt( ElectronFakeRateError(pt, eta, phi, TH2DAsymErr::kSysErrHigh)*
                      ElectronFakeRateError(pt, eta, phi, TH2DAsymErr::kSysErrHigh) +
                      ElectronFakeRateError(pt, eta, phi, TH2DAsymErr::kStatErrHigh)*
                      ElectronFakeRateError(pt, eta, phi, TH2DAsymErr::kStatErrHigh)
    );
}

//--------------------------------------------------------------------------------------------------
Double_t FakeRate::MuonFakeRate(Double_t pt, Double_t eta, Double_t phi)
{
  // Calculate the muon fake rate given pt, eta, and phi
  Double_t prob = 0.0;

  if (fIsInit) {
    if (fUse2DFakeRate) {
      if (fUseFitFunction) {
        cerr << "Error: Using 2D Fake Rates with Fit Function is not currently supported.\n";      
      } else {
        if (fMuonFakeRateHist_PtEta) {
          Int_t ptbin = fMuonFakeRateHist_PtEta->GetXaxis()->FindFixBin(pt);
          Int_t etabin = fMuonFakeRateHist_PtEta->GetYaxis()->FindFixBin(eta);    
          prob = fMuonFakeRateHist_PtEta->GetBinContent(ptbin,etabin);
        } else {
          Fatal("ElectronFakeRate","Error: fMuonFakeRateHist_PtEta was not loaded properly.");
        }
      }
    } else {
      if (fUseFitFunction) {
        if (fMuonFakeRateFit_Pt) {
          prob = fMuonFakeRateFit_Pt->Eval(pt);
        } else {
          Fatal("ElectronFakeRate","Error: fMuonFakeRateFit_Pt was not loaded properly.");
        }
      } else {
        if (fMuonFakeRateHist_Pt) {
          Int_t ptbin = fMuonFakeRateHist_Pt->GetXaxis()->FindFixBin(pt);
          prob = fMuonFakeRateHist_Pt->GetBinContent(ptbin);
        } else {
          Fatal("ElectronFakeRate","Error: fMuonFakeRateHist_Pt was not loaded properly.");
        }
      }
    }
  } else {
    Fatal("ElectronFakeRate","Error: FakeRate was not properly initialized.");
  }
  if (prob > 1) {
    cerr << "Fake Rate is larger than 1.0" << endl;
  }
  return prob;
}

//--------------------------------------------------------------------------------------------------
Double_t FakeRate::MuonFakeRateError(Double_t pt, Double_t eta, Double_t phi, 
                                     TH2DAsymErr::EErrType errorType)
{
  // Calculate the muon fake rate error given pt, eta, and phi.

  Double_t error = 0.0;

  if (fIsInit) {
    if (fUse2DFakeRate) {
      if (fUseFitFunction) {
        cerr << "Error: Using 2D Fake Rates with Fit Function is not currently supported.\n";      
      } else {
        if (fMuonFakeRateHist_PtEta) {
          return fMuonFakeRateHist_PtEta->GetError(pt, eta, errorType);                 
        } else {
          Fatal("ElectronFakeRate","Error: fMuonFakeRateHist_PtEta_sysError was not loaded properly.");
        }
      }
    } 
  } else {
    Fatal("ElectronFakeRate","Error: FakeRate was not properly initialized.");
  }
  return error;
}

//--------------------------------------------------------------------------------------------------
Double_t FakeRate::MuonFakeRateStatErrorLow(Double_t pt, Double_t eta, Double_t phi)
{
  // Calculate the muon fake rate lower statistical error given pt, eta, and phi.

  return MuonFakeRateError(pt, eta, phi, TH2DAsymErr::kStatErrLow);
}

//--------------------------------------------------------------------------------------------------
Double_t FakeRate::MuonFakeRateStatErrorHigh(Double_t pt, Double_t eta, Double_t phi)
{
  // Calculate the muon fake rate upper statistical error given pt, eta, and phi.

  return MuonFakeRateError(pt, eta, phi, TH2DAsymErr::kStatErrHigh);
}

//--------------------------------------------------------------------------------------------------
Double_t FakeRate::MuonFakeRateSysErrorLow(Double_t pt, Double_t eta, Double_t phi)
{
  // Calculate the muon fake rate lower systematic error given pt, eta, and phi.

  return MuonFakeRateError(pt, eta, phi, TH2DAsymErr::kSysErrLow);
}

//--------------------------------------------------------------------------------------------------
Double_t FakeRate::MuonFakeRateSysErrorHigh(Double_t pt, Double_t eta, Double_t phi)
{
  // Calculate the muon fake rate upper systematic error given pt, eta, and phi.

  return MuonFakeRateError(pt, eta, phi, TH2DAsymErr::kSysErrHigh);
}

//--------------------------------------------------------------------------------------------------
Double_t FakeRate::MuonFakeRateErrorLow(Double_t pt, Double_t eta, Double_t phi)
{
  // Calculate the muon fake rate total lower error given pt, eta, and phi.

  return TMath::Sqrt( MuonFakeRateError(pt, eta, phi, mithep::TH2DAsymErr::kSysErrLow)*
                      MuonFakeRateError(pt, eta, phi, mithep::TH2DAsymErr::kSysErrLow) +
                      MuonFakeRateError(pt, eta, phi, mithep::TH2DAsymErr::kStatErrLow)*
                      MuonFakeRateError(pt, eta, phi, mithep::TH2DAsymErr::kStatErrLow)
    );
}

//--------------------------------------------------------------------------------------------------
Double_t FakeRate::MuonFakeRateErrorHigh(Double_t pt, Double_t eta, Double_t phi)
{
  // Calculate the muon fake rate total upper error given pt, eta, and phi.

  return TMath::Sqrt( MuonFakeRateError(pt, eta, phi, TH2DAsymErr::kSysErrHigh)*
                      MuonFakeRateError(pt, eta, phi, TH2DAsymErr::kSysErrHigh) +
                      MuonFakeRateError(pt, eta, phi, TH2DAsymErr::kStatErrHigh)*
                      MuonFakeRateError(pt, eta, phi, TH2DAsymErr::kStatErrHigh)
    );
}
