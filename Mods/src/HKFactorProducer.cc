// $Id: HKFactorProducer.cc,v 1.10 2011/03/18 13:46:57 ceballos Exp $

#include "MitPhysics/Mods/interface/HKFactorProducer.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include <TTree.h>
#include <TFile.h>

using namespace mithep;

ClassImp(mithep::HKFactorProducer)

//--------------------------------------------------------------------------------------------------
HKFactorProducer::HKFactorProducer(const char *name, const char *title) : 
  BaseMod(name,title),
  fProcessID(102),
  fInputFileName("setme"),
  fMCBosonsName(ModNames::gkMCBosonsName),
  fMCEvInfoName(Names::gkMCEvtInfoBrn),
  fIsData(kFALSE),
  fMakePDFNtuple(kFALSE),
  fPt_histo(0),
  fMCEventInfo(0),
  fOutputFile(0),
  fOutputName("ntuple.root")
{
  // Constructor
}

//--------------------------------------------------------------------------------------------------
HKFactorProducer::~HKFactorProducer()
{
  // Destructor
  delete fPt_histo;
}

//--------------------------------------------------------------------------------------------------
void HKFactorProducer::Process()
{
  // Process entries of the tree.

  Double_t theWeight = 1.0;

  if(fIsData == kFALSE){
    // get the bosons
    MCParticleCol *GenBosons = GetObjThisEvt<MCParticleCol>(fMCBosonsName);

    LoadBranch(fMCEvInfoName);

    // only accept the exact process id
    if (fProcessID == fMCEventInfo->ProcessId()) { 

      Double_t ptH = -1.0;
      for (UInt_t i=0; i<GenBosons->GetEntries(); ++i) {
	if(GenBosons->At(i)->PdgId() == MCParticle::kH) {
    	  ptH = GenBosons->At(i)->Pt();
    	  break;
	}
      }

      if(ptH >= 0) {
	// calculate bin size
	Double_t binsize = (fPt_histo->GetXaxis()->GetXmax()-fPt_histo->GetXaxis()->GetXmin())/fPt_histo->GetNbinsX();
	Int_t bin = 0;
	// underflow protection: use underflow entry
	if(ptH >= fPt_histo->GetXaxis()->GetXmin()){
	  bin = Int_t((ptH-fPt_histo->GetXaxis()->GetXmin())/binsize) + 1;
	}
	// overflow protection: use overflow entry
	if(bin > fPt_histo->GetNbinsX()) bin=fPt_histo->GetNbinsX()+1;
	theWeight = fPt_histo->GetBinContent(bin);

	if (0) {
          cout << "Bin Size: " << binsize << ", Higgs Pt: " << ptH
               << ", Bin: "<< bin  << ", KFactor: "<< fPt_histo->GetBinContent(bin) << endl;
	}

	if (GetFillHist()) {
          hDHKFactor[0]->Fill(0.5);
          hDHKFactor[1]->Fill(TMath::Min(ptH,499.999));
          hDHKFactor[2]->Fill(TMath::Min(ptH,499.999),theWeight);
          hDHKFactor[3]->Fill(TMath::Max(TMath::Min(theWeight,3.999),-3.999));
	}
      }
    } else if (fProcessID == 999){
      theWeight = fMCEventInfo->Weight();
      if (GetFillHist()) hDHKFactor[3]->Fill(TMath::Max(TMath::Min(theWeight,3.999),-3.999));
      theWeight = theWeight/fabs(theWeight);
      if (GetFillHist()) hDHKFactor[0]->Fill(0.5,theWeight);
    }
    else {
      if (GetFillHist()) hDHKFactor[0]->Fill(0.5,1.0);
    }
    // process id distribution
    if (GetFillHist()) hDHKFactor[4]->Fill(TMath::Min((Double_t)fMCEventInfo->ProcessId(),999.499));

    if (fMakePDFNtuple == kTRUE){
      fTreeVariables[0] = fMCEventInfo->Weight();
      fTreeVariables[1] = fMCEventInfo->Scale();
      fTreeVariables[2] = fMCEventInfo->Id1();
      fTreeVariables[3] = fMCEventInfo->X1();
      fTreeVariables[4] = fMCEventInfo->Pdf1();
      fTreeVariables[5] = fMCEventInfo->Id2();
      fTreeVariables[6] = fMCEventInfo->X2();
      fTreeVariables[7] = fMCEventInfo->Pdf2();
      fTree->Fill();
    }
  }

  TParameter<Double_t> *NNLOWeight = new TParameter<Double_t>("NNLOWeight", theWeight);
  AddObjThisEvt(NNLOWeight);
}

//--------------------------------------------------------------------------------------------------
void HKFactorProducer::SlaveBegin()
{
  // Book branch and histograms if wanted.

  if(fIsData == kFALSE){
    ReqBranch(fMCEvInfoName, fMCEventInfo);
  }

  if (!fPt_histo) {
    Info("SlaveBegin", "Using %s as input data file", fInputFileName.Data());
    fPt_histo = new  HWWKfactorList("KFactorList", fInputFileName);
  }

  if (GetFillHist()) {
    char sb[1024];
    sprintf(sb,"hDHKFactor_%d", 0);  hDHKFactor[0]  = new TH1D(sb,sb,1,0,1); 
    sprintf(sb,"hDHKFactor_%d", 1);  hDHKFactor[1]  = new TH1D(sb,sb,500,0.,500.); 
    sprintf(sb,"hDHKFactor_%d", 2);  hDHKFactor[2]  = new TH1D(sb,sb,500,0.,500.); 
    sprintf(sb,"hDHKFactor_%d", 3);  hDHKFactor[3]  = new TH1D(sb,sb,400,-4.0,4.0); 
    sprintf(sb,"hDHKFactor_%d", 4);  hDHKFactor[4]  = new TH1D(sb,sb,1000,-0.5,999.5); 
    for(Int_t i=0; i<5; i++) AddOutput(hDHKFactor[i]);
  }

  //***********************************************************************************************
  //Create Ntuple Tree  
  //***********************************************************************************************
  if (fMakePDFNtuple == kTRUE){
    printf("... init PDF ntuple ...\n");
    fOutputFile = new TFile(fOutputName, "RECREATE");
    fTree = new TTree("PDFTree", "PDFTree");
    const char* TreeFormat;
    TreeFormat = "weight/F:lq:lid1:lx1:lpdf1:lid2:lx2:lpdf2";
    fTree->Branch("H", &fTreeVariables,TreeFormat);    
  }
}

//--------------------------------------------------------------------------------------------------
void HKFactorProducer::SlaveTerminate(){
  fOutputFile->cd();
  fTree->Write();
  fOutputFile->Close();
}
