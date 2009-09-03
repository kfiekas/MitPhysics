// $Id: HKFactorProducer.cc,v 1.3 2009/06/15 15:00:21 loizides Exp $

#include "MitPhysics/Mods/interface/HKFactorProducer.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>

using namespace mithep;

ClassImp(mithep::HKFactorProducer)

//--------------------------------------------------------------------------------------------------
HKFactorProducer::HKFactorProducer(const char *name, const char *title) : 
  BaseMod(name,title),
  fProcessID(102),
  fInputFileName("setme"),
  fMCBosonsName(ModNames::gkMCBosonsName),
  fMCEvInfoName(Names::gkMCEvtInfoBrn),
  fPt_histo(0),
  fMCEventInfo(0)
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

  // get the bosons
  MCParticleCol *GenBosons = GetObjThisEvt<MCParticleCol>(fMCBosonsName);

  LoadBranch(fMCEvInfoName);

  Double_t theWeight = 1.0;

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
      Double_t binsize = fPt_histo->GetXaxis()->GetXmax()/fPt_histo->GetNbinsX();
      // get bin
      Int_t bin = Int_t((ptH/binsize)) + 1;
      // overflow protection: use last entry
      if(bin > fPt_histo->GetNbinsX()) bin=fPt_histo->GetNbinsX();
      theWeight = fPt_histo->GetBinContent(bin);

      if (0) {
        cout << "Bin Size: " << binsize << ", Higgs Pt: " << ptH
             << ", Bin: "<< bin  << ", KFactor: "<< fPt_histo->GetBinContent(bin) << endl;
      }

      if (GetFillHist()) {
        hDHKFactor[0]->Fill(0.5);
        hDHKFactor[1]->Fill(TMath::Min(ptH,499.999));
        hDHKFactor[2]->Fill(TMath::Min(ptH,499.999),theWeight);
        hDHKFactor[3]->Fill(TMath::Min(theWeight,9.999));
      }
    }
  }
  // process id distribution
  if (GetFillHist()) hDHKFactor[4]->Fill(TMath::Min(fMCEventInfo->ProcessId(),999.499));

  TParameter<Double_t> *NNLOWeight = new TParameter<Double_t>("NNLOWeight", theWeight);
  AddObjThisEvt(NNLOWeight);
}

//--------------------------------------------------------------------------------------------------
void HKFactorProducer::SlaveBegin()
{
  // Book branch and histograms if wanted.

  ReqBranch(fMCEvInfoName, fMCEventInfo);

  if (!fPt_histo) {
    Info("SlaveBegin", "Using %s as input data file", fInputFileName.Data());
    fPt_histo = new  HWWKfactorList("KFactorList", fInputFileName);
  }

  if (GetFillHist()) {
    char sb[1024];
    sprintf(sb,"hDHKFactor_%d", 0);  hDHKFactor[0]  = new TH1D(sb,sb,1,0,1); 
    sprintf(sb,"hDHKFactor_%d", 1);  hDHKFactor[1]  = new TH1D(sb,sb,500,0.,500.); 
    sprintf(sb,"hDHKFactor_%d", 2);  hDHKFactor[2]  = new TH1D(sb,sb,500,0.,500.); 
    sprintf(sb,"hDHKFactor_%d", 3);  hDHKFactor[3]  = new TH1D(sb,sb,400,0.,4.); 
    sprintf(sb,"hDHKFactor_%d", 4);  hDHKFactor[4]  = new TH1D(sb,sb,1000,-0.5,999.5); 
    for(Int_t i=0; i<5; i++) AddOutput(hDHKFactor[i]);
  }
}
