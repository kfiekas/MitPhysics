// $Id: HKFactorProducer.cc,v 1.31 2009/04/04 09:40:35 ceballos Exp $

#include "MitPhysics/Mods/interface/HKFactorProducer.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>

using namespace mithep;

ClassImp(mithep::HKFactorProducer)

//--------------------------------------------------------------------------------------------------
HKFactorProducer::HKFactorProducer(const char *name, const char *title) : 
  BaseMod(name,title),
  fDebug(kFALSE),
  fFillHist(kTRUE),
  fPt_histo(0),
  fProcessID(102),
  fInputFilename("dummy.root"),
  fMCBosonsName(ModNames::gkMCBosonsName),
  fMCEventInfoName("MCEventInfo"),
  fMCEventInfo(0)
{
}

HKFactorProducer::~HKFactorProducer()
{
  delete fPt_histo;
}

//--------------------------------------------------------------------------------------------------
void HKFactorProducer::Process()
{
  // Process entries of the tree.
  if (!fPt_histo) {
    cout << "HKFactorProducer: using " << fInputFilename.c_str() << " as input data file" << endl;
    fPt_histo = new  HWWKfactorList("KFactorList", fInputFilename.c_str() );
  }

  // these arrays will be filled in the loop of particles
  ObjArray<MCParticle> *GenBosons = dynamic_cast<ObjArray<MCParticle>* > 
                                      (FindObjThisEvt(fMCBosonsName.Data()));

  LoadBranch(fMCEventInfoName);
  Double_t theWeight = 1.0;
  if (fProcessID == fMCEventInfo->ProcessId()) { // Only accept the exact process id
    Double_t ptH = -1.0;
    for (UInt_t i=0; i<GenBosons->GetEntries(); ++i) {
      if(GenBosons->At(i)->PdgId() == MCParticle::kH) {
    	ptH = GenBosons->At(i)->Pt();
    	break;
      }
    }
  
    if(ptH >= 0) { // This should always be true, but just in case
      //calculate bin size
      Double_t binsize = fPt_histo->GetXaxis()->GetXmax()/fPt_histo->GetNbinsX();
      // which bin?
      Int_t bin = Int_t((ptH/binsize)) + 1;
      // overflow protection: use last entry
      if(bin > fPt_histo->GetNbinsX()) bin=fPt_histo->GetNbinsX();
      theWeight = fPt_histo->GetBinContent(bin);

      if (fDebug) {
        cout << "Bin Size: " << binsize;
        cout << ", Higgs Pt: " << ptH;
        cout << ", Bin: "<< bin;
        cout << ", KFactor: "<< fPt_histo->GetBinContent(bin) << endl;
      }
      if (fFillHist) {
        hDHKFactor[0]->Fill(0.5);
        hDHKFactor[1]->Fill(TMath::Min(ptH,499.999));
        hDHKFactor[2]->Fill(TMath::Min(ptH,499.999),theWeight);
        hDHKFactor[3]->Fill(TMath::Min(theWeight,9.999));
      }
    }
  }

  TParameter<Double_t> *NNLOWeight = new TParameter<Double_t>("NNLOWeight", theWeight);
  AddObjThisEvt(NNLOWeight);

}

//--------------------------------------------------------------------------------------------------
void HKFactorProducer::SlaveBegin()
{
  // Book branch and histograms if wanted.

  ReqBranch(fMCEventInfoName, fMCEventInfo);

  // fill histograms
  if (fFillHist) {
    char sb[1024];
    sprintf(sb,"hDHKFactor_%d", 0);  hDHKFactor[0]  = new TH1D(sb,sb,1,0,1); 
    sprintf(sb,"hDHKFactor_%d", 1);  hDHKFactor[1]  = new TH1D(sb,sb,500,0.,500.); 
    sprintf(sb,"hDHKFactor_%d", 2);  hDHKFactor[2]  = new TH1D(sb,sb,500,0.,500.); 
    sprintf(sb,"hDHKFactor_%d", 3);  hDHKFactor[3]  = new TH1D(sb,sb,400,0.,4.); 
    for(Int_t i=0; i<4; i++) AddOutput(hDHKFactor[i]);
  }
}
