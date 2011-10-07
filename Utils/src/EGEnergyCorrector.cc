// $Id: EGEnergyCorrector.cc,v 1.1 2011/09/08 15:51:24 bendavid Exp $

#include "MitPhysics/Utils/interface/EGEnergyCorrector.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/GBRForest.h"
#include "MitAna/DataTree/interface/StableData.h"
#include <TFile.h>
#include <TRandom3.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"


ClassImp(mithep::EGEnergyCorrector)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
EGEnergyCorrector::EGEnergyCorrector() :
fReadereb(0),
fReaderebvariance(0),
fReaderee(0),
fReadereevariance(0),
fMethodname("BDTG method"),
fIsInitialized(kFALSE),
fIsMC(kFALSE),
fVals(0)
{
  // Constructor.
}


//--------------------------------------------------------------------------------------------------
EGEnergyCorrector::~EGEnergyCorrector()
{
  
  if (fVals) delete [] fVals;
  if (fReadereb) delete fReadereb;
  if (fReaderebvariance) delete fReaderebvariance;  
  if (fReaderee) delete fReaderee;
  if (fReadereevariance) delete fReadereevariance;
}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrector::Initialize(Bool_t ismc, TString phfixstring, TString phfixfile, TString regweights) {
    fIsInitialized = kTRUE;
    fIsMC = ismc;
    fPhFix.initialise(std::string(phfixstring),std::string(phfixfile));

    if (fVals) delete [] fVals;
    if (fReadereb) delete fReadereb;
    if (fReaderebvariance) delete fReaderebvariance;  
    if (fReaderee) delete fReaderee;
    if (fReadereevariance) delete fReadereevariance;    
    
    fVals = new Float_t[18];
    
    TFile *fgbr = new TFile(regweights,"READ");
    fReadereb = (GBRForest*)fgbr->Get("EBCorrection");
    fReaderebvariance = (GBRForest*)fgbr->Get("EBUncertainty");  
    fReaderee = (GBRForest*)fgbr->Get("EECorrection");
    fReadereevariance = (GBRForest*)fgbr->Get("EEUncertainty");      
    fgbr->Close();

}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrector::CorrectEnergyWithError(Photon *p) {
  
  std::pair<double,double> correction = CorrectedEnergyWithError(p);
  //printf("photon: e = %5f, eerr = %5f, sceta = %5f\n",p->E(),p->EnergyErr(),p->SCluster()->Eta());
  //printf("gbr   : e = %5f, eerr = %5f\n",correction.first,correction.second);
  FourVectorM mom = p->Mom();
  double scale = correction.first/mom.E();
  p->SetMom(scale*mom.X(), scale*mom.Y(), scale*mom.Z(), scale*mom.E());
  p->SetEnergyErr(correction.second);
  p->SetEnergyErrSmeared(correction.second);
  
  //printf("initial pt = %5f, regression pt = %5f\n",mom.Pt(),p->Pt());
  
  return;

}

//--------------------------------------------------------------------------------------------------
std::pair<double,double> EGEnergyCorrector::CorrectedEnergyWithError(const Photon *p) {
  
  const SuperCluster *s = p->SCluster();
 // const BasicCluster *b = s->Seed();  
  fPhFix.setup(p->E(),s->Eta(),s->Phi(),p->R9()); 
  
  double scpse = s->PreshowerEnergy();
  
  Bool_t isbarrel = (s->AbsEta()<1.5);

  if (isbarrel) {
    fVals[0]  = s->RawEnergy();
    fVals[1]  = p->R9();
    fVals[2]  = s->Eta();
    fVals[3]  = s->Phi();
    fVals[4]  = p->E55()/s->RawEnergy();
    fVals[5]  = fPhFix.etaC();
    fVals[6]  = fPhFix.etaS();
    fVals[7]  = fPhFix.etaM();
    fVals[8]  = fPhFix.phiC();
    fVals[9]  = fPhFix.phiS();
    fVals[10] = fPhFix.phiM();    
    fVals[11] = p->HadOverEm();
    fVals[12] = s->EtaWidth();
    fVals[13] = s->PhiWidth();
    fVals[14] = p->CoviEtaiEta();
  }
  else {
    fVals[0]  = s->RawEnergy();
    fVals[1]  = p->R9();
    fVals[2]  = s->Eta();
    fVals[3]  = s->Phi();
    fVals[4]  = p->E55()/s->RawEnergy();
    fVals[5]  = scpse/s->RawEnergy();
    fVals[6]  = fPhFix.xZ();
    fVals[7]  = fPhFix.xC();
    fVals[8]  = fPhFix.xS();
    fVals[9]  = fPhFix.xM();
    fVals[10] = fPhFix.yZ();
    fVals[11] = fPhFix.yC();
    fVals[12] = fPhFix.yS();
    fVals[13] = fPhFix.yM();
    fVals[14] = p->HadOverEm();
    fVals[15] = s->EtaWidth();
    fVals[16] = s->PhiWidth();
    fVals[17] = p->CoviEtaiEta();    
  }
    
  const Double_t varscale = 1.253;
  Double_t den;
  const GBRForest *reader;
  const GBRForest *readervar;
  if (isbarrel) {
    den = s->RawEnergy();
    reader = fReadereb;
    readervar = fReaderebvariance;
  }
  else {
    den = s->RawEnergy() + scpse;
    reader = fReaderee;
    readervar = fReadereevariance;
  }
  
  Double_t ecor = reader->GetResponse(fVals)*den;
  Double_t ecorerr = readervar->GetResponse(fVals)*den*varscale;
  
  return std::pair<double,double>(ecor,ecorerr);
}