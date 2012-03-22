// $Id: EGEnergyCorrector.cc,v 1.5 2011/12/11 00:03:05 bendavid Exp $

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
    
    fVals = new Float_t[73];
    
    TFile *fgbr = new TFile(regweights,"READ");
    fReadereb = (GBRForest*)fgbr->Get("EBCorrection");
    fReaderebvariance = (GBRForest*)fgbr->Get("EBUncertainty");  
    fReaderee = (GBRForest*)fgbr->Get("EECorrection");
    fReadereevariance = (GBRForest*)fgbr->Get("EEUncertainty");      
    fgbr->Close();

}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrector::CorrectEnergyWithError(Photon *p, const VertexCol *vtxs, UInt_t version) {
  
  std::pair<double,double> correction;
  if (version == 1) {
    correction = CorrectedEnergyWithError(p);
  }
  else if (version == 2) {
    correction = CorrectedEnergyWithErrorV2(p,vtxs);
  }

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

//--------------------------------------------------------------------------------------------------
std::pair<double,double> EGEnergyCorrector::CorrectedEnergyWithErrorV2(const Photon *p, const VertexCol *vtxs) {
  
  const SuperCluster *s = p->SCluster();
  const BasicCluster *b = s->Seed();
  
  const BasicCluster *b2 = 0;
  Double_t ebcmax = -99.;
  for (UInt_t i=0; i<s->ClusterSize(); ++i) {
    const BasicCluster *bc = s->Cluster(i);
    if (bc->Energy() > ebcmax && bc !=b) {
      b2 = bc;
      ebcmax = bc->Energy();
    }
  }

  const BasicCluster *bclast = 0;
  Double_t ebcmin = 1e6;
  for (UInt_t i=0; i<s->ClusterSize(); ++i) {
    const BasicCluster *bc = s->Cluster(i);
    if (bc->Energy() < ebcmin && bc !=b) {
      bclast = bc;
      ebcmin = bc->Energy();
    }
  }

  const BasicCluster *bclast2 = 0;
  ebcmin = 1e6;
  for (UInt_t i=0; i<s->ClusterSize(); ++i) {
    const BasicCluster *bc = s->Cluster(i);
    if (bc->Energy() < ebcmin && bc !=b && bc!=bclast) {
      bclast2 = bc;
      ebcmin = bc->Energy();
    }
  }

  
  double scpse = s->PreshowerEnergy();
  
  Bool_t isbarrel = (s->AbsEta()<1.5);
  Bool_t hasbc2 = b2 && b2->Energy()>0.;
  Bool_t hasbclast = bclast && bclast->Energy()>0.;
  Bool_t hasbclast2 = bclast2 && bclast2->Energy()>0.;
  
  
  if (isbarrel) {
    fVals[0]  = s->RawEnergy();
    fVals[1]  = p->R9();
    fVals[2]  = s->Eta();
    fVals[3]  = s->Phi();
    fVals[4]  = p->E55()/s->RawEnergy();   
    fVals[5] = p->HadOverEm();
    fVals[6] = s->EtaWidth();
    fVals[7] = s->PhiWidth();
    
    fVals[8] = b->Eta()-s->Eta();
    fVals[9] = MathUtils::DeltaPhi(b,s);
    fVals[10] = b->Energy()/s->RawEnergy();
    fVals[11] = b->E3x3()/b->Energy();
    fVals[12] = b->E5x5()/b->Energy();
    fVals[13] = TMath::Sqrt(b->CoviEtaiEta());
    fVals[14] = TMath::Sqrt(b->CoviPhiiPhi());
    fVals[15] = b->CoviEtaiPhi();
    fVals[16] = b->EMax()/b->Energy();
    fVals[17] = log(b->E2nd()/b->EMax());
    fVals[18] = log(b->ETop()/b->EMax());
    fVals[19] = log(b->EBottom()/b->EMax());
    fVals[20] = log(b->ELeft()/b->EMax());
    fVals[21] = log(b->ERight()/b->EMax());
    fVals[22] = (b->ETop()-b->EBottom())/(b->ETop()+b->EBottom());
    fVals[23] = (b->ELeft()-b->ERight())/(b->ELeft()+b->ERight());

    fVals[24] = hasbc2 ? (b2->Eta()-s->Eta()) : 0.;
    fVals[25] = hasbc2 ? MathUtils::DeltaPhi(b2,s) : 0.;
    fVals[26] = hasbc2 ? b2->Energy()/s->RawEnergy() : 0.;
    fVals[27] = hasbc2 ? b2->E3x3()/b2->Energy() : 0.;
    fVals[28] = hasbc2 ? b2->E5x5()/b2->Energy() : 0.;
    fVals[29] = hasbc2 ? TMath::Sqrt(b2->CoviEtaiEta()) : 0.;
    fVals[30] = hasbc2 ? TMath::Sqrt(b2->CoviPhiiPhi()) : 0.;
    //fVals[31] = hasbc2 ? b2->CoviEtaiPhi() : 0.;
    fVals[31] = hasbc2 ? b->CoviEtaiPhi() : 0.;
    fVals[32] = hasbc2 ? b2->EMax()/b2->Energy() : 0.;
    fVals[33] = hasbc2 ? log(b2->E2nd()/b2->EMax()) : 0.;
    fVals[34] = hasbc2 ? log(b2->ETop()/b2->EMax()) : 0.;
    fVals[35] = hasbc2 ? log(b2->EBottom()/b2->EMax()) : 0.;
    fVals[36] = hasbc2 ? log(b2->ELeft()/b2->EMax()) : 0.;
    fVals[37] = hasbc2 ? log(b2->ERight()/b2->EMax()) : 0.;
    fVals[38] = hasbc2 ? (b2->ETop()-b2->EBottom())/(b2->ETop()+b2->EBottom()) : 0.;
    fVals[39] = hasbc2 ? (b2->ELeft()-b2->ERight())/(b2->ELeft()+b2->ERight()) : 0.;

    fVals[40] = hasbclast ? (bclast->Eta()-s->Eta()) : 0.;
    fVals[41] = hasbclast ? MathUtils::DeltaPhi(bclast,s) : 0.;
    fVals[42] = hasbclast ? bclast->Energy()/s->RawEnergy() : 0.;
    fVals[43] = hasbclast ? bclast->E3x3()/bclast->Energy() : 0.;
    fVals[44] = hasbclast ? bclast->E5x5()/bclast->Energy() : 0.;
    fVals[45] = hasbclast ? TMath::Sqrt(bclast->CoviEtaiEta()) : 0.;
    fVals[46] = hasbclast ? TMath::Sqrt(bclast->CoviPhiiPhi()) : 0.;
    fVals[47] = hasbclast ? bclast->CoviEtaiPhi() : 0.;    

    fVals[48] = hasbclast2 ? (bclast2->Eta()-s->Eta()) : 0.;
    fVals[49] = hasbclast2 ? MathUtils::DeltaPhi(bclast2,s) : 0.;
    fVals[50] = hasbclast2 ? bclast2->Energy()/s->RawEnergy() : 0.;
    fVals[51] = hasbclast2 ? bclast2->E3x3()/bclast2->Energy() : 0.;
    fVals[52] = hasbclast2 ? bclast2->E5x5()/bclast2->Energy() : 0.;
    fVals[53] = hasbclast2 ? TMath::Sqrt(bclast2->CoviEtaiEta()) : 0.;
    fVals[54] = hasbclast2 ? TMath::Sqrt(bclast2->CoviPhiiPhi()) : 0.;
    fVals[55] = hasbclast2 ? bclast2->CoviEtaiPhi() : 0.;        
    
    fVals[56] = b->IEta();
    fVals[57] = b->IPhi();
    fVals[58] = b->IEta()%5;
    fVals[59] = b->IPhi()%2;
    fVals[60] = (TMath::Abs(b->IEta())<=25)*(b->IEta()%25) + (TMath::Abs(b->IEta())>25)*((b->IEta()-25*TMath::Abs(b->IEta())/b->IEta())%20);
    fVals[61] = b->IPhi()%20;
    fVals[62] = b->EtaCry();
    fVals[63] = b->PhiCry();

    fVals[64] = hasbc2 ? b2->IEta() : 0.;
    fVals[65] = hasbc2 ? b2->IPhi() : 0.;
    fVals[66] = hasbc2 ? b2->IEta()%5 : 0.;
    fVals[67] = hasbc2 ? b2->IPhi()%2 : 0.;
    fVals[68] = hasbc2 ? (TMath::Abs(b2->IEta())<=25)*(b2->IEta()%25) + (TMath::Abs(b2->IEta())>25)*((b2->IEta()-25*TMath::Abs(b2->IEta())/b2->IEta())%20) : 0.;
    fVals[69] = hasbc2 ? b2->IPhi()%20 : 0.;
    fVals[70] = hasbc2 ? b2->EtaCry() : 0.;
    fVals[71] = hasbc2 ? b2->PhiCry() : 0.;    
    
    fVals[72] = vtxs->GetEntries();
    
  }
  else {
    fVals[0]  = s->RawEnergy();
    fVals[1]  = p->R9();
    fVals[2]  = s->Eta();
    fVals[3]  = s->Phi();
    fVals[4]  = p->E55()/s->RawEnergy();
    fVals[5] = s->EtaWidth();
    fVals[6] = s->PhiWidth();
    fVals[7] = vtxs->GetEntries();
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