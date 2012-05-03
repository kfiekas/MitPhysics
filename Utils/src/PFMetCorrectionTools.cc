// $Id: PFMetCorrectionTools.cc,adapted from RedNtpTree.cc from Chiara 2012/04/30  Heng $

#include "MitPhysics/Utils/interface/PFMetCorrectionTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/Photon.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include <TFile.h>
#include <TRandom3.h>
// #include "CondFormats/EgammaObjects/interface/GBRForest.h"
// #include "Cintex/Cintex.h"
#include <TLorentzVector.h>
using namespace mithep;
ClassImp(mithep::PFMetCorrectionTools)



//--------------------------------------------------------------------------------------------------
PFMetCorrectionTools::PFMetCorrectionTools()  
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------

// pfjet resolutions. taken from AN-2010-371
Double_t PFMetCorrectionTools::ErrEt( Double_t Et, Double_t Eta) {
  
  Double_t InvPerr2;

  Double_t N = 0.;
  Double_t m = 0.;
  Double_t S = 0.;
  Double_t C = 0.;

  if(fabs(Eta) < 0.5 ) {
    N = 3.96859;
    S = 0.18348;
    C = 0.;
    m = 0.62627;
  } else if( fabs(Eta) < 1. ) {
    N = 3.55226;
    S = 0.24026;
    C = 0.;
    m = 0.52571;
  } else if( fabs(Eta) < 1.5 ) {
    N = 4.54826;
    S = 0.22652;
    C = 0.;
    m = 0.58963;
  } else if( fabs(Eta) < 2. ) {
    N = 4.62622;
    S = 0.23664;
    C = 0.;
    m = 0.48738;
  } else if( fabs(Eta) < 3. ) {
    N = 2.53324;
    S = 0.34306;
    C = 0.;
    m = 0.28662;
  } else if( fabs(Eta) < 5. ) {
    N = 2.95397;
    S = 0.11619;
    C = 0.;
    m = 0.96086;
  }

  // this is the absolute resolution (squared), not sigma(pt)/pt
  // so have to multiply by pt^2, thats why m+1 instead of m-1
  InvPerr2 =  (N * fabs(N) ) + (S * S) * pow(Et, m+1) + (C * C) * Et * Et ;


  return sqrt(InvPerr2)/Et;

}

void PFMetCorrectionTools::correctMet(Met *met, const Photon *phHard, const Photon *phSoft,  Bool_t smearing, Bool_t scale, const PFJetCol *fPFJet, const GenJetCol *fGenJet, const JetCol *fcorrJet) {
  //fPFJet is the AKt5PFJets 
  //fGenJet is the AKT5GenJets
  //fcorrJets is the corrected jets collections

  TLorentzVector jetSumSmeared;
  jetSumSmeared.SetXYZT(0.,0.,0.,0);
  
  TLorentzVector jetSumUnsmeared;
  jetSumUnsmeared.SetXYZT(0.,0.,0.,0);

   if( !fPFJet || !fcorrJet) {
     return;
   }

  if(fPFJet->GetEntries() != fcorrJet->GetEntries()) return;

  // associating reco - gen met                                                                                                            
  for( unsigned int i=0; i<fPFJet->GetEntries(); ++i){
    const Jet *recojet=fPFJet->At(i);
    const Jet *corrjet=fcorrJet->At(i);
    
    if( !recojet || !corrjet )  return;
       
    if (! phSoft || !phHard ) {
      return;
    }
    
    // remove identified photons
    if ((phHard && MathUtils::DeltaR(recojet->RawMom(),phHard->Mom())<0.5) ||(phSoft && MathUtils::DeltaR(recojet->RawMom(),phSoft->Mom())<0.5)) continue; 
    // smearing via association with genjets
    int match        = -999;
    Double_t DRmin = 999.;     
    if (fGenJet){
      for(unsigned int j=0; j< fGenJet->GetEntries(); ++j){
	const GenJet *genjet=fGenJet->At(j);
	if(!genjet) {
	  std::cout<<" genjet not there..."<<std::endl;
	  continue;
	}
	
	Double_t DR = MathUtils::DeltaR(recojet->RawMom(),genjet->Mom());
	
	Double_t expres = ErrEt(corrjet->RawMom().Pt(),recojet->RawMom().Eta());  
	
	if(DR < DRmin && (corrjet->RawMom().Pt()-genjet->Mom().Pt())/corrjet->RawMom().Pt() < 5. * expres) {
	  match = (int) j;
	  DRmin = DR;
	}
      }
      if( match > -1 )
	if(DRmin > 0.1 + 0.3 * exp(-0.05*(fGenJet->At(match)->Mom().Pt()-10)))  match = -999;
    }
    else match=-999;
    
    // smearing for non-associated jets, using expected resolutions
    float smear = -999.;
    if (fabs(recojet->RawMom().Eta())<=1.1)                               smear = 1.06177;
    if (fabs(recojet->RawMom().Eta())<=1.7 && fabs(recojet->RawMom().Eta())>1.1) smear = 1.08352;
    if (fabs(recojet->RawMom().Eta())<=2.3 && fabs(recojet->RawMom().Eta())>1.7) smear = 1.02911;
    if (fabs(recojet->RawMom().Eta())>2.3)                                smear = 1.15288;
    
    
    
    Double_t shift=0;
    if( match>-1 && fGenJet )
      shift = (smear-1) * (corrjet->RawMom().Pt() - fGenJet->At(match)->Mom().Pt())/corrjet->RawMom().Pt();
    else {
      Double_t expres = ErrEt(recojet->RawMom().Pt(),recojet->RawMom().Eta());
      Double_t relsmear = expres * sqrt(smear*smear-1);
      gRandom->SetSeed(1);
      shift = gRandom->Gaus(0.,relsmear);
    }
    
    float ptSmeared  = recojet->RawMom().Pt();
    float eneSmeared = recojet->RawMom().Pt();
    
    if(smearing && shift>-1 && shift < 2) {
      ptSmeared  *= 1 + shift;
      eneSmeared *= 1 + shift;
    }
    
    // JEC scaling to correct for residual jet corrections
    if(scale) {
      Double_t factor=1;
      if(TMath::Abs(recojet->RawMom().Eta())<1.5) factor = 1.015;
      else if(TMath::Abs(recojet->RawMom().Eta())<3) factor = 1.04;
      else factor = 1.15;
      ptSmeared  *= factor;
      eneSmeared *= factor;
    }
    
    TLorentzVector thisJetSmeared;
    thisJetSmeared.SetPtEtaPhiE(ptSmeared,recojet->RawMom().Eta(),recojet->RawMom().Phi(),eneSmeared);
    
    TLorentzVector thisJetUnsmeared;
    
    thisJetUnsmeared.SetPtEtaPhiE(recojet->RawMom().Pt(),recojet->RawMom().Eta(),recojet->RawMom().Phi(),recojet->RawMom().E());
    
    if (recojet->RawMom().Pt()>10 && TMath::Abs(recojet->RawMom().Eta())<4.7) {
      jetSumSmeared   += thisJetSmeared;
      jetSumUnsmeared += thisJetUnsmeared;
    }
    
  }

  //  TLorentzVector correctedMet;
  //  correctedMet = uncormet + jetSumUnsmeared - jetSumSmeared;
  //  Double_t px=met.Pt()*cos(met.Phi())+jetSumUnsmeared.Px()-jetSumSmeared.Px();
  //  Double_t py=met.Pt()*sin(met.Phi())+jetSumUnsmeared.Py()-jetSumSmeared.Py();
  Double_t px=met->Px()+jetSumUnsmeared.Px()-jetSumSmeared.Px();
  Double_t py=met->Py()+jetSumUnsmeared.Py()-jetSumSmeared.Py();
  //  Double_t e = sqrt(px*px+py*py);
  //  met.SetPxPyPzE(px,px,0,e);
  met->SetMex(px);
  met->SetMey(py);
  return;
}

void PFMetCorrectionTools::shiftMet(Met *uncormet, Bool_t fIsData) {

  //  TLorentzVector correctedMet;
  Double_t spfMet=uncormet->SumEt();
  // correction for METx, METy bias
  Double_t px=0;
  Double_t py=0;
  // data
  if(fIsData){
    px =uncormet->Px() -0.00563109*spfMet+0.959742;
    py =uncormet->Py() +0.00586162*spfMet-0.540137;
  // MC
  }else{
    px = uncormet->Px()-0.00069992*spfMet+0.430059;
    py = uncormet->Py()+0.00262869*spfMet+0.210784;
  }
  //  e = sqrt(px*px+py*py);
  
  //  Met.SetPxPyPzE(px,py,0,e);
  uncormet->SetMex(px);
  uncormet->SetMey(py);
  return;
}

