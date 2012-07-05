// $Id: PFMetCorrectionTools.cc,v 1.5 2012/06/21 16:06:08 mtouch Exp $

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

void PFMetCorrectionTools::correctMet(Met *met, const Photon *phHard, const Photon *phSoft,  Bool_t smearing, Bool_t scale, const PFJetCol *fPFJet, const GenJetCol *fGenJet, const JetCol *fcorrJet, UInt_t evt ) {
  //fPFJet is the AKt5PFJets 
  //fGenJet is the AKT5GenJets
  //fcorrJets is the corrected jets collections

  TLorentzVector jetSumSmeared;
  jetSumSmeared.SetXYZT(0.,0.,0.,0);
  
  TLorentzVector jetSumUnsmeared;
  jetSumUnsmeared.SetXYZT(0.,0.,0.,0);
  TRandom3 *jSmearRan= new TRandom3(evt);

  
  if( !fPFJet || !fcorrJet) {
    return;
  }
  // printf("im here 1 \n");
  //printf("there are %d akt5pfjets and    %d corrJets  in event %d \n",fPFJet->GetEntries(),fcorrJet->GetEntries(),evt);
  // printf(" ind         uncorEta      uncorPt          factor     ptSmeared\n");

  if(fPFJet->GetEntries() != fcorrJet->GetEntries()) return;
  //  Int_t jetColIdx =0;
  //  Int_t uncorrJetColIdx=0;
  // associating reco - gen met                                                                                                            
  for( UInt_t i=0; i<fPFJet->GetEntries(); ++i){
    const Jet *recojet=fPFJet->At(i);
    float ptSmeared  = recojet->Mom().Pt();
    float eneSmeared = recojet->Mom().E();    
    if ( (phHard && MathUtils::DeltaR(recojet->Mom(),phHard->Mom())<0.5 )|| (phSoft && MathUtils::DeltaR(recojet->Mom(),phSoft->Mom())<0.5)){
      //         std::cout<<"removed jet  "<<recojet->Pt()<< "  "<<recojet->Eta()<<std::endl;
      continue;
    }

    // if (recojet->Mom().Pt()>1.) uncorrJetColIdx += 1;

//     if (corrjet->Pt() < 5.0) continue;
//     jetColIdx += 1;

    if (smearing){
      const Jet *corrjet=fcorrJet->At(i);
      if( !recojet || !corrjet )  return;
      
      if (! phSoft || !phHard ) {
	return;
      }


      //  printf(" jet %i pt is %f       corr jet pt is  %f   \n",i,recojet->Mom().Pt(),corrjet->Mom().Pt());
      // smearing via association with genjets
      Int_t match        = -999;
      Double_t DRmin = 999.;     
      if (fGenJet){
	for(unsigned int j=0; j< fGenJet->GetEntries(); ++j){
	  const GenJet *genjet=fGenJet->At(j);
	  if(!genjet) {
	    std::cout<<" genjet not there..."<<std::endl;
	    continue;
	  }
	  
	  Double_t DR = MathUtils::DeltaR(recojet->Mom(),genjet->Mom());
	  
	  Double_t expres = ErrEt(corrjet->Mom().Pt(),recojet->Mom().Eta());  
	  
	  if(DR < DRmin && (corrjet->Mom().Pt()-genjet->Mom().Pt())/corrjet->Mom().Pt() < 5. * expres) {
	    match = (int) j;
	    DRmin = DR;
	  }
	}
	if( match > -1 && DRmin > 0.1 + 0.3 * exp(-0.05*(fGenJet->At(match)->Mom().Pt()-10)))  match = -999;
      }
      else match=-999;
      
      // smearing for non-associated jets, using expected resolutions
      float smear = -999.;      if (fabs(recojet->Mom().Eta())<=1.1)                               smear = 1.06177;
      if (fabs(recojet->Mom().Eta())<=1.7 && fabs(recojet->Mom().Eta())>1.1) smear = 1.08352;
      if (fabs(recojet->Mom().Eta())<=2.3 && fabs(recojet->Mom().Eta())>1.7) smear = 1.02911;
      if (fabs(recojet->Mom().Eta())>2.3)                                smear = 1.15288;      
      
      
      Double_t shift=0;
      if( match>-1 && fGenJet )
	shift = (smear-1) * (corrjet->Mom().Pt() - fGenJet->At(match)->Mom().Pt())/corrjet->Mom().Pt();
      else {
	Double_t expres = ErrEt(recojet->Mom().Pt(),recojet->Mom().Eta());
	Double_t relsmear = expres * sqrt(smear*smear-1);
	jSmearRan->SetSeed(evt+(Int_t)(recojet->Mom().Eta()*1000));
	shift = jSmearRan->Gaus(0.,relsmear);
      }
      
           
      if(smearing && shift>-1 && shift < 2) {
	ptSmeared  *= 1 + shift;
	eneSmeared *= 1 + shift;
      }
//       if (match>-1 && match < (Double_t)(fGenJet->GetEntries()))
// 	printf("%i       %f    %f    %f      %f       %f        %f        %f          %f         %f \n",i, recojet->Mom().Eta(),recojet->Mom().Pt(),corrjet->Eta(),corrjet->Pt(),fGenJet->At(match)->Mom().Pt(),fGenJet->At(match)->Mom().Eta(),smear,shift,ptSmeared);
//       else
//         printf("%i       %f    %f    %f      %f       %f        %f        %f          %f         %f \n",i, recojet->Mom().Eta(),recojet->Mom().Pt(),corrjet->Eta(),corrjet->Pt(),0.,0.,smear,shift,ptSmeared);
    }
    //   printf("jet %i has ptsmeared = %f  \n",i,ptSmeared);
   
    // JEC scaling to correct for residual jet corrections
    if(scale) {
      Double_t factor=1;
      //     printf("im here 2 \n");
      //printf("%i       %f    %f    %f      %f \n",i, recojet->Mom().Eta(),recojet->Mom().Pt(),corrjet->Eta(),corrjet->Pt(),fGenJet->At(match)->Mom().Pt(),fGenJet->At(match)->Mom().Eta());
      if(TMath::Abs(recojet->Mom().Eta())<1.5) factor = 1.015;
      else if(TMath::Abs(recojet->Mom().Eta())<3) factor = 1.04;
      else factor = 1.15;
      ptSmeared  *= factor;
      eneSmeared *= factor;
      //printf("%i       %f      %f       %f       %f \n",i, recojet->Mom().Eta(),recojet->Mom().Pt(),factor,ptSmeared);
    }
    //printf("ptSmeared is %d  \n",ptSmeared);
    TLorentzVector thisJetSmeared;
    thisJetSmeared.SetPtEtaPhiE(ptSmeared,recojet->Mom().Eta(),recojet->Mom().Phi(),eneSmeared);
    
    TLorentzVector thisJetUnsmeared;
    thisJetUnsmeared.SetPtEtaPhiE(recojet->Mom().Pt(),recojet->Mom().Eta(),recojet->Mom().Phi(),recojet->Mom().E());
    
    if (recojet->Mom().Pt()>10 && TMath::Abs(recojet->Mom().Eta())<4.7) {
      jetSumSmeared   += thisJetSmeared;
      jetSumUnsmeared += thisJetUnsmeared;
    }
    
  }
  //  printf(" %d pfakt5 jets pass the cut and %d corrected jet pass the cut \n",uncorrJetColIdx,jetColIdx);
  Double_t px=met->Px()+jetSumUnsmeared.Px()-jetSumSmeared.Px();
  Double_t py=met->Py()+jetSumUnsmeared.Py()-jetSumSmeared.Py();

  met->SetMex(px);
  met->SetMey(py);
  return;
}

void PFMetCorrectionTools::shiftMet(Met *uncormet, Bool_t fIsData, Double_t spfMet) {
  //  Double_t spfMet=uncormet->SumEt();
  // correction for METx, METy bias
  Double_t px=0;
  Double_t py=0;
  // data
  if(fIsData){
    //    px =uncormet->Px() -0.00563109*spfMet+0.959742;
    //    py =uncormet->Py() +0.00586162*spfMet-0.540137; 
    px = uncormet->Px()-0.006239*spfMet+0.662; // change shift factor for 2012 correction
    py = uncormet->Py()+0.004613*spfMet-0.673;
    // MC
  }else{  
    px =uncormet->Px() +0.00135*spfMet-0.021;
    py =uncormet->Py() +0.00371*spfMet-0.826;
    
  }
  //  e = sqrt(px*px+py*py);
  
  //  Met.SetPxPyPzE(px,py,0,e);
  uncormet->SetMex(px);
  uncormet->SetMey(py);
  return;
}

