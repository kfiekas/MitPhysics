#include "MitPhysics/Utils/interface/JetTools.h"

ClassImp(mithep::JetTools)

using namespace mithep;
 
JetTools::JetTools()
{
   // Constructor
}

JetTools::~JetTools() 
{
  // Destructor.
}

//Remember to remove the signal from particles before inputting into the function
Double_t JetTools::NJettiness(const ParticleOArr *particles, const JetOArr *jets, bool UseQ, Double_t Y){
  if(particles->GetEntries() <= 0) return 0.0;

  Double_t fval = 0.0;
  Double_t fvalpart;

  for(int i=0;i<int(particles->GetEntries());i++){
    fvalpart = (particles->At(i)->Pt()) * TMath::Exp(-TMath::Abs(particles->At(i)->Eta()-Y)); 

    for(int j=0;j<int(jets->GetEntries());j++){
      fvalpart = TMath::Min(fvalpart,(jets->At(j)->Pt()) * 
                 (2 * TMath::CosH(TMath::Abs(jets->At(j)->Eta()-particles->At(i)->Eta()))
		- 2 * TMath::Cos(MathUtils::DeltaPhi(jets->At(j)->Phi(),particles->At(i)->Phi()))));
    }
    fval = fval + fvalpart;
  }

  if(UseQ == kTRUE) fval = fval / particles->At(0)->Pt();

  return fval;
}

Double_t JetTools::NJettiness(const TrackOArr *tracks, const JetOArr *jets, bool UseQ, Double_t Y){
  if(tracks->GetEntries() <= 0) return 0.0;

  Double_t fval = 0.0;
  Double_t fvalpart;

  for(int i=0;i<int(tracks->GetEntries());i++){
    fvalpart = (tracks->At(i)->Pt()) * TMath::Exp(-TMath::Abs(tracks->At(i)->Eta()-Y));     

    for(int j=0;j<int(jets->GetEntries());j++){
      fvalpart = TMath::Min(fvalpart,(jets->At(j)->Pt()) * 
                 (2 * TMath::CosH(TMath::Abs(jets->At(j)->Eta()-tracks->At(i)->Eta()))
		- 2 * TMath::Cos(MathUtils::DeltaPhi(jets->At(j)->Phi(),tracks->At(i)->Phi()))));
    }
    fval = fval + fvalpart;
  }

  if(UseQ == kTRUE) fval = fval / tracks->At(0)->Pt();
  
  return fval;
}

Double_t JetTools::NJettiness(const JetOArr *jetsS, const JetOArr *jets, bool UseQ, Double_t Y){
  if(jetsS->GetEntries() <= 0) return 0.0;

  Double_t fval = 0.0;
  Double_t fvalpart;

  for(int i=0;i<int(jetsS->GetEntries());i++){
    fvalpart = (jetsS->At(i)->Pt()) * TMath::Exp(-TMath::Abs(jetsS->At(i)->Eta()-Y));     

    for(int j=0;j<int(jets->GetEntries());j++){
      fvalpart = TMath::Min(fvalpart,(jets->At(j)->Pt()) * 
                 (2 * TMath::CosH(TMath::Abs(jets->At(j)->Eta()-jetsS->At(i)->Eta()))
		- 2 * TMath::Cos(MathUtils::DeltaPhi(jets->At(j)->Phi(),jetsS->At(i)->Phi()))));
    }
    fval = fval + fvalpart;
  }

  if(UseQ == kTRUE) fval = fval / jetsS->At(0)->Pt();
  
  return fval;
}

Double_t JetTools::NJettiness(const CaloTowerOArr *calos, const JetOArr *jets, bool UseQ, Double_t Y){
  if(calos->GetEntries() <= 0) return 0.0;

  Double_t fval = 0.0;
  Double_t fvalpart;

  for(int i=0;i<int(calos->GetEntries());i++){
    fvalpart = (calos->At(i)->Pt()) * TMath::Exp(-TMath::Abs(calos->At(i)->Eta()-Y));     

    for(int j=0;j<int(jets->GetEntries());j++){
      fvalpart = TMath::Min(fvalpart,(jets->At(j)->Pt()) * 
                 (2 * TMath::CosH(TMath::Abs(jets->At(j)->Eta()-calos->At(i)->Eta()))
		- 2 * TMath::Cos(MathUtils::DeltaPhi(jets->At(j)->Phi(),calos->At(i)->Phi()))));
    }
    fval = fval + fvalpart;
  }

  if(UseQ == kTRUE) fval = fval / calos->At(0)->Pt();
  
  return fval;
}

//M_r
Double_t JetTools::M_r(const ParticleOArr *particles){

  if(particles->GetEntries() < 2) return -999.;

  Double_t E0  = particles->At(0)->E();
  Double_t E1  = particles->At(1)->E();
  Double_t Pz0 = particles->At(0)->Pz();
  Double_t Pz1 = particles->At(1)->Pz();

  Double_t den =  TMath::Power(Pz0-Pz1, 2) - TMath::Power(E0-E1,2);
  if(den <= 0) return -100.;

  return 2.0*TMath::Sqrt(TMath::Power(E0*Pz1 - E1*Pz0, 2)/den);
}

//Beta_r
Double_t JetTools::Beta_r(const ParticleOArr *particles){

  if(particles->GetEntries() < 2) return -999.;

  Double_t E0  = particles->At(0)->E();
  Double_t E1  = particles->At(1)->E();
  Double_t Pz0 = particles->At(0)->Pz();
  Double_t Pz1 = particles->At(1)->Pz();

  return (E0-E1)/(Pz0-Pz1);
}

//M_r_t
Double_t JetTools::M_r_t(const ParticleOArr *particles, const Met *met){

  if(particles->GetEntries() < 2) return -999.;

  Double_t Pt0    = particles->At(0)->Pt();
  Double_t Pt1    = particles->At(1)->Pt();
  Double_t etmiss = met->Pt();

  Double_t Px0    = particles->At(0)->Px();
  Double_t Px1    = particles->At(1)->Px();
  Double_t metx   = met->Px();
  Double_t Py0    = particles->At(0)->Py();
  Double_t Py1    = particles->At(1)->Py();
  Double_t mety   = met->Py();

  return TMath::Sqrt(0.5*etmiss*(Pt0 + Pt1) - 0.5*(metx*(Px0 + Px1) + mety*(Py0 + Py1)));
}

//Razor
Double_t JetTools::Razor(const ParticleOArr *particles, const Met *met){
  if(particles->GetEntries() < 2) return -999.;

  Double_t mr  = M_r(particles);
  Double_t mrt = M_r_t(particles,met);
  
  if(mr != 0) return mrt/mr;
  
  return -999.;
}

//Cosine Omega
Double_t JetTools::CosineOmega(const ParticleOArr *particles){
  if(particles->GetEntries() < 2) return -999.;

  TLorentzVector v_L1(particles->At(0)->Px(),particles->At(0)->Py(),particles->At(0)->Pz(),particles->At(0)->E());
  TLorentzVector v_L2(particles->At(1)->Px(),particles->At(1)->Py(),particles->At(1)->Pz(),particles->At(1)->E());

  Double_t beta = (v_L1.P()-v_L2.P())/(v_L1.Pz()-v_L2.Pz());

  TVector3 B;
  B.SetXYZ(0.0,0.0,-1.0*beta);

  v_L1.Boost(B);
  v_L2.Boost(B);

  Double_t cosomega = v_L1.Vect().Dot(v_L2.Vect())/(v_L1.P()*v_L2.P());

  return cosomega;
}

//Transverse Higgs mass
Double_t JetTools::MtHiggs(const CompositeParticle *dilepton, const Met *met, int nsel){
  double mtHiggs = -999.0;
  double enell,enenn,enex,eney,mll,mnu;
  
  if     (nsel == 0){ // Use of Mt mass and mnu == mll
    enell = TMath::Sqrt(dilepton->Pt()*dilepton->Pt() + dilepton->Mt()*dilepton->Mt());
    enenn = TMath::Sqrt(met->Pt() *met->Pt()  + dilepton->Mt()*dilepton->Mt());
    enex  = dilepton->Px() + met->Px();
    eney  = dilepton->Py() + met->Py();
    mll   = dilepton->Mass();
    mnu   = mll;
  }
  else if(nsel == 1){ // Use of Mt mass and mnu == 0
    enell = TMath::Sqrt(dilepton->Pt()*dilepton->Pt() + dilepton->Mt()*dilepton->Mt());
    enenn = TMath::Sqrt(met->Pt() *met->Pt()  + 0.0*0.0);
    enex  = dilepton->Px() + met->Px();
    eney  = dilepton->Py() + met->Py();
    mll   = dilepton->Mass();
    mnu   = 0.0;
  }
  else if(nsel == 2){ // Use of M mass and mnu == mll
    enell = TMath::Sqrt(dilepton->Pt()*dilepton->Pt() + dilepton->Mass()*dilepton->Mass());
    enenn = TMath::Sqrt(met->Pt() *met->Pt()  + dilepton->Mass()*dilepton->Mass());
    enex  = dilepton->Px() + met->Px();
    eney  = dilepton->Py() + met->Py();
    mll   = dilepton->Mass();
    mnu   = mll;
  }
  else if(nsel == 3){ // Use of M mass and mnu == 0
    enell = TMath::Sqrt(dilepton->Pt()*dilepton->Pt() + dilepton->Mass()*dilepton->Mass());
    enenn = TMath::Sqrt(met->Pt() *met->Pt()  + 0.0*0.0);
    enex  = dilepton->Px() + met->Px();
    eney  = dilepton->Py() + met->Py();
    mll   = dilepton->Mass();
    mnu   = 0.0;
  }
  else {
    return -999.;
  }

  mtHiggs = mll*mll + mnu*mnu + 2.0*(enell*enenn - enex*enex - eney*eney);
  if(mtHiggs <= 0) mtHiggs = 0.0;
  else             mtHiggs = TMath::Sqrt(mtHiggs);

  return mtHiggs;
}
