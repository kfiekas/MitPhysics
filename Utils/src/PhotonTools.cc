// $Id: PhotonTools.cc,v 1.1 2011/04/12 22:14:21 bendavid Exp $

#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include <TFile.h>

ClassImp(mithep::PhotonTools)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
PhotonTools::PhotonTools()  
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
Bool_t PhotonTools::PassConversionId(const Photon *p, const DecayParticle *c) {

  if (!c) return kTRUE;
  
  ThreeVector dirconvsc = ThreeVector(p->SCluster()->Point()) - c->Position();
  Double_t deta = c->Eta()-dirconvsc.Eta();
  Double_t dphi = MathUtils::DeltaPhi(c->Phi(),dirconvsc.Phi());
  Double_t eoverp = p->SCluster()->Energy()/c->P();
  
  if (p->IsEB() && eoverp>2.0) return kFALSE;
  if (p->IsEE() && eoverp>3.0) return kFALSE;
  
  if (p->IsEE() && TMath::Abs(deta)>0.01) return kFALSE;
  if (p->IsEE() && TMath::Abs(dphi)>0.01) return kFALSE;

  return kTRUE;
    
}

//--------------------------------------------------------------------------------------------------
Bool_t PhotonTools::PassElectronVeto(const Photon *p, const ElectronCol *els) {

  Bool_t pass = kTRUE;
  for (UInt_t i=0; i<els->GetEntries(); ++i) {
    const Electron *e = els->At(i);
    if (e->SCluster()==p->SCluster() && e->GsfTrk()->NExpectedHitsInner()==0) {
      pass = kFALSE;
    }
  }
  
  return pass;
}

//--------------------------------------------------------------------------------------------------
Bool_t PhotonTools::PassElectronVetoConvRecovery(const Photon *p, const ElectronCol *els, const DecayParticleCol *conversions, const BaseVertex *v) {

  Bool_t pass = kTRUE;
  for (UInt_t i=0; i<els->GetEntries(); ++i) {
    const Electron *e = els->At(i);
    if (e->SCluster()==p->SCluster() && e->GsfTrk()->NExpectedHitsInner()==0 && ElectronTools::PassConversionFilter(e, conversions, 
                                                         v, 0, 1e-6, 2.0, kFALSE, kFALSE) ) {
      pass = kFALSE;
    }
  }
  
  return pass;
}

//--------------------------------------------------------------------------------------------------
Bool_t PhotonTools::PassTriggerMatching(const Photon *p, const TriggerObjectCol *trigobjs)
{
  
  for (UInt_t i=0; i<trigobjs->GetEntries(); ++i) {
    const TriggerObject *trigobj = trigobjs->At(i);
    if (trigobj->TriggerType()==TriggerObject::TriggerCluster || trigobj->TriggerType()==TriggerObject::TriggerElectron || trigobj->TriggerType()==TriggerObject::TriggerPhoton) {
      if (MathUtils::DeltaR(p->SCluster(),trigobj)<0.3) {
        return kTRUE;
      }
    }
  }
  
  return kFALSE;
  
  
}

//--------------------------------------------------------------------------------------------------
const DecayParticle *PhotonTools::MatchedConversion(const Photon *p, const DecayParticleCol *conversions, 
                                               const BaseVertex *vtx, Int_t nWrongHitsMax, Double_t probMin,
                                               Double_t lxyMin, Double_t dRMin) {
  
  const DecayParticle *match = 0;
  Double_t drsmallest = 999.;
  for (UInt_t i=0; i<conversions->GetEntries(); ++i) {
    const DecayParticle *c = conversions->At(i);
    ThreeVector dirconvsc = ThreeVector(p->SCluster()->Point()) - c->Position();
    Double_t dr = MathUtils::DeltaR(*c,dirconvsc);
    if (dr<dRMin && dr<drsmallest && c->Prob()>probMin && c->LxyCorrected(vtx)>lxyMin) {
      Int_t nhb1 = dynamic_cast<const StableData*>(c->DaughterDat(0))->NHitsBeforeVtx();
      Int_t nhb2 = dynamic_cast<const StableData*>(c->DaughterDat(1))->NHitsBeforeVtx();
      if (TMath::Max(nhb1,nhb2)<=nWrongHitsMax) {
        drsmallest = dr;
        match = c;
      }
    }
    
  }
  
  return match;
  
}

//--------------------------------------------------------------------------------------------------
const DecayParticle *PhotonTools::MatchedConversion(const Track *t, const DecayParticleCol *conversions, 
                                               const BaseVertex *vtx, Int_t nWrongHitsMax, Double_t probMin,
                                               Double_t lxyMin) {
  
  for (UInt_t i=0; i<conversions->GetEntries(); ++i) {
    const DecayParticle *c = conversions->At(i);
    if (c->Prob()>probMin && c->LxyCorrected(vtx)>lxyMin) {
      Int_t nhb1 = dynamic_cast<const StableData*>(c->DaughterDat(0))->NHitsBeforeVtx();
      Int_t nhb2 = dynamic_cast<const StableData*>(c->DaughterDat(1))->NHitsBeforeVtx();
      if (TMath::Max(nhb1,nhb2)<=nWrongHitsMax) {
        const Track *ct1 = dynamic_cast<const ChargedParticle*>(c->Daughter(0))->Trk();
        const Track *ct2 = dynamic_cast<const ChargedParticle*>(c->Daughter(1))->Trk();
        if (t==ct1 || t==ct2) return c;
      }
    }
    
  }
  
  return 0;
  
}

PhotonTools::DiphotonR9EtaCats PhotonTools::DiphotonR9EtaCat(const Photon *p1, const Photon *p2) {
  
  if (p1->IsEB() && p2->IsEB()) {
    if (p1->R9()>0.93 && p2->R9()>0.93) return kCat1;
    else return kCat2;
    
  }
  else {
    if (p1->R9()>0.93 && p2->R9()>0.93) return kCat3;
    else return kCat4; 
  }
  
}

PhotonTools::DiphotonR9EtaConversionCats PhotonTools::DiphotonR9EtaConversionCat(const Photon *p1, const Photon *p2, const DecayParticleCol *conversions, const BaseVertex *v) {
  
  const DecayParticle *conv1 = MatchedConversion(p1, conversions, v);
  const DecayParticle *conv2 = MatchedConversion(p2, conversions, v);
    
  if (p1->IsEB() && p2->IsEB()) {
    if (p1->R9()>0.93 && p2->R9()>0.93) return kNewCat1;
    else if (conv1||conv2) return kNewCat2;
    else return kNewCat3;
    
  }
  else {
    if (p1->R9()>0.93 && p2->R9()>0.93) return kNewCat4;
    else if (conv1||conv2) return kNewCat5;
    else return kNewCat6;
  }
  
}