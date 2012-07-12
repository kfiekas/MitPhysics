#include "MitAna/DataTree/interface/PileupEnergyDensityFwd.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/Track.h"
#include "MitAna/DataTree/interface/Vertex.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#include "MitPhysics/Skim/interface/BaseH4lSkim.h"

using namespace mithep;

ClassImp(mithep::BaseH4lSkim)

BaseH4lSkim::BaseH4lSkim(const char *name, const char *title):
  BaseMod         (name,title),
  fTotal          (0),
  fSelected       (0),
  fMuonName       (Names::gkMuonBrn),
  fElectronName   (Names::gkElectronBrn),
  fPrimVtxName    (Names::gkPVBrn),
  fPfCandidateName(Names::gkPFCandidatesBrn),
  fTracksName     (Names::gkTrackBrn),
  fMuons          (0),
  fElectrons      (0),
  fPrimVerts      (0),
  fPfCandidates   (0),
  fTracks         (0),
  fBestPv         (0)
{
  // Constructor
}

//----------------------------------------------------------------------------------------
BaseH4lSkim::~BaseH4lSkim()
{
  // Destructor
}

//----------------------------------------------------------------------------------------
bool BaseH4lSkim::SetBestPv()
{
  //
  // Get primary vertices Assumes primary vertices are ordered by sum-pT^2 (as should be in CMSSW)
  // NOTE: if no PV is found from fitting tracks, the beamspot is used
  //
  bool hasGoodPv = false;

  // reset vertex choice (cool we are not responsible for the memory)
  fBestPv = 0;
  
  for (UInt_t i=0; i<fPrimVerts->GetEntries(); ++i) {
    const Vertex *pv = fPrimVerts->At(i);

    // Select best Pv for corrected d0; if no Pv passing cuts, the first Pv in the collection will
    // be used
    if(!pv->IsValid())
      continue;
    if (pv->Ndof()           < 4)
      continue;
    if (fabs(pv->Z())        > 24)
      continue;
    if (pv->Position().Rho() > 2)
      continue;

    hasGoodPv = kTRUE;
    fBestPv   = pv;
    break;
  }

  if (!fBestPv)
    fBestPv = fPrimVerts->At(0);

  return hasGoodPv;
}

//--------------------------------------------------------------------------------------------------
Float_t BaseH4lSkim::ComputePfMuonIso(const Muon *muon, const Double_t dRMax)
{
  const Double_t dRMin    = 0;
  const Double_t neuPtMin = 1.0;
  const Double_t dzMax    = 0.1;

  Double_t zLepton = (muon->BestTrk()) ? muon->BestTrk()->DzCorrected(*fBestPv) : 0.0;

  Float_t    iso = 0;
  for (UInt_t ipf=0; ipf<fPfCandidates->GetEntries(); ipf++) {
    const PFCandidate *pfcand = fPfCandidates->At(ipf);

    if (!pfcand->HasTrk() && (pfcand->Pt()<=neuPtMin))
      continue;  // pT cut on neutral particles

    // exclude THE muon
    if (pfcand->TrackerTrk() && muon->TrackerTrk() &&
        (pfcand->TrackerTrk() == muon->TrackerTrk()))
      continue;

    // dz cut
    Double_t dz = (pfcand->BestTrk()) ? fabs(pfcand->BestTrk()->DzCorrected(*fBestPv)-zLepton) : 0;
    if (dz >= dzMax)
      continue;

    // check iso cone
    Double_t dr = MathUtils::DeltaR(muon->Mom(), pfcand->Mom());
    if (dr<dRMax && dr>=dRMin)
      iso += pfcand->Pt();
  }

  return iso;
}
//----------------------------------------------------------------------------------------
Float_t BaseH4lSkim::ComputePfElecIso(const Electron *electron, const Double_t dRMax)
{
  const Double_t dRMin    = 0;
  const Double_t neuPtMin = 1.0;
  const Double_t dzMax    = 0.1;

  Double_t zLepton = (electron->BestTrk()) ? electron->BestTrk()->DzCorrected(*fBestPv) : 0.0;

  Float_t iso = 0;
  for (UInt_t ipf=0; ipf<fPfCandidates->GetEntries(); ipf++) {
    const PFCandidate *pfcand = fPfCandidates->At(ipf);

    if (!pfcand->HasTrk() && (pfcand->Pt()<=neuPtMin))
      continue;  // pT cut on neutral particles

    // dz cut
    Double_t dz = (pfcand->BestTrk()) ? fabs(pfcand->BestTrk()->DzCorrected(*fBestPv)-zLepton) : 0;
    if (dz >= dzMax)
      continue;

    // remove THE electron
    if (pfcand->TrackerTrk() && electron->TrackerTrk() &&
	(pfcand->TrackerTrk()==electron->TrackerTrk()))
      continue;
    if (pfcand->GsfTrk()     && electron->GsfTrk()     && 
	(pfcand->GsfTrk()==electron->GsfTrk()))
      continue;

    // check iso cone
    Double_t dr = MathUtils::DeltaR(electron->Mom(), pfcand->Mom());
    if (dr<dRMax && dr>=dRMin) {
      // eta-strip veto for photons
      if ((pfcand->PFType() == PFCandidate::eGamma) &&
	  fabs(electron->Eta() - pfcand->Eta()) < 0.025)
	continue;

      // Inner cone (one tower = dR < 0.07) veto for non-photon neutrals
      if (!pfcand->HasTrk() && (pfcand->PFType() == PFCandidate::eNeutralHadron) &&
         (MathUtils::DeltaR(electron->Mom(), pfcand->Mom()) < 0.07))
	continue;

      iso += pfcand->Pt();
    }
  }

  return iso;
}

//----------------------------------------------------------------------------------------
bool BaseH4lSkim::PassMuonPreselNoIp(const Muon *mu)
{
  if (mu->Pt() < 5)
    return false;
  if (fabs(mu->Eta()) > 2.4)
    return false;
  if (!(mu->IsTrackerMuon() || mu->IsGlobalMuon()))
    return false;

  return true;
}

//----------------------------------------------------------------------------------------
bool BaseH4lSkim::PassElecPreselNoIp(const Electron *ele)
{
  if (ele->Pt() < 5)
    return false;
  if (fabs(ele->Eta()) >= 2.5)
    return false;
  if (ele->CorrectedNExpectedHitsInner() > 1)
    return false;

  return true;
}
//----------------------------------------------------------------------------
unsigned BaseH4lSkim::makePFnoPUArray()
{
  for(UInt_t i = 0; i < fPfCandidates->GetEntries(); i++) {
    const PFCandidate *pf = fPfCandidates->At(i);
    assert(pf);

    if(pf->PFType() == PFCandidate::eHadron) {
      if(pf->HasTrackerTrk() && 
         fPrimVerts->At(0)->HasTrack(pf->TrackerTrk()) &&
         fPrimVerts->At(0)->TrackWeight(pf->TrackerTrk()) > 0) {

	pfNoPileUpflag.push_back(1);

      } else { 

        Bool_t vertexFound = kFALSE;
        const Vertex *closestVtx = 0;
        Double_t dzmin = 10000;
	
	for(UInt_t j = 0; j < fPrimVerts->GetEntries(); j++) {
	  const Vertex *vtx = fPrimVerts->At(j);
	  assert(vtx);
	  
	  if(pf->HasTrackerTrk() && 
	     vtx->HasTrack(pf->TrackerTrk()) &&
	     vtx->TrackWeight(pf->TrackerTrk()) > 0) {
	    vertexFound = kTRUE;
	    closestVtx = vtx;
	    break;
	  }
	  Double_t dz = fabs(pf->SourceVertex().Z() - vtx->Z());
	  if(dz < dzmin) {
	    closestVtx = vtx;
	    dzmin = dz;
	  }
	}

	if(vertexFound || closestVtx != fPrimVerts->At(0)) {
	  pfNoPileUpflag.push_back(0);
	} else {
	  pfNoPileUpflag.push_back(1);
	}
      }
    } else {
      pfNoPileUpflag.push_back(1);
    }
  }

  return pfNoPileUpflag.size();
}
//----------------------------------------------------------------------------------------
bool BaseH4lSkim::muon2012CutBasedIDTight(const mithep::Muon *mu)
{
  Float_t charged_iso = 0;
  
  for(unsigned k=0; k<fPfCandidates->GetEntries(); ++k) {
    const mithep::PFCandidate *pf = (mithep::PFCandidate*)((*fPfCandidates)[k]);
    Double_t dr = mithep::MathUtils::DeltaR(mu->Phi(),mu->Eta(), pf->Phi(), pf->Eta());
    if (dr > 0.4) continue;
    if (pf->Charge() != 0 && (pf->HasTrackerTrk()||pf->HasGsfTrk()) ) {
      if (abs(pf->PFType()) == mithep::PFCandidate::eElectron || abs(pf->PFType()) == mithep::PFCandidate::eMuon) continue;
      charged_iso += pf->Pt();
    }
  }
  double iso_sum = charged_iso;

  if(!mu->IsGlobalMuon())
    return false;
  if(mu->GlobalTrk()->RChi2() > 10)
    return false;
  if(mu->NValidHits() == 0 )
    return false;
  if(mu->NMatches() <= 1)
    return false;
  if(fabs(mu->TrackerTrk()->D0Corrected(*fBestPv)) > 0.2)
    return false;
  if(fabs(mu->TrackerTrk()->DzCorrected(*fBestPv)) > 0.5)
    return false;
  if(mu->BestTrk()->NPixelHits() == 0)
    return false;
  if(mu->NTrkLayersHit() <= 5)
    return false; 
  if(iso_sum/mu->Pt() > 0.12 )
    return false;

  return true;
}
//----------------------------------------------------------------------------------------
bool BaseH4lSkim::muon2012CutBasedIDTightForTagProbe(const mithep::Muon *mu)
{
  
  Float_t sumChargedHadronPt = 0;
  Float_t sumNeutralHadronPt = 0;
  Float_t sumPhotonPt = 0;
  Float_t sumPUPt = 0;

  //use delta beta correction isolation
  for(unsigned k=0; k<fPfCandidates->GetEntries(); ++k) {
    const mithep::PFCandidate *pf = (mithep::PFCandidate*)((*fPfCandidates)[k]);
    Double_t dr = mithep::MathUtils::DeltaR(mu->Phi(),mu->Eta(), pf->Phi(), pf->Eta());
    if(dr > 0.4) continue;
    if(pf->TrackerTrk() && mu->TrackerTrk() && pf->TrackerTrk() == mu->TrackerTrk()) continue;

    if(pfNoPileUpflag[k]){ //pf no pileup
      if(abs(pf->PFType()) == mithep::PFCandidate::eHadron)
	sumChargedHadronPt += pf->Pt();
      else if (pf->Pt() > 0.5 && abs(pf->PFType()) == mithep::PFCandidate::eNeutralHadron)
	sumNeutralHadronPt += pf->Pt();
      else if (pf->Pt() > 0.5 && abs(pf->PFType()) == mithep::PFCandidate::eGamma)
	sumPhotonPt += pf->Pt();
    }
    else {//pf pile up
      if(pf->Charge() != 0)
	sumPUPt += pf->Pt();
    }
  
  }
  double iso_sum =  sumChargedHadronPt+ max(0.,sumNeutralHadronPt+sumPhotonPt-0.5*sumPUPt);

  if(!mu->IsGlobalMuon())
    return false;
  if(mu->GlobalTrk()->RChi2() > 10)
    return false;
  if(mu->NValidHits() == 0 )
    return false;
  if(mu->NMatches() <= 1)
    return false;
  if(fabs(mu->TrackerTrk()->D0Corrected(*fBestPv)) > 0.2)
    return false;
  if(fabs(mu->TrackerTrk()->DzCorrected(*fBestPv)) > 0.5)
    return false;
  if(mu->BestTrk()->NPixelHits() == 0)
    return false;
  if(mu->NTrkLayersHit() <= 5)
    return false; 
  if(iso_sum/mu->Pt() > 0.12 )
    return false;
  return true;
}
//----------------------------------------------------------------------------------------
bool BaseH4lSkim::electron2012CutBasedIDMedium(const mithep::Electron *ele)
{

  Float_t charged_iso = 0;

  for(unsigned k=0; k<fPfCandidates->GetEntries(); ++k) {
    const mithep::PFCandidate *pf = (mithep::PFCandidate*)((*fPfCandidates)[k]);
    Double_t dr = mithep::MathUtils::DeltaR(ele->Phi(),ele->Eta(), pf->Phi(), pf->Eta());
    if (dr > 0.3) continue;
    if (pf->Charge() != 0 && (pf->HasTrackerTrk()||pf->HasGsfTrk()) ) {
      if (abs(pf->PFType()) == mithep::PFCandidate::eElectron || abs(pf->PFType()) == mithep::PFCandidate::eMuon) continue;
      if (fabs(ele->SCluster()->Eta()) > 1.479 && dr < 0.015) continue;
      charged_iso += pf->Pt();
    }
  }

  double iso_sum = charged_iso;

  if(fabs(ele->BestTrk()->D0Corrected(*fBestPv)) > 0.02)
    return false;	
  if(fabs(ele->BestTrk()->DzCorrected(*fBestPv)) > 0.1)
    return false;
  if(fabs(1 - ele->ESuperClusterOverP())/(ele->SCluster()->Et()*TMath::CosH(ele->SCluster()->Eta())) > 0.05)
    return false;
  if(ele->CorrectedNExpectedHitsInner() > 1)
    return false;

  // barrel
  if(fabs(ele->Eta()) < 1.5 ) {
    if(fabs(ele->DeltaEtaSuperClusterTrackAtVtx()) >  0.004)
      return false;
    if(fabs(ele->DeltaPhiSuperClusterTrackAtVtx()) > 0.06)
      return false;
    if(ele->CoviEtaiEta() > 0.01)
      return false;
    if(ele->HadronicOverEm() > 0.12)
      return false;
    if(iso_sum/ele->Pt() > 0.15)
      return false;	
  } else if(fabs(ele->Eta()) < 2.5){
    // endcap
    if(fabs(ele->DeltaEtaSuperClusterTrackAtVtx()) > 0.007)
      return false;	 
    if(fabs(ele->DeltaPhiSuperClusterTrackAtVtx()) > 0.03)
      return false;
    if(ele->CoviEtaiEta() > 0.03)
      return false;	
    if(ele->HadronicOverEm() > 0.10)
      return false;
    if(ele->Pt() > 20 && iso_sum/ele->Pt() > 0.15)
      return false;
    if(ele->Pt() < 20 && iso_sum/ele->Pt() > 0.10)
      return false;
  } else {
    return false;
  }

  return true;
}
//----------------------------------------------------------------------------------------
bool BaseH4lSkim::electron2012CutBasedIDMediumForTagProbe(const mithep::Electron *ele)
{
  
  Float_t charged_iso = 0;
  
  for(unsigned k=0; k<fPfCandidates->GetEntries(); ++k) {
    const mithep::PFCandidate *pf = (mithep::PFCandidate*)((*fPfCandidates)[k]);
    if( !(pfNoPileUpflag[k]) ) continue; 
    Double_t dr = mithep::MathUtils::DeltaR(ele->Phi(),ele->Eta(), pf->Phi(), pf->Eta());
    if (dr > 0.3) continue;
    if (pf->Charge() != 0 && (pf->HasTrackerTrk()||pf->HasGsfTrk()) ) {
      if (abs(pf->PFType()) == mithep::PFCandidate::eElectron || abs(pf->PFType()) == mithep::PFCandidate::eMuon) continue;
      if (fabs(ele->SCluster()->Eta()) > 1.479 && dr < 0.015) continue;
      charged_iso += pf->Pt();
    }
  }

  double iso_sum = charged_iso;

  if(fabs(ele->BestTrk()->D0Corrected(*fBestPv)) > 0.02)
    return false;	
  if(fabs(ele->BestTrk()->DzCorrected(*fBestPv)) > 0.1)
    return false;
  if(fabs(1 - ele->ESuperClusterOverP())/(ele->SCluster()->Et()*TMath::CosH(ele->SCluster()->Eta())) > 0.05)
    return false;
  if(ele->CorrectedNExpectedHitsInner() > 1)
    return false;

  // barrel
  if(fabs(ele->Eta()) < 1.5 ) {
    if(fabs(ele->DeltaEtaSuperClusterTrackAtVtx()) >  0.004)
      return false;
    if(fabs(ele->DeltaPhiSuperClusterTrackAtVtx()) > 0.06)
      return false;
    if(ele->CoviEtaiEta() > 0.01)
      return false;
    if(ele->HadronicOverEm() > 0.12)
      return false;
    if(iso_sum/ele->Pt() > 0.15)
      return false;	
  } else if(fabs(ele->Eta()) < 2.5){
    // endcap
    if(fabs(ele->DeltaEtaSuperClusterTrackAtVtx()) > 0.007)
      return false;	 
    if(fabs(ele->DeltaPhiSuperClusterTrackAtVtx()) > 0.03)
      return false;
    if(ele->CoviEtaiEta() > 0.03)
      return false;	
    if(ele->HadronicOverEm() > 0.10)
      return false;
    if(ele->Pt() > 20 && iso_sum/ele->Pt() > 0.15)
      return false;
    if(ele->Pt() < 20 && iso_sum/ele->Pt() > 0.10)
      return false;
  } else {
    return false;
  }

  return true;
}
//--------------------------------------------------------------------------------------------------
void BaseH4lSkim::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void BaseH4lSkim::SlaveBegin()
{
}

//--------------------------------------------------------------------------------------------------
void BaseH4lSkim::SlaveTerminate()
{
  cout << " BaseH4lSkim::SlaveTerminate - selected events: " <<  fSelected << endl
       << " BaseH4lSkim::SlaveTerminate - total events:    " <<  fTotal << endl; 
}

//--------------------------------------------------------------------------------------------------
void BaseH4lSkim::Terminate()
{
}

//--------------------------------------------------------------------------------------------------
void BaseH4lSkim::BeginRun()
{
}

//--------------------------------------------------------------------------------------------------
void BaseH4lSkim::EndRun()
{
}

//--------------------------------------------------------------------------------------------------
void BaseH4lSkim::Process()
{
}
