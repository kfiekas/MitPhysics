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
  fTrigMask       (0),
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
    if (pv->NTracksFit()     < 9)
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
  //if (fabs(mu->Ip3dPVSignificance()) >= 100)
  //  return false;
  if (mu->Pt() < 5)
    return false;
  if (fabs(mu->Eta()) > 2.4)
    return false;
  if (!(mu->IsTrackerMuon() || mu->IsGlobalMuon()))
    return false;
  if (mu->IsoR03SumPt()/mu->Pt() >= 0.7)
    return false;

  return true;
}

//----------------------------------------------------------------------------------------
bool BaseH4lSkim::PassElecPreselNoIp(const Electron *ele)
{
  // if (fabs(ele->Ip3dPVSignificance()) >= 100)
  //   return false;
  if (ele->Pt() < 5)
    return false;
  if (fabs(ele->Eta()) >= 2.5)
    return false;
  if (ele->CorrectedNExpectedHitsInner() > 1)
    return false;
  if (ele->TrackIsolationDr03()/ele->Pt() >= 0.7)
    return false;

  return true;
}

//----------------------------------------------------------------------------------------
bool BaseH4lSkim::PassWwMuonSel(const Muon *mu)
{
  if (mu->Pt() < 5)
    return false;
  if (fabs(mu->Eta()) > 2.4)
    return false;
  if (mu->BestTrk()->PtErr()/mu->Pt() > 0.1)
    return false;
  if (fabs(mu->BestTrk()->DzCorrected(*fBestPv)) > 0.1)
    return false;
  Bool_t isGlobal  = (mu->IsGlobalMuon()) && (mu->BestTrk()->RChi2() < 10) &&
                      (mu->NMatches() > 1) && (mu->NValidHits() > 0);
  Bool_t isTracker = (mu->IsTrackerMuon()) &&
                      (mu->Quality().Quality(MuonQuality::TMLastStationTight));
  if (!isGlobal && !isTracker)
    return false;
  int ntrkhits = (mu->HasTrackerTrk()) ? mu->TrackerTrk()->NHits() : 0;
  if (ntrkhits < 10)
    return false;
  if (mu->BestTrk()->NPixelHits() < 1)
    return false;
  if (fabs(mu->BestTrk()->D0Corrected(*fBestPv)) > 0.02)
    return false;

  // note: this isn't really ww muon isolation
  float pfiso = ComputePfMuonIso(mu,0.3);

  if (pfiso > 0.2)
    return false;

  return true;
}
//----------------------------------------------------------------------------------------
bool BaseH4lSkim::PassVbtfWp80(const Electron *ele)
{
  if (ele->SCluster()->Et() < 5)
    return false;
  if (fabs(ele->SCluster()->Eta()) > 2.5)
    return false;

  float deta      = fabs(ele->DeltaEtaSuperClusterTrackAtVtx());
  float dphi      = fabs(ele->DeltaPhiSuperClusterTrackAtVtx());
  float sigieie   = ele->CoviEtaiEta();
  float HoE       = ele->HadronicOverEm();
  float missHits  = ele->CorrectedNExpectedHitsInner();

  float maxDeta, maxDphi, maxHoE, maxSigieie, maxMissHits;

  if (ele->IsEB()) {
    maxDeta = 0.004;
    maxDphi = 0.06 ;
    maxHoE = 0.04;
    maxSigieie = 0.01;
    maxMissHits = 0;
  } else {
    maxDeta = 0.007;
    maxDphi = 0.03;
    maxHoE = 0.025;
    maxSigieie = 0.03;
    maxMissHits = 0;
  }

  if (deta     > maxDeta)
    return false;
  if (dphi     > maxDphi)
    return false;
  if (sigieie  > maxSigieie)
    return false;
  if (HoE      > maxHoE)
    return false;
  if (missHits > maxMissHits)
    return false;

  // note: I'm not absolutely certain this isolation is correct
  // float det_iso = ele->ele->TrackIsolationDr04 + ele->
  cout << " WARNING: no isolation in here yet..." << endl;

  return true;
}
//--------------------------------------------------------------------------------------------------
bool BaseH4lSkim::PassElecTagSel(const Electron *ele)
{
  bool passID=true, passIso=true;

  double scet   = ele->SCluster()->Et();
  double sceta  = ele->SCluster()->Eta();

  if (scet < 5)
    passID = false;
  if (fabs(sceta) > 2.5)
    passID = false;

  if (fabs(ele->BestTrk()->D0Corrected(*fBestPv)) > 0.02) // cut on impact parameter
    passID = false;
  if (fabs(ele->BestTrk()->DzCorrected(*fBestPv)) > 0.1)  // cut on impact parameter in z
    passID = false;

  // conversion rejection
  double misshits = ele->CorrectedNExpectedHitsInner();
  if (misshits > 0)
    passID = false;

  // not implemented for bambu, so leave it out for the moment...
  //if (ele->isConv)
  //  passID = false;

  // barrel/endcap dependent requirments
  double pfiso  = ComputePfElecIso(ele,0.4);
  double sieie  = ele->CoviEtaiEta();
  double dphiin = ele->DeltaPhiSuperClusterTrackAtVtx();
  double detain = ele->DeltaPhiSuperClusterTrackAtVtx();
  double hoe    = ele->HadronicOverEm();

  if (fabs(sceta)<1.479) {
    // barrel
    if (pfiso        > 0.13*scet)
      passIso = false;
    if (sieie        > 0.01)
      passID = false;
    if (fabs(dphiin) > 0.06)
      passID = false;
    if (fabs(detain) > 0.004)
      passID = false;
    if (hoe          > 0.04)
      passID = false;
  }
  else {
    // endcap
    if (pfiso        > 0.09*(scet))
      passIso = false;
    if (sieie        > 0.03)
      passID = false;
    if (fabs(dphiin) > 0.03)
      passID = false;
    if (fabs(detain) > 0.007)
      passID = false;
    if (hoe          > 0.10)
      passID = false;
  }

  return (passID && passIso);
}

//--------------------------------------------------------------------------------------------------
//bool BaseH4lSkim::PassKsTagprobeMuonId(const TMuon *muon, const Double_t rho)
//{
//  if (muon->nTkLayers  < 6)        return kFALSE;
//  if (muon->nPixHits   < 1)        return kFALSE;
//  if (fabs(muon->d0)   > 0.2)      return kFALSE;
//  if (fabs(muon->dz)   > 0.5)      return kFALSE;
//  if (muon->muNchi2    > 10)       return kFALSE;
//  if (muon->nMatch     < 2)        return kFALSE;
//  if (muon->nValidHits < 1)        return kFALSE;
//  if (!(muon->typeBits & kGlobal)) return kFALSE;
//  if (!(muon->typeBits & kPfMuon)) return kFALSE;
//  Double_t iso = muon->pfChIso04 +
//              TMath::Max(muon->pfNeuIso04 + muon->pfGamIso04 - 0.5*(muon->puIso04),Double_t(0));
//  if (iso > 0.12*(muon->pt)) return kFALSE;
//  return kTRUE;
//}

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
