#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/Track.h"
#include "MitAna/DataTree/interface/Vertex.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/TriggerTable.h"
#include "MitAna/DataTree/interface/TriggerObjectsTable.h"
#include "MitAna/DataTree/interface/TriggerName.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/DataFormats/interface/Vect4M.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"

#include "MitPhysics/Skim/interface/H4lSkim.h"

using namespace mithep;

ClassImp(mithep::H4lSkim)
//--------------------------------------------------------------------------------------------------
H4lSkim::H4lSkim(const char *name, const char *title):
  BaseMod         (name,title),
  fSkipIfHLTFail  (kFALSE),
  flpt_gt_5       (0), 
  flpt_gt_10      (0), 
  fselected       (0), 
  fMuonName       (Names::gkMuonBrn),
  fElectronName   (Names::gkElectronBrn),
  fPrimVtxName    (Names::gkPVBrn),
  fPFCandidateName(Names::gkPFCandidatesBrn),
  fMuons          (0),
  fElectrons      (0),
  fPrimVerts      (0),
  fPFCandidates   (0)
{
  // Constructor
}

//--------------------------------------------------------------------------------------------------
H4lSkim::~H4lSkim()
{
  // Destructor
}	

//--------------------------------------------------------------------------------------------------      
void H4lSkim::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void H4lSkim::SlaveBegin()
{
  ReqBranch(fMuonName,       fMuons);
  ReqBranch(fElectronName,   fElectrons);
  ReqBranch(fPrimVtxName,    fPrimVerts);
  ReqBranch(fPFCandidateName,fPFCandidates);  
}

//--------------------------------------------------------------------------------------------------
void H4lSkim::SlaveTerminate()
{
}

//--------------------------------------------------------------------------------------------------
void H4lSkim::Terminate()
{
  cout << "selected events : " <<  fselected << endl; 
}

//--------------------------------------------------------------------------------------------------
void H4lSkim::BeginRun()
{
}

//--------------------------------------------------------------------------------------------------
void H4lSkim::EndRun()
{
}

//--------------------------------------------------------------------------------------------------
void H4lSkim::Process()
{
  // Load branches
  LoadBranch(fMuonName);
  LoadBranch(fElectronName);
  LoadBranch(fPrimVtxName);
  LoadBranch(fPFCandidateName);

  // trigger
  //
  // Get HLT info. Trigger objects can be matched by name to the corresponding trigger that passed.
  // note: TriggerName::Id() is bambu numbering scheme, fTriggerIdsv[itrig] is kevin's
  ULong64_t trigbits=0;
  if (HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);
    for (UInt_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if (!trigname)
	continue;
      // if event passed this trigger, set the trigger's bit in trigbits
      if (fTrigMask->At(trigname->Id()))
	trigbits |= fTriggerIdsv[itrig];
    }  
  }
  if (fSkipIfHLTFail && (trigbits==0))
    return;

  // Get primary vertices
  // Assumes primary vertices are ordered by sum-pT^2 (as should be in CMSSW)
  // NOTE: if no PV is found from fitting tracks, the beamspot is used
  const Vertex *bestPV = 0;
  Bool_t hasGoodPV = kFALSE;  
  for (UInt_t i=0; i<fPrimVerts->GetEntries(); ++i) {
    const Vertex *pv = fPrimVerts->At(i);
    
    // Select best PV for corrected d0; if no PV passing cuts, the first PV in the collection will
    // be used
    //if (!pv->IsValid()) continue;
    if (pv->NTracksFit()     < 9)
      continue;
    if (pv->Ndof()	     < 4)
      continue;
    if (fabs(pv->Z())	     > 24)
      continue;
    if (pv->Position().Rho() > 2)
      continue;    
    hasGoodPV = kTRUE;
    
    if (!bestPV)
      bestPV = pv;
  }
  if (!bestPV)
    bestPV = fPrimVerts->At(0);

  fVertex.SetPosition(bestPV->X(),bestPV->Y(),bestPV->Z());
  fVertex.SetErrors(bestPV->XErr(),bestPV->YErr(),bestPV->ZErr());

  flpt_gt_5=0;
  flpt_gt_10=0;

  // loop through muons
  for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i); 
    if (!mu->IsTrackerMuon())
      continue; 
    if (fabs(mu->Eta()) > 2.5)
      continue; 
    if (mu->Pt() >5)
      flpt_gt_5++; 
    if (mu->Pt() >10)
      flpt_gt_10++; 
  }
  //
  // loop through electrons.
  //
  ElectronTools eleTools;        // helper class for electron ID decisions
  for (UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    const Electron *ele = fElectrons->At(i);  
    if (fabs(ele->Eta()) > 2.5)
      continue;
    if (!eleTools.PassSpikeRemovalFilter(ele))
      continue;  // spike cleaning
    if (ele->Pt() >5)
      flpt_gt_5++;
    if (ele->Pt() >10)
      flpt_gt_10++;
  }

  if ( !(flpt_gt_5>=4 && flpt_gt_10>=2)) SkipEvent();
  else fselected++;
}

//--------------------------------------------------------------------------------------------------
ULong64_t H4lSkim::MatchHLT(const Double_t eta, const Double_t phi)
{
  ULong64_t bits = 0;
  
  const Double_t hltMatchR = 0.2;
  
  if (HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);
    for (UInt_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if (!trigname)
	continue;
      
      const TList *list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data());
      if (!list)
	continue;

      TIter iter(list->MakeIterator());
      const TriggerObject *to = dynamic_cast<const TriggerObject*>(iter.Next()); 
   
      while(to) {             
        if (to->IsHLT()) {          
	  
	  if (fTriggerObjNames1v[itrig].Length()>0 &&
	      fTriggerObjNames1v[itrig].CompareTo(to->ModuleName())==0) {
	    Bool_t match = kTRUE;
	    if (to->Pt() < fTriggerObjMinPt1v[itrig])
	      match=kFALSE;  // minimum pT threshold on trigger object
	    if (MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR)
	      match=kFALSE;  // eta-phi matching
	    if (match)
	      bits |= fTriggerObjIds1v[itrig];
	  }
	  
	  if (fTriggerObjNames2v[itrig].Length()>0 &&
	      fTriggerObjNames2v[itrig].CompareTo(to->ModuleName())==0) {
	    Bool_t match = kTRUE;
	    if (to->Pt() < fTriggerObjMinPt2v[itrig])
	      match=kFALSE;  // minimum pT threshold on trigger object
	    if (MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR)
	      match=kFALSE;  // eta-phi matching
	    if (match)
	      bits |= fTriggerObjIds2v[itrig];
	  }
	  
	  if (fTriggerObjNames1v[itrig].Length()==0 && fTriggerObjNames2v[itrig].Length()==0) {
	    Bool_t match = kTRUE;
	    if (to->Pt() < fTriggerObjMinPt1v[itrig])
	      match=kFALSE;  // minimum pT threshold on trigger object
	    if (MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR)
	      match=kFALSE;  // eta-phi matching
	    if (match)
	      bits |= fTriggerObjIds1v[itrig];
          }
        } 
        to = dynamic_cast<const TriggerObject*>(iter.Next());
      }    
    }
  }
  
  return bits;
}

ULong64_t H4lSkim::MatchHLT(const Double_t pt, const Double_t eta, const Double_t phi)
{
  ULong64_t bits = 0;
  
  const Double_t hltMatchR = 0.2;
  const Double_t hltMatchPtFrac = 1;
  
  if (HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);
    for (UInt_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if (!trigname) continue;
      
      const TList *list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data());
      if (!list) continue;
      TIter iter(list->MakeIterator());
      const TriggerObject *to = dynamic_cast<const TriggerObject*>(iter.Next()); 
    
      while(to) {         
        if (to->IsHLT()) {
	  if (fTriggerObjNames1v[itrig].Length()>0 &&
	      fTriggerObjNames1v[itrig].CompareTo(to->ModuleName())==0) {
	    Bool_t match = kTRUE;
	    if (to->Pt() < fTriggerObjMinPt1v[itrig])
	      match=kFALSE;  // minimum pT threshold on trigger object
	    if (MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR)
	      match=kFALSE;  // eta-phi matching
	    if (fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt()))
              match=kFALSE;  // pT matching
	    if (match)
	      bits |= fTriggerObjIds1v[itrig];
	  }
	  
	  if (fTriggerObjNames2v[itrig].Length()>0 &&
	      fTriggerObjNames2v[itrig].CompareTo(to->ModuleName())==0) {
	    Bool_t match = kTRUE;
	    if (to->Pt() < fTriggerObjMinPt2v[itrig])
	      match=kFALSE;  // minimum pT threshold on trigger object
	    if (MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR)
	      match=kFALSE;  // eta-phi matching
	    if (fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt()))
              match=kFALSE;  // pT matching
	    if (match)
	      bits |= fTriggerObjIds2v[itrig];
	  }
	  
	  if (fTriggerObjNames1v[itrig].Length()==0 && fTriggerObjNames2v[itrig].Length()==0) {
	    Bool_t match = kTRUE;
	    if (to->Pt() < fTriggerObjMinPt1v[itrig])
	      match=kFALSE;  // minimum pT threshold on trigger object
	    if (MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR)
	      match=kFALSE;  // eta-phi matching
	    if (fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt()))
              match=kFALSE;  // pT matching
	    if (match)
	      bits |= fTriggerObjIds1v[itrig];
          }
        }
        to = dynamic_cast<const TriggerObject*>(iter.Next());
      }    
    }
  }
  
  return bits;
}
      
//----------------------------------------------------------------------------------------
Bool_t H4lSkim::isMuFO(const Muon *mu) { 

  // Use tracker track when available
  const Track *muTrk=0;
  if (mu->HasTrackerTrk())
    muTrk = mu->TrackerTrk();
  else if (mu->HasGlobalTrk())
    muTrk = mu->GlobalTrk();
  else if (mu->HasStandaloneTrk())
    muTrk = mu->StandaloneTrk();
  assert(muTrk);                  
  
  Float_t dz = muTrk->DzCorrected(fVertex);
  if (fabs(dz) > 0.2)
    return kFALSE;   // needed to remove cosmics in HF sample

  // HF seems to not want tkquality, SS does
  if (muTrk->Pt()>20) { 
    if (!(mu->IsGlobalMuon() && mu->NSegments()>0))
      return kFALSE;
  }
  else {
    if (!(mu->IsGlobalMuon()  && mu->NSegments()>0) &&
       !(mu->IsTrackerMuon() && (mu->Quality().QualityMask().Mask() &
				 MuonQuality::TMOneStationLoose)     ))
      return kFALSE;
  }

  return kTRUE;
}
//--------------------------------------------------------------------------------------------------
Bool_t H4lSkim::isLooseEleFO(const Electron *ele)
{
  if (fabs(ele->BestTrk()->DzCorrected(fVertex)) > 0.1)     return kFALSE;

  float isocut=-1;
  float dIso = 0.7;
  if (ele->IsEB()) {
    if (ele->Pt() < 20) isocut = 0.3  + dIso;
    else               isocut = 0.45 + dIso;
  } else {
    if (ele->Pt() < 20) isocut = 0.3 + dIso;
    else               isocut = 0.5 + dIso;
  }

  float pfIso04 = computePFEleIso(ele,0.4);
  if (pfIso04/ele->Pt() > isocut) return kFALSE;

  float hecut=-1;
  if (fabs(ele->SCluster()->Eta())<1.479) hecut = 0.12;
  else                                   hecut = 0.10;

  if (ele->HadronicOverEm() > hecut)      return kFALSE;
    
  return kTRUE;
}
//----------------------------------------------------------------------------------------
Float_t H4lSkim::computePFEleIso(const Electron *electron, const Double_t dRMax)
{
  const Double_t dRMin    = 0;
  const Double_t neuPtMin = 1.0;
  const Double_t dzMax    = 0.1;
    
  Double_t zLepton = (electron->BestTrk()) ? electron->BestTrk()->DzCorrected(fVertex) : 0.0;
  
  Float_t iso = 0;
  for (UInt_t ipf=0; ipf<fPFCandidates->GetEntries(); ipf++) {
    const PFCandidate *pfcand = fPFCandidates->At(ipf);
    
    if (!pfcand->HasTrk() && (pfcand->Pt()<=neuPtMin))
      continue;  // pT cut on neutral particles
    
    // dz cut
    Double_t dz = (pfcand->BestTrk()) ? fabs(pfcand->BestTrk()->DzCorrected(fVertex) - zLepton) : 0;
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
      if ((pfcand->PFType() == PFCandidate::eGamma) && fabs(electron->Eta()-pfcand->Eta()) < 0.025)
	continue;
      
      // Inner cone (one tower = dR < 0.07) veto for non-photon neutrals
      if (!pfcand->HasTrk() && (pfcand->PFType() == PFCandidate::eNeutralHadron) && 
	  (MathUtils::DeltaR(electron->Mom(), pfcand->Mom()) < 0.07)) continue;
      
      iso += pfcand->Pt();
    }
  }
  
  return iso;
}
