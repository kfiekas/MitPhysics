// $Id: MuonIDMod.cc,v 1.28 2010/05/12 19:06:53 ceballos Exp $

#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::MuonIDMod)

//--------------------------------------------------------------------------------------------------
  MuonIDMod::MuonIDMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fMuonBranchName(Names::gkMuonBrn),
  fCleanMuonsName(ModNames::gkCleanMuonsName),  
  fVertexName("PrimaryVertexes"),
  fMuonIDType("Loose"),
  fMuonIsoType("TrackCaloSliding"),  
  fMuonClassType("Global"),  
  fTrackIsolationCut(3.0),
  fCaloIsolationCut(3.0),
  fCombIsolationCut(5.0),
  fMuonPtMin(10),
  fApplyD0Cut(kTRUE),
  fD0Cut(0.025),
  fReverseIsoCut(kFALSE),
  fReverseD0Cut(kFALSE),
  fMuIDType(kIdUndef),
  fMuIsoType(kIsoUndef),
  fMuClassType(kClassUndef),
  fMuons(0),
  fVertices(0),
  fMuonTools(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void MuonIDMod::Process()
{
  // Process entries of the tree. 

  LoadEventObject(fMuonBranchName, fMuons);
  if (fApplyD0Cut) {
    LoadEventObject(fVertexName,     fVertices);
  }

  MuonOArr *CleanMuons = new MuonOArr;
  CleanMuons->SetName(fCleanMuonsName);

  for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i);

    Bool_t pass = kFALSE;
    Double_t pt = 0; // make sure pt is taken from the correct track!
    switch (fMuClassType) {
      case kAll:
        pass = kTRUE;
        if (mu->HasTrk())
          pt = mu->Pt();
        break;
      case kGlobal:
        pass = mu->HasGlobalTrk() && mu->IsTrackerMuon() &&
	       mu->Quality().Quality(MuonQuality::TrackerMuonArbitrated);
        if (pass) 
          pt = mu->GlobalTrk()->Pt();
        break;
      case kSta:
        pass = mu->HasStandaloneTrk();
        if (pass) 
          pt = mu->StandaloneTrk()->Pt();
        break;
      case kTrackerMuon:
        pass = mu->HasTrackerTrk() && mu->IsTrackerMuon() &&
	       mu->Quality().Quality(MuonQuality::TrackerMuonArbitrated);
        if (pass) 
          pt = mu->TrackerTrk()->Pt();
        break;
      case kCaloMuon:
        pass = mu->HasTrackerTrk() && mu->IsCaloMuon();
        if (pass) 
          pt = mu->TrackerTrk()->Pt();
        break;
      case kTrackerBased:
        pass = mu->HasTrackerTrk();
        if (pass) 
          pt = mu->TrackerTrk()->Pt();
        break;
      default:
        break;
    }

    if (!pass)
      continue;

    if (pt <= fMuonPtMin) 
      continue;

    Bool_t idpass = kFALSE;
    switch (fMuIDType) {
      case kLoose:
        idpass = mu->Quality().Quality(MuonQuality::TMOneStationLoose) &&
                 mu->Quality().Quality(MuonQuality::TM2DCompatibilityLoose) &&
		 mu->BestTrk()->NHits() > 10 &&
		 mu->BestTrk()->Chi2()/mu->BestTrk()->Ndof() < 10 &&
		 mu->NSegments() > 0;
        break;
      case kTight:
        idpass = mu->Quality().Quality(MuonQuality::TMOneStationTight) &&
                 mu->Quality().Quality(MuonQuality::TM2DCompatibilityTight) &&
		 mu->BestTrk()->NHits() > 10 &&
		 mu->BestTrk()->Chi2()/mu->BestTrk()->Ndof() < 10 &&
		 mu->NSegments() > 0;
        break;
      case kMinimal:
        idpass = mu->BestTrk()->NHits() > 10 &&
		 mu->BestTrk()->Chi2()/mu->BestTrk()->Ndof() < 10 &&
		 mu->NSegments() > 0;
        break;
      case kNoId:
        idpass = kTRUE;
        break;
      default:
        break;
    }

    if (!idpass)
      continue;

    Bool_t isocut = kFALSE;
    switch (fMuIsoType) {
      case kTrackCalo:
        isocut = (mu->IsoR03SumPt() < fTrackIsolationCut) &&
          (mu->IsoR03EmEt() + mu->IsoR03HadEt() < fCaloIsolationCut);
        break;
      case kTrackCaloCombined:
        isocut = (1.0 * mu->IsoR03SumPt() + 1.0 * mu->IsoR03EmEt() + 
                  1.0 * mu->IsoR03HadEt() < fCombIsolationCut);
        break;
      case kTrackCaloSliding:
        { 
          Double_t totalIso = 1.0 * mu->IsoR03SumPt() + 
                              1.0 * mu->IsoR03EmEt() + 
                              1.0 * mu->IsoR03HadEt();
          if (totalIso < (mu->Pt()*0.15) )
            isocut = kTRUE;

	  if     (fReverseIsoCut == kTRUE &&
    	          isocut == kFALSE && totalIso < 10)
	    isocut = kTRUE;
          else if(fReverseIsoCut == kTRUE)
	    isocut = kFALSE;
	}
        break;
      case kNoIso:
        isocut = kTRUE;
        break;
      case kCustomIso:
      default:
        break;
    }

    if (isocut == kFALSE)
      continue;

    if (fApplyD0Cut) {
      Bool_t d0cut = kFALSE;
      const Track *mt = mu->BestTrk();
      if (!mt)
        continue;
      Double_t d0_real = 1e30;
      for(UInt_t i0 = 0; i0 < fVertices->GetEntries(); i0++) {
        Double_t pD0 = mt->D0Corrected(*fVertices->At(i0));
        if(TMath::Abs(pD0) < TMath::Abs(d0_real)) 
          d0_real = TMath::Abs(pD0);
      }
      if(d0_real < fD0Cut) d0cut = kTRUE;

      if     (fReverseD0Cut == kTRUE &&
              d0cut == kFALSE && d0_real < 0.05)
	d0cut = kTRUE;
      else if(fReverseD0Cut == kTRUE)
	d0cut = kFALSE;

      if (d0cut == kFALSE)
        continue;
    }

    // add good muon
    CleanMuons->Add(mu);
  }

  // sort according to pt
  CleanMuons->Sort();

  // add objects for other modules to use
  AddObjThisEvt(CleanMuons);  
}

//--------------------------------------------------------------------------------------------------
void MuonIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the muon collection branch.

  ReqEventObject(fMuonBranchName, fMuons, kTRUE);

  if (fApplyD0Cut)
    ReqEventObject(fVertexName, fVertices, kTRUE);

  fMuonTools = new MuonTools;

  if (fMuonIDType.CompareTo("Tight") == 0) 
    fMuIDType = kTight;
  else if (fMuonIDType.CompareTo("Loose") == 0) 
    fMuIDType = kLoose;
  else if (fMuonIDType.CompareTo("Minimal") == 0) 
    fMuIDType = kMinimal;
  else if (fMuonIDType.CompareTo("NoId") == 0) 
    fMuIDType = kNoId;
  else if (fMuonIDType.CompareTo("Custom") == 0) {
    fMuIDType = kCustomId;
    SendError(kWarning, "SlaveBegin",
              "Custom muon identification is not yet implemented.");
  } else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified muon identification %s is not defined.",
              fMuonIDType.Data());
    return;
  }

  if (fMuonIsoType.CompareTo("TrackCalo") == 0)
    fMuIsoType = kTrackCalo;
  else if (fMuonIsoType.CompareTo("TrackCaloCombined") == 0)
    fMuIsoType = kTrackCaloCombined;
  else if (fMuonIsoType.CompareTo("TrackCaloSliding") == 0)
    fMuIsoType = kTrackCaloSliding;
  else if (fMuonIsoType.CompareTo("NoIso") == 0)
    fMuIsoType = kNoIso;
  else if (fMuonIsoType.CompareTo("Custom") == 0) {
    fMuIsoType = kCustomIso;
    SendError(kWarning, "SlaveBegin",
              "Custom muon isolation is not yet implemented.");
  } else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified muon isolation %s is not defined.",
              fMuonIsoType.Data());
    return;
  }

  if (fMuonClassType.CompareTo("All") == 0) 
    fMuClassType = kAll;
  else if (fMuonClassType.CompareTo("Global") == 0) 
    fMuClassType = kGlobal;
  else if (fMuonClassType.CompareTo("Standalone") == 0) 
    fMuClassType = kSta;
  else if (fMuonClassType.CompareTo("TrackerMuon") == 0) 
    fMuClassType = kTrackerMuon;
  else if (fMuonClassType.CompareTo("CaloMuon") == 0) 
    fMuClassType = kCaloMuon;
  else if (fMuonClassType.CompareTo("TrackerBased") == 0) 
    fMuClassType = kTrackerBased;
  else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified muon class %s is not defined.",
              fMuonClassType.Data());
    return;
  }
}
