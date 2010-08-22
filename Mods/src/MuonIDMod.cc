// $Id: MuonIDMod.cc,v 1.32 2010/08/19 14:37:17 ceballos Exp $

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
  fMuonIDType("Minimal"),
  fMuonIsoType("TrackCaloSliding"),  
  fMuonClassType("Global"),  
  fTrackIsolationCut(3.0),
  fCaloIsolationCut(3.0),
  fCombIsolationCut(5.0),
  fMuonPtMin(10),
  fApplyD0Cut(kTRUE),
  fD0Cut(0.020),
  fEtaCut(2.4),
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
    Double_t pt = 0;  // make sure pt is taken from the correct track!
    Double_t eta = 0; // make sure eta is taken from the correct track!
    switch (fMuClassType) {
      case kAll:
        pass = kTRUE;
        if (mu->HasTrk()) {
          pt  = mu->Pt();
	  eta = TMath::Abs(mu->Eta());
	}
        break;
      case kGlobal:
        pass = mu->HasGlobalTrk() && mu->IsTrackerMuon();
        if (pass && mu->TrackerTrk()) {
          pt  = mu->TrackerTrk()->Pt();
	  eta = TMath::Abs(mu->TrackerTrk()->Eta());
        }
	else {
          pt  = mu->Pt();
	  eta = TMath::Abs(mu->Eta());
	}
	break;
      case kSta:
        pass = mu->HasStandaloneTrk();
        if (pass) {
          pt  = mu->StandaloneTrk()->Pt();
          eta = TMath::Abs(mu->StandaloneTrk()->Eta());
	}
        break;
      case kTrackerMuon:
        pass = mu->HasTrackerTrk() && mu->IsTrackerMuon() &&
	       mu->Quality().Quality(MuonQuality::TrackerMuonArbitrated);
        if (pass) {
          pt  = mu->TrackerTrk()->Pt();
          eta = TMath::Abs(mu->TrackerTrk()->Eta());
	}
        break;
      case kCaloMuon:
        pass = mu->HasTrackerTrk() && mu->IsCaloMuon();
        if (pass) {
          pt  = mu->TrackerTrk()->Pt();
          eta = TMath::Abs(mu->TrackerTrk()->Eta());
	}
        break;
      case kTrackerBased:
        pass = mu->HasTrackerTrk();
        if (pass) {
          pt  = mu->TrackerTrk()->Pt();
          eta = TMath::Abs(mu->TrackerTrk()->Eta());
	}
        break;
      default:
        break;
    }

    if (!pass)
      continue;

    if (pt <= fMuonPtMin) 
      continue;

    if (eta >= fEtaCut) 
      continue;

    Bool_t idpass = kFALSE;
    switch (fMuIDType) {
      case kLoose:
        idpass = mu->BestTrk() != 0 &&
	         mu->Quality().Quality(MuonQuality::TMOneStationLoose) &&
                 mu->Quality().Quality(MuonQuality::TM2DCompatibilityLoose) &&
		 mu->BestTrk()->NHits() > 10 &&
		 mu->BestTrk()->Chi2()/mu->BestTrk()->Ndof() < 10 &&
		 //mu->NSegments() > 0 &&
		 mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight);
        break;
      case kTight:
        idpass = mu->BestTrk() !=  0 &&
	         mu->Quality().Quality(MuonQuality::TMOneStationTight) &&
                 mu->Quality().Quality(MuonQuality::TM2DCompatibilityTight) &&
		 mu->BestTrk()->NHits() > 10 &&
		 mu->BestTrk()->Chi2()/mu->BestTrk()->Ndof() < 10 &&
		 //mu->NSegments() > 0 &&
		 mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight);
        break;
      case kMinimal:
        idpass = mu->BestTrk() != 0 &&
	         mu->BestTrk()->NHits() > 10 &&
		 mu->BestTrk()->Chi2()/mu->BestTrk()->Ndof() < 10 &&
		 //mu->NSegments() > 0 &&
		 mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight);
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
      Bool_t passD0cut = MuonTools::PassD0Cut(mu, fVertices, fD0Cut, fReverseD0Cut);
      if (!passD0cut)
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
