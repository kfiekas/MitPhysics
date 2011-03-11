// $Id: MuonIDMod.cc,v 1.40 2011/03/07 12:45:52 ceballos Exp $

#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::MuonIDMod)

//--------------------------------------------------------------------------------------------------
  MuonIDMod::MuonIDMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fMuonBranchName(Names::gkMuonBrn),
  fCleanMuonsName(ModNames::gkCleanMuonsName),  
  fNonIsolatedMuonsName("random"),  
  fNonIsolatedElectronsName("random"),  
  fVertexName(ModNames::gkGoodVertexesName),
  fTrackName(Names::gkTrackBrn),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fMuonIDType("WWMuId"),
  fMuonIsoType("TrackCaloSliding"),
  fMuonClassType("Global"),  
  fTrackIsolationCut(3.0),
  fCaloIsolationCut(3.0),
  fCombIsolationCut(0.15),
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
  fTracks(0),
  fPFCandidates(0),
  fNonIsolatedMuons(0),
  fNonIsolatedElectrons(0),
  fPileupEnergyDensityName(Names::gkPileupEnergyDensityBrn),
  fPileupEnergyDensity(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void MuonIDMod::Process()
{
  // Process entries of the tree. 

  if(fMuIsoType != kPFIsoNoL) {
    LoadEventObject(fMuonBranchName, fMuons);
  }
  else {
    fMuons = GetObjThisEvt<MuonOArr>(fMuonBranchName);
  }
  LoadEventObject(fTrackName, fTracks);
  LoadEventObject(fPFCandidatesName, fPFCandidates);
  if(fMuIsoType == kTrackCaloSliding) {
    LoadEventObject(fPileupEnergyDensityName, fPileupEnergyDensity);
  }
  MuonOArr *CleanMuons = new MuonOArr;
  CleanMuons->SetName(fCleanMuonsName);

  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);

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

    Double_t RChi2 = 0.0;
    if     (mu->HasGlobalTrk()) {
      RChi2 = mu->GlobalTrk()->Chi2()/mu->GlobalTrk()->Ndof();
    }
    else if(mu->BestTrk() != 0){
      RChi2 = mu->BestTrk()->Chi2()/mu->BestTrk()->Ndof();
    }
    Bool_t idpass = kFALSE;
    switch (fMuIDType) {
      case kWMuId:
        idpass = mu->BestTrk() != 0 &&
	         mu->BestTrk()->NHits() > 10 &&
		 RChi2 < 10.0 &&
		(mu->NSegments() > 1 || mu->NMatches() > 1) &&
		 mu->BestTrk()->NPixelHits() > 0 &&
		 mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight);
        break;
      case kZMuId:
        idpass = mu->BestTrk() != 0 &&
	         mu->BestTrk()->NHits() > 10 &&
		(mu->NSegments() > 1 || mu->NMatches() > 1) &&
		 mu->BestTrk()->NPixelHits() > 0 &&
		 mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight);
        break;
      case kLoose:
        idpass = mu->BestTrk() != 0 &&
	         mu->Quality().Quality(MuonQuality::TMOneStationLoose) &&
                 mu->Quality().Quality(MuonQuality::TM2DCompatibilityLoose) &&
		 mu->BestTrk()->NHits() > 10 &&
		 RChi2 < 10.0 &&
		 mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight);
        break;
      case kTight:
        idpass = mu->BestTrk() !=  0 &&
	         mu->Quality().Quality(MuonQuality::TMOneStationTight) &&
                 mu->Quality().Quality(MuonQuality::TM2DCompatibilityTight) &&
		 mu->BestTrk()->NHits() > 10 &&
		 RChi2 < 10.0 &&
		 mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight);
        break;
      case kWWMuId:
        idpass = mu->BestTrk() != 0 &&
	         mu->BestTrk()->NHits() > 10 &&
		 RChi2 < 10.0 &&
		(mu->NSegments() > 1 || mu->NMatches() > 1) &&
		 mu->BestTrk()->NPixelHits() > 0 &&
		 mu->Quality().Quality(MuonQuality::GlobalMuonPromptTight) &&
		 mu->BestTrk()->PtErr()/mu->BestTrk()->Pt() < 0.1;
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
        isocut = (1.0 * mu->IsoR03SumPt() +
	          1.0 * mu->IsoR03EmEt()  + 
		  1.0 * mu->IsoR03HadEt() < fCombIsolationCut);
        break;
      case kTrackCaloSliding:
        { 
          //Double_t beta = IsolationTools::BetaM(fTracks, mu, fVertices->At(0), 0.0, 0.2, 0.3, 0.02); 
          //if(beta == 0) beta = 1.0;
          const PileupEnergyDensity *rho =  fPileupEnergyDensity->At(0);
          Double_t totalIso =  mu->IsoR03SumPt() + TMath::Max(mu->IsoR03EmEt() + mu->IsoR03HadEt() - rho->Rho() * TMath::Pi() * 0.3 * 0.3, 0.0);
          if (totalIso < (mu->Pt()*fCombIsolationCut) )
            isocut = kTRUE;

	  if     (fReverseIsoCut == kTRUE &&
    	          isocut == kFALSE && totalIso < 10)
	    isocut = kTRUE;
          else if(fReverseIsoCut == kTRUE)
	    isocut = kFALSE;
	}
        break;
      case kTrackCaloSlidingNoCorrection:
        { 
          Double_t totalIso =  1.0 * mu->IsoR03SumPt() + 
                               1.0 * mu->IsoR03EmEt()  + 
                               1.0 * mu->IsoR03HadEt();
          if (totalIso < (mu->Pt()*fCombIsolationCut) )
            isocut = kTRUE;
	}
        break;
      case kPFIso:
        {
          Double_t beta = IsolationTools::BetaM(fTracks, mu, fVertices->At(0), 0.0, 0.2, 0.3, 0.02); 
          if(beta == 0) beta = 1.0;
          Double_t totalIso =  IsolationTools::PFMuonIsolation(mu, fPFCandidates, fVertices->At(0), 0.2, 0.5, 0.3, 0.02, 0, beta, fNonIsolatedMuons, fNonIsolatedElectrons);
          if (totalIso < (mu->Pt()*fCombIsolationCut) )
            isocut = kTRUE;
	}
        break;
      case kPFIsoNoL:
        {
          fNonIsolatedMuons     = GetObjThisEvt<MuonCol>(fNonIsolatedMuonsName);
          fNonIsolatedElectrons = GetObjThisEvt<ElectronCol>(fNonIsolatedElectronsName);

          Double_t beta = IsolationTools::BetaM(fTracks, mu, fVertices->At(0), 0.0, 0.2, 0.3, 0.02); 
          if(beta == 0) beta = 1.0;
          Double_t totalIso =  IsolationTools::PFMuonIsolation(mu, fPFCandidates, fVertices->At(0), 0.2, 0.5, 0.3, 0.02, 3, beta, fNonIsolatedMuons, fNonIsolatedElectrons);
          if (totalIso < (mu->Pt()*fCombIsolationCut) )
            isocut = kTRUE;
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
      Bool_t passD0cut = MuonTools::PassD0Cut(mu, fVertices, fD0Cut);
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

   // In this case we cannot have a branch
  if (fMuonIsoType.CompareTo("PFIsoNoL") != 0) {
    ReqEventObject(fMuonBranchName, fMuons, kTRUE);
  }
  ReqEventObject(fTrackName, fTracks, kTRUE);
  ReqEventObject(fPFCandidatesName, fPFCandidates, kTRUE);
  if (fMuonIsoType.CompareTo("TrackCaloSliding") == 0) {
    ReqEventObject(fPileupEnergyDensityName, fPileupEnergyDensity, kTRUE);
  }

  if (fMuonIDType.CompareTo("WMuId") == 0) 
    fMuIDType = kWMuId;
  else if (fMuonIDType.CompareTo("ZMuId") == 0) 
    fMuIDType = kZMuId;
  else if (fMuonIDType.CompareTo("Tight") == 0) 
    fMuIDType = kTight;
  else if (fMuonIDType.CompareTo("Loose") == 0) 
    fMuIDType = kLoose;
  else if (fMuonIDType.CompareTo("WWMuId") == 0) 
    fMuIDType = kWWMuId;
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
  else if (fMuonIsoType.CompareTo("TrackCaloSlidingNoCorrection") == 0)
    fMuIsoType = kTrackCaloSlidingNoCorrection;
  else if (fMuonIsoType.CompareTo("PFIso") == 0)
    fMuIsoType = kPFIso;
  else if (fMuonIsoType.CompareTo("PFIsoNoL") == 0)
    fMuIsoType = kPFIsoNoL;
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
