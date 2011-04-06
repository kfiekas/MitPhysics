// $Id: PhotonIDMod.cc,v 1.17 2011/04/06 18:03:48 fabstoec Exp $

#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"

using namespace mithep;

ClassImp(mithep::PhotonIDMod)

//--------------------------------------------------------------------------------------------------
PhotonIDMod::PhotonIDMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPhotonBranchName  (Names::gkPhotonBrn),
  fGoodPhotonsName   (ModNames::gkGoodPhotonsName),
  fTrackBranchName   (Names::gkTrackBrn),
  fBeamspotBranchName(Names::gkBeamSpotBrn),
  fPileUpDenName     (Names::gkPileupEnergyDensityBrn),
  fPhotonIDType("Custom"),
  fPhotonIsoType("Custom"),
  fPhotonPtMin(15.0),
  fHadOverEmMax(0.02),
  fApplySpikeRemoval(kFALSE),
  fApplyPixelSeed(kTRUE),
  fPhotonR9Min(0.5),
  fPhIdType(kIdUndef),
  fPhIsoType(kIsoUndef),
  fFiduciality(kTRUE),
  fEtaWidthEB(0.01),
  fEtaWidthEE(0.028),
  fAbsEtaMax(2.5),
  fApplyR9Min(kFALSE),
  fEffAreaEcal(0.139),
  fEffAreaHcal(0.056),
  fEffAreaTrack(0.296),
  fPhotons(0),
  fTracks(0),
  fBeamspots(0),
  fPileUpDen(0)
  
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void PhotonIDMod::Process()
{
  // Process entries of the tree. 

  LoadEventObject(fPhotonBranchName,   fPhotons);
  LoadEventObject(fTrackBranchName,    fTracks);
  LoadEventObject(fBeamspotBranchName, fBeamspots);
  LoadEventObject(fPileUpDenName,      fPileUpDen);

  const BaseVertex *bsp = NULL;
  if(fBeamspots->GetEntries() > 0)
    bsp = dynamic_cast<const BaseVertex*>(fBeamspots->At(0));

  Double_t _tRho = -1.;
  if(fPileUpDen->GetEntries() > 0)
    _tRho = (Double_t) fPileUpDen->At(0)->Rho();
  
  PhotonOArr *GoodPhotons = new PhotonOArr;
  GoodPhotons->SetName(fGoodPhotonsName);

  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {    
    const Photon *ph = fPhotons->At(i);        

    if (ph->Pt() <= fPhotonPtMin) 
      continue;


    
    Bool_t passSpikeRemovalFilter = kTRUE;

    if (ph->SCluster() && ph->SCluster()->Seed()) {
      if(ph->SCluster()->Seed()->Energy() > 5.0 && 
         ph->SCluster()->Seed()->EMax() / ph->SCluster()->Seed()->E3x3() > 0.95
        ) {
        passSpikeRemovalFilter = kFALSE;
      }
    }

    // For Now Only use the EMax/E3x3 prescription.
    //if(ph->SCluster()->Seed()->Energy() > 5.0 && 
    //   (1 - (ph->SCluster()->Seed()->E1x3() + ph->SCluster()->Seed()->E3x1() - 2*ph->SCluster()->Seed()->EMax())) > 0.95
    //  ) {
    //  passSpikeRemovalFilter = kFALSE;
    //}
    if (fApplySpikeRemoval && !passSpikeRemovalFilter) continue;
    


    if (ph->HadOverEm() >= fHadOverEmMax) 
      continue;

    if (fApplyPixelSeed == kTRUE &&
        ph->HasPixelSeed() == kTRUE) 
      continue;

    Bool_t idcut = kFALSE;
    switch (fPhIdType) {
      case kTight:
        idcut = ph->IsTightPhoton();
        break;
      case kLoose:
        idcut = ph->IsLoosePhoton();
       break;
      case kLooseEM:
        idcut = ph->IsLooseEM();
        break;
      case kCustomId:
        idcut = kTRUE;
      default:
        break;
    }

    if (!idcut) 
      continue;

    Bool_t isocut = kFALSE;
    switch (fPhIsoType) {      
      case kNoIso:
        isocut = kTRUE;
        break;
      case kCombinedIso:
        {
          Double_t totalIso = ph->HollowConeTrkIsoDr04()+
                              ph->EcalRecHitIsoDr04() +
                              ph->HcalTowerSumEtDr04();
          if (totalIso/ph->Pt() < 0.25)
            isocut = kTRUE;
        }
        break;
      case kCustomIso:
        {
          if ( ph->HollowConeTrkIsoDr04() < (1.5 + 0.001*ph->Pt()) && ph->EcalRecHitIsoDr04()<(2.0+0.006*ph->Pt()) && ph->HcalTowerSumEtDr04()<(2.0+0.0025*ph->Pt()) )
            isocut = kTRUE;
        }
	break;
	
    case kMITPUCorrected:
      {
	// compute the PU corrections only if Rho is available
	// ... otherwise (_tRho = -1.0) it's the std isolation
	isocut = kTRUE;
	Double_t EcalCorrISO =   ph->EcalRecHitIsoDr04();
	if(_tRho > -0.5 ) EcalCorrISO -= _tRho * fEffAreaEcal;
	if ( EcalCorrISO > (2.0+0.006*ph->Pt()) ) isocut = kFALSE; 
	if ( isocut ) {
	  Double_t HcalCorrISO = ph->HcalTowerSumEtDr04(); 
	  if(_tRho > -0.5 ) HcalCorrISO -= _tRho * fEffAreaHcal;
	  if ( HcalCorrISO > (2.0+0.0025*ph->Pt()) ) isocut = kFALSE;
	}
	if ( isocut ) {
	  Double_t TrackCorrISO = IsolationTools::TrackIsolationNoPV(ph, bsp, 0.4, 0.04, 0.0, 0.015, 0.1, TrackQuality::highPurity, fTracks);
	  if(_tRho > -0.5 )
	    TrackCorrISO -= _tRho * fEffAreaTrack;
	  if ( TrackCorrISO > (1.5 + 0.001*ph->Pt()) ) isocut = kFALSE;
	}
	break;
      }
    default:
      break;
    }
    
    if (!isocut) 
      continue;
    
    if ( fApplyR9Min && ph->R9() <= fPhotonR9Min) 
      continue;

    if (fFiduciality == kTRUE &&
        ph->IsEB() == kFALSE && ph->IsEE() == kFALSE) 
      continue;

    if ((ph->IsEB() == kTRUE && ph->CoviEtaiEta() >= fEtaWidthEB) ||
        (ph->IsEE() == kTRUE && ph->CoviEtaiEta() >= fEtaWidthEE))
      continue;

    if (ph->AbsEta() >= fAbsEtaMax) 
      continue;

    // add good electron
    GoodPhotons->Add(fPhotons->At(i));
  }

  // sort according to pt
  GoodPhotons->Sort();

  // add to event for other modules to use
  AddObjThisEvt(GoodPhotons);  
}

//--------------------------------------------------------------------------------------------------
void PhotonIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the photon collection branch.

  ReqEventObject(fPhotonBranchName,   fPhotons,   kTRUE);
  ReqEventObject(fTrackBranchName,    fTracks,    kTRUE);
  ReqEventObject(fBeamspotBranchName, fBeamspots, kTRUE);
  ReqEventObject(fPileUpDenName,      fPileUpDen, kTRUE);

  if (fPhotonIDType.CompareTo("Tight") == 0) 
    fPhIdType = kTight;
  else if (fPhotonIDType.CompareTo("Loose") == 0) 
    fPhIdType = kLoose;
  else if (fPhotonIDType.CompareTo("LooseEM") == 0) 
    fPhIdType = kLooseEM;
  else if (fPhotonIDType.CompareTo("Custom") == 0) {
    fPhIdType = kCustomId;
  } else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified photon identification %s is not defined.",
              fPhotonIDType.Data());
    return;
  }

  if (fPhotonIsoType.CompareTo("NoIso") == 0 )
    fPhIsoType = kNoIso;
  else if (fPhotonIsoType.CompareTo("CombinedIso") == 0 )
    fPhIsoType = kCombinedIso;
  else if (fPhotonIsoType.CompareTo("Custom") == 0 )
    fPhIsoType = kCustomIso;
  else if (fPhotonIsoType.CompareTo("MITPUCorrected") == 0 )
    fPhIsoType = kMITPUCorrected;
  else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified photon isolation %s is not defined.",
              fPhotonIsoType.Data());
    return;
  }
}
