
// $Id: PhotonIDMod.cc,v 1.29 2011/12/19 23:45:00 bendavid Exp $

#include "TDataMember.h"
#include "TTree.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"

#include <TSystem.h>

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
  fConversionName    ("MergedConversions"),
  fElectronName      ("Electrons"),
  fGoodElectronName  ("Electrons"),
  fPVName            (Names::gkPVBeamSpotBrn),
  // MC specific stuff...
  fMCParticleName    (Names::gkMCPartBrn),
  fPileUpName        (Names::gkPileupInfoBrn),

  fPhotonIDType      ("Custom"),
  fPhotonIsoType     ("Custom"),
  fPhotonPtMin       (15.0),
  fHadOverEmMax      (0.02),
  fApplySpikeRemoval (kFALSE),
  fApplyPixelSeed    (kTRUE),
  fApplyElectronVeto (kFALSE),
  fInvertElectronVeto(kFALSE),
  fApplyElectronVetoConvRecovery(kFALSE),
  fApplyConversionId (kFALSE),
  fApplyTriggerMatching(kFALSE),
  fPhotonR9Min       (0.5),
  fPhIdType          (kIdUndef),
  fPhIsoType         (kIsoUndef),
  fFiduciality       (kTRUE),
  fEtaWidthEB        (0.01),
  fEtaWidthEE        (0.028),
  fAbsEtaMax         (999.99),
  fApplyR9Min        (kFALSE),

  fEffAreaEcalEE     (0.089),
  fEffAreaHcalEE     (0.156),
  fEffAreaTrackEE    (0.261),

  fEffAreaEcalEB     (0.183),
  fEffAreaHcalEB     (0.062),
  fEffAreaTrackEB    (0.306),

  fPhotons(0),
  fTracks(0),
  fBeamspots(0),
  fPileUpDen(0),
  fConversions(0),
  fElectrons(0),
  fGoodElectrons(0),
  fPV(0),
  fMCParticles       (0),
  fPileUp            (0),

  // MVA ID Stuff
  fbdtCutBarrel      (0.0744), //cuts give the same effiiciency (relative to preselection) with cic
  fbdtCutEndcap      (0.0959), //cuts give the same effiiciency (relative to preselection) with cic  
  fVariableType      (10), //please use 4 which is the correct type
  fEndcapWeights      (gSystem->Getenv("CMSSW_BASE")+TString("/src/MitPhysics/data/TMVAClassificationPhotonID_Endcap_PassPreSel_Variable_10_BDTnCuts2000_BDT.weights.xml")),
  fBarrelWeights      (gSystem->Getenv("CMSSW_BASE")+TString("/src/MitPhysics/data/TMVAClassificationPhotonID_Barrel_PassPreSel_Variable_10_BDTnCuts2000_BDT.weights.xml")),


  fDoMCR9Scaling     (kFALSE),
  fMCR9ScaleEB       (1.0),
  fMCR9ScaleEE       (1.0),
  fDoMCSigIEtaIEtaScaling(kFALSE),
  fDoMCWidthScaling(kFALSE),
  fDoMCErrScaling     (kFALSE),
  fMCErrScaleEB       (1.0),
  fMCErrScaleEE       (1.0),

  fPhotonsFromBranch(true),
  fPVFromBranch      (true),
  fGoodElectronsFromBranch (kTRUE),
  fIsData(false)
  
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void PhotonIDMod::Process()
{
  // Process entries of the tree. 

  LoadEventObject(fPhotonBranchName,   fPhotons);
  
  const BaseVertex *bsp = NULL;
  Double_t _tRho = -1.;
  const TriggerObjectCol *trigObjs = 0;  
  if (fPhotons->GetEntries()>0) {
    LoadEventObject(fTrackBranchName,    fTracks);
    LoadEventObject(fBeamspotBranchName, fBeamspots);
    LoadEventObject(fPileUpDenName,      fPileUpDen);
    LoadEventObject(fConversionName,     fConversions);
    LoadEventObject(fElectronName,       fElectrons);
    LoadEventObject(fGoodElectronName,   fGoodElectrons);
    LoadEventObject(fPVName,             fPV);
    
    if(fBeamspots->GetEntries() > 0)
      bsp = fBeamspots->At(0);

    if(fPileUpDen->GetEntries() > 0)
      _tRho = (Double_t) fPileUpDen->At(0)->RhoRandomLowEta();

    
    //get trigger object collection if trigger matching is enabled
    if (fApplyTriggerMatching) {
      trigObjs = GetHLTObjects(fTrigObjectsName);
    }    
  
  }
  
  PhotonOArr *GoodPhotons = new PhotonOArr;
  GoodPhotons->SetName(fGoodPhotonsName);
  GoodPhotons->SetOwner(kTRUE);

  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {    
    // need to cpoy the photon in order to be able to scale R9 etc.
    Photon *ph = new Photon(*fPhotons->At(i));        

    
    if (fFiduciality == kTRUE &&
        (ph->SCluster()->AbsEta()>=2.5 || (ph->SCluster()->AbsEta()>=1.4442 && ph->SCluster()->AbsEta()<=1.566) ) ) 
      continue;
    
    if (fInvertElectronVeto && PhotonTools::PassElectronVeto(ph,fGoodElectrons) && false) {
      continue;
    }
    
    // -----------------------------------------------------------------------------------
    // Do all the scaling ....
    if (fDoMCErrScaling && !fIsData) {
      if (ph->SCluster()->AbsEta()<1.5) PhotonTools::ScalePhotonError(ph,fMCErrScaleEB);
      else PhotonTools::ScalePhotonError(ph,fMCErrScaleEE);
    }

    if (fDoMCR9Scaling && !fIsData) {
      if (ph->SCluster()->AbsEta()<1.5) PhotonTools::ScalePhotonR9(ph,fMCR9ScaleEB);
      else PhotonTools::ScalePhotonR9(ph,fMCR9ScaleEE);
    }

    if (fDoMCSigIEtaIEtaScaling && !fIsData) {
      if (ph->SCluster()->AbsEta()<1.5) ph->SetCoviEtaiEta(0.87*ph->CoviEtaiEta() + 0.0011);
      else ph->SetCoviEtaiEta(0.99*ph->CoviEtaiEta());
    }

    if (fDoMCWidthScaling && !fIsData) {
      ph->SetEtaWidth(0.99*ph->EtaWidth());
      ph->SetPhiWidth(0.99*ph->PhiWidth());
    }

    

    // ---------------------------------------------------------------------
    // check if we use the CiC Selection. If yes, bypass all the below...
    if(fPhIdType == kBaseLineCiC) {
      if( PhotonTools::PassCiCSelection(ph, fPV->At(0), fTracks, fElectrons, fPV, _tRho, fPhotonPtMin, fApplyElectronVeto) )
	GoodPhotons->AddOwned(ph);
      continue; // go to next Photons
    }
    // ---------------------------------------------------------------------

    //loose photon preselection for subsequent mva
    if(fPhIdType == kMITPhSelection ) {
      if( ph->Pt()>fPhotonPtMin && PhotonTools::PassSinglePhotonPresel(ph,fElectrons,fConversions,bsp,fTracks,fPV->At(0),_tRho,fApplyElectronVeto,fInvertElectronVeto) ) {
        GoodPhotons->AddOwned(ph);
      }
      continue;
    }

    // add MingMings MVA ID on single Photon level
    if(fPhIdType == kMITMVAId ) {
      if( ph->Pt()>fPhotonPtMin && PhotonTools::PassSinglePhotonPresel(ph,fElectrons,fConversions,bsp,fTracks,fPV->At(0),_tRho,fApplyElectronVeto) && fTool.PassMVASelection(ph, fPV->At(0) ,fTracks, fPV, _tRho ,fbdtCutBarrel,fbdtCutEndcap, fElectrons, fApplyElectronVeto) ) {
	GoodPhotons->AddOwned(ph);
      }
      continue;
    }  // go to next Photon

    if (ph->Pt() <= fPhotonPtMin) 
      continue;    // add good electron

    Bool_t isbarrel = ph->SCluster()->AbsEta()<1.5;
    
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

    if (fApplyElectronVeto && !PhotonTools::PassElectronVeto(ph,fElectrons) ) continue;

    if (fApplyElectronVetoConvRecovery && !PhotonTools::PassElectronVetoConvRecovery(ph,fElectrons,fConversions,bsp) ) continue;

    if (fApplyConversionId && !PhotonTools::PassConversionId(ph,PhotonTools::MatchedConversion(ph,fConversions,bsp))) continue;

    if (fApplyTriggerMatching && !PhotonTools::PassTriggerMatching(ph,trigObjs)) continue;

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
	Double_t fEffAreaEcal = fEffAreaEcalEB;
	Double_t fEffAreaHcal = fEffAreaHcalEB;
	Double_t fEffAreaTrack = fEffAreaTrackEB;

	if( !isbarrel ) {
	  fEffAreaEcal = fEffAreaEcalEE;
	  fEffAreaHcal = fEffAreaHcalEE;
	  fEffAreaTrack = fEffAreaTrackEE;
	}	  

	//std::cout<<"-----------------------------------"<<std::endl;
	//std::cout<<" MIT phooton Et = "<<ph->Et()<<std::endl;		  
	Double_t EcalCorrISO =   ph->EcalRecHitIsoDr04();
	//std::cout<<" ecaliso4    = "<<EcalCorrISO<<std::endl;
	if(_tRho > -0.5 ) EcalCorrISO -= _tRho * fEffAreaEcal;
	//std::cout<<" ecaliso4Corr = "<<EcalCorrISO<<"  ("<<2.0+0.006*ph->Pt()<<")"<<std::endl;
	if ( EcalCorrISO > (2.0+0.006*ph->Pt()) ) isocut = kFALSE; 
	if ( isocut || true ) {
	  Double_t HcalCorrISO = ph->HcalTowerSumEtDr04(); 
	  //std::cout<<" hcaliso4     = "<<HcalCorrISO<<std::endl;
	  if(_tRho > -0.5 ) HcalCorrISO -= _tRho * fEffAreaHcal;
	  //std::cout<<" hcaliso4Corr = "<<HcalCorrISO<<"  ("<<2.0+0.0025*ph->Pt()<<")"<<std::endl;
	  if ( HcalCorrISO > (2.0+0.0025*ph->Pt()) ) isocut = kFALSE;
	}
	if ( isocut || true ) {
	  Double_t TrackCorrISO = IsolationTools::TrackIsolationNoPV(ph, bsp, 0.4, 0.04, 0.0, 0.015, 0.1, TrackQuality::highPurity, fTracks);
	  //std::cout<<" TrackIso       = "<<TrackCorrISO<<std::endl;	  
	  if(_tRho > -0.5 )
	    TrackCorrISO -= _tRho * fEffAreaTrack;
	  //std::cout<<" TrackIsoCorr   = "<<TrackCorrISO<<"  ("<<1.5 + 0.001*ph->Pt()<<")"<<std::endl;	  
	  if ( TrackCorrISO > (1.5 + 0.001*ph->Pt()) ) isocut = kFALSE;
	  //std::cout<<"      passes ? "<<isocut<<std::endl;
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

    if ((isbarrel && ph->CoviEtaiEta() >= fEtaWidthEB) ||
        (!isbarrel && ph->CoviEtaiEta() >= fEtaWidthEE))
      continue;

    if (ph->AbsEta() >= fAbsEtaMax) 
      continue;
      
    
    // add good electron
    GoodPhotons->AddOwned(ph);
  }

  // sort according to pt
  GoodPhotons->Sort();

  // add to event for other modules to use
  AddObjThisEvt(GoodPhotons);  
  
  return; 
  
}

//--------------------------------------------------------------------------------------------------
void PhotonIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the photon collection branch.

  ReqEventObject(fPhotonBranchName,   fPhotons,     fPhotonsFromBranch);
  ReqEventObject(fTrackBranchName,    fTracks,      kTRUE);
  ReqEventObject(fBeamspotBranchName, fBeamspots,   kTRUE);
  ReqEventObject(fConversionName,     fConversions, kTRUE);
  ReqEventObject(fElectronName,       fElectrons,   kTRUE);
  ReqEventObject(fGoodElectronName,       fGoodElectrons,   fGoodElectronsFromBranch);  
  ReqEventObject(fPVName, fPV, fPVFromBranch);
  ReqEventObject(fPileUpDenName,      fPileUpDen, kTRUE);

  if (!fIsData) {
    ReqBranch(fPileUpName,            fPileUp);
    ReqBranch(fMCParticleName,        fMCParticles);
  }   
  
  if (fPhotonIDType.CompareTo("Tight") == 0) 
    fPhIdType = kTight;
  else if (fPhotonIDType.CompareTo("Loose") == 0) 
    fPhIdType = kLoose;
  else if (fPhotonIDType.CompareTo("LooseEM") == 0) 
    fPhIdType = kLooseEM;
  else if (fPhotonIDType.CompareTo("Custom") == 0)
    fPhIdType = kCustomId;
  else if (fPhotonIDType.CompareTo("BaseLineCiC") == 0) {
    fPhIdType = kBaseLineCiC;
    fPhotonIsoType = "NoIso";
  }
  else if (fPhotonIDType.CompareTo("MITMVAId") == 0) {
    fPhIdType = kMITMVAId;
    fPhotonIsoType = "NoIso";
    fTool.InitializeMVA(fVariableType,fEndcapWeights,fBarrelWeights);
  }
  else if (fPhotonIDType.CompareTo("MITSelection") == 0)  {
    fPhIdType = kMITPhSelection;
    fPhotonIsoType = "NoIso";    
  }
  else {
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
