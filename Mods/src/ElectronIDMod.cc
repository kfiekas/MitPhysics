// $Id: ElectronIDMod.cc,v 1.87 2011/05/12 21:31:02 sixie Exp $

#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::ElectronIDMod)

//--------------------------------------------------------------------------------------------------
ElectronIDMod::ElectronIDMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fElectronBranchName(Names::gkElectronBrn),
  fConversionBranchName(Names::gkMvfConversionBrn),
  fGoodElectronsName(ModNames::gkGoodElectronsName),  
  fNonIsolatedMuonsName("random"),  
  fNonIsolatedElectronsName("random"),  
  fVertexName(ModNames::gkGoodVertexesName),
  fBeamSpotName(Names::gkBeamSpotBrn),
  fTrackName(Names::gkTrackBrn),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fElectronIDType("CustomTight"),
  fElectronIsoType("PFIso"),
  fTrigObjectsName("HLTModTrigObjs"),
  fElectronPtMin(10),
  fElectronEtMin(0.0),  
  fElectronEtaMax(2.5),
  fIDLikelihoodCut(0.75),
  fTrackIsolationCut(5.0),
  fCaloIsolationCut(5.0),
  fEcalJuraIsoCut(5.0),
  fHcalIsolationCut(5.0),
  fCombIsolationCut(-1.0),
  fApplyConvFilterType1(kTRUE),
  fApplyConvFilterType2(kFALSE),
  fNWrongHitsMax(0),
  fNExpectedHitsInnerCut(999),
  fCombinedIdCut(kFALSE),
  fApplySpikeRemoval(kTRUE),
  fApplyD0Cut(kTRUE),
  fApplyDZCut(kTRUE),
  fChargeFilter(kTRUE),
  fD0Cut(0.020),
  fDZCut(0.20),
  fWhichVertex(-1),
  fApplyTriggerMatching(kFALSE),
  fApplyEcalSeeded(kFALSE),
  fApplyCombinedIso(kTRUE),
  fApplyEcalFiducial(kFALSE),
  fElectronsFromBranch(kTRUE),
  fElIdType(ElectronTools::kIdUndef),
  fElIsoType(ElectronTools::kIsoUndef),
  fElectrons(0),
  fConversions(0),
  fVertices(0),
  fBeamSpot(0),
  fTracks(0),
  fPFCandidates(0),
  fNonIsolatedMuons(0),
  fNonIsolatedElectrons(0),
  fLH(0),
  fPileupEnergyDensityName(Names::gkPileupEnergyDensityBrn),
  fPileupEnergyDensity(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronIDMod::Likelihood(const Electron *ele) const
{
  LikelihoodMeasurements measurements;
  measurements.pt = ele->Pt();
  measurements.subdet = (fabs(ele->Eta())<1.479) ? 0 : 1;
  measurements.deltaPhi = TMath::Abs(ele->DeltaPhiSuperClusterTrackAtVtx());
  measurements.deltaEta = TMath::Abs(ele->DeltaEtaSuperClusterTrackAtVtx());
  measurements.eSeedClusterOverPout = ele->ESeedClusterOverPout();
  measurements.eSuperClusterOverP = ele->ESuperClusterOverP();
  measurements.hadronicOverEm = ele->HadronicOverEm();
  measurements.sigmaIEtaIEta = ele->CoviEtaiEta();
  measurements.sigmaIPhiIPhi = TMath::Sqrt(ele->SCluster()->Seed()->CoviPhiiPhi());
  measurements.fBrem = ele->FBrem();
  measurements.nBremClusters = ele->NumberOfClusters() - 1;
  double likelihood = fLH->result(measurements);

  if(likelihood > fIDLikelihoodCut) return kTRUE;
  return kFALSE;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronIDMod::PassIDCut(const Electron *ele, ElectronTools::EElIdType idType) const
{

  Bool_t idcut = kFALSE;
  switch (idType) {
    case ElectronTools::kTight:
      idcut = ele->PassTightID();
      break;
    case ElectronTools::kLoose:
      idcut = ele->PassLooseID();
      break;
    case ElectronTools::kLikelihood:
      idcut = Likelihood(ele);
      break;
    case ElectronTools::kNoId:
      idcut = kTRUE;
      break;
    case ElectronTools::kCustomIdLoose:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kCustomIdLoose);
      break;
    case ElectronTools::kCustomIdTight:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kCustomIdTight);
      break;
    case ElectronTools::kVBTFWorkingPointFakeableId:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPointFakeableId);
      break;
    case ElectronTools::kVBTFWorkingPoint95Id:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPoint95Id);
      break;
    case ElectronTools::kVBTFWorkingPoint90Id:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPoint90Id);
      break;
    case ElectronTools::kVBTFWorkingPoint85Id:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPoint85Id);
      break;
    case ElectronTools::kVBTFWorkingPoint80Id:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPoint80Id);
      break;
    case ElectronTools::kVBTFWorkingPointLowPtId:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPointLowPtId);
      break;
    case ElectronTools::kVBTFWorkingPoint70Id:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPoint70Id);
      break;
    default:
      break;
  }
  
  return idcut;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronIDMod::PassIsolationCut(const Electron *ele, ElectronTools::EElIsoType isoType,
                                       const TrackCol *tracks, const Vertex *vertex, 
				       const Double_t rho) const
{

  Bool_t isocut = kFALSE;
  switch (isoType) {
    case ElectronTools::kTrackCalo:
      isocut = (ele->TrackIsolationDr03() < fTrackIsolationCut) &&
        (ele->CaloIsolation() < fCaloIsolationCut);
      break;
    case ElectronTools::kTrackJura:
      isocut = (ele->TrackIsolationDr03() < ele->Pt()*fTrackIsolationCut) &&
               (ele->EcalRecHitIsoDr03()  < ele->Pt()*fEcalJuraIsoCut) &&
               (ele->HcalTowerSumEtDr03() < ele->Pt()*fHcalIsolationCut);
      break;
    case ElectronTools::kTrackJuraCombined:
      isocut = (ele->TrackIsolationDr03() + ele->EcalRecHitIsoDr03() 
                - 1.5 < fCombIsolationCut);
      break;
    case ElectronTools::kTrackJuraSliding:
    {
      Double_t totalIso = ele->TrackIsolationDr03() + TMath::Max(ele->EcalRecHitIsoDr03() + ele->HcalTowerSumEtDr03() - rho * TMath::Pi() * 0.3 * 0.3, 0.0);
      if(ele->SCluster()->AbsEta() < 1.479) totalIso = ele->TrackIsolationDr03() + TMath::Max(TMath::Max(ele->EcalRecHitIsoDr03() - 1.0, 0.0) + ele->HcalTowerSumEtDr03() - rho * TMath::Pi() * 0.3 * 0.3, 0.0);
      if (totalIso < (ele->Pt()*fCombIsolationCut) )
        isocut = kTRUE;
    }
    break;
    case ElectronTools::kTrackJuraSlidingNoCorrection:
    {
      Double_t totalIso = ele->TrackIsolationDr03() + (ele->EcalRecHitIsoDr03() + ele->HcalTowerSumEtDr03());
      if(ele->SCluster()->AbsEta() < 1.479) totalIso = ele->TrackIsolationDr03() + (TMath::Max(ele->EcalRecHitIsoDr03() - 1.0, 0.0) + ele->HcalTowerSumEtDr03());
      if (totalIso < (ele->Pt()*fCombIsolationCut) )
        isocut = kTRUE;
    }
    break;
    case ElectronTools::kPFIso:
    {
      Double_t pfIsoCutValue = 9999;
      if(fCombIsolationCut > 0){
        pfIsoCutValue = fCombIsolationCut;
      } else {
        if (fabs(ele->SCluster()->Eta()) < 1.479) {
          if (ele->Pt() > 20) {
            pfIsoCutValue = 0.18;
          } else {
            pfIsoCutValue = 0.14;
          }
        } else {
          pfIsoCutValue = 0.10;
        }
      }
      Double_t totalIso = IsolationTools::PFElectronIsolation(ele, fPFCandidates, vertex, 0.2, 1.0, 0.4, 0.0);
      if (totalIso < (ele->Pt()*pfIsoCutValue) )
        isocut = kTRUE;     
    }
    break;
    case ElectronTools::kPFIsoNoL:
    {
      Double_t beta = IsolationTools::BetaE(tracks, ele, vertex, 0.0, 0.2, 0.3, 0.02); 
      if(beta == 0) beta = 1.0;
      Double_t totalIso = IsolationTools::PFElectronIsolation(ele, fPFCandidates, vertex, fNonIsolatedMuons, fNonIsolatedElectrons, 0.1, 1.0, 0.4, 0.0, 3, beta);
      if (totalIso < (ele->Pt()*fCombIsolationCut) )
        isocut = kTRUE;
    }
    break;
    case ElectronTools::kVBTFWorkingPoint95Iso:
      isocut = ElectronTools::PassCustomIso(ele, ElectronTools::kVBTFWorkingPoint95Iso, fApplyCombinedIso);
      break;
    case ElectronTools::kVBTFWorkingPoint90Iso:
      isocut = ElectronTools::PassCustomIso(ele, ElectronTools::kVBTFWorkingPoint90Iso, fApplyCombinedIso);
      break;
    case ElectronTools::kVBTFWorkingPoint85Iso:
      isocut = ElectronTools::PassCustomIso(ele, ElectronTools::kVBTFWorkingPoint85Iso, fApplyCombinedIso);
      break;
    case ElectronTools::kVBTFWorkingPoint80Iso:
      isocut = ElectronTools::PassCustomIso(ele, ElectronTools::kVBTFWorkingPoint80Iso, fApplyCombinedIso);
      break;
    case ElectronTools::kVBTFWorkingPoint70Iso:
      isocut = ElectronTools::PassCustomIso(ele, ElectronTools::kVBTFWorkingPoint70Iso, fApplyCombinedIso);
      break;
    case ElectronTools::kNoIso:
      isocut = kTRUE;
      break;
    case ElectronTools::kCustomIso:
    default:
      break;
  }

  return isocut;
}


//--------------------------------------------------------------------------------------------------
void ElectronIDMod::Process()
{
  // Process entries of the tree. 

  if(fElIsoType != ElectronTools::kPFIsoNoL) {
    LoadEventObject(fElectronBranchName, fElectrons);
  }
  else {
    fElectrons    = GetObjThisEvt<ElectronOArr>(fElectronBranchName);
    fNonIsolatedMuons	  = GetObjThisEvt<MuonCol>(fNonIsolatedMuonsName);
    fNonIsolatedElectrons = GetObjThisEvt<ElectronCol>(fNonIsolatedElectronsName);
  }
  LoadEventObject(fBeamSpotName, fBeamSpot);
  LoadEventObject(fTrackName, fTracks);
  LoadEventObject(fPFCandidatesName, fPFCandidates);
  if(fElIsoType == ElectronTools::kTrackJuraSliding) {
    LoadEventObject(fPileupEnergyDensityName, fPileupEnergyDensity);
  }
  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);

  //get trigger object collection if trigger matching is enabled
  const TriggerObjectCol *trigObjs = 0;
  if (fApplyTriggerMatching) {
    trigObjs = GetHLTObjects(fTrigObjectsName);
  }
  
  ElectronOArr *GoodElectrons = new ElectronOArr;
  GoodElectrons->SetName(fGoodElectronsName);

  for (UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    const Electron *e = fElectrons->At(i);        

    if (e->SCluster() == 0) 
      continue;
    
    if (e->Pt() < fElectronPtMin) 
      continue;
    
    if (e->SCluster()->Et() < fElectronEtMin)
      continue;    
    
    if (e->AbsEta() > fElectronEtaMax) 
      continue;
    
    if (fApplyEcalFiducial && ( (e->SCluster()->AbsEta()>1.4442 && e->SCluster()->AbsEta()<1.5666) || e->SCluster()->AbsEta()>2.5 )) {
      continue;
    }
    
    if (fApplyEcalSeeded && !e->IsEcalDriven()) {
      continue;
    }    
    
    //apply trigger matching
    Bool_t matchTrigger = fApplyTriggerMatching && ElectronTools::PassTriggerMatching(e,trigObjs);
    if (fApplyTriggerMatching && !matchTrigger)
      continue;
    
    //apply ECAL spike removal    
    Bool_t spikecut = ElectronTools::PassSpikeRemovalFilter(e);
    if (fApplySpikeRemoval && !spikecut)
      continue;

    //apply id cut
    Bool_t idcut = PassIDCut(e, fElIdType);
    if (!idcut) 
      continue;

    //apply Isolation Cut
    Double_t Rho = 0.0;
    if(fElIsoType == ElectronTools::kTrackJuraSliding) {
     Rho = fPileupEnergyDensity->At(0)->Rho();
    }
    Bool_t isocut = PassIsolationCut(e, fElIsoType, fTracks, fVertices->At(0), Rho);
    if (!isocut)
      continue;

    // apply conversion filters
    Bool_t passConvVetoType1 = kFALSE;
    if (fApplyConvFilterType1) {
      LoadEventObject(fConversionBranchName, fConversions);
      passConvVetoType1 = ElectronTools::PassConversionFilter(e, fConversions, 
                                                         fBeamSpot->At(0), 0, 1e-6, 2.0, kTRUE, kFALSE);      
    }
    else {
      passConvVetoType1 = kTRUE;
    }

    if (passConvVetoType1 == kFALSE) continue;

    Bool_t passConvVetoType2 = kFALSE;
    if (fApplyConvFilterType2) {
      passConvVetoType2 = TMath::Abs(e->ConvPartnerDCotTheta()) >= 0.02 || 
                          TMath::Abs(e->ConvPartnerDist())      >= 0.02;
    }
    else {
      passConvVetoType2 = kTRUE;
    }
    
    if (passConvVetoType2 == kFALSE) continue;

    // apply NExpectedHitsInner Cut
    if(fNExpectedHitsInnerCut < 999 && 
       e->CorrectedNExpectedHitsInner() > fNExpectedHitsInnerCut) continue;

    // apply d0 cut
    if (fApplyD0Cut) {
      Bool_t passD0cut = kTRUE;
      if(fWhichVertex >= -1) passD0cut = ElectronTools::PassD0Cut(e, fVertices, fD0Cut, fWhichVertex);
      else                   passD0cut = ElectronTools::PassD0Cut(e, fBeamSpot, fD0Cut);
      if (!passD0cut)
        continue;
    }

    // apply dz cut
    if (fApplyDZCut) {
      Bool_t passDZcut = ElectronTools::PassDZCut(e, fVertices, fDZCut, fWhichVertex);
      if (!passDZcut)
        continue;
    }

    // apply charge filter
    if(fChargeFilter == kTRUE) {
      Bool_t passChargeFilter = ElectronTools::PassChargeFilter(e);
      if (!passChargeFilter) continue;
    }

    // apply full combined id, using Tight cuts
    if(fCombinedIdCut == kTRUE) {
      fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
      LoadEventObject(fConversionBranchName, fConversions);
      Int_t result = ElectronTools::PassTightId(e, *&fVertices, fConversions, 2);
      if(result != 15) continue;
    }

    // add good electron
    GoodElectrons->Add(e);
  }

  // sort according to pt
  GoodElectrons->Sort();

  // add to event for other modules to use
  AddObjThisEvt(GoodElectrons);  
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the electron collection branch.

   // In this case we cannot have a branch
  if (fElectronIsoType.CompareTo("PFIsoNoL") != 0 ) {
    ReqEventObject(fElectronBranchName, fElectrons,fElectronsFromBranch);
  }
  ReqEventObject(fBeamSpotName, fBeamSpot, kTRUE);
  ReqEventObject(fTrackName, fTracks, kTRUE);
  ReqEventObject(fPFCandidatesName, fPFCandidates, kTRUE);
  if (fElectronIsoType.CompareTo("TrackJuraSliding") == 0 ) {
    ReqEventObject(fPileupEnergyDensityName, fPileupEnergyDensity, kTRUE);
  }

  if(fCombinedIdCut == kTRUE) {
    fElectronIDType  	  = "NoId";
    fElectronIsoType 	  = "NoIso";
    fApplyConvFilterType1 = kFALSE;
    fApplyConvFilterType2 = kFALSE;
    fApplyD0Cut           = kFALSE;
    fApplyDZCut           = kFALSE;
  }

  if (fApplyConvFilterType1 || fCombinedIdCut == kTRUE)
    ReqEventObject(fConversionBranchName, fConversions, kTRUE);

  Setup();
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMod::Setup()
{
  // Set all options properly before execution.

  if (fElectronIDType.CompareTo("Tight") == 0) 
    fElIdType = ElectronTools::kTight;
  else if (fElectronIDType.CompareTo("Loose") == 0) 
    fElIdType = ElectronTools::kLoose;
  else if (fElectronIDType.CompareTo("Likelihood") == 0) 
    fElIdType = ElectronTools::kLikelihood;
  else if (fElectronIDType.CompareTo("NoId") == 0) 
    fElIdType = ElectronTools::kNoId;
  else if (fElectronIDType.CompareTo("ZeeId") == 0) 
    fElIdType = ElectronTools::kZeeId;
  else if (fElectronIDType.CompareTo("CustomLoose") == 0) 
    fElIdType = ElectronTools::kCustomIdLoose;
  else if (fElectronIDType.CompareTo("CustomTight") == 0) 
    fElIdType = ElectronTools::kCustomIdTight;
  else if (fElectronIDType.CompareTo("VBTFWorkingPointFakeableId") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPointFakeableId;
  else if (fElectronIDType.CompareTo("VBTFWorkingPoint95Id") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPoint95Id;
  else if (fElectronIDType.CompareTo("VBTFWorkingPoint90Id") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPoint90Id;
  else if (fElectronIDType.CompareTo("VBTFWorkingPoint80Id") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPoint80Id;
  else if (fElectronIDType.CompareTo("VBTFWorkingPointLowPtId") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPointLowPtId;
  else if (fElectronIDType.CompareTo("VBTFWorkingPoint85Id") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPoint85Id;
  else if (fElectronIDType.CompareTo("VBTFWorkingPoint70Id") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPoint70Id;
  
   else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified electron identification %s is not defined.",
              fElectronIDType.Data());
    return;
  }

  if (fElectronIsoType.CompareTo("TrackCalo") == 0 )
    fElIsoType = ElectronTools::kTrackCalo;
  else if (fElectronIsoType.CompareTo("TrackJura") == 0) 
    fElIsoType = ElectronTools::kTrackJura;
  else if(fElectronIsoType.CompareTo("TrackJuraCombined") == 0)
    fElIsoType = ElectronTools::kTrackJuraCombined;
  else if(fElectronIsoType.CompareTo("TrackJuraSliding") == 0)
    fElIsoType = ElectronTools::kTrackJuraSliding;
  else if(fElectronIsoType.CompareTo("TrackJuraSlidingNoCorrection") == 0)
    fElIsoType = ElectronTools::kTrackJuraSlidingNoCorrection;
  else if (fElectronIsoType.CompareTo("PFIso") == 0 )
    fElIsoType = ElectronTools::kPFIso;
  else if (fElectronIsoType.CompareTo("PFIsoNoL") == 0 )
    fElIsoType = ElectronTools::kPFIsoNoL;
  else if (fElectronIsoType.CompareTo("NoIso") == 0 )
    fElIsoType = ElectronTools::kNoIso;
  else if (fElectronIsoType.CompareTo("ZeeIso") == 0 )
    fElIsoType = ElectronTools::kZeeIso;
  else if (fElectronIsoType.CompareTo("VBTFWorkingPoint95Iso") == 0 )
    fElIsoType = ElectronTools::kVBTFWorkingPoint95Iso;
  else if (fElectronIsoType.CompareTo("VBTFWorkingPoint90Iso") == 0 )
    fElIsoType = ElectronTools::kVBTFWorkingPoint90Iso;
  else if (fElectronIsoType.CompareTo("VBTFWorkingPoint85Iso") == 0 )
    fElIsoType = ElectronTools::kVBTFWorkingPoint85Iso;
  else if (fElectronIsoType.CompareTo("VBTFWorkingPoint80Iso") == 0 )
    fElIsoType = ElectronTools::kVBTFWorkingPoint80Iso;
  else if (fElectronIsoType.CompareTo("VBTFWorkingPoint70Iso") == 0 )
    fElIsoType = ElectronTools::kVBTFWorkingPoint70Iso;
  else if (fElectronIsoType.CompareTo("Custom") == 0 ) {
    fElIsoType = ElectronTools::kCustomIso;
    SendError(kWarning, "SlaveBegin",
              "Custom electron isolation is not yet implemented.");
  } else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified electron isolation %s is not defined.",
              fElectronIsoType.Data());
    return;
  }


}

//--------------------------------------------------------------------------------------------------
void ElectronIDMod::Terminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}
