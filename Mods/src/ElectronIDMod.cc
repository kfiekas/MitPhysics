// $Id: ElectronIDMod.cc,v 1.105 2011/09/25 18:14:03 sixie Exp $

#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

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
  fIDLikelihoodCut(-999.0),
  fTrackIsolationCut(5.0),
  fCaloIsolationCut(5.0),
  fEcalJuraIsoCut(5.0),
  fHcalIsolationCut(5.0),
  fCombIsolationCut(0.1),
  fCombRelativeIsolationCut(0.10),
  fPFIsolationCut(-1.0),
  fApplyConvFilterType1(kTRUE),
  fApplyConvFilterType2(kFALSE),
  fNWrongHitsMax(0),
  fNExpectedHitsInnerCut(999),
  fInvertNExpectedHitsInnerCut(kFALSE),
  fCombinedIdCut(kFALSE),
  fApplySpikeRemoval(kTRUE),
  fApplyD0Cut(kTRUE),
  fApplyDZCut(kTRUE),
  fChargeFilter(kTRUE),
  fD0Cut(0.020),
  fDZCut(0.10),
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
  fIntRadius(0.0),
  fNonIsolatedMuons(0),
  fNonIsolatedElectrons(0),
  fLH(0),
  fPileupEnergyDensityName(Names::gkPileupEnergyDensityBrn),
  fPileupEnergyDensity(0),
  fElectronIDMVA(0),
  fElectronMVAWeights_Subdet0Pt10To20(""),
  fElectronMVAWeights_Subdet1Pt10To20(""),
  fElectronMVAWeights_Subdet2Pt10To20(""),
  fElectronMVAWeights_Subdet0Pt20ToInf(""),
  fElectronMVAWeights_Subdet1Pt20ToInf(""),
  fElectronMVAWeights_Subdet2Pt20ToInf("")
  
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronIDMod::PassLikelihoodID(const Electron *ele) const
{

  Double_t LikelihoodValue = ElectronTools::Likelihood(fLH, ele);

  double likCut = fIDLikelihoodCut;
  if(likCut > -900){
    if(ele->Pt() > 20){
      if(ele->SCluster()->AbsEta() < 1.479){
        if(ele->NumberOfClusters() - 1 == 0) likCut = 3.5;
	else                                 likCut = 4.0;
      }
      else  {                                
        if(ele->NumberOfClusters() - 1 == 0) likCut = 4.0;
	else                                 likCut = 4.0;
      }
    }
    else {
      if(ele->SCluster()->AbsEta() < 1.479){
        if(ele->NumberOfClusters() - 1 == 0) likCut =  4.0;
	else                                 likCut =  4.5;
      }
      else  {                                
        if(ele->NumberOfClusters() - 1 == 0) likCut =  4.0;
	else                                 likCut =  4.0;
      }
    }
  }
  if (LikelihoodValue > likCut) return kTRUE;
  return kFALSE;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronIDMod::PassMVAID(const Electron *el, ElectronTools::EElIdType idType, 
                                const Vertex *vertex) const
{
  Double_t MVAValue = fElectronIDMVA->MVAValue(el, vertex);
 
  Int_t subdet = 0;
  if (el->SCluster()->AbsEta() < 1.0) subdet = 0;
  else if (el->SCluster()->AbsEta() < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (el->Pt() > 20.0) ptBin = 1;
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t MVACut = -9999;
  if (MVABin == 0) MVACut = 0.362;
  if (MVABin == 1) MVACut = 0.467;
  if (MVABin == 2) MVACut = 0.320; 
  if (MVABin == 3) MVACut = 0.968;
  if (MVABin == 4) MVACut = 0.961;
  if (MVABin == 5) MVACut = 0.929;
 
  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronIDMod::PassIDCut(const Electron *ele, ElectronTools::EElIdType idType, 
                                const Vertex *vertex) const
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
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPointFakeableId) &&
              PassLikelihoodID(ele);
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
    case ElectronTools::kMVA_BDTG_V3:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPointFakeableId) &&
        PassMVAID(ele, ElectronTools::kMVA_BDTG_V3, vertex);
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
    case ElectronTools::kCombinedRelativeConeAreaCorrected:
    {
      Double_t totalIso = ele->TrackIsolationDr03() + ele->EcalRecHitIsoDr03() + ele->HcalTowerSumEtDr03() - rho * TMath::Pi() * 0.3 * 0.3;
      if (totalIso < (ele->Pt()*fCombRelativeIsolationCut) )
        isocut = kTRUE;
    }
    break;
    case ElectronTools::kPFIso:
    {
      Double_t pfIsoCutValue = 9999;
      if(fPFIsolationCut > 0){
        pfIsoCutValue = fPFIsolationCut;
      } else {
        if (ele->SCluster()->AbsEta() < 1.479) {
          if (ele->Pt() > 20) {
            pfIsoCutValue = 0.13;
          } else {
            pfIsoCutValue = 0.13;
          }
        } else {
          if (ele->Pt() > 20) {
            pfIsoCutValue = 0.09;
          } else {
            pfIsoCutValue = 0.09;
          }
	}
      }
      Double_t totalIso = IsolationTools::PFElectronIsolation(ele, fPFCandidates, vertex, 0.1, 1.0, 0.4, fIntRadius);
      if (totalIso < (ele->Pt()*pfIsoCutValue) )
        isocut = kTRUE;
    }
    break;
    case ElectronTools::kPFIsoNoL:
    {
      Double_t beta = IsolationTools::BetaE(tracks, ele, vertex, 0.0, 0.2, 0.3, 0.02); 
      if(beta == 0) beta = 1.0;
      Double_t totalIso = IsolationTools::PFElectronIsolation(ele, fPFCandidates, vertex, fNonIsolatedMuons, fNonIsolatedElectrons, 0.1, 1.0, 0.4, 0.0, 3, beta);
      if (totalIso < (ele->Pt()*fPFIsolationCut) )
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
  if(fElIsoType == ElectronTools::kTrackJuraSliding || 
     fElIsoType == ElectronTools::kCombinedRelativeConeAreaCorrected) {
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

    //apply Isolation Cut
    Double_t Rho = 0.0;
    if( fElIsoType == ElectronTools::kTrackJuraSliding
        || fElIsoType == ElectronTools::kCombinedRelativeConeAreaCorrected ) {
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
    if(fInvertNExpectedHitsInnerCut == kFALSE && fNExpectedHitsInnerCut < 999 && 
       e->CorrectedNExpectedHitsInner() > fNExpectedHitsInnerCut) continue;

    // apply NExpectedHitsInner inverted Cut
    if(fInvertNExpectedHitsInnerCut == kTRUE && fNExpectedHitsInnerCut < 999 && 
       e->CorrectedNExpectedHitsInner() <= fNExpectedHitsInnerCut) continue;

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

    //apply id cut
    Bool_t idcut = PassIDCut(e, fElIdType, fVertices->At(0));
    if (!idcut) 
      continue;

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
  if (fElectronIsoType.CompareTo("TrackJuraSliding") == 0 
    || fElectronIsoType.CompareTo("CombinedRelativeConeAreaCorrected") == 0 ) {
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
  else if (fElectronIDType.CompareTo("Likelihood") == 0) {
    if (!fLH) { cout << "Error: Likelihood not initialized.\n"; assert(0); }
    fElIdType = ElectronTools::kLikelihood;
  } else if (fElectronIDType.CompareTo("NoId") == 0) 
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
  else if (fElectronIDType.CompareTo("MVA_BDTG_V3") == 0) {
    fElIdType = ElectronTools::kMVA_BDTG_V3; 
    if (!fLH) { cout << "Error: Likelihood not initialized.\n"; assert(0);}
  } else {
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
  else if(fElectronIsoType.CompareTo("CombinedRelativeConeAreaCorrected") == 0)
    fElIsoType = ElectronTools::kCombinedRelativeConeAreaCorrected;
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


  //If we use MVA ID, need to load MVA weights
  if (fElIdType == ElectronTools::kMVA_BDTG_V3) {
    fElectronIDMVA = new ElectronIDMVA();
    fElectronIDMVA->Initialize("BDTG method",
                               fElectronMVAWeights_Subdet0Pt10To20,
                               fElectronMVAWeights_Subdet1Pt10To20,
                               fElectronMVAWeights_Subdet2Pt10To20,
                               fElectronMVAWeights_Subdet0Pt20ToInf,
                               fElectronMVAWeights_Subdet1Pt20ToInf,
                               fElectronMVAWeights_Subdet2Pt20ToInf,
                               fLH);
  }

}

//--------------------------------------------------------------------------------------------------
void ElectronIDMod::Terminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}
