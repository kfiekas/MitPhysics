// $Id: ElectronIDMod.cc,v 1.67 2010/10/09 20:15:33 bendavid Exp $

#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
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
  fVertexName(ModNames::gkGoodVertexesName),
  fElectronIDType("CustomTight"),
  fElectronIsoType("TrackJuraSliding"),
  fTrigObjectsName("HLTModTrigObjs"),
  fElectronPtMin(10),
  fElectronEtMin(0.0),  
  fElectronEtaMax(2.5),
  fIDLikelihoodCut(0.9),
  fTrackIsolationCut(5.0),
  fCaloIsolationCut(5.0),
  fEcalJuraIsoCut(5.0),
  fHcalIsolationCut(5.0),
  fCombIsolationCut(5.0),
  fApplyConvFilterType1(kTRUE),
  fApplyConvFilterType2(kFALSE),
  fWrongHitsRequirement(kTRUE),
  fNExpectedHitsInnerCut(999),
  fCombinedIdCut(kFALSE),
  fApplySpikeRemoval(kTRUE),
  fApplyD0Cut(kTRUE),
  fChargeFilter(kTRUE),
  fD0Cut(0.020),
  fReverseIsoCut(kFALSE),
  fReverseD0Cut(kFALSE),
  fApplyTriggerMatching(kFALSE),
  fApplyEcalSeeded(kFALSE),
  fApplyCombinedIso(kTRUE),
  fApplyEcalFiducial(kFALSE),
  fElIdType(ElectronTools::kIdUndef),
  fElIsoType(ElectronTools::kIsoUndef),
  fElectrons(0),
  fConversions(0),
  fVertices(0)
{
  // Constructor.
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
      idcut = (ele->IDLikelihood() > fIDLikelihoodCut);
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
    case ElectronTools::kVBTFWorkingPoint70Id:
      idcut = ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPoint70Id);
      break;
    default:
      break;
  }
  
  return idcut;
}


//--------------------------------------------------------------------------------------------------
Bool_t ElectronIDMod::PassIsolationCut(const Electron *ele, ElectronTools::EElIsoType isoType) const
{

  Bool_t isocut = kFALSE;
  switch (isoType) {
    case ElectronTools::kTrackCalo:
      isocut = (ele->TrackIsolationDr03() < fTrackIsolationCut) &&
        (ele->CaloIsolation() < fCaloIsolationCut);
      break;
    case ElectronTools::kTrackJura:
      isocut = (ele->TrackIsolationDr03() < fTrackIsolationCut) &&
        (ele->EcalRecHitIsoDr03() < fEcalJuraIsoCut) &&
        (ele->HcalTowerSumEtDr03() < fHcalIsolationCut);
      break;
    case ElectronTools::kTrackJuraCombined:
      isocut = (ele->TrackIsolationDr03() + ele->EcalRecHitIsoDr03() 
                - 1.5 < fCombIsolationCut);
      break;
    case ElectronTools::kTrackJuraSliding:
    {
      Double_t totalIso = ele->TrackIsolationDr03() + 
                          TMath::Max(ele->EcalRecHitIsoDr03() - 1.0, 0.0)  +
			  ele->HcalTowerSumEtDr03();
      if (totalIso < (ele->Pt()*0.10) )
        isocut = kTRUE;
      
      if     (fReverseIsoCut == kTRUE &&
              isocut == kFALSE && totalIso < 10)
        isocut = kTRUE;
      else if(fReverseIsoCut == kTRUE)
        isocut = kFALSE;
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

  LoadEventObject(fElectronBranchName, fElectrons);

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
    Bool_t isocut = PassIsolationCut(e, fElIsoType);
    if (!isocut)
      continue;

    // apply conversion filters
    Bool_t passConvVetoType1 = kFALSE;
    if (fApplyConvFilterType1) {
      LoadEventObject(fConversionBranchName, fConversions);
      passConvVetoType1 = ElectronTools::PassConversionFilter(e, fConversions, 
                                                         fWrongHitsRequirement);      
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
       e->BestTrk()->NExpectedHitsInner() > fNExpectedHitsInnerCut) continue;

    // apply d0 cut
    if (fApplyD0Cut) {
      fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
      Bool_t passD0cut = ElectronTools::PassD0Cut(e, fVertices, fD0Cut, fReverseD0Cut);
      if (!passD0cut)
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

  ReqEventObject(fElectronBranchName, fElectrons, kTRUE);

  if(fCombinedIdCut == kTRUE) {
    fElectronIDType  	  = "NoId";
    fElectronIsoType 	  = "NoIso";
    fApplyConvFilterType1 = kFALSE;
    fApplyConvFilterType2 = kFALSE;
    fApplyD0Cut           = kFALSE;
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
  else if (fElectronIDType.CompareTo("VBTFWorkingPoint95Id") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPoint95Id;
  else if (fElectronIDType.CompareTo("VBTFWorkingPoint90Id") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPoint90Id;
  else if (fElectronIDType.CompareTo("VBTFWorkingPoint80Id") == 0) 
    fElIdType = ElectronTools::kVBTFWorkingPoint80Id;
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


