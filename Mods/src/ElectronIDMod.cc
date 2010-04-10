// $Id: ElectronIDMod.cc,v 1.57 2010/04/10 19:38:45 sixie Exp $

#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
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
  fVertexName(string("PrimaryVertexes").c_str()),
  fElectronIDType("CustomTight"),
  fElectronIsoType("TrackJuraSliding"),
  fElectronPtMin(10),
  fIDLikelihoodCut(0.9),
  fTrackIsolationCut(5.0),
  fCaloIsolationCut(5.0),
  fEcalJuraIsoCut(5.0),
  fHcalIsolationCut(5.0),
  fCombIsolationCut(5.0),
  fApplyConvFilter(kTRUE),
  fWrongHitsRequirement(kTRUE),
  fApplySpikeRemoval(kTRUE),
  fApplyD0Cut(kTRUE),
  fChargeFilter(kTRUE),
  fD0Cut(0.25),
  fReverseIsoCut(kFALSE),
  fReverseD0Cut(kFALSE),
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
        (ele->EcalRecHitIsoDr04() < fEcalJuraIsoCut) &&
        (ele->HcalIsolation() < fHcalIsolationCut);
      break;
    case ElectronTools::kTrackJuraCombined:
      isocut = (ele->TrackIsolationDr03() + ele->EcalRecHitIsoDr04() 
                - 1.5 < fCombIsolationCut);
      break;
    case ElectronTools::kTrackJuraSliding:
    {
      Double_t totalIso = ele->TrackIsolationDr03() + ele->EcalRecHitIsoDr04() - 1.5;
      if (totalIso < (ele->Pt()-10.0)*4.5/20.0 ||
          totalIso <= 0)
        isocut = kTRUE;
      
      if     (fReverseIsoCut == kTRUE &&
              isocut == kFALSE && totalIso < 10)
        isocut = kTRUE;
      else if(fReverseIsoCut == kTRUE)
        isocut = kFALSE;
    }
    break;
    case ElectronTools::kVBTFWorkingPoint95Iso:
      isocut = ElectronTools::PassCustomIso(ele, ElectronTools::kVBTFWorkingPoint95Iso);
      break;
    case ElectronTools::kVBTFWorkingPoint90Iso:
      isocut = ElectronTools::PassCustomIso(ele, ElectronTools::kVBTFWorkingPoint90Iso);
      break;
    case ElectronTools::kVBTFWorkingPoint80Iso:
      isocut = ElectronTools::PassCustomIso(ele, ElectronTools::kVBTFWorkingPoint80Iso);
      break;
    case ElectronTools::kVBTFWorkingPoint70Iso:
      isocut = ElectronTools::PassCustomIso(ele, ElectronTools::kVBTFWorkingPoint70Iso);
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

  ElectronOArr *GoodElectrons = new ElectronOArr;
  GoodElectrons->SetName(fGoodElectronsName);

  for (UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    const Electron *e = fElectrons->At(i);        

    if (e->Pt() <= fElectronPtMin) 
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

    // apply conversion filter
    Bool_t passConvVeto = kFALSE;
    if (fApplyConvFilter) {
      LoadEventObject(fConversionBranchName, fConversions);
      passConvVeto = ElectronTools::PassConversionFilter(e, fConversions, 
                                                         fWrongHitsRequirement);      
    } else {
      passConvVeto = kTRUE;
    }
    if (passConvVeto == kFALSE) continue;
    
    // apply d0 cut
    if (fApplyD0Cut) {
      LoadEventObject(fVertexName, fVertices);
      Bool_t passD0cut = ElectronTools::PassD0Cut(e, fVertices, fD0Cut, fReverseD0Cut);
      if (!passD0cut)
        continue;
    }

    //apply charge filter
    if(fChargeFilter == kTRUE) {
      Bool_t passChargeFilter = ElectronTools::PassChargeFilter(e);
      if (!passChargeFilter) continue;
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

  if (fApplyConvFilter)
    ReqEventObject(fConversionBranchName, fConversions, kTRUE);

  if (fApplyD0Cut)
    ReqEventObject(fVertexName, fVertices, kTRUE);

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


