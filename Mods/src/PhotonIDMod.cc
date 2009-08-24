// $Id: PhotonIDMod.cc,v 1.9 2009/08/21 15:51:53 ceballos Exp $

#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::PhotonIDMod)

//--------------------------------------------------------------------------------------------------
PhotonIDMod::PhotonIDMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPhotonBranchName(Names::gkPhotonBrn),
  fGoodPhotonsName(ModNames::gkGoodPhotonsName),  
  fPhotonIDType("Custom"),
  fPhotonIsoType("CombinedIso"),
  fPhotonPtMin(15.0),
  fHadOverEmMax(0.05),
  fApplyPixelSeed(kTRUE),
  fPhotonR9Min(0.5),
  fPhIdType(kIdUndef),
  fPhIsoType(kIsoUndef),
  fFiduciality(kTRUE),
  fEtaWidthEB(0.013),
  fEtaWidthEE(0.031),
  fPhotons(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void PhotonIDMod::Process()
{
  // Process entries of the tree. 

  LoadEventObject(fPhotonBranchName, fPhotons);

  PhotonOArr *GoodPhotons = new PhotonOArr;
  GoodPhotons->SetName(fGoodPhotonsName);

  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {    
    const Photon *ph = fPhotons->At(i);        

    if (ph->Pt() <= fPhotonPtMin) 
      continue;
    
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
        Double_t totalIso = ph->HollowConeTrkIsoDr04()+
                            ph->EcalRecHitIsoDr04() +
                            ph->HcalTowerSumEtDr04();
        if (totalIso/ph->Pt() < 0.25)
          isocut = kTRUE;
        break;
      case kCustomIso:
      default:
        break;
    }

    if (!isocut) 
      continue;

    if (ph->R9() <= fPhotonR9Min) 
      continue;

    if (fFiduciality == kTRUE &&
        ph->IsEB() == kFALSE && ph->IsEE() == kFALSE) 
      continue;

    if ((ph->IsEB() == kTRUE && ph->SCluster()->EtaWidth() >= fEtaWidthEB) ||
        (ph->IsEE() == kTRUE && ph->SCluster()->EtaWidth() >= fEtaWidthEE))
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

  ReqEventObject(fPhotonBranchName, fPhotons, kTRUE);

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
  else if (fPhotonIsoType.CompareTo("Custom") == 0 ) {
    fPhIsoType = kCustomIso;
    SendError(kWarning, "SlaveBegin",
              "Custom photon isolation is not yet implemented.");
  } else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified photon isolation %s is not defined.",
              fPhotonIsoType.Data());
    return;
  }
}
