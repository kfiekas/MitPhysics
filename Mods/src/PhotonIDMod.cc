// $Id: PhotonIDMod.cc,v 1.8 2008/11/28 11:35:29 ceballos Exp $

#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::PhotonIDMod)

//--------------------------------------------------------------------------------------------------
PhotonIDMod::PhotonIDMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPhotonBranchName(Names::gkPhotonBrn),
  fGoodPhotonsName(ModNames::gkGoodPhotonsName),  
  fPhotonIDType("Tight"),
  fPhotonIsoType("NoIso"),
  fPhotonPtMin(15.0),
  fPhotons(0),
  fPhIdType(kIdUndef),
  fPhIsoType(kIsoUndef)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void PhotonIDMod::Process()
{
  // Process entries of the tree. 

  LoadBranch(fPhotonBranchName);

  PhotonOArr *GoodPhotons = new PhotonOArr;
  GoodPhotons->SetName(fGoodPhotonsName);

  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {    
    const Photon *ph = fPhotons->At(i);        

    if (ph->Pt() <= fPhotonPtMin) 
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
      case kCustomIso:
      default:
        break;
    }

    if (!isocut) 
      continue;

    // add good electron
    GoodPhotons->Add(fPhotons->At(i));
  }

  // add to event for other modules to use
  AddObjThisEvt(GoodPhotons);  
}

//--------------------------------------------------------------------------------------------------
void PhotonIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the electron collection branch.

  ReqBranch(fPhotonBranchName, fPhotons);

  if (fPhotonIDType.CompareTo("Tight") == 0) 
    fPhIdType = kTight;
  else if (fPhotonIDType.CompareTo("Loose") == 0) 
    fPhIdType = kLoose;
  else if (fPhotonIDType.CompareTo("LooseEM") == 0) 
    fPhIdType = kLooseEM;
  else if (fPhotonIDType.CompareTo("Custom") == 0) {
    fPhIdType = kCustomId;
    SendError(kWarning, "SlaveBegin",
              "Custom electron identification is not yet implemented.");
  } else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified electron identification %s is not defined.",
              fPhotonIDType.Data());
    return;
  }

  if (fPhotonIsoType.CompareTo("NoIso") == 0 )
    fPhIsoType = kNoIso;
  else if (fPhotonIsoType.CompareTo("Custom") == 0 ) {
    fPhIsoType = kCustomIso;
    SendError(kWarning, "SlaveBegin",
              "Custom electron isolation is not yet implemented.");
  } else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified electron isolation %s is not defined.",
              fPhotonIsoType.Data());
    return;
  }
}
