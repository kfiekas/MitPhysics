// $Id: PhotonCleaningMod.cc,v 1.1 2008/11/29 18:43:25 sixie Exp $

#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::PhotonCleaningMod)

//--------------------------------------------------------------------------------------------------
PhotonCleaningMod::PhotonCleaningMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCleanElectronsName(ModNames::gkCleanElectronsName),        
  fGoodPhotonsName(ModNames::gkGoodPhotonsName),        
  fCleanPhotonsName(ModNames::gkCleanPhotonsName),
  fMinDeltaRToElectron(0.3)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void PhotonCleaningMod::Process()
{
  // Process entries of the tree.

  // get input collections
  ElectronOArr *CleanElectrons = GetObjThisEvt<ElectronOArr>(fCleanElectronsName);
  PhotonOArr *GoodPhotons = GetObjThisEvt<PhotonOArr>(fGoodPhotonsName);

  // create output collection
  PhotonOArr *CleanPhotons = new PhotonOArr;
  CleanPhotons->SetName(fCleanPhotonsName);

  // Remove any photon that overlaps in eta, phi with an isolated electron.
  for (UInt_t i=0; i<GoodPhotons->GetEntries(); ++i) {
    const Photon *ph = GoodPhotons->At(i);        

    Bool_t isElectronOverlap =  false;
     
    //Check for overlap with an electron
    if (CleanElectrons) {
      for (UInt_t j=0; j<CleanElectrons->Entries(); j++) {
        Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->Mom(),ph->Mom());  
        if (deltaR < fMinDeltaRToElectron) {
          isElectronOverlap = kTRUE;
          break;	 	 
        }      
      }
    }

    if (isElectronOverlap)
      continue;

    CleanPhotons->Add(ph);     
  }

  // sort according to pt
  CleanPhotons->Sort();

  // add to event for other modules to use
  AddObjThisEvt(CleanPhotons);
}
