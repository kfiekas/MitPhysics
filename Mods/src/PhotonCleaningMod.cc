// $Id: PhotonCleaningMod.cc,v 1.7 2009/11/02 13:39:23 ceballos Exp $

#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
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
  const PhotonCol   *GoodPhotons    = GetObjThisEvt<PhotonCol>(fGoodPhotonsName);
  const ElectronCol *CleanElectrons = GetObjThisEvt<ElectronCol>(fCleanElectronsName);

  // create output collection
  PhotonOArr *CleanPhotons = new PhotonOArr;
  CleanPhotons->SetName(fCleanPhotonsName);

  // remove any photon that overlaps in eta, phi with an isolated electron.
  UInt_t n = GoodPhotons->GetEntries();
  for (UInt_t i=0; i<n; ++i) {
    const Photon *ph = GoodPhotons->At(i);        

    Bool_t isElectronOverlap =  false;
     
    // check for overlap with an electron
    if (CleanElectrons) {
      UInt_t n2 = CleanElectrons->GetEntries();
      for (UInt_t j=0; j<n2; j++) {
        Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->SCluster()->Phi(), 
                                            CleanElectrons->At(j)->SCluster()->Eta(),
                                             ph->SCluster()->Phi(),ph->SCluster()->Eta());  
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
