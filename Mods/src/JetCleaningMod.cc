// $Id: JetCleaningMod.cc,v 1.5 2008/11/29 18:45:43 sixie Exp $

#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::JetCleaningMod)

//--------------------------------------------------------------------------------------------------
JetCleaningMod::JetCleaningMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCleanElectronsName(ModNames::gkCleanElectronsName),        
  fCleanPhotonsName(ModNames::gkCleanPhotonsName),        
  fGoodJetsName(ModNames::gkGoodJetsName),        
  fCleanJetsName(ModNames::gkCleanJetsName),
  fMinDeltaRToElectron(0.3),
  fMinDeltaRToPhoton(0.3)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void JetCleaningMod::Process()
{
  // Process entries of the tree.

  // get input collections
  ElectronOArr *CleanElectrons = GetObjThisEvt<ElectronOArr>(fCleanElectronsName);
  PhotonOArr *CleanPhotons = GetObjThisEvt<PhotonOArr>(fCleanPhotonsName);
  JetOArr *GoodJets = GetObjThisEvt<JetOArr>(fGoodJetsName);

  // create output collection
  JetOArr *CleanJets = new JetOArr;
  CleanJets->SetName(fCleanJetsName);

  // Remove any jet that overlaps in eta, phi with an isolated electron.    
  for (UInt_t i=0; i<GoodJets->GetEntries(); ++i) {
    const Jet *jet = GoodJets->At(i);        

    Bool_t isElectronOverlap = kFALSE;
    Bool_t isPhotonOverlap = kFALSE;

    //Check for overlap with an electron
    if (CleanElectrons) {
      for (UInt_t j=0; j<CleanElectrons->Entries(); j++) {
        Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->Mom(),jet->Mom());  
        if (deltaR < fMinDeltaRToElectron) {
          isElectronOverlap = kTRUE;
          break;	 	 
        }      
      }
    }

    if (isElectronOverlap) continue;

    //Check for overlap with a photon
    if (CleanPhotons) {
      for (UInt_t j=0; j<CleanPhotons->Entries(); j++) {
        Double_t deltaR = MathUtils::DeltaR(CleanPhotons->At(j)->Mom(),jet->Mom());  
        if (deltaR < fMinDeltaRToPhoton) {
          isPhotonOverlap = kTRUE;
          break;	 	 
        }      
      }
    }

    if (isPhotonOverlap) continue;

    CleanJets->Add(jet);     
  }

  // add to event for other modules to use
  AddObjThisEvt(CleanJets);
}
