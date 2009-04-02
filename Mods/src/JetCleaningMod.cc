// $Id: JetCleaningMod.cc,v 1.9 2009/03/12 16:00:46 bendavid Exp $

#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::JetCleaningMod)

//--------------------------------------------------------------------------------------------------
JetCleaningMod::JetCleaningMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCleanElectronsName(ModNames::gkCleanElectronsName),        
  fCleanMuonsName(ModNames::gkCleanMuonsName),        
  fCleanPhotonsName(ModNames::gkCleanPhotonsName),        
  fGoodJetsName(ModNames::gkGoodJetsName),        
  fCleanJetsName(ModNames::gkCleanJetsName),
  fMinDeltaRToElectron(0.3),
  fMinDeltaRToMuon(0.3),
  fMinDeltaRToPhoton(0.3)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void JetCleaningMod::Process()
{
  // Process entries of the tree.

  // get input collections
  const JetCol  *GoodJets       = GetObjThisEvt<JetCol>(fGoodJetsName);
  const ElectronCol *CleanElectrons = 0;
  if (!fCleanElectronsName.IsNull())
    CleanElectrons = GetObjThisEvt<ElectronCol>(fCleanElectronsName);
  const MuonCol *CleanMuons = 0;
  if (!fCleanMuonsName.IsNull())
    CleanMuons = GetObjThisEvt<MuonCol>(fCleanMuonsName);
  const PhotonCol   *CleanPhotons   = 0;
  if (!fCleanPhotonsName.IsNull())
  CleanPhotons    = GetObjThisEvt<PhotonCol>(fCleanPhotonsName);

  // create output collection
  JetOArr *CleanJets = new JetOArr;
  CleanJets->SetName(fCleanJetsName);

  // remove any jet that overlaps in eta, phi with an isolated electron.    
  for (UInt_t i=0; i<GoodJets->GetEntries(); ++i) {
    const Jet *jet = GoodJets->At(i);        

    // check for overlap with an Electron
    Bool_t isElectronOverlap = kFALSE;
    if (CleanElectrons) {
      UInt_t n = CleanElectrons->GetEntries();
      for (UInt_t j=0; j<n; ++j) {
        Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->Mom(),jet->Mom());  
        if (deltaR < fMinDeltaRToElectron) {
          isElectronOverlap = kTRUE;
          break;	 	 
        }      
      }
    }

    if (isElectronOverlap) continue;

    // check for overlap with an Muon
    Bool_t isMuonOverlap = kFALSE;
    if (CleanMuons) {
      UInt_t n = CleanMuons->GetEntries();
      for (UInt_t j=0; j<n; ++j) {
        Double_t deltaR = MathUtils::DeltaR(CleanMuons->At(j)->Mom(),jet->Mom());  
        if (deltaR < fMinDeltaRToMuon) {
          isMuonOverlap = kTRUE;
          break;	 	 
        }      
      }
    }

    if (isMuonOverlap) continue;

    // check for overlap with a photon
    Bool_t isPhotonOverlap = kFALSE;
    if (CleanPhotons) {
      UInt_t n = CleanPhotons->GetEntries();
      for (UInt_t j=0; j<n; ++j) {
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

  // sort according to pt
  CleanJets->Sort();

  // add to event for other modules to use
  AddObjThisEvt(CleanJets);
}
