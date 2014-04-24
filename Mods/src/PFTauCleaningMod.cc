#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/PFTauCleaningMod.h"

using namespace mithep;

ClassImp(mithep::PFTauCleaningMod)

//--------------------------------------------------------------------------------------------------
PFTauCleaningMod::PFTauCleaningMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCleanElectronsName(ModNames::gkCleanElectronsName),        
  fCleanMuonsName(ModNames::gkCleanMuonsName),        
  fGoodPFTausName(ModNames::gkGoodPFTausName),        
  fCleanPFTausName(ModNames::gkCleanPFTausName),
  fMinDeltaRToElectron(0.3),
  fMinDeltaRToMuon(0.3)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void PFTauCleaningMod::Process()
{
  // Process entries of the tree.

  // get input collections
  const PFTauCol *GoodPFTaus = GetObjThisEvt<PFTauCol>(fGoodPFTausName);
  const ElectronCol *CleanElectrons = GetObjThisEvt<ElectronCol>(fCleanElectronsName);
  const MuonCol *CleanMuons = GetObjThisEvt<MuonCol>(fCleanMuonsName);

  // create output collection
  PFTauOArr *CleanPFTaus = new PFTauOArr;
  CleanPFTaus->SetName(fCleanPFTausName);

  // remove any Tau that overlaps in eta, phi with an isolated electron.
  UInt_t n = GoodPFTaus->GetEntries();
  for (UInt_t i=0; i<n; ++i) {
    const PFTau *tau = GoodPFTaus->At(i);        

    Bool_t isElectronOverlap =  false;
     
    // check for overlap with an electron
    if (CleanElectrons) {
      UInt_t n1 = CleanElectrons->GetEntries();
      for (UInt_t j=0; j<n1; j++) {
        Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->Mom(),
	                                    tau->Mom());  
        if (deltaR < fMinDeltaRToElectron) {
          isElectronOverlap = kTRUE;
          break;	 	 
        }      
      }
    }

    if (isElectronOverlap)
      continue;

    Bool_t isMuonOverlap =  false;
     
    // check for overlap with an Muon
    if (CleanMuons) {
      UInt_t n2 = CleanMuons->GetEntries();
      for (UInt_t j=0; j<n2; j++) {
        Double_t deltaR = MathUtils::DeltaR(CleanMuons->At(j)->Mom(),
	                                    tau->Mom());  
        if (deltaR < fMinDeltaRToMuon) {
          isMuonOverlap = kTRUE;
          break;	 	 
        }      
      }
    }

    if (isMuonOverlap)
      continue;

    CleanPFTaus->Add(tau);     
  }

  // sort according to pt
  CleanPFTaus->Sort();

  // add to event for other modules to use
  AddObjThisEvt(CleanPFTaus);
}
