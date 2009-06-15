// $Id: TauCleaningMod.cc,v 1.3 2009/04/09 08:45:49 loizides Exp $

#include "MitPhysics/Mods/interface/TauCleaningMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/CaloTauCol.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::TauCleaningMod)

//--------------------------------------------------------------------------------------------------
TauCleaningMod::TauCleaningMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCleanElectronsName(ModNames::gkCleanElectronsName),        
  fCleanMuonsName(ModNames::gkCleanMuonsName),        
  fGoodTausName(ModNames::gkGoodTausName),        
  fCleanCaloTausName(ModNames::gkCleanTausName),
  fMinDeltaRToElectron(0.3),
  fMinDeltaRToMuon(0.3)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void TauCleaningMod::Process()
{
  // Process entries of the tree.

  // get input collections
  const CaloTauCol  *GoodTaus = GetObjThisEvt<CaloTauCol>(fGoodTausName);
  const ElectronCol *CleanElectrons = GetObjThisEvt<ElectronCol>(fCleanElectronsName);
  const MuonCol *CleanMuons = GetObjThisEvt<MuonCol>(fCleanMuonsName);

  // create output collection
  CaloTauOArr *CleanCaloTaus = new CaloTauOArr;
  CleanCaloTaus->SetName(fCleanCaloTausName);

  // remove any Tau that overlaps in eta, phi with an isolated electron.
  UInt_t n = GoodTaus->GetEntries();
  for (UInt_t i=0; i<n; ++i) {
    const CaloTau *tau = GoodTaus->At(i);        

    Bool_t isElectronOverlap =  false;
     
    // check for overlap with an electron
    if (CleanElectrons) {
      UInt_t n1 = CleanElectrons->GetEntries();
      for (UInt_t j=0; j<n1; j++) {
        Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->Mom(),
	                                    tau->SourceCaloJet()->Mom());  
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
	                                    tau->SourceCaloJet()->Mom());  
        if (deltaR < fMinDeltaRToMuon) {
          isMuonOverlap = kTRUE;
          break;	 	 
        }      
      }
    }

    if (isMuonOverlap)
      continue;

    CleanCaloTaus->Add(tau);     
  }

  // sort according to pt
  CleanCaloTaus->Sort();

  // add to event for other modules to use
  AddObjThisEvt(CleanCaloTaus);
}
