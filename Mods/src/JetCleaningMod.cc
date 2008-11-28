// $Id: JetCleaningMod.cc,v 1.3 2008/11/27 16:30:27 loizides Exp $

#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::JetCleaningMod)

//--------------------------------------------------------------------------------------------------
JetCleaningMod::JetCleaningMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCleanElectronsName(ModNames::gkCleanElectronsName),        
  fGoodJetsName(ModNames::gkGoodJetsName),        
  fCleanJetsName(ModNames::gkCleanJetsName)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void JetCleaningMod::Process()
{
  // Process entries of the tree.

  // get input collections
  ElectronOArr *CleanElectrons = GetObjThisEvt<ElectronOArr>(fCleanElectronsName);
  JetOArr *GoodJets = GetObjThisEvt<JetOArr>(fGoodJetsName);

  // create output collection
  JetOArr *CleanJets = new JetOArr;
  CleanJets->SetName(fCleanJetsName);

  // Remove any jet that overlaps in eta, phi with an isolated electron.
  for (UInt_t i=0; i<GoodJets->GetEntries(); ++i) {
    const Jet *jet = GoodJets->At(i);        

    Bool_t isElectronOverlap =  false;
     
    //Check for overlap with an electron
    for (UInt_t j=0; j<CleanElectrons->Entries(); j++) {
      Double_t deltaR = MathUtils::DeltaR(CleanElectrons->At(j)->Mom(),jet->Mom());  
      if (deltaR < 0.3) {
	isElectronOverlap = kTRUE;
	break;	 	 
      }      
    }

    if (isElectronOverlap)
      continue;

    CleanJets->Add(jet);     
  }

  // add to event for other modules to use
  AddObjThisEvt(CleanJets);
}
