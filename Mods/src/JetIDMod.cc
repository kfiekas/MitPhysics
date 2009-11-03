// $Id: JetIDMod.cc,v 1.18 2009/11/03 08:39:19 ceballos Exp $

#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"

using namespace mithep;

ClassImp(mithep::JetIDMod)

//--------------------------------------------------------------------------------------------------
JetIDMod::JetIDMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fJetsName(ModNames::gkPubJetsName),
  fGoodJetsName(ModNames::gkGoodJetsName),  
  fUseJetCorrection(kTRUE),
  fJetPtCut(35.0),
  fJetEtaMaxCut(5.0),
  fJetEEMFractionMinCut(0.01)
{
  // Constructor.
}

//-------------------------------------------------------------------------------------------------
void JetIDMod::Process()
{
  // Process entries of the tree. 

  const JetCol *inJets = GetObjThisEvt<JetCol>(fJetsName);
  if (!inJets) {
    SendError(kAbortModule, "Process", 
              "Pointer to input jet collection %s is null.",
              fJetsName.Data());
    return;
  }

  JetOArr *GoodJets = new JetOArr; 
  GoodJets->SetName(fGoodJetsName);

  // loop over jets
  for (UInt_t i=0; i<inJets->GetEntries(); ++i) {
    const Jet *jet = inJets->At(i);

    if (jet->AbsEta() > fJetEtaMaxCut) 
      continue;
    
    Double_t jetpt;
    if (fUseJetCorrection)
      jetpt = jet->Pt();
    else
      jetpt = jet->RawMom().Pt();

    if (jetpt < fJetPtCut)
      continue;
    
    Bool_t passEEMFractionMinCut = kTRUE;
    if(fJetEEMFractionMinCut > 0){
      const CaloJet *caloJet = dynamic_cast<const CaloJet*>(jet->MakeCopy());
      // The 2.6 value is hardcoded, no reason to change that value in CMS
      if (caloJet && caloJet->AbsEta() < 2.6 &&
          caloJet->EnergyFractionEm() < fJetEEMFractionMinCut)
        passEEMFractionMinCut = kFALSE;
      delete caloJet;
    }
    if(passEEMFractionMinCut == kFALSE)
      continue;

    // add good jet to collection
    GoodJets->Add(jet);             
  }

  // sort according to pt
  GoodJets->Sort();
  
  // add to event for other modules to use
  AddObjThisEvt(GoodJets);  
}

