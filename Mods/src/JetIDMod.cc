// $Id: JetIDMod.cc,v 1.7 2008/12/10 11:44:33 loizides Exp $

#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::JetIDMod)

//--------------------------------------------------------------------------------------------------
JetIDMod::JetIDMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fJetBranchName(Names::gkCaloJetBrn),
  fGoodJetsName(ModNames::gkGoodJetsName),  
  fUseJetCorrection(kTRUE),
  fJetEtCut(35.0),
  fJets(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void JetIDMod::Process()
{
  // Process entries of the tree. 

  LoadBranch(fJetBranchName);

  JetOArr *GoodJets = new JetOArr; 
  GoodJets->SetName(fGoodJetsName);

  // loop over jets
  for (UInt_t i=0; i<fJets->GetEntries(); ++i) {
    const Jet *jet = fJets->At(i);

    if (jet->AbsEta() > 5.0) 
      continue;
    
    Double_t jetet = jet->Et();
    if (fUseJetCorrection)
      jetet *= jet->L2RelativeCorrectionScale() * jet->L3AbsoluteCorrectionScale();

    if (jetet < fJetEtCut)
      continue;
    
    // add good jet to collection
    GoodJets->Add(jet);             
  }

  // sort according to pt
  GoodJets->Sort();
  
  // add to event for other modules to use
  AddObjThisEvt(GoodJets);  
}

//--------------------------------------------------------------------------------------------------
void JetIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the jet collection branch.

  ReqBranch(fJetBranchName, fJets);
}
