// $Id: JetCorrectionMod.cc,v 1.17 2009/06/15 15:00:21 loizides Exp $

#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::JetCorrectionMod)

//--------------------------------------------------------------------------------------------------
JetCorrectionMod::JetCorrectionMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fJetsName(ModNames::gkPubJetsName),
  fCorrectedJetsName(ModNames::gkCorrectedJetsName),  
  fJetCorrector(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
JetCorrectionMod::~JetCorrectionMod()
{
  if (fJetCorrector) {
    delete fJetCorrector;
    fJetCorrector = 0;
  }
} 

//--------------------------------------------------------------------------------------------------
void JetCorrectionMod::SlaveBegin()
{
  //initialize jet corrector for L2 and L3 corrections with the configured tag
  std::string levels = "L2:L3";
  fJetCorrector = new CombinedJetCorrector(levels,std::string(fCorrectionTag));
}

//--------------------------------------------------------------------------------------------------
void JetCorrectionMod::Process()
{
  // Process entries of the tree. 

  const JetCol *inJets = GetObjThisEvt<JetCol>(fJetsName);
  if (!inJets) {
    SendError(kAbortModule, "Process", 
              "Pointer to input jet collection %s is null.",
              fJetsName.Data());
    return;
  }

  JetOArr *CorrectedJets = new JetOArr;
  CorrectedJets->SetOwner(kTRUE);
  CorrectedJets->SetName(fCorrectedJetsName);

  std::vector<double> corrections;

  // loop over jets
  for (UInt_t i=0; i<inJets->GetEntries(); ++i) {
    const Jet *inJet = inJets->At(i);

    //copy input jet, using special function to copy full derived class
    Jet *jet = inJet->MakeCopy();

    //cache uncorrected momentum
    const FourVectorM rawMom = jet->RawMom();
    
    //compute correction factors
    corrections = fJetCorrector->getSubCorrections(rawMom.Pt(),rawMom.Eta(),rawMom.E());
    Double_t l2Factor = corrections.at(0);
    Double_t l3Factor = corrections.at(1)/l2Factor;

    //set and enable correction factors in the output jet
    jet->SetL2RelativeCorrectionScale(l2Factor);
    jet->SetL3AbsoluteCorrectionScale(l3Factor);     
    jet->EnableCorrection(mithep::Jet::L2);
    jet->EnableCorrection(mithep::Jet::L3);  

    if (0) {
      printf("In L2 = %5f, Out L2 = %5f\n",inJet->L2RelativeCorrectionScale(),jet->L2RelativeCorrectionScale());
      printf("In L3 = %5f, Out L3 = %5f\n",inJet->L3AbsoluteCorrectionScale(),jet->L3AbsoluteCorrectionScale());
      printf("In RawPt = %5f, Out RawPt = %5f\n",inJet->RawMom().Pt(),jet->RawMom().Pt());
      printf("In Pt = %5f, Out Pt = %5f\n",inJet->Pt(),jet->Pt());
    }

    // add corrected jet to collection
    CorrectedJets->AddOwned(jet);             
  }

  // sort according to pt
  CorrectedJets->Sort();
  
  // add to event for other modules to use
  AddObjThisEvt(CorrectedJets);
}

