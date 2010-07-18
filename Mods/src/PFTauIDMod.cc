// $Id: PFTauIDMod.cc,v 1.9 2010/07/06 08:33:16 sixie Exp $

#include "MitPhysics/Mods/interface/PFTauIDMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataTree/interface/PFTauCol.h"

using namespace mithep;

ClassImp(mithep::PFTauIDMod)

//------------------------------------------------------------------------------
PFTauIDMod::PFTauIDMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPFTausName(Names::gkPFTauBrn),
  fGoodPFTausName(ModNames::gkGoodPFTausName),
  fPtMin(15.0),
  fPtLeadChargedHadronPFCandMin(5.0),
  fIsoChargedHadronPtSumMax(2.0),
  fIsoGammaEtSumMax(3.0),
  fSignalMassMin(0.13),
  fSignalMassMax(2.00),
  fPFTaus(0)
{
  // Constructor.
}

//------------------------------------------------------------------------------
void PFTauIDMod::Process()
{
  // Process entries of the tree. 

  LoadBranch(fPFTausName);

  PFTauOArr *GoodTaus = new PFTauOArr;
  GoodTaus->SetName(fGoodPFTausName);

  for (UInt_t i=0; i<fPFTaus->GetEntries(); ++i) {    
    const PFTau *tau = fPFTaus->At(i);        

    if (tau->NSignalPFCands() == 0)
      continue;

    CompositeParticle tauSystem;
    CompositeParticle tauChargedSystem;
    UInt_t nTrk = 0;
    for (UInt_t j=0; j<tau->NSignalPFCands(); ++j) {
      tauSystem.AddDaughter(tau->SignalPFCand(j));
      if(tau->SignalPFCand(j)->Charge() != 0){
        nTrk++;
        tauChargedSystem.AddDaughter(tau->SignalPFCand(j));
      }
    }

    if (tauSystem.Pt() <= fPtMin)
      continue;
    
    if (nTrk != 1 && nTrk != 3)
      continue;

    if(!tau->LeadChargedHadronPFCand())
      continue;

    if(tau->LeadChargedHadronPFCand()->Pt() <= fPtLeadChargedHadronPFCandMin)
      continue;

    if (TMath::Abs(tau->Charge()) != 1)
      continue;

    if (tau->IsoChargedHadronPtSum() >= fIsoChargedHadronPtSumMax)
      continue;

    if (tau->IsoGammaEtSum() >= fIsoGammaEtSumMax)
      continue;

    if (tauChargedSystem.Mass() <= fSignalMassMin || tauChargedSystem.Mass() >= fSignalMassMax)
      continue;

    // add good tau to output collection
    GoodTaus->Add(tau);
  }

  // sort according to pt
  GoodTaus->Sort();

  // add to event for other modules to use
  AddObjThisEvt(GoodTaus);  
}

//------------------------------------------------------------------------------
void PFTauIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the Tau collection branch.

  ReqBranch(fPFTausName, fPFTaus);
}
