// $Id: TauIDMod.cc,v 1.1 2009/04/08 10:11:44 ceballos Exp $

#include "MitPhysics/Mods/interface/TauIDMod.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::TauIDMod)

//------------------------------------------------------------------------------
TauIDMod::TauIDMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCaloTausName(Names::gkCaloTauBrn),
  fCaloTaus(0),
  fGoodTausName(ModNames::gkGoodTausName),
  fTauPtMin(10.0),
  fTauJetPtMin(20.0),
  fNSignalTracksMax(3),
  fNIsoTracksMax(0),
  fSignalTracksMassMax(2.0),
  fIsoTrackPtSumMax(5.0)
{
  // Constructor.
}

//------------------------------------------------------------------------------
void TauIDMod::Process()
{
  // Process entries of the tree. 

  LoadBranch(fCaloTausName);

  CaloTauOArr *GoodTaus = new CaloTauOArr;
  GoodTaus->SetName(fGoodTausName);

  for (UInt_t i=0; i<fCaloTaus->GetEntries(); ++i) {    
    const CaloTau *tau = fCaloTaus->At(i);        

    if (tau->Pt() <= fTauPtMin)
      continue;
    
    if (tau->SourceCaloJet()->Pt() <= fTauJetPtMin)
      continue;
    
    if (tau->NSignalTracks() == 0 || tau->NSignalTracks() > fNSignalTracksMax)
      continue;
    
    if (tau->NIsoTracks() > fNIsoTracksMax)
      continue;

    if (tau->SignalTracksMass() > fSignalTracksMassMax)
      continue;

    if (tau->IsoTrackPtSum() > fIsoTrackPtSumMax)
      continue;

    // Always apply this requirement
    if (TMath::Abs(tau->Charge()) != 1)
      continue;

    // add good electron
    GoodTaus->Add(tau);
  }

  // sort according to pt
  GoodTaus->Sort();

  // add to event for other modules to use
  AddObjThisEvt(GoodTaus);  
}

//------------------------------------------------------------------------------
void TauIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the Tau collection branch.

  ReqBranch(fCaloTausName, fCaloTaus);

}
