#include "MitAna/DataTree/interface/CaloTauCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/TauIDMod.h"

using namespace mithep;

ClassImp(mithep::TauIDMod)

//------------------------------------------------------------------------------
TauIDMod::TauIDMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCaloTausName(Names::gkCaloTauBrn),
  fGoodTausName(ModNames::gkGoodTausName),
  fPtMin(10.0),
  fJetPtMin(20.0),
  fLeadTrackPtMin(6.0),
  fNSignalTracksMax(3),
  fNIsoTracksMax(0),
  fSignalTracksMassMax(2.0),
  fIsoTrackPtSumMax(5.0),
  fEnergyFractionEmMax(0.95),
  fHCalEtOverPtMin(0.1),
  fCaloTaus(0)
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

    if (tau->Pt() <= fPtMin)
      continue;
    
    if (!tau->SourceCaloJet() || tau->SourceCaloJet()->Pt() <= fJetPtMin)
      continue;
    
    if (tau->NSignalTracks() == 0 || tau->NSignalTracks() > fNSignalTracksMax)
      continue;
    
    if (!tau->LeadTrack() || tau->LeadTrack()->Pt() < fLeadTrackPtMin)
      continue;
    
    if (tau->NIsoTracks() > fNIsoTracksMax)
      continue;

    if (tau->SignalTracksMass() > fSignalTracksMassMax)
      continue;

    if (tau->IsoTrackPtSum() > fIsoTrackPtSumMax)
      continue;

    if (tau->SourceCaloJet()->EnergyFractionEm() > fEnergyFractionEmMax)
      continue;

    if (tau->LeadTrack3x3HCalEt()/tau->LeadTrack()->Pt() < fHCalEtOverPtMin)
      continue;

    // always apply these requirements
    if (TMath::Abs(tau->Charge()) != 1)
      continue;

    if (TMath::Abs(TMath::Abs(tau->Eta())-1.53) < 0.02)
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
void TauIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the Tau collection branch.

  ReqBranch(fCaloTausName, fCaloTaus);
}
