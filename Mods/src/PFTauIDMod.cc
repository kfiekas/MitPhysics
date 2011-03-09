// $Id: PFTauIDMod.cc,v 1.3 2011/02/17 14:08:49 ceballos Exp $

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
  fIsHPSSel(kFALSE),
  fHPSIso("loose"),
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

    if(fIsHPSSel == kFALSE) {
      if (tau->NSignalPFCands() == 0)
	continue;

      CompositeParticle tauSystem;
      CompositeParticle tauChargedSystem;
      UInt_t nTrk = 0;
      for (UInt_t j=0; j<tau->NSignalPFCands(); ++j) {
	tauSystem.AddDaughter(tau->SignalPFCand(j));
	if(tau->SignalPFCand(j) != 0 &&
	   tau->SignalPFCand(j)->Charge() != 0){
	  nTrk++;
	  tauChargedSystem.AddDaughter(tau->SignalPFCand(j));
	}
      }

      if (nTrk != 1 && nTrk != 3)
	continue;

      if (TMath::Abs(tau->Charge()) - 1.0 > 0.0001)
	continue;

      if (tauSystem.Pt() <= fPtMin)
	continue;

      if(!tau->LeadChargedHadronPFCand())
	continue;

      if(tau->LeadChargedHadronPFCand()->Pt() <= fPtLeadChargedHadronPFCandMin)
	continue;

      if (tau->IsoChargedHadronPtSum() >= fIsoChargedHadronPtSumMax)
	continue;

      if (tau->IsoGammaEtSum() >= fIsoGammaEtSumMax)
	continue;

      if (tauChargedSystem.Mass() <= fSignalMassMin || tauChargedSystem.Mass() >= fSignalMassMax)
	continue;
    }
    else { // if we're doing hps selection:
      if(tau->Pt() <= fPtMin)
	continue;
      if(tau->DiscriminationByDecayModeFinding() < 0.5)
        continue;
      if(tau->DiscriminationAgainstElectron() < 0.5)
        continue;
      if(tau->DiscriminationAgainstMuon() < 0.5)
        continue;
      // "loose" should be default, but others can be used
      if(fHPSIso.Contains("loose",TString::kIgnoreCase)) {
	if(tau->DiscriminationByLooseIsolation() < 0.5)
	  continue;
      }
      else if(fHPSIso.Contains("med",TString::kIgnoreCase)) {
	if(tau->DiscriminationByMediumIsolation() < 0.5)
	  continue;
      }
      else if(fHPSIso.Contains("tight",TString::kIgnoreCase)) {
	if(tau->DiscriminationByTightIsolation() < 0.5)
	  continue;
      }
      else {
	SendError(kWarning,"Process","ERROR: HPS Isolation not properly defined!");
      }
    }

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
  if(fIsHPSSel) {
    fPFTausName = TString(Names::gkHPSTauBrn);
  }
  else
    fPFTausName = TString(Names::gkPFTauBrn);

  ReqBranch(fPFTausName, fPFTaus);
}
