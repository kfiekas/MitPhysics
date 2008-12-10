// $Id: JetIDMod.cc,v 1.6 2008/11/28 13:07:38 ceballos Exp $

#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::MergeLeptonsMod)

//--------------------------------------------------------------------------------------------------
mithep::MergeLeptonsMod::MergeLeptonsMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fElName(ModNames::gkCleanElectronsName),
  fMuName(ModNames::gkCleanMuonsName),
  fMergedName(ModNames::gkMergedLeptonsName),
  fElIn(0),
  fMuIn(0),
  fColOut(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void mithep::MergeLeptonsMod::Process()
{
  // Merge the two input collections and publish merged collection. 

  fElIn = GetObjThisEvt<ElectronCol>(fElName);
  fMuIn = GetObjThisEvt<MuonCol>(fElName);

  UInt_t nents = 0;
  if (fElIn) 
    nents += fElIn->GetEntries();
  if (fMuIn) 
    nents += fMuIn->GetEntries();

  fColOut = new mithep::ParticleOArr(nents, GetMergedName());

  if (fElIn)
    fColOut->Add(fElIn);
  if (fMuIn)
    fColOut->Add(fMuIn);

  // sort according to pt
  fColOut->Sort();

  // add to event for other modules to use
  AddObjThisEvt(fColOut);
}
