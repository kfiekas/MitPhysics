// $Id: MergeLeptonsMod.cc,v 1.3 2009/06/15 15:00:21 loizides Exp $

#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "TH1D.h"

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
void mithep::MergeLeptonsMod::BeginRun()
{
  
}

//--------------------------------------------------------------------------------------------------
void mithep::MergeLeptonsMod::Process()
{
  // Merge the two input collections and publish merged collection. 

  fElIn = GetObjThisEvt<ElectronCol>(fElName);
  fMuIn = GetObjThisEvt<MuonCol>(fMuName);
  if(fElIn){;
  fRecoWElectrons->Fill(fElIn->GetEntries());
  }
  if(fMuIn){
  fRecoWMuons->Fill(fMuIn->GetEntries());
  }
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

//--------------------------------------------------------------------------------------------------
void mithep::MergeLeptonsMod::SlaveBegin()
{
	AddTH1(fRecoWMuons,"fRecoWMuons","Number of Reconstructed Muons;N_{muons};#",10,-0.5,9.5);
	AddTH1(fRecoWElectrons,"fRecoWElectrons","Number of Reconstructed Electrons;N_{electrons};#",10,-0.5,9.5);
}
//--------------------------------------------------------------------------------------------------
void mithep::MergeLeptonsMod::SlaveTerminate()
{

}

