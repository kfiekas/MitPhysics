// $Id: runSingleLeptonSkim.C,v 1.1 2008/12/10 17:31:35 loizides Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#include "MitAna/PhysicsMod//interface/PlotKineMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/SelMods/interface/GenericSelMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void runSingleLeptonSkim(const char *files, UInt_t nTestEvs=0)
{
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  const char *muInput  = Names::gkMuonBrn;
  const char *elInput  = Names::gkElectronBrn;
  const Double_t ptMin = 20;

  MuonIDMod *muId = new MuonIDMod;
  muId->SetInputName(muInput);
  muId->SetIDType("Loose");
  muId->SetIsoType("TrackCalo");
  muId->SetTrackIsoCut(10.0);
  muId->SetCaloIsoCut(100.0);
  muId->SetPtMin(ptMin);

  ElectronIDMod *elId = new ElectronIDMod;
  elId->SetInputName(elInput);
  elId->SetIDType("NoId");
  elId->SetIsoType("TrackCalo");
  elId->SetTrackIsoCut(10.0);
  elId->SetCaloIsoCut(100.0);
  elId->SetPtMin(ptMin);

  ElectronCleaningMod *elCl = new ElectronCleaningMod;
  elCl->SetGoodElectronsName(elId->GetOutputName());
  elCl->SetCleanMuonsName(muId->GetOutputName());
  elId->Add(elCl);

  MergeLeptonsMod *merger = new MergeLeptonsMod;
  merger->SetMuonsName(muId->GetOutputName());
  merger->SetElectronsName(elCl->GetOutputName());

  GenericSelMod<Particle> *selMod = new GenericSelMod<Particle>;
  selMod->SetPtMin(ptMin);
  selMod->SetColName(merger->GetOutputName());
  merger->Add(selMod);

  OutputMod *outMod = new OutputMod;
  outMod->Keep("*");
  outMod->SetFileName("single_lepton_skim");
  outMod->SetPathName(".");
  selMod->Add(outMod);

  if (nTestEvs) {
    PlotKineMod<Muon> *plotMus = new PlotKineMod<Muon>("PlotMus");
    plotMus->SetPtMin(ptMin);
    plotMus->SetInputName(muId->GetInputName());
    muId->Add(plotMus);

    PlotKineMod<Electron> *plotEls = new PlotKineMod<Electron>("PlotEls");
    plotEls->SetPtMin(ptMin);
    plotEls->SetInputName(elId->GetInputName());
    elCl->Add(plotEls);

    PlotKineMod<Particle> *plotMerged = new PlotKineMod<Particle>("PlotMerged");
    plotMerged->SetLoadBranch(0);
    plotMerged->SetInputName(merger->GetOutputName());
    merger->Add(plotMerged);
  }

  // set up analysis
  Analysis *ana = new Analysis;
  ana->AddSuperModule(muId);
  ana->AddSuperModule(elId);
  ana->AddSuperModule(merger);

  if (nTestEvs)
    ana->SetProcessNEvents(nTestEvs);
  ana->AddFile(files);
  ana->Run(nTestEvs>0);
}
