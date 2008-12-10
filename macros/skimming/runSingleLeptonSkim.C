// $Id: skimtest.C,v 1.1 2008/12/04 14:14:59 loizides Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#include "MitAna/PhysicsMod//interface/PlotKineMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
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

  const Double_t ptmin = 20;

  MuonIDMod *muId = new MuonIDMod;  
  muId->SetIDType("Loose");
  muId->SetIsoType("TrackCalo");
  muId->SetTrackIsoCut(3.0);
  muId->SetCaloIsoCut(3.0);
  muId->SetPtMin(ptmin);

  ElectronIDMod *elId = new ElectronIDMod;
  elId->SetIDType("Loose");
  elId->SetIsoType("TrackCalo");
  elId->SetTrackIsoCut(3.0); //def is 5?
  elId->SetCaloIsoCut(3.0);  //def is 5?
  elId->SetPtMin(ptmin);

  ElectronCleaningMod *elCl = new ElectronCleaningMod;
  elCl->SetGoodElectronsName(elId->GetOutputName());
  elCl->SetCleanMuonsName(muId->GetCleanName());
  elId->Add(elCl);

  MergeLeptonsMod *merger = new MergeLeptonsMod;
  merger->SetMuonsName(muId->GetOutputName());
  merger->SetElectronsName(elCl->GetOutputName());

  GenericSelMod<Particle> *selmod = new GenericSelMod<Particle>;
  selmod->SetPtMin(ptmin);
  selmod->SetColName(merger->GetOutputName());
  merger->Add(selmod);

  OutputMod *omod = new OutputMod;
  omod->Keep("*");
  omod->SetFileName("single_lepton_skim");
  omod->SetPathName(".");
  selmod->Add(omod);

  if (nTestEvs) {
    PlotKineMod<Muon> *plotmus = new PlotKineMod<Muon>("PlotMuons");
    plotmus->SetPtMin(ptmin);
    plotmus->SetInputName(muId->GetInputName());
    muId->Add(plotmus);

    PlotKineMod<Electron> *plotels = new PlotKineMod<Electron>("PlotEls");
    plotels->SetPtMin(ptmin);
    plotels->SetInputName(elId->GetInputName());
    elCl->Add(plotels);
    
    PlotKineMod<Particle> *plotmerged = new PlotKineMod<Particle>("PlotMergedLeptons");
    plotmerged->SetLoadBranch(0);
    plotmerged->SetInputName(merger->GetOutputName());
    merger->Add(plotmerged);
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
