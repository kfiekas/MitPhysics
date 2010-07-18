// $Id: doubleLepton.C,v 1.3 2010/05/12 19:05:51 ceballos Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#include "MitAna/PhysicsMod/interface/PlotKineMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/SelMods/interface/GenericSelMod.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void doubleLepton(const char *catalogDir   = "/home/sixie/catalog",
                  const char *book	   = "cern/filler/014",
                  const char *dataset	   = "run2010a-eg-rrjun14",
                  const char *fileset	   = "0003",
                  int nsel = 0, int NEvents = 10000)
{
  TString skimName("doubleLepton");
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  //------------------------------------------------------------------------------------------------
  // organize selection
  //------------------------------------------------------------------------------------------------
  HLTMod *hltmod = new HLTMod;
  hltmod->SetPrintTable(kFALSE);
  hltmod->AddTrigger("HLT_Mu3");
  hltmod->AddTrigger("HLT_Mu5");
  hltmod->AddTrigger("HLT_Mu9");
  hltmod->AddTrigger("HLT_Mu0_Track0_Jpsi");
  hltmod->AddTrigger("HLT_Mu3_Track0_Jpsi");
  hltmod->AddTrigger("HLT_Mu5_Track0_Jpsi");
  hltmod->AddTrigger("HLT_Ele10_LW_L1R");
  hltmod->AddTrigger("HLT_Ele10_LW_EleId_L1R");
  hltmod->AddTrigger("HLT_Ele15_LW_L1R");
  hltmod->AddTrigger("HLT_Ele20_LW_L1R");
  hltmod->AddTrigger("HLT_Photon10_L1R");
  hltmod->AddTrigger("HLT_Photon15_L1R");
  hltmod->AddTrigger("HLT_Photon20_L1R");
  hltmod->SetTrigObjsName("myhltobjs");

  const char     *muInput  = Names::gkMuonBrn;
  const char     *elInput  = Names::gkElectronBrn;

  MuonIDMod *muId = new MuonIDMod;  
  muId->SetInputName (muInput);
  muId->SetPtMin     (3.0);
  muId->SetApplyD0Cut(kTRUE);
  muId->SetD0Cut     (0.2);
  muId->SetClassType ("All");
  muId->SetIDType    ("NoId");
  muId->SetIsoType   ("NoIso");

  ElectronIDMod *elId = new ElectronIDMod;
  elId->SetInputName                 (elInput);
  elId->SetPtMin                     (10.0);
  elId->SetApplyConversionFilterType1(kFALSE);
  elId->SetApplyConversionFilterType2(kFALSE);
  elId->SetChargeFilter              (kFALSE);
  elId->SetApplySpikeRemoval         (kFALSE);
  elId->SetApplyD0Cut                (kTRUE);
  elId->SetD0Cut                     (0.2);
  elId->SetNExpectedHitsInnerCut     (999);
  elId->SetIDType                    ("NoId");
  elId->SetIsoType                   ("NoIso");

  MergeLeptonsMod *merger = new MergeLeptonsMod;
  merger->SetMuonsName    (muId->GetOutputName());
  merger->SetElectronsName(elId->GetOutputName());

  GenericSelMod<Particle> *selMod = new GenericSelMod<Particle>;
  selMod->SetPtMin(3.0);
  selMod->SetMinCounts(1);
  selMod->SetColName(merger->GetOutputName());
  merger->Add(selMod);

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  OutputMod *outMod = new OutputMod;
  outMod->Keep("*");
  selMod->Add(outMod);
  TString rootFile = "";
  rootFile += TString("skims/") + dataset + TString("/");
  rootFile += skimName;
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  printf("\nRoot output: %s\n\n",rootFile.Data());  
  outMod->SetFileName(rootFile);
  outMod->SetPathName(".");

  //------------------------------------------------------------------------------------------------
  // set up analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->AddSuperModule(hltmod);
  hltmod->Add(muId);
  muId->Add(elId);
  elId->Add(merger);
  if (NEvents>0)
    ana->SetProcessNEvents(NEvents);

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  printf("\nRely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Fileset: %s <-\n\n",book,dataset,fileset);
  Catalog *c = new Catalog(catalogDir);
  Dataset *d = c->FindDataset(book,dataset,fileset);
  ana->AddDataset(d);

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(kFALSE);
}
