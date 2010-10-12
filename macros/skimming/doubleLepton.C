// $Id: doubleLepton.C,v 1.4 2010/07/18 21:16:18 ceballos Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/SelMods/interface/DilepSelMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void doubleLepton(const char *catalogDir   = "/home/ceballos/catalog",
                  const char *book	   = "/cern/filler/014a",
                  const char *dataset	   = "p10-wjets-mg-v26",
                  const char *fileset	   = "0003",
                  int nsel = 0, int NEvents = 1000)
{
  TString skimName("doubleLepton");
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  //------------------------------------------------------------------------------------------------
  // organize selection
  //------------------------------------------------------------------------------------------------
  const char     *muInput  = Names::gkMuonBrn;
  const char     *elInput  = Names::gkElectronBrn;

  MuonIDMod *muId = new MuonIDMod;  
  muId->SetInputName (muInput);
  muId->SetPtMin     (10.0);
  muId->SetApplyD0Cut(kFALSE);
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
  elId->SetApplyD0Cut                (kFALSE);
  elId->SetNExpectedHitsInnerCut     (999);
  elId->SetIDType                    ("NoId");
  elId->SetIsoType                   ("NoIso");

  ElectronCleaningMod *elCl = new ElectronCleaningMod;
  elCl->SetGoodElectronsName(elId->GetOutputName());
  elCl->SetCleanMuonsName   (muId->GetOutputName());
  elId->Add(elCl);

  MergeLeptonsMod *merger = new MergeLeptonsMod;
  merger->SetMuonsName    (muId->GetOutputName());
  merger->SetElectronsName(elCl->GetOutputName());

  DilepSelMod *selMod = new DilepSelMod;
  selMod->SetCleanLeptonsName(merger->GetOutputName());
  selMod->SetMinPt(10.0);
  selMod->SetMinDilMass(2.0);
  selMod->SetFillHist(kFALSE);
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
  ana->AddSuperModule(muId);
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
