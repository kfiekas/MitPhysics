// $Id: doubleLepton.C,v 1.3 2009/04/30 12:13:37 ceballos Exp $

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
#endif

//--------------------------------------------------------------------------------------------------
void doubleLepton(const char *fileset    = "0000",
		  const char *dataset    = "p10-wjets-mg-v26",
		  const char *book       = "cern/filler/013",
		  const char *catalogDir = "/home/mitprod/catalog",
		  Int_t       nEvents    = -1)
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
  const Double_t  ptMinMax  = 10;
  const Double_t  ptMin     = 10;

  MuonIDMod *muId = new MuonIDMod;  
  muId->SetInputName (muInput);
  muId->SetPtMin     (ptMin);
  muId->SetApplyD0Cut(kFALSE);
  muId->SetClassType ("Global");
  muId->SetIDType    ("NoId");
  muId->SetIsoType   ("NoIso");

  ElectronIDMod *elId = new ElectronIDMod;
  elId->SetInputName            (elInput);
  elId->SetPtMin                (ptMin);
  elId->SetApplyConversionFilter(kFALSE);
  elId->SetApplySpikeRemoval    (kFALSE);
  elId->SetApplyD0Cut           (kFALSE);
  elId->SetIDType               ("NoId");
  elId->SetIsoType              ("NoIso");

  MergeLeptonsMod *merger = new MergeLeptonsMod;
  merger->SetMuonsName    (muId->GetOutputName());
  merger->SetElectronsName(elId->GetOutputName());

  GenericSelMod<Particle> *selMod = new GenericSelMod<Particle>;
  selMod->SetMinMaxPt(ptMinMax);
  selMod->SetPtMin(ptMin);
  selMod->SetMinCounts(2);
  selMod->SetColName(merger->GetOutputName());
  merger->Add(selMod);

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  OutputMod *outMod = new OutputMod;
  outMod->Keep("*");
  selMod->Add(outMod);
  TString rootFile = skimName;
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
  if (nEvents>0)
    ana->SetProcessNEvents(nEvents);

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
