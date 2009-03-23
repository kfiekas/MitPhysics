// $Id: jetPlusIsoTrack.C,v 1.1 2009/03/22 09:04:13 loizides Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#include "MitAna/PhysicsMod//interface/PlotKineMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitPhysics/SelMods/interface/JetPlusIsoTrackSelMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void jetPlusIsoTrack(const char *fileset    = "",
                     const char *dataset    = "s8-wm-id9",
                     const char *book       = "mit/filler/006",
                     const char *catalogDir = "/home/mitprod/catalog",
                     Int_t       nEvents    = -1)
{
  TString skimName("jetPlusIsoTrack");
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  //------------------------------------------------------------------------------------------------
  // organize selection
  //------------------------------------------------------------------------------------------------
  const char     *jetInput   = Names::gkSC5JetBrn;
  const char     *gsfTracks  = "GsfTracks";
  const Double_t  jetPtMin   = 30;
  const Double_t  trackPtMin = 10;

  JetIDMod *jetId = new JetIDMod;  
  jetId->SetInputName  (jetInput);
  jetId->SetUseCorrection(kFALSE);
  jetId->SetPtCut      (jetPtMin); 

  JetPlusIsoTrackSelMod *selMod = new JetPlusIsoTrackSelMod;
  selMod->SetTrackPtMin(trackPtMin);
  selMod->SetJetColName(jetId->GetOutputName());
  selMod->SetTrackerTrackColName(Names::gkTrackBrn);
  selMod->SetGsfTrackColName(gsfTracks);

  //------------------------------------------------------------------------------------------------
  // link modules together
  //------------------------------------------------------------------------------------------------
  
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
  ana->AddSuperModule(jetId);
  ana->AddSuperModule(selMod);
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
