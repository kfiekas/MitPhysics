// $Id: jetPlusIsoTrack.C,v 1.1 2009/03/13 13:02:49 sixie Exp $
//root -l $CMSSW_BASE/src/MitPhysics/macros/skimming/jetPlusIsoTrack.C+\(\"0000\",\"s8-qcd_250_500-mg-id9\",\"mit/filler/006\",\"/home/mitprod/catalog\",999999999\)
//root -l $CMSSW_BASE/src/MitPhysics/macros/skimming/jetPlusIsoTrack.C+\(\"0000\",\"s8-qcdem_30_80-id9\",\"mit/filler/006\",\"/home/mitprod/catalog\",999999999\)
//root -l $CMSSW_BASE/src/MitPhysics/macros/skimming/jetPlusIsoTrack.C+\(\"0000\",\"s8-qcdbc_30_80-id9\",\"mit/filler/006\",\"/home/mitprod/catalog\",999999999\)

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
                     int         nEvents    = 999999999)
{
  TString skimName("jetPlusIsoTrack");
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  //------------------------------------------------------------------------------------------------
  // organize selection
  //------------------------------------------------------------------------------------------------
  const char     *jetInput  = Names::gkSC5JetBrn;
  const char     *gsfTracks = "GsfTracks";
  const Double_t  jetPtMin    = 30;
  const Double_t  trackPtMin    = 10;

  JetIDMod *jetId = new JetIDMod;  
  jetId->SetInputName  (jetInput);
  jetId->SetUseCorrection(kFALSE);
  jetId->SetPtCut      (jetPtMin); 

  JetPlusIsoTrackSelMod *selMod = new JetPlusIsoTrackSelMod;
  //selMod->SetJetPtMin(jetPtMin); //don't use this for now. rely on jetID pt cut.
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
  ana->Run(true);
}
