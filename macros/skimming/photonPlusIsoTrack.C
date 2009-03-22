// $Id: photonPlusIsoTrack.C,v 1.1 2009/03/13 13:02:49 sixie Exp $
//root -l $CMSSW_BASE/src/MitPhysics/macros/skimming/photonPlusIsoTrack.C+\(\"0000\",\"s8-pj80-mg-id9\",\"mit/filler/006\",\"/home/mitprod/catalog\",999999999\)

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
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/SelMods/interface/PhotonPlusIsoTrackSelMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void photonPlusIsoTrack(const char *fileset    = "",
                        const char *dataset    = "s8-wm-id9",
                        const char *book       = "mit/filler/006",
                        const char *catalogDir = "/home/mitprod/catalog",
                        int         nEvents    = 999999999)
{
  TString skimName("photonPlusIsoTrack");
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  //------------------------------------------------------------------------------------------------
  // organize selection
  //------------------------------------------------------------------------------------------------
  const char     *muInput   = Names::gkMuonBrn;
  const char     *elInput   = Names::gkElectronBrn;
  const char     *phInput   = Names::gkPhotonBrn;
  const char     *gsfTracks = "GsfTracks";
  const Double_t  ptMin     = 10;

  MuonIDMod *muId = new MuonIDMod;  
  muId->SetInputName  (muInput);
  muId->SetIDType     ("Loose");
  muId->SetIsoType    ("TrackCalo");
  muId->SetTrackIsoCut( 10.0);
  muId->SetCaloIsoCut (100.0);
  muId->SetPtMin      (ptMin);

  ElectronIDMod *elId = new ElectronIDMod;
  elId->SetInputName  (elInput);
  elId->SetIDType     ("NoId");
  elId->SetIsoType    ("TrackCalo");
  elId->SetTrackIsoCut( 10.0);
  elId->SetCaloIsoCut (100.0);
  elId->SetPtMin      (ptMin);

  PhotonIDMod *photonId = new PhotonIDMod;  
  photonId->SetInputName  (phInput);
  photonId->SetIDType     ("Loose");
  photonId->SetIsoType    ("CombinedIso");
  photonId->SetHadOverEmMax(0.1);
  photonId->SetPtMin      (ptMin);
  photonId->SetApplyPixelSeed(kFALSE);

  ElectronCleaningMod *elCl = new ElectronCleaningMod;
  elCl->SetGoodElectronsName(elId->GetOutputName());
  elCl->SetCleanMuonsName   (muId->GetOutputName());
  elId->Add(elCl);

  PhotonCleaningMod *photonCl = new PhotonCleaningMod;
  photonCl->SetCleanElectronsName(elCl->GetOutputName());
  photonCl->SetGoodPhotonsName(photonId->GetOutputName());
  photonId->Add(photonCl);

  PhotonPlusIsoTrackSelMod *selMod = new PhotonPlusIsoTrackSelMod;
  selMod->SetPhotonPtMin(ptMin);
  selMod->SetTrackPtMin(ptMin);
  selMod->SetPhotonColName(photonCl->GetOutputName());
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
  ana->AddSuperModule(muId);
  ana->AddSuperModule(elId);
  ana->AddSuperModule(photonId);
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
