//root -l -q -b $CMSSW_BASE/src/MitHiggs/macros/runMacros/runHwwExampleAnalysis.C+\(\"0000\",\"noskim\",\"s8-h190ww2l-gf-mc3\",\"mit/filler/011\",\"/home/mitprod/catalog\",\"HwwExampleAnalysis\",1000,1\)

// $Id: runPhysicsExample.C,v 1.1 2010/03/12 13:52:01 bendavid Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <exception>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/PDFProducerMod.h"
#include "MitPhysics/Mods/interface/HKFactorProducer.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/CaloMetCorrectionMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/TauIDMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/TauCleaningMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/MetCol.h" 
#include "MitAna/DataTree/interface/CaloMetCol.h"
#include "MitAna/PhysicsMod/interface/FullExampleMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void executePhysicsExample(const char *fileset  = "",
                             const char *skim         = "",
                             const char *dataset    = "",
                             const char *book       = "",
                             const char *catalogDir = "",
                             const char *outputName = "",
                             int   sampleID         = -1,
                             int   nEvents          = -1)
{
  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  //------------------------------------------------------------------------------------------------
  // set up information
  //------------------------------------------------------------------------------------------------
  Bool_t isFastSim = kFALSE;
  if(sampleID >= 100) isFastSim = kTRUE;

  bool useHLTE29 = true;

  //------------------------------------------------------------------------------------------------
  // generator information
  //------------------------------------------------------------------------------------------------
  GeneratorMod *generatorMod = new GeneratorMod;

  //------------------------------------------------------------------------------------------------
  // HLT information
  //------------------------------------------------------------------------------------------------
  HLTMod *hltmod = new HLTMod;
  if(useHLTE29 == false) {
    hltmod->AddTrigger("HLT_Ele10_SW_L1R");
    hltmod->AddTrigger("HLT_Ele15_SW_LooseTrackIso_L1R");
    hltmod->AddTrigger("HLT_Ele15_SW_EleId_L1R");
    hltmod->AddTrigger("HLT_Ele15_LW_L1R");
    hltmod->AddTrigger("HLT_Ele15_SC10_LW_L1R");
    hltmod->AddTrigger("HLT_Ele20_SW_L1R");
    hltmod->AddTrigger("HLT_IsoMu9");
  } else {
    hltmod->AddTrigger("HLT_Ele10_LW_EleId_L1R");
    hltmod->AddTrigger("HLT_Mu9");
    hltmod->SetBitsName("HLTBits_E29");
  }
  hltmod->SetTrigObjsName("myhltobjs");


  //------------------------------------------------------------------------------------------------
  // publisher Mod
  //------------------------------------------------------------------------------------------------
  PublisherMod<CaloJet,Jet> *pubJet = new PublisherMod<CaloJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5Jets");
  pubJet->SetOutputName("PubAKt5Jets");

  PublisherMod<Met,Met> *pubMet = new PublisherMod<Met,Met>("MetPub");
  pubMet->SetInputName("TCMet");
  pubMet->SetOutputName("PubTCMet");

  //------------------------------------------------------------------------------------------------
  // Apply Jet/Met Corrections
  //------------------------------------------------------------------------------------------------
  JetCorrectionMod *jetCorr = new JetCorrectionMod;  
  jetCorr->SetCorrectionTag("Summer09_L2Relative_AK5Calo:Summer09_L3Absolute_AK5Calo");
  jetCorr->SetInputName(pubJet->GetOutputName());

  CaloMetCorrectionMod *metCorr = new CaloMetCorrectionMod;
  metCorr->SetInputName(pubMet->GetOutputName());
  metCorr->SetCorrectedJetsName(jetCorr->GetOutputName());


  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  MuonIDMod           *muonID           = new MuonIDMod;  
  ElectronIDMod       *electronID       = new ElectronIDMod;
  electronID->SetIDType(TString("CustomTight"));
  PhotonIDMod         *photonID       = new PhotonIDMod;
  photonID->SetIDType(TString("Custom"));
  TauIDMod *tauID = new TauIDMod;
  JetIDMod            *jetID            = new JetIDMod;
  jetID->SetInputName(jetCorr->GetOutputName());
  jetID->SetUseCorrection(kTRUE); 
  jetID->SetPtCut(30.0);
  jetID->SetEtaMaxCut(5.0);
  jetID->SetOutputName(ModNames::gkGoodJetsName);
  ElectronCleaningMod *electronCleaning = new ElectronCleaningMod;
  PhotonCleaningMod   *photonCleaning   = new PhotonCleaningMod;
  TauCleaningMod      *tauCleaning      = new TauCleaningMod;
  JetCleaningMod      *jetCleaning      = new JetCleaningMod;
  jetCleaning->SetGoodJetsName(ModNames::gkGoodJetsName);
  jetCleaning->SetCleanJetsName(ModNames::gkCleanJetsName);

  //------------------------------------------------------------------------------------------------
  // merge modules
  //------------------------------------------------------------------------------------------------
  MergeLeptonsMod *mergeLeptonsMod = new MergeLeptonsMod;
  mergeLeptonsMod->SetMuonsName(muonID->GetOutputName());
  mergeLeptonsMod->SetElectronsName(electronCleaning->GetOutputName());

  //------------------------------------------------------------------------------------------------
  // analyses modules
  //------------------------------------------------------------------------------------------------
  FullExampleMod *analysisMod = new FullExampleMod;
  analysisMod->SetMuonName(muonID->GetOutputName());
  analysisMod->SetMuonsFromBranch(kFALSE);
  analysisMod->SetElectronName(electronID->GetOutputName());
  analysisMod->SetElectronsFromBranch(kFALSE);


  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  generatorMod->Add(muonID);
  muonID->Add(electronID);
  electronID->Add(photonID);
  photonID->Add(tauID);
  tauID->Add(pubJet);
  pubJet->Add(pubMet);
  pubMet->Add(jetCorr);
  jetCorr->Add(metCorr);
  metCorr->Add(jetID);
  jetID->Add(electronCleaning);
  electronCleaning->Add(photonCleaning);
  photonCleaning->Add(tauCleaning);
  tauCleaning->Add(jetCleaning);
  jetCleaning->Add(mergeLeptonsMod);
  mergeLeptonsMod->Add(analysisMod);

  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseHLT(kFALSE);
  ana->SetKeepHierarchy(kTRUE);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);
  ana->SetSuperModule(generatorMod);

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  printf("\nRely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n\n",book,dataset,skim,fileset);
  Catalog *c = new Catalog(catalogDir);
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(book,dataset,fileset);
  else 
    d = c->FindDataset(book,skimdataset.Data(),fileset);
  ana->AddDataset(d);
  //ana->AddFile("rfio:/castor/cern.ch/user/p/paus/filler/011/s09-ttbar-mc3/s09-ttbar-mc3_000_10.root");

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  printf("\nRoot output: %s\n\n",rootFile.Data());  
  ana->SetOutputName(rootFile.Data());

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());  

  return;
}

//--------------------------------------------------------------------------------------------------
void runPhysicsExample(const char *fileset      = "",
                         const char *skim         = "noskim",
                         const char *dataset      = "s09-ttbar-mc3",
                         const char *book         = "cern/filler/011",
                         const char *catalogDir   = "/home/mitprod/catalog",
                         const char *outputName   = "PhysicsExample",
                         int         nEvents      = 10000, 
                         int         runTypeIndex = -1)
{
  TString outfileName = TString(outputName);
  executePhysicsExample(fileset,skim,dataset,book,catalogDir,outfileName,runTypeIndex,nEvents);
}
