//root -l -q -b $CMSSW_BASE/src/MitHiggs/macros/runMacros/runHwwExampleAnalysis.C+\(\"0000\",\"noskim\",\"s8-h190ww2l-gf-mc3\",\"mit/filler/011\",\"/home/mitprod/catalog\",\"HwwExampleAnalysis\",1000,1\)

// $Id: runPhysicsExample.C,v 1.3 2010/05/10 16:17:02 bendavid Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <exception>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
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
void runPhysicsExample(const char *catalogDir = "/home/mitprod/catalog",
		       const char *book       = "cern/filler/013",
                       const char *dataset    = "p10-ggwwll-v26",
                       const char *fileset    = "0000",
                       const char *skim       = "noskim",
                       const char *outputName = "histo",
                       int   sampleID	      = -1,
                       int   nEvents	      = 1000)
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
  Bool_t useHLTE29      = kFALSE;
  Bool_t applyISRFilter = kFALSE;
  Bool_t applyMllGenCut = kFALSE;
  Bool_t isData         = kFALSE;

//   RunLumiSelectionMod *runLumiSelectionMod = new RunLumiSelectionMod;
//   runLumiSelectionMod->SetAcceptMC(kTRUE);
//   runLumiSelectionMod->AddJSONFile("Cert_132440-133511_StreamExpress_Commissioning10-Express_DQM_JSON.txt");

  if(sampleID > 1000) isData = kTRUE;

  //------------------------------------------------------------------------------------------------
  // generator information
  //------------------------------------------------------------------------------------------------
  GeneratorMod *generatorMod = new GeneratorMod;
  generatorMod->SetPrintDebug(kFALSE);
  generatorMod->SetPtLeptonMin(0.0);
  generatorMod->SetEtaLeptonMax(2.7);
  generatorMod->SetPtPhotonMin(15.0);
  generatorMod->SetEtaPhotonMax(2.7);
  generatorMod->SetPtRadPhotonMin(10.0);
  generatorMod->SetEtaRadPhotonMax(2.7);
  generatorMod->SetIsData(isData);
  generatorMod->SetFillHist(!isData);
  if(applyMllGenCut == kTRUE){
    generatorMod->SetPdgIdCut(23);
    generatorMod->SetMassMaxCut(50.);
  }
  generatorMod->SetApplyISRFilter(applyISRFilter);

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
    hltmod->AddTrigger("HLT_Mu9");
    hltmod->AddTrigger("HLT_Ele10_LW_EleId_L1R");
    hltmod->AddTrigger("HLT_Ele15_SW_EleId_L1R");
  } else {
    hltmod->AddTrigger("HLT_Mu9");
    hltmod->AddTrigger("HLT_Ele10_LW_EleId_L1R");
    hltmod->AddTrigger("HLT_Ele15_SW_EleId_L1R");
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

  PublisherMod<CaloMet> *pubCaloMet = new PublisherMod<CaloMet>;
  pubCaloMet->SetName("CaloMetPub");
  pubCaloMet->SetInputName("CorMuonMet");
  pubCaloMet->SetOutputName("pubCaloMet");

  //------------------------------------------------------------------------------------------------
  // Apply Jet Corrections
  //------------------------------------------------------------------------------------------------
  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  jetCorr->AddCorrectionFromRelease("CondFormats/JetMETObjects/data/Summer09_7TeV_ReReco332_L2Relative_AK5Calo.txt"); 
  jetCorr->AddCorrectionFromRelease("CondFormats/JetMETObjects/data/Summer09_7TeV_ReReco332_L3Absolute_AK5Calo.txt");  
  jetCorr->SetInputName(pubJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");

  //------------------------------------------------------------------------------------------------
  // Apply Met Corrections
  //------------------------------------------------------------------------------------------------
  CaloMetCorrectionMod *metCaloCorr = new CaloMetCorrectionMod;
  metCaloCorr->SetInputName(pubCaloMet->GetOutputName());
  metCaloCorr->SetCorrectedJetsName(jetCorr->GetOutputName());
  metCaloCorr->SetOutputName("pubCaloCorrectedMet");

  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  MuonIDMod           *muonID        = new MuonIDMod;  
  ElectronIDMod       *electronID    = new ElectronIDMod;
  PhotonIDMod         *photonID      = new PhotonIDMod;
  TauIDMod            *tauID         = new TauIDMod;
  JetIDMod            *jetID         = new JetIDMod;
  jetID->SetInputName(jetCorr->GetOutputName());
  jetID->SetPtCut(20.0);
  jetID->SetEtaMaxCut(5.0);
  jetID->SetJetEEMFractionMinCut(0.01);
  jetID->SetOutputName("GoodJets");

  ElectronCleaningMod *electronCleaning = new ElectronCleaningMod;
  PhotonCleaningMod   *photonCleaning   = new PhotonCleaningMod;
  TauCleaningMod      *tauCleaning      = new TauCleaningMod;
  JetCleaningMod      *jetCleaning      = new JetCleaningMod;
  jetCleaning->SetGoodJetsName("GoodJets");
  jetCleaning->SetCleanJetsName("CleanJets");

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
  generatorMod->Add(hltmod);
  hltmod->Add(muonID);
  muonID->Add(electronID);
  electronID->Add(photonID);
  photonID->Add(tauID);
  tauID->Add(pubJet);
  pubJet->Add(pubMet); 
  pubMet->Add(pubCaloMet); 
  pubCaloMet->Add(jetCorr);
  jetCorr->Add(metCaloCorr);
  metCaloCorr->Add(jetID);
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
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kFALSE);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);
  ana->SetSuperModule(generatorMod);
  ana->SetPrintScale(100);
  if(useHLTE29 == true){
    ana->SetHLTTreeName("HLT_E29");
    ana->SetHLTObjsName("HLTObjects_E29");
  }

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
  //ana->AddFile("root://castorcms//castor/cern.ch/user/p/paus/filler/011/s09-ttbar-7-mc3/*.root");
  //ana->AddFile("/build/bendavid/XX-MITDATASET-XX_000-valskim-Run132605.root");

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

  ana->SetCacheSize(64*1024*1024);
  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());  

  return;
}
