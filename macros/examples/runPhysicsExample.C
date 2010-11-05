//root -l -q -b $CMSSW_BASE/src/MitHiggs/macros/runMacros/runHwwExampleAnalysis.C+\(\"0000\",\"noskim\",\"s8-h190ww2l-gf-mc3\",\"mit/filler/011\",\"/home/mitprod/catalog\",\"HwwExampleAnalysis\",1000,1\)

// $Id: runPhysicsExample.C,v 1.9 2010/10/29 16:19:23 ceballos Exp $

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
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
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
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/MetCol.h" 
#include "MitAna/DataTree/interface/CaloMetCol.h"
#include "MitPhysics/SelMods/interface/HwwExampleAnalysisMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void runPhysicsExample(const char *catalogDir = "/home/mitprod/catalog",
		       const char *book	      = "cern/filefi/015",
                       const char *dataset    = "f10-ww2l-z2-v12",
                       const char *fileset    = "0000",
                       const char *skim       = "noskim",
                       const char *outputName = "histo",
                       int   sampleID	      = -1,
                       int   nEvents	      = 1000000)
{
  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 3;

  //------------------------------------------------------------------------------------------------
  // set up information
  //------------------------------------------------------------------------------------------------
  Bool_t applyISRFilter = kFALSE;
  Bool_t applyMllGenCut = kFALSE;
  Bool_t isData         = kFALSE;
  Bool_t isElData       = kFALSE;
  int processId         = -999999999; // use 999 for MCatNLO MC sample, 102 for H->WW
  TString fInputFilenameKF = "/home/ceballos/releases/CMSSW_3_8_5/src/MitPhysics/data/HWW_KFactors_160_10TeV.dat";

  if(sampleID > 1000) isData   = kTRUE;
  if(sampleID > 2000) isElData = kTRUE;

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
    generatorMod->SetMassMinCut( 0.);
    generatorMod->SetMassMaxCut(50.);
  }
  generatorMod->SetApplyISRFilter(applyISRFilter);

  HKFactorProducer *hKFactorProducer = new HKFactorProducer;
  hKFactorProducer->SetProcessID(processId);
  hKFactorProducer->SetInputFilename(fInputFilenameKF);
  hKFactorProducer->SetIsData(isData);
  hKFactorProducer->SetFillHist(!isData);

  //------------------------------------------------------------------------------------------------
  // Run RunLumiSelectionMod
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSelectionMod = new RunLumiSelectionMod;
  runLumiSelectionMod->SetAcceptMC(!isData);    
  runLumiSelectionMod->AddJSONFile("/home/ceballos/releases/CMSSW_3_8_5/src/json/Cert_TopNov5_Merged_135821-149442_allPVT.txt");

  //------------------------------------------------------------------------------------------------
  // PV filter selection
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof(4);
  goodPVFilterMod->SetMaxAbsZ(24.0);
  goodPVFilterMod->SetMaxRho(2.0);

  //------------------------------------------------------------------------------------------------
  // HLT information
  //------------------------------------------------------------------------------------------------
  HLTMod *hltmod = new HLTMod;

  if     (isData == kFALSE){
    hltmod->AddTrigger("HLT_Mu9");
    hltmod->AddTrigger("!HLT_Mu9");
  }
  else if(isElData == kFALSE){
    hltmod->AddTrigger("HLT_Mu9",136033,147116);
    hltmod->AddTrigger("HLT_Mu9&HLT_Ele10_LW_L1R",136033,139980);
    hltmod->AddTrigger("HLT_Mu9&HLT_Ele15_SW_L1R",140058,141882);
    hltmod->AddTrigger("HLT_Mu9&HLT_Ele15_SW_CaloEleId_L1R",141956,144114);
    hltmod->AddTrigger("HLT_Mu9&HLT_Ele17_SW_CaloEleId_L1R",146428,147116);

    hltmod->AddTrigger("HLT_Mu15_v1",147196,999999);
    hltmod->AddTrigger("HLT_Mu15_v1&HLT_Ele17_SW_TightEleId_L1R",147196,148058);
    hltmod->AddTrigger("HLT_Mu15_v1&HLT_Ele17_SW_TighterEleIdIsol_L1R_v2",148819,149064);
    hltmod->AddTrigger("HLT_Mu15_v1&HLT_Ele17_SW_TighterEleIdIsol_L1R_v3",149181,999999);
  }
  else {
    hltmod->AddTrigger("!HLT_Mu9&HLT_Ele10_LW_L1R",136033,139980);
    hltmod->AddTrigger("!HLT_Mu9&HLT_Ele15_SW_L1R",140058,141882);
    hltmod->AddTrigger("!HLT_Mu9&HLT_Ele15_SW_CaloEleId_L1R",141956,144114);
    hltmod->AddTrigger("!HLT_Mu9&HLT_Ele17_SW_CaloEleId_L1R",146428,147116);

    hltmod->AddTrigger("!HLT_Mu15_v1&HLT_Ele17_SW_TightEleId_L1R",147196,148058);
    hltmod->AddTrigger("!HLT_Mu15_v1&HLT_Ele17_SW_TighterEleIdIsol_L1R_v2",148819,149064);
    hltmod->AddTrigger("!HLT_Mu15_v1&HLT_Ele17_SW_TighterEleIdIsol_L1R_v3",149181,999999);
  }
  hltmod->SetTrigObjsName("myhltobjs");

  //------------------------------------------------------------------------------------------------
  // publisher Mod
  //------------------------------------------------------------------------------------------------
  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");

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
  jetCorr->AddCorrectionFromRelease("CondFormats/JetMETObjects/data/Spring10_L2Relative_AK5PF.txt"); 
  jetCorr->AddCorrectionFromRelease("CondFormats/JetMETObjects/data/Spring10_L3Absolute_AK5PF.txt");  
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
  muonID->SetClassType("Global");
  muonID->SetIDType("Minimal");
  muonID->SetIsoType("TrackCaloSliding");
  muonID->SetApplyD0Cut(kTRUE);

  ElectronIDMod       *electronID    = new ElectronIDMod;
  electronID->SetIDType("VBTFWorkingPoint80Id");
  electronID->SetIsoType("TrackJuraSliding");
  electronID->SetApplyConversionFilterType1(kFALSE);
  electronID->SetApplyConversionFilterType2(kTRUE);
  electronID->SetChargeFilter(kFALSE);
  electronID->SetApplyD0Cut(kTRUE);
  electronID->SetNExpectedHitsInnerCut(0);

  PhotonIDMod         *photonID      = new PhotonIDMod;
  TauIDMod            *tauID         = new TauIDMod;
  JetIDMod            *jetID         = new JetIDMod;
  jetID->SetInputName(jetCorr->GetOutputName());
  jetID->SetPtCut(25.0);
  jetID->SetEtaMaxCut(5.0);
  jetID->SetJetEEMFractionMinCut(0.0);
  jetID->SetOutputName("GoodJets");

  ElectronCleaningMod *electronCleaning = new ElectronCleaningMod;
  PhotonCleaningMod   *photonCleaning   = new PhotonCleaningMod;
  TauCleaningMod      *tauCleaning      = new TauCleaningMod;
  JetCleaningMod      *jetCleaning      = new JetCleaningMod;
  jetCleaning->SetGoodJetsName("GoodJets");
  jetCleaning->SetCleanJetsName("CleanJets");

  JetIDMod            *jetIDNoPtCut     = new JetIDMod;
  jetIDNoPtCut->SetInputName(jetCorr->GetOutputName());
  jetIDNoPtCut->SetPtCut(0.0);
  jetIDNoPtCut->SetEtaMaxCut(5.0);
  jetIDNoPtCut->SetJetEEMFractionMinCut(0.0);
  jetIDNoPtCut->SetOutputName("GoodJetsNoPtCut");

  JetCleaningMod      *jetCleaningNoPtCut = new JetCleaningMod;
  jetCleaningNoPtCut->SetGoodJetsName("GoodJetsNoPtCut");
  jetCleaningNoPtCut->SetCleanJetsName("CleanJetsNoPtCut");

  //------------------------------------------------------------------------------------------------
  // merge modules
  //------------------------------------------------------------------------------------------------
  MergeLeptonsMod *mergeLeptonsMod = new MergeLeptonsMod;
  mergeLeptonsMod->SetMuonsName(muonID->GetOutputName());
  mergeLeptonsMod->SetElectronsName(electronCleaning->GetOutputName());

  //------------------------------------------------------------------------------------------------
  // analyses modules
  //------------------------------------------------------------------------------------------------
  HwwExampleAnalysisMod *analysisMod = new HwwExampleAnalysisMod;
  analysisMod->SetMetName(pubMet->GetOutputName());
  analysisMod->SetCleanJetsName(jetCleaning->GetOutputName());
  analysisMod->SetCleanJetsNoPtCutName(jetCleaningNoPtCut->GetOutputName());

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  generatorMod->Add(hKFactorProducer);
  hKFactorProducer->Add(runLumiSelectionMod);
  runLumiSelectionMod->Add(goodPVFilterMod);
  goodPVFilterMod->Add(hltmod);
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
  jetCleaning->Add(jetIDNoPtCut);
  jetIDNoPtCut->Add(jetCleaningNoPtCut);
  jetCleaningNoPtCut->Add(mergeLeptonsMod);
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
