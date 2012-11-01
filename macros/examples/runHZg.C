//root -l -q -b /home/sixie/CMSSW_5_3_2_patch4/src/MitPhysics/macros/examples/runHZg.C+\(\"0000\",\"noskim\",\"r12a-del-pr-v1\",\"cern/filefi/028\",\"/afs/cern.ch/user/s/sixie/catalog\",\"AllNtupler\",-1,-1\)

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/PDFProducerMod.h"
#include "MitPhysics/Mods/interface/HKFactorProducer.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/CaloMetCorrectionMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PFTauIDMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/PFTauCleaningMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/PartonFlavorHistoryMod.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/MetCol.h" 
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/CaloMetCol.h"
#include "MitPhysics/SelMods/interface/GenericSelMod.h"
#include "MitPhysics/Mods/interface/LeptonPairPhotonTreeWriter.h"
#include "MitPhysics/Mods/interface/PhotonMvaMod.h"
#include "MitPhysics/Mods/interface/MVASystematicsMod.h"

#endif

using namespace mithep;

//==================================================================================================
/*
 * Triggers of interest
 *

  kHLT_Mu9                                   = 0x0000001,
  kHLT_Mu11                                  = 0x0000002,
  kHLT_Mu13_v1                               = 0x0000004,
  kHLT_Mu15_v1                               = 0x0000008,
  kHLT_DoubleMu5_v1                          = 0x0000010,
  kHLT_Jet15U                                = 0x0000020,
  kHLT_Jet30U                                = 0x0000040,
  kHLT_Jet50U                                = 0x0000080,
  kHLT_Photon10_L1R                          = 0x0000100,
  kHLT_Photon10_Cleaned_L1R                  = 0x0000200,
  kHLT_Photon15_L1R                          = 0x0000400,
  kHLT_Photon15_Cleaned_L1R                  = 0x0000800,
  kHLT_Photon20_L1R                          = 0x0001000,
  kHLT_Photon20_Cleaned_L1R                  = 0x0002000,
  kHLT_Photon30_Cleaned_L1R                  = 0x0004000,
  kHLT_Ele15_SW_L1R                          = 0x0008000,
  kHLT_Ele15_LW_L1R                          = 0x0010000,
  kHLT_Ele15_SW_CaloEleId_L1R                = 0x0020000,
  kHLT_Ele17_SW_L1R                          = 0x0040000,
  kHLT_Ele17_SW_CaloEleId_L1R                = 0x0080000,       
  kHLT_Ele17_SW_TightEleId_L1R               = 0x0100000,
  kHLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1 = 0x0200000,
  kHLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1  = 0x0400000,
  kHLT_Ele17_SW_TighterEleIdIsol_L1R_v2      = 0x0800000,
  kHLT_DoubleEle10_SW_L1R                    = 0x1000000,
  kHLT_DoubleEle15_SW_L1R_v1                 = 0x2000000

  */
    
//==================================================================================================
/*
 * Run on a BAMBU fileset
 *
 * Example usage:
 *   root -l -q -b runZeeNtupler.C+\(\"0000\",\"p10-zee-v26\",\"cern/filler/014a\",\"/home/ceballos/catalog\",1,0,1,-1,0,1\)
 *
 * Output file name has standard format: <dataset>_<fileset>_ntuple.root
 *
 */
void runHZg( 
  const char *fileset  = "",
  const char *skim         = "noskim",
  const char *dataset    = "s8-ttbar-id9",
  const char *book       = "mit/filler/006",
  const char *catalogDir = "/home/mitprod/catalog",
  const char *outputName = "HwwHiggsNtupleMaker",
  int   nEvents          = -1,
  int   sampleID         = -1
  )
{
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 3;

  //******************************************************************
  //Set up Options
  //******************************************************************
  bool usePDFProducer     = false;
  string pdfSetName = "";

  bool isData             = false;
  bool isDataDMuon        = false;
  bool isDataDElectron    = false;

  int fDecay = sampleID; 
  if (fDecay < 0) {
    isData = true;
    if (fDecay == -1) isDataDMuon = true;
    if (fDecay == -2) isDataDElectron = true;
  }

  cout << "Summarize Run Options: " << "fDecay == " << fDecay << " "
       << endl;



  //******************************************************************
  //Modules
  //******************************************************************

  // Generator info
  GeneratorMod *GeneratorMod1 = new GeneratorMod;
  GeneratorMod1->SetPrintDebug(kFALSE);
  GeneratorMod1->SetPtLeptonMin(0.0);
  GeneratorMod1->SetEtaLeptonMax(2.7);
  GeneratorMod1->SetPtPhotonMin(0.0);
  GeneratorMod1->SetEtaPhotonMax(2.7);
  GeneratorMod1->SetPtRadPhotonMin(5.0);
  GeneratorMod1->SetEtaRadPhotonMax(2.7);
  GeneratorMod1->SetIsData(isData);
  GeneratorMod1->SetFillHist(!isData);


  // HLT info

  HLTMod *hltmod_dmu = new HLTMod;
  hltmod_dmu->SetAbortIfNotAccepted(kTRUE);
  hltmod_dmu->SetTrigObjsName("myhltobjs_dmu");
  if (!isData) {
    hltmod_dmu->AddTrigger("HLT_Mu17_Mu8_v*");
  }
  if(isData == true && isDataDMuon == true) {
    hltmod_dmu->AddTrigger("HLT_DoubleMu7_v1",150000,161176);
    hltmod_dmu->AddTrigger("HLT_DoubleMu7_v1",161179,163261);
    hltmod_dmu->AddTrigger("HLT_DoubleMu7_v2",163262,164237);
    hltmod_dmu->AddTrigger("HLT_Mu13_Mu8_v2" ,165085,165888);
    hltmod_dmu->AddTrigger("HLT_Mu13_Mu8_v2" ,165900,167043);
    hltmod_dmu->AddTrigger("HLT_Mu13_Mu8_v4" ,167044,170053);
    hltmod_dmu->AddTrigger("HLT_Mu13_Mu8_v6" ,170054,173198);
    hltmod_dmu->AddTrigger("HLT_Mu13_Mu8_v7" ,173199,178380);
    hltmod_dmu->AddTrigger("HLT_Mu17_Mu8_v10" ,178381,179889);
    hltmod_dmu->AddTrigger("HLT_Mu17_TkMu8_v3" ,178381,179889);
    hltmod_dmu->AddTrigger("HLT_Mu17_Mu8_v11"  ,179890,180000);
    hltmod_dmu->AddTrigger("HLT_Mu17_TkMu8_v4" ,179890,180000);

    hltmod_dmu->AddTrigger("HLT_Mu17_Mu8_v16"  ,190456,190738);
    hltmod_dmu->AddTrigger("HLT_Mu17_TkMu8_v9" ,190456,190738);
    hltmod_dmu->AddTrigger("HLT_Mu17_Mu8_v16"  ,190739,191419);
    hltmod_dmu->AddTrigger("HLT_Mu17_TkMu8_v9" ,190739,191419);
    hltmod_dmu->AddTrigger("HLT_Mu17_Mu8_v16"  ,191420,193686);
    hltmod_dmu->AddTrigger("HLT_Mu17_TkMu8_v9" ,191420,193686);
    hltmod_dmu->AddTrigger("HLT_Mu17_Mu8_v17"  ,193687,196045);
    hltmod_dmu->AddTrigger("HLT_Mu17_TkMu8_v10",193687,196045);
    hltmod_dmu->AddTrigger("HLT_Mu17_Mu8_v18"  ,196046,197669);
    hltmod_dmu->AddTrigger("HLT_Mu17_TkMu8_v11",196046,197669);
    hltmod_dmu->AddTrigger("HLT_Mu17_Mu8_v19"  ,197770,199631);
    hltmod_dmu->AddTrigger("HLT_Mu17_TkMu8_v12",197770,199631);
    hltmod_dmu->AddTrigger("HLT_Mu17_Mu8_v21"  ,199632,999999);
    hltmod_dmu->AddTrigger("HLT_Mu17_TkMu8_v13",199632,999999);
  }

  HLTMod *hltmod_del = new HLTMod;
  hltmod_del->SetAbortIfNotAccepted(kTRUE);
  hltmod_del->SetTrigObjsName("myhltobjs_del");
  if (!isData) {
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
  }
  if(isData == true && isDataDElectron == true) {
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",150000,161176);
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",161179,163261);
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",163262,164237);
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4",165085,165888);
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5",165900,166967);
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6",166968,170053);
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",170054,170759);
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",170760,173198);
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",173199,178380);
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9",178381,179889);
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10",179890,180000);

    hltmod_del->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15" ,190456,190738);
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16" ,190739,191419);
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17" ,191420,193686);
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17",193687,196045);
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17",196046,197669);
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18",197770,199631);
    hltmod_del->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18",199632,999999);
  }



  //------------------------------------------------------------------------------------------------
  // Run RunLumiSelectionMod
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSelection = new RunLumiSelectionMod;      
  runLumiSelection->SetAcceptMC(!isData);
  runLumiSelection->SetAcceptAll(kTRUE);
  runLumiSelection->SetAbortIfNotAccepted(kFALSE);

  //------------------------------------------------------------------------------------------------
  // PV filter selection
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof(4);
  goodPVFilterMod->SetMaxAbsZ(24.0);
  goodPVFilterMod->SetMaxRho(2.0);
  goodPVFilterMod->SetAbortIfNotAccepted(kTRUE);
  goodPVFilterMod->SetIsMC(!isData);


  //------------------------------------------------------------------------------------------------
  // Photon ID
  //------------------------------------------------------------------------------------------------
  PhotonIDMod *myPhId = new PhotonIDMod;
  myPhId -> SetPtMin    (10.);
  myPhId -> SetAbsEtaMax   (2.5);
  myPhId -> SetGoodElectronsFromBranch( true );
  myPhId -> SetPhotonsFromBranch( true);
  myPhId -> SetOutputName("outputPhotons");
  myPhId -> SetIDType("TrivialSelection");
  myPhId -> DoMCSmear( true );
  myPhId -> SetMCSmearFactors(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02);
  myPhId -> SetIsData(isData);

  PhotonIDMod *myPhIdNoSmear = new PhotonIDMod;
  myPhIdNoSmear -> SetPtMin    (10.);
  myPhIdNoSmear -> SetAbsEtaMax   (2.5);
  myPhIdNoSmear -> SetGoodElectronsFromBranch( true );
  myPhIdNoSmear -> SetPhotonsFromBranch(true);
  myPhIdNoSmear -> SetOutputName("outputPhotonsNoSmear");
  myPhIdNoSmear -> SetIDType("TrivialSelection");
  myPhIdNoSmear -> DoMCSmear( false );
  myPhIdNoSmear -> SetMCSmearFactors(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02);
  myPhIdNoSmear -> SetIsData(isData);

  //------------------------------------------------------------------------------------------------
  // Photon Regression 
  //------------------------------------------------------------------------------------------------
  PhotonMvaMod *photreg = new PhotonMvaMod;
  photreg->SetRegressionVersion(3);
  photreg->SetRegressionWeights("/afs/cern.ch/user/s/sixie/CMSSW_analysis/src/MitPhysics/data/gbrv3ph_52x.root");
  photreg->SetOutputName("GoodPhotonsRegr");
  photreg->SetApplyShowerRescaling(kTRUE);
  photreg->SetIsData(isData);
//   photreg->SetPtMin(10.0);

  //------------------------------------------------------------------------------------------------
  // MVA systematics
  //------------------------------------------------------------------------------------------------
  MVASystematicsMod *sysMod = new MVASystematicsMod;
  sysMod->SetMCR9Scale(1.0035, 1.0035);  
  sysMod->SetIsData(isData);


  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile =  TString("./") + TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  printf("\nRoot output: %s\n\n",rootFile.Data());  

  //------------------------------------------------------------------------------------------------
  //
  // HZG Ntupler
  //
  //------------------------------------------------------------------------------------------------

  LeptonPairPhotonTreeWriter *eegtree = new LeptonPairPhotonTreeWriter;
  eegtree ->SetInputPhotonsName(myPhIdNoSmear->GetOutputName());
  eegtree ->SetPhotonsFromBranch( false );
  eegtree ->SetGoodElectronsFromBranch( true );
  eegtree ->SetIsData(isData);
  eegtree ->SetGoodMuonsFromBranch( true ); 
  eegtree ->SetYear(2012); 
  eegtree ->SetVerbose(false); 
  eegtree ->SetDoMuonChannel(false); 
  eegtree->SetTupleName("HZeegEvents");

  LeptonPairPhotonTreeWriter *mmgtree = new LeptonPairPhotonTreeWriter;
  mmgtree ->SetInputPhotonsName(myPhIdNoSmear->GetOutputName());
  mmgtree ->SetPhotonsFromBranch( false );
  mmgtree ->SetGoodElectronsFromBranch( true );
  mmgtree ->SetIsData(isData);
  mmgtree ->SetGoodMuonsFromBranch( true ); 
  mmgtree ->SetYear(2012); 
  mmgtree ->SetVerbose(false); 
  mmgtree ->SetDoElectronChannel(false); 
  mmgtree->SetTupleName("HZmmgEvents");




  TString rootFileHZZ4lNtuple = TString("./");
  rootFileHZZ4lNtuple += TString(outputName);
  rootFileHZZ4lNtuple += TString("_HZgNtuple_") + TString(dataset) + TString("_") + TString(skim); 
  rootFileHZZ4lNtuple += TString("_") + TString(fileset);
  rootFileHZZ4lNtuple += TString(".root");


  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // Chain modules together
  GeneratorMod1->Add(runLumiSelection);

  if (TString(dataset).Contains("-h")) {    
    runLumiSelection->Add(sysMod);
  }

  runLumiSelection->Add(goodPVFilterMod);
  goodPVFilterMod->Add(myPhIdNoSmear);

  if(!isData || isDataDMuon) {
    myPhIdNoSmear->Add(hltmod_dmu);
    hltmod_dmu->Add(mmgtree);   
  } 
  if(!isData || isDataDElectron) {
    myPhIdNoSmear->Add(hltmod_del);
    hltmod_del->Add(eegtree);   
  } 


  TFile::SetReadaheadSize(128*1024*1024);

  //------------------------------------------------------------------------------------------------
  //
  // setup analysis object
  //
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kFALSE);
  if(nEvents >= 0) 
    ana->SetProcessNEvents(nEvents);

  ana->AddSuperModule(GeneratorMod1);
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

//   ana->AddDataset(d);

  //sync files
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_10_1_A67.root");
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_1_1_1JI.root");
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_11_1_Lpt.root");
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_12_1_MZg.root");
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_13_1_FX4.root");
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_14_1_2kc.root");
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_15_1_tUn.root");
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_16_1_kaX.root");
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_17_1_k18.root");
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_2_1_hQb.root");
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_3_1_OMp.root");
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_4_1_WGK.root");
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_5_1_Goy.root");
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_6_1_CqI.root");
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_7_1_xxV.root");
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_8_1_h4N.root");
  ana->AddFile("/temp/khahn/bambu/s12-h125zg-gf-v9/s12-h125zg-gf-v9_000_9_1_MNL.root");



  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  ana->SetOutputName(rootFile.Data());
  ana->SetCacheSize(0);


  //
  // run analysis after successful initialisation
  //
  ana->Run(!gROOT->IsBatch());
}

