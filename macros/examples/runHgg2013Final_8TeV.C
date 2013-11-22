// $Id: runHgg2013Final_8TeV.C,v 1.3 2013/11/20 18:33:07 mingyang Exp $
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TProfile.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitAna/PhysicsMod/interface/MCProcessSelectionMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PhotonPairSelector.h"
#include "MitPhysics/Mods/interface/PhotonTreeWriter.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/PhotonMvaMod.h"
#include "MitPhysics/Mods/interface/MVASystematicsMod.h"
#include "MitPhysics/Mods/interface/SeparatePileUpMod.h"

#endif

//--------------------------------------------------------------------------------------------------
void runHgg2013Final_8TeV(const char *fileset    = "0000",
			  const char *skim       = "noskim",
			  //const char *dataset = "f11--h120gg-gf-v14b-pu",
			  //const char *dataset = "r11a-pho-j16-v1",   
			  //const char *dataset = "meridiani2012-diphoj-v9",
			  //const char *dataset = "s12-h125gg-gf-v9",
			  //const char *dataset = "s12-pj40-2em-v9",
			  //const char *dataset = "s12-diphoj-3-v9",
			  //const char *dataset = "s12-zllm50-v9",
			  //const char *dataset = "f11--h121gg-gf-v14b-pu",
			  //const char *dataset = "r12a-pho-pr-v1",
			  //const char *dataset = "s12-pj40-2em-v9",
			  //const char *dataset = "s12-h125gg-vbf-v7a", 
			  //const char *dataset = "s12-h125gg-vh-v7a",   
			  //const char *dataset = "s12-h125gg-gf-v7a", 
			  //const char *dataset = "r12a-pho-j13-v1", 
			  //const char *dataset = "r12b-dph-j13-v1", 
			  //const char *dataset = "s12-zllm50-v7a", 
			  //const char *dataset = "r12a-pho-j22-v1", 
			  //const char *dataset = "r12b-dph-j22-v1", 
			  //const char *dataset = "r12c-dph-j22-v1", 
			  //const char *dataset = "r12d-dph-j22-v1", 
			  //const char *dataset = "s12-h124gg-gf-v7n",
			  //const char *dataset = "s12-diphoj-m60-v7n",
			  //const char *dataset = "s12-pjm80-2em-v7n",
			  //const char *dataset = "s12-h120gg-vh-v7n",
			  //const char *dataset = "r12a-pho-j22-v1",
			  const char *dataset = "s12-h120gg-vh-v7n",
			  //const char *dataset = "s12-h145gg-vh-v7n",
			  //const char *dataset = "s12-h120gg-gf-v7n",
			  const char *book       = "t2mit/filefi/030",
			  //const char *book       = "local/filefi/029",
			  //const char *book       = "t2mit/filefi/031",
			  const char *catalogDir = "/home/cmsprod/catalog",
			  const char *outputName = "hgg",
			  int         nEvents    = 100)
{
  //------------------------------------------------------------------------------------------------
  // some parameters get passed through the environment
  //------------------------------------------------------------------------------------------------
  char json[1024], overlap[1024];
  float overlapCut = -1;
  
  if (gSystem->Getenv("MIT_PROD_JSON"))
    sprintf(json,   "%s",gSystem->Getenv("MIT_PROD_JSON"));
  else {
    sprintf(json, "%s", "~");
    //printf(" JSON file was not properly defined. EXIT!\n");
    //return;
  } 
  
  TString jsonFile = TString("/home/mingyang/cms/json/") + TString(json);
  //TString jsonFile = TString("/home/mingyang/cms/json/") + TString("Cert_136033-149442_7TeV_Dec22ReReco_Collisions10_JSON_v4.txt");
  Bool_t  isData   = ( (jsonFile.CompareTo("/home/mingyang/cms/json/~") != 0) );
  
  if (gSystem->Getenv("MIT_PROD_OVERLAP")) {
    sprintf(overlap,"%s",gSystem->Getenv("MIT_PROD_OVERLAP"));
    if (EOF == sscanf(overlap,"%f",&overlapCut)) {
      printf(" Overlap was not properly defined. EXIT!\n");
      return;
    }
  }
  else {
     sprintf(overlap,"%s", "-1.0");
    //printf(" OVERLAP file was not properly defined. EXIT!\n");
    //return;
  } 

  printf("\n Initialization worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace mithep;
  gDebugMask  = Debug::kGeneral;
  gDebugLevel = 3;

  //------------------------------------------------------------------------------------------------
  // set up information
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSel = new RunLumiSelectionMod;
  runLumiSel->SetAcceptMC(kTRUE);                          // Monte Carlo events are always accepted

  MCProcessSelectionMod *mcselmod = new MCProcessSelectionMod;
  
  MVASystematicsMod *sysMod = new MVASystematicsMod;
  sysMod->SetMCR9Scale(1.0035, 1.0035);  
  sysMod->SetIsData(isData);
  
  // only select on run- and lumisection numbers when valid json file present
  if ((jsonFile.CompareTo("/home/mingyang/cms/json/~") != 0) &&
      (jsonFile.CompareTo("/home/mingyang/cms/json/-") != 0)   ) {
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  if ((jsonFile.CompareTo("/home/mingyang/cms/json/-") == 0)   ) {
    printf("\n WARNING -- Looking at data without JSON file: always accept.\n\n");
    runLumiSel->SetAbortIfNotAccepted(kFALSE);   // accept all events if there is no valid JSON file
  }

  printf("\n Run lumi worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // HLT information
  //------------------------------------------------------------------------------------------------
  HLTMod *hltModM = new HLTMod("HLTModM");
  hltModM->AddTrigger("HLT_Mu9");
  hltModM->AddTrigger("HLT_Mu11");
  hltModM->AddTrigger("HLT_Mu15_v1");
  hltModM->SetTrigObjsName("MyHltMuonObjs");
  hltModM->SetAbortIfNotAccepted(kFALSE);

  HLTMod *hltModE = new HLTMod("HLTModE");
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",150000,161176);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",161179,163261);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",163262,164237);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4",165085,165888);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5",165900,166967);
  hltModE->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6",166968,170053);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",170054,170759);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",170760,173198);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",173199,178380);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9",178381,179889);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10",179890,189999);
  hltModE->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",190000,999999);
  hltModE->SetTrigObjsName("MyHltElecObjs");
  hltModE->SetAbortIfNotAccepted(isData);

  HLTMod *hltModES = new HLTMod("HLTModES");
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",160000,189999);
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",160000,189999);
  hltModES->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",160000,189999);
  hltModES->AddTrigger("HLT_Ele25_WP80_PFMT40_v*",160000,189999);
  hltModES->AddTrigger("HLT_Ele27_WP80_PFMT50_v*",160000,189999);
  hltModES->AddTrigger("HLT_Ele27_WP80_v*",190000,999999);
  hltModES->SetTrigObjsName("MyHltElecSObjs");
  hltModES->SetAbortIfNotAccepted(isData);

  HLTMod *hltModP = new HLTMod("HLTModP");
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v*",160000,161176);
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v*",161216,165633);
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v*",161216,165633);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v*",161216,165633);
  hltModP->AddTrigger("HLT_Photon20_R9Id_Photon18_R9Id_v*",161216,165633);

  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_R9Id_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",165970,173198);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",165970,173198);
  
  hltModP->AddTrigger("HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v*",165970,166967);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",167039,173198);
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9Id_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_CaloIdXL_IsoXL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon26_R9Id_Photon18_R9Id_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",173236,178380);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",173236,178380);
  
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9IdT_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon26_R9IdT_Photon18_CaloIdXL_IsoXL_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon26_R9IdT_Photon18_R9IdT_Mass60_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v*",178420,189999);
  hltModP->AddTrigger("HLT_Photon36_R9Id_Photon22_R9Id_v*",178420,189999);

  //hltModP->AddTrigger("HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v*",190000,999999); 
  //hltModP->AddTrigger("HLT_Photon36_R9Id85_Photon22_R9Id85_v*",190000,999999); 

  hltModP->AddTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass*",190000,999999); 
  hltModP->AddTrigger("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v*",190000,999999); 
    
  hltModP->SetTrigObjsName("MyHltPhotObjs");
  hltModP->SetAbortIfNotAccepted(isData);
  //------------------------------------------------------------------------------------------------
  // split pfcandidates to PFPU and PFnoPU
  //------------------------------------------------------------------------------------------------
  SeparatePileUpMod* SepPUMod = new SeparatePileUpMod;
  //  SepPUMod->SetUseAllVerteces(kFALSE);
  // SepPUMod->SetVertexName("OutVtxCiC");
  SepPUMod->SetPFNoPileUpName("pfnopileupcands");
  SepPUMod->SetPFPileUpName("pfpileupcands");
  SepPUMod->SetCheckClosestZVertex(kTRUE);
   
  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof         (4.0);
  goodPVFilterMod->SetMaxAbsZ         (24.0);
  goodPVFilterMod->SetMaxRho          (2.0);
  goodPVFilterMod->SetAbortIfNotAccepted(kFALSE);
  goodPVFilterMod->SetIsMC(!isData);

  GoodPVFilterMod *goodPVFilterModE = new GoodPVFilterMod("GoodPVFilterModE");
  goodPVFilterModE->SetOutputName("GoodVertexesE");
  goodPVFilterModE->SetMinVertexNTracks(0);
  goodPVFilterModE->SetMinNDof         (4.0);
  goodPVFilterModE->SetMaxAbsZ         (24.0);
  goodPVFilterModE->SetMaxRho          (2.0);  
  goodPVFilterModE->SetIsMC(!isData);
  
  GoodPVFilterMod *goodPVFilterModES = new GoodPVFilterMod("GoodPVFilterModES");
  goodPVFilterModES->SetOutputName("GoodVertexesES");
  goodPVFilterModES->SetMinVertexNTracks(0);
  goodPVFilterModES->SetMinNDof         (4.0);
  goodPVFilterModES->SetMaxAbsZ         (24.0);
  goodPVFilterModES->SetMaxRho          (2.0);    
  

  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  //-----------------------------------
  // Lepton Selection 
  //-----------------------------------
  ElectronIDMod* eleIdMod = new ElectronIDMod;
  eleIdMod -> SetPtMin(20);  
  eleIdMod -> SetEtaMax(2.5);
  eleIdMod -> SetApplyEcalFiducial(true);
  eleIdMod -> SetIDType("Hgg_LeptonTag_2012IdHCP");  
  eleIdMod -> SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
  eleIdMod -> SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
  eleIdMod -> SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
  eleIdMod -> SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
  eleIdMod -> SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
  eleIdMod -> SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml")))); 
  eleIdMod -> SetWhichVertex(-1);
  eleIdMod -> SetD0Cut(0.02);
  eleIdMod -> SetDZCut(0.2); //h  
  eleIdMod -> SetIsoType("PFIso_HggLeptonTag2012HCP"); //h
  eleIdMod -> SetOutputName("HggLeptonTagElectrons");
  eleIdMod -> SetRhoType(RhoUtilities::CMS_RHO_RHOKT6PFJETS);
  eleIdMod -> SetPFNoPileUpName("pfnopileupcands");
  eleIdMod -> SetInvertNExpectedHitsInnerCut(kFALSE);
  eleIdMod -> SetNExpectedHitsInnerCut(1);   
  eleIdMod -> SetApplyConversionFilterType1(kTRUE);
  eleIdMod -> SetPVName(Names::gkPVBeamSpotBrn);   
  
  ElectronIDMod* softEleIdMod = new ElectronIDMod;
  softEleIdMod -> SetPtMin(10);  
  softEleIdMod -> SetEtaMax(2.5);
  softEleIdMod -> SetApplyEcalFiducial(true);
  softEleIdMod -> SetIDType("Hgg_LeptonTag_2012IdHCP");  
  softEleIdMod -> SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
  softEleIdMod -> SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
  softEleIdMod -> SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
  softEleIdMod -> SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
  softEleIdMod -> SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
  softEleIdMod -> SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml")))); 
  softEleIdMod -> SetWhichVertex(-1);
  softEleIdMod -> SetD0Cut(0.02);
  softEleIdMod -> SetDZCut(0.2); //h  
  softEleIdMod -> SetIsoType("PFIso_HggLeptonTag2012HCP"); //h
  softEleIdMod -> SetOutputName("HggLeptonTagSoftElectrons");
  softEleIdMod -> SetRhoType(RhoUtilities::CMS_RHO_RHOKT6PFJETS);
  softEleIdMod -> SetPFNoPileUpName("pfnopileupcands");
  softEleIdMod -> SetInvertNExpectedHitsInnerCut(kFALSE);
  softEleIdMod -> SetNExpectedHitsInnerCut(1);   
  softEleIdMod -> SetApplyConversionFilterType1(kTRUE);
  softEleIdMod -> SetPVName(Names::gkPVBeamSpotBrn);

 
  MuonIDMod* muonIdMod = new MuonIDMod;
  // base kinematics
  muonIdMod -> SetPtMin(20.);
  muonIdMod -> SetEtaCut(2.4);
  // base ID
  //muonIdMod -> SetIDType("Tight");
  muonIdMod -> SetIDType("muonPOG2012CutBasedIDTight");
  muonIdMod -> SetWhichVertex(-1); // this is a 'hack'.. but hopefully good enough...
  muonIdMod -> SetD0Cut(0.2);
  muonIdMod -> SetDZCut(0.5);
  muonIdMod -> SetIsoType("PFIsoBetaPUCorrected"); //h
  muonIdMod -> SetPFIsoCut(0.2); //h
  muonIdMod -> SetOutputName("HggLeptonTagMuons");
  muonIdMod -> SetPFNoPileUpName("pfnopileupcands");
  muonIdMod -> SetPFPileUpName("pfpileupcands");
  muonIdMod -> SetPVName(Names::gkPVBeamSpotBrn); 

  MuonIDMod* softMuonIdMod = new MuonIDMod;
  // base kinematics
  softMuonIdMod -> SetPtMin(10.);
  softMuonIdMod -> SetEtaCut(2.4);
  // base ID
  //softMuonIdMod -> SetIDType("Tight");
  softMuonIdMod -> SetIDType("muonPOG2012CutBasedIDTight");
  softMuonIdMod -> SetWhichVertex(-1); // this is a 'hack'.. but hopefully good enough...
  softMuonIdMod -> SetD0Cut(0.2);
  softMuonIdMod -> SetDZCut(0.5);
  softMuonIdMod -> SetIsoType("PFIsoBetaPUCorrected"); //h
  softMuonIdMod -> SetPFIsoCut(0.2); //h
  softMuonIdMod -> SetOutputName("HggLeptonTagSoftMuons");
  softMuonIdMod -> SetPFNoPileUpName("pfnopileupcands");
  softMuonIdMod -> SetPFPileUpName("pfpileupcands");
  softMuonIdMod -> SetPVName(Names::gkPVBeamSpotBrn); 
  
  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");
  
  PublisherMod<PFJet,Jet> *pubJetOpen = new PublisherMod<PFJet,Jet>("JetPubOpen");
  pubJetOpen->SetInputName("AKt5PFJets");
  pubJetOpen->SetOutputName("PubAKt5PFJetsOpen");  
  
  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  /*if(isData){ 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/GR_P_V42_AN3_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/GR_P_V42_AN3_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/GR_P_V42_AN3_L3Absolute_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/GR_P_V42_AN3_L2L3Residual_AK5PF.txt")).Data())); 
  }
  else {
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START53_V15_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START53_V15_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/START53_V15_L3Absolute_AK5PF.txt")).Data())); 
    }*/
  /*if(isData){ 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L3Absolute_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L2L3Residual_AK5PF.txt")).Data())); 
  }
  else {
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_MC_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_MC_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_MC_L3Absolute_AK5PF.txt")).Data())); 
    }*/
  if(isData){ 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L3Absolute_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt")).Data())); 
  }
  else {
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L3Absolute_AK5PF.txt")).Data())); 
  }
  jetCorr->SetInputName(pubJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");    
  
  
  Bool_t excludedoubleprompt = kFALSE;
  if (TString(dataset).Contains("-pj")) {
    mcselmod->ExcludeProcess(18);
    mcselmod->ExcludeProcess(114);
    excludedoubleprompt = kTRUE;
  }
  
  if (TString(dataset).Contains("-qcd")) {
    excludedoubleprompt = kTRUE;
  }
  
  Bool_t is25 = kFALSE;
  if (TString(book).Contains("025")) is25 = kTRUE;
  
  PhotonMvaMod *photreg = new PhotonMvaMod;
  //photreg->SetRegressionVersion(3);
  //photreg->SetRegressionWeights(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/gbrv3ph_52x.root")).Data()));
  //photreg->SetOutputName("GoodPhotonsRegr");
  //photreg->SetApplyShowerRescaling(kTRUE);
  //photreg->SetIsData(isData);
  photreg->SetRegressionVersion(5);
  photreg->SetRegressionWeights(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/regweights_v5_forest_ph.root")).Data()));
  photreg->SetOutputName("GoodPhotonsRegr");
  //photreg->SetApplyShowerRescaling(kTRUE);
  photreg->SetIsData(isData);

  
  PhotonPairSelector         *photcic = new PhotonPairSelector("PhotonPairSelectorCiC");
  photcic->SetOutputName("GoodPhotonsCIC");
  photcic->SetOutputVtxName("OutVtxCiC");        
  photcic->SetPhotonSelType("CiCPFSelection");
  photcic->SetVertexSelType("CiCMVA2012Selection");
  photcic->SetUseSingleLegConversions(kFALSE);
  photcic->DoMCSmear(kTRUE);
  photcic->DoDataEneCorr(kFALSE);
  photcic->SetPhotonsFromBranch(kFALSE);
  photcic->SetInputPhotonsName(photreg->GetOutputName());
  
  //------------------------------------------2012 HCP--------------------------------------------------------------
  photcic->SetMCSmearFactors2012HCP(0.0075,0.0075,0.0086,0.0086,0.0122,0.0188,0.0163,0.0198,0.0186,0.0192);
  photcic->SetMCSmearFactors2012HCPMVA(0.0075,0.0075,0.0086,0.0086,0.0122,0.0188,0.0163,0.0198,0.0186,0.0192);
  photcic->UseSpecialSmearForDPMVA(true);

  photcic->AddEnCorrPerRun2012HCP(190645,190781,0.9894,0.9894,0.9922,0.9922,0.9876,0.9982,0.9856,0.9919,0.9812,0.9862);
  photcic->AddEnCorrPerRun2012HCP(190782,191042,0.9961,0.9961,0.9989,0.9989,0.9909,1.0014,0.9854,0.9917,0.976,0.981);
  photcic->AddEnCorrPerRun2012HCP(191043,191720,0.9902,0.9902,0.9931,0.9931,0.9858,0.9963,0.988,0.9944,0.9795,0.9845);
  photcic->AddEnCorrPerRun2012HCP(191721,193833,0.9894,0.9894,0.9922,0.9922,0.9876,0.9982,0.9853,0.9917,0.9791,0.9841);
  photcic->AddEnCorrPerRun2012HCP(193834,194116,0.99,0.99,0.9929,0.9929,0.9865,0.997,0.9842,0.9906,0.9793,0.9843);
  photcic->AddEnCorrPerRun2012HCP(194117,194427,0.9907,0.9907,0.9935,0.9935,0.987,0.9975,0.9855,0.9918,0.9805,0.9855);
  photcic->AddEnCorrPerRun2012HCP(194428,194618,0.9901,0.9901,0.9929,0.9929,0.9867,0.9973,0.9841,0.9905,0.9786,0.9836);
  photcic->AddEnCorrPerRun2012HCP(194619,194789,0.9904,0.9904,0.9932,0.9932,0.9874,0.9979,0.9876,0.994,0.9788,0.9838);
  photcic->AddEnCorrPerRun2012HCP(194790,195111,0.991,0.991,0.9938,0.9938,0.9887,0.9992,0.99,0.9963,0.9803,0.9853);
  photcic->AddEnCorrPerRun2012HCP(195112,195377,0.9912,0.9912,0.994,0.994,0.9871,0.9976,0.9885,0.9948,0.9817,0.9866);
  photcic->AddEnCorrPerRun2012HCP(195378,195398,0.9903,0.9903,0.9931,0.9931,0.9863,0.9968,0.9867,0.993,0.9801,0.9851);
  photcic->AddEnCorrPerRun2012HCP(195399,195657,0.9908,0.9908,0.9936,0.9936,0.9888,0.9993,0.9866,0.9929,0.9806,0.9856);
  photcic->AddEnCorrPerRun2012HCP(195658,195918,0.9914,0.9914,0.9942,0.9942,0.9878,0.9983,0.9898,0.9962,0.9797,0.9846);
  photcic->AddEnCorrPerRun2012HCP(195919,196198,0.9908,0.9908,0.9936,0.9936,0.9874,0.998,0.9878,0.9942,0.9796,0.9846);
  photcic->AddEnCorrPerRun2012HCP(196199,196356,0.9915,0.9915,0.9943,0.9943,0.9877,0.9983,0.988,0.9944,0.9781,0.9831);
  photcic->AddEnCorrPerRun2012HCP(196357,198115,0.991,0.991,0.9938,0.9938,0.9871,0.9977,0.9859,0.9923,0.9791,0.9841);
  photcic->AddEnCorrPerRun2012HCP(198116,198940,0.9906,0.9906,0.9934,0.9934,0.9864,0.997,0.9877,0.994,0.9818,0.9868);
  photcic->AddEnCorrPerRun2012HCP(198941,199317,0.9907,0.9907,0.9936,0.9936,0.9865,0.997,0.9856,0.992,0.9805,0.9855);
  photcic->AddEnCorrPerRun2012HCP(199318,199428,0.9904,0.9904,0.9933,0.9933,0.9868,0.9973,0.9857,0.9921,0.9808,0.9858);
  photcic->AddEnCorrPerRun2012HCP(199429,199697,0.9907,0.9907,0.9935,0.9935,0.987,0.9976,0.9853,0.9916,0.9803,0.9853);
  photcic->AddEnCorrPerRun2012HCP(199698,199832,0.991,0.991,0.9938,0.9938,0.9879,0.9985,0.9845,0.9909,0.9803,0.9853);
  photcic->AddEnCorrPerRun2012HCP(199833,199960,0.9912,0.9912,0.994,0.994,0.9877,0.9983,0.9875,0.9938,0.9823,0.9873);
  photcic->AddEnCorrPerRun2012HCP(199961,200151,0.9913,0.9913,0.9942,0.9942,0.9865,0.997,0.986,0.9924,0.9817,0.9867);
  photcic->AddEnCorrPerRun2012HCP(200152,200490,0.9912,0.9912,0.994,0.994,0.9879,0.9984,0.9878,0.9942,0.9809,0.9858);
  photcic->AddEnCorrPerRun2012HCP(200491,200991,0.9919,0.9919,0.9947,0.9947,0.9866,0.9971,0.9861,0.9924,0.9814,0.9864);
  photcic->AddEnCorrPerRun2012HCP(200992,201201,0.9909,0.9909,0.9937,0.9937,0.9864,0.9969,0.985,0.9914,0.9833,0.9883);
  photcic->AddEnCorrPerRun2012HCP(201202,201624,0.9915,0.9915,0.9943,0.9943,0.9886,0.9992,0.9857,0.992,0.9827,0.9876);
  photcic->AddEnCorrPerRun2012HCP(201625,201707,0.9917,0.9917,0.9945,0.9945,0.9878,0.9983,0.9856,0.9919,0.982,0.987);
  photcic->AddEnCorrPerRun2012HCP(201708,202059,0.9916,0.9916,0.9944,0.9944,0.9869,0.9974,0.9866,0.9929,0.9833,0.9883);
  photcic->AddEnCorrPerRun2012HCP(202060,202204,0.9919,0.9919,0.9947,0.9947,0.9871,0.9977,0.985,0.9914,0.9827,0.9877);
  photcic->AddEnCorrPerRun2012HCP(202205,202332,0.9923,0.9923,0.9951,0.9951,0.9883,0.9989,0.9887,0.995,0.9825,0.9874);
  photcic->AddEnCorrPerRun2012HCP(202333,202972,0.9921,0.9921,0.9949,0.9949,0.9888,0.9994,0.9857,0.992,0.9812,0.9861);
  photcic->AddEnCorrPerRun2012HCP(202973,203002,0.9916,0.9916,0.9944,0.9944,0.9857,0.9962,0.986,0.9924,0.9816,0.9865);
  photcic->AddEnCorrPerRun2012HCP(203003,203852,0.993,0.993,0.9958,0.9958,0.9831,0.9936,0.9705,0.9769,0.9728,0.9778);
  photcic->AddEnCorrPerRun2012HCP(203853,204099,0.9906,0.9906,0.9935,0.9935,0.9886,0.9991,0.9851,0.9914,0.9777,0.9827);
  photcic->AddEnCorrPerRun2012HCP(204100,204562,0.991,0.991,0.9939,0.9939,0.9894,0.9999,0.9868,0.9931,0.9814,0.9864);
  photcic->AddEnCorrPerRun2012HCP(204563,205085,0.991,0.991,0.9938,0.9938,0.9901,1.0006,0.9871,0.9934,0.9789,0.9839);
  photcic->AddEnCorrPerRun2012HCP(205086,205310,0.9911,0.9911,0.9939,0.9939,0.9886,0.9991,0.9868,0.9932,0.9782,0.9832);
  photcic->AddEnCorrPerRun2012HCP(205311,205617,0.9909,0.9909,0.9938,0.9938,0.989,0.9995,0.9852,0.9916,0.9766,0.9816);
  photcic->AddEnCorrPerRun2012HCP(205618,205825,0.9914,0.9914,0.9942,0.9942,0.9877,0.9982,0.9866,0.993,0.9792,0.9842);
  photcic->AddEnCorrPerRun2012HCP(205826,206207,0.9921,0.9921,0.9949,0.9949,0.9908,1.0012,0.9882,0.9945,0.9816,0.9866);
  photcic->AddEnCorrPerRun2012HCP(206208,206389,0.9918,0.9918,0.9946,0.9946,0.989,0.9995,0.9842,0.9906,0.9794,0.9844);
  photcic->AddEnCorrPerRun2012HCP(206390,206483,0.9917,0.9917,0.9945,0.9945,0.989,0.9995,0.9894,0.9957,0.9824,0.9873);
  photcic->AddEnCorrPerRun2012HCP(206484,206597,0.9917,0.9917,0.9945,0.9945,0.9892,0.9997,0.9886,0.9949,0.9791,0.9841);
  photcic->AddEnCorrPerRun2012HCP(206598,206896,0.9913,0.9913,0.9941,0.9941,0.9893,0.9998,0.9846,0.991,0.9789,0.9839);
  photcic->AddEnCorrPerRun2012HCP(206897,207220,0.9923,0.9923,0.9952,0.9952,0.9892,0.9997,0.9872,0.9935,0.9813,0.9863);
  photcic->AddEnCorrPerRun2012HCP(207221,207315,0.9921,0.9921,0.9949,0.9949,0.9901,1.0006,0.987,0.9934,0.9798,0.9848);
  photcic->AddEnCorrPerRun2012HCP(207316,207489,0.9922,0.9922,0.995,0.995,0.9896,1.0001,0.9881,0.9944,0.9791,0.9841);
  photcic->AddEnCorrPerRun2012HCP(207490,207919,0.9922,0.9922,0.9951,0.9951,0.9901,1.0006,0.9866,0.993,0.9769,0.9819);
  photcic->AddEnCorrPerRun2012HCP(207920,208351,0.9919,0.9919,0.9947,0.9947,0.9895,1,0.9867,0.993,0.9813,0.9863);
  photcic->AddEnCorrPerRun2012HCP(208352,208686,0.9924,0.9924,0.9952,0.9952,0.9906,1.0011,0.9871,0.9934,0.9809,0.9859);
  
  

  //-----------------------------------------------------------------------------------------------------------------
  
  //  photcic->SetDoMCR9Scaling(kTRUE);
  //  photcic->SetMCR9Scale(1.0035, 1.0035);
  photcic->SetDoShowerShapeScaling(kFALSE); 
  photcic->SetShowerShapeType("2012ShowerShape");
  //photcic->SetDoMCErrScaling(kTRUE);
  //photcic->SetMCErrScale(1.07, 1.045); 
  //photcic->SetMCErrScale(1, 1); //ming:scale(sigE/E)
  photcic->SetJetsName(jetCorr->GetOutputName());    
  //photcic->SetRescaledBeamspotWidth(5.0);
  photcic->SetIsData(isData);
  photcic->SetApplyLeptonTag(kTRUE);
  photcic->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  photcic->SetLeptonTagMuonsName("HggLeptonTagMuons");  
  photcic->Set2012HCP(kTRUE);
  photcic->SetIdMVAType("2013FinalIdMVA_8TeV");
  
  PhotonPairSelector         *photcicnoeleveto = new PhotonPairSelector("PhotonPairSelectorCiCInvertEleVeto");
  photcicnoeleveto->SetOutputName("GoodPhotonsCICNoEleVeto");
  photcicnoeleveto->SetOutputVtxName("OutVtxCiCInvertEleVeto");      
  photcicnoeleveto->SetPhotonSelType("CiCPFSelection");
  photcicnoeleveto->SetVertexSelType("CiCMVA2012Selection");
  photcicnoeleveto->SetUseSingleLegConversions(kFALSE);
  photcicnoeleveto->DoMCSmear(kTRUE);
  photcicnoeleveto->DoDataEneCorr(kFALSE);
  photcicnoeleveto->SetPhotonsFromBranch(kFALSE);
  photcicnoeleveto->SetInputPhotonsName(photreg->GetOutputName());
  
  //------------------------------------------2012 HCP--------------------------------------------------------------
  photcicnoeleveto->SetMCSmearFactors2012HCP(0.0075,0.0075,0.0086,0.0086,0.0122,0.0188,0.0163,0.0198,0.0186,0.0192);
  photcicnoeleveto->SetMCSmearFactors2012HCPMVA(0.0075,0.0075,0.0086,0.0086,0.0122,0.0188,0.0163,0.0198,0.0186,0.0192);
  photcicnoeleveto->UseSpecialSmearForDPMVA(true);

  photcicnoeleveto->AddEnCorrPerRun2012HCP(190645,190781,0.9894,0.9894,0.9922,0.9922,0.9876,0.9982,0.9856,0.9919,0.9812,0.9862);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(190782,191042,0.9961,0.9961,0.9989,0.9989,0.9909,1.0014,0.9854,0.9917,0.976,0.981);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(191043,191720,0.9902,0.9902,0.9931,0.9931,0.9858,0.9963,0.988,0.9944,0.9795,0.9845);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(191721,193833,0.9894,0.9894,0.9922,0.9922,0.9876,0.9982,0.9853,0.9917,0.9791,0.9841);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(193834,194116,0.99,0.99,0.9929,0.9929,0.9865,0.997,0.9842,0.9906,0.9793,0.9843);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(194117,194427,0.9907,0.9907,0.9935,0.9935,0.987,0.9975,0.9855,0.9918,0.9805,0.9855);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(194428,194618,0.9901,0.9901,0.9929,0.9929,0.9867,0.9973,0.9841,0.9905,0.9786,0.9836);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(194619,194789,0.9904,0.9904,0.9932,0.9932,0.9874,0.9979,0.9876,0.994,0.9788,0.9838);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(194790,195111,0.991,0.991,0.9938,0.9938,0.9887,0.9992,0.99,0.9963,0.9803,0.9853);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(195112,195377,0.9912,0.9912,0.994,0.994,0.9871,0.9976,0.9885,0.9948,0.9817,0.9866);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(195378,195398,0.9903,0.9903,0.9931,0.9931,0.9863,0.9968,0.9867,0.993,0.9801,0.9851);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(195399,195657,0.9908,0.9908,0.9936,0.9936,0.9888,0.9993,0.9866,0.9929,0.9806,0.9856);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(195658,195918,0.9914,0.9914,0.9942,0.9942,0.9878,0.9983,0.9898,0.9962,0.9797,0.9846);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(195919,196198,0.9908,0.9908,0.9936,0.9936,0.9874,0.998,0.9878,0.9942,0.9796,0.9846);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(196199,196356,0.9915,0.9915,0.9943,0.9943,0.9877,0.9983,0.988,0.9944,0.9781,0.9831);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(196357,198115,0.991,0.991,0.9938,0.9938,0.9871,0.9977,0.9859,0.9923,0.9791,0.9841);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(198116,198940,0.9906,0.9906,0.9934,0.9934,0.9864,0.997,0.9877,0.994,0.9818,0.9868);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(198941,199317,0.9907,0.9907,0.9936,0.9936,0.9865,0.997,0.9856,0.992,0.9805,0.9855);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(199318,199428,0.9904,0.9904,0.9933,0.9933,0.9868,0.9973,0.9857,0.9921,0.9808,0.9858);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(199429,199697,0.9907,0.9907,0.9935,0.9935,0.987,0.9976,0.9853,0.9916,0.9803,0.9853);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(199698,199832,0.991,0.991,0.9938,0.9938,0.9879,0.9985,0.9845,0.9909,0.9803,0.9853);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(199833,199960,0.9912,0.9912,0.994,0.994,0.9877,0.9983,0.9875,0.9938,0.9823,0.9873);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(199961,200151,0.9913,0.9913,0.9942,0.9942,0.9865,0.997,0.986,0.9924,0.9817,0.9867);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(200152,200490,0.9912,0.9912,0.994,0.994,0.9879,0.9984,0.9878,0.9942,0.9809,0.9858);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(200491,200991,0.9919,0.9919,0.9947,0.9947,0.9866,0.9971,0.9861,0.9924,0.9814,0.9864);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(200992,201201,0.9909,0.9909,0.9937,0.9937,0.9864,0.9969,0.985,0.9914,0.9833,0.9883);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(201202,201624,0.9915,0.9915,0.9943,0.9943,0.9886,0.9992,0.9857,0.992,0.9827,0.9876);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(201625,201707,0.9917,0.9917,0.9945,0.9945,0.9878,0.9983,0.9856,0.9919,0.982,0.987);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(201708,202059,0.9916,0.9916,0.9944,0.9944,0.9869,0.9974,0.9866,0.9929,0.9833,0.9883);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(202060,202204,0.9919,0.9919,0.9947,0.9947,0.9871,0.9977,0.985,0.9914,0.9827,0.9877);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(202205,202332,0.9923,0.9923,0.9951,0.9951,0.9883,0.9989,0.9887,0.995,0.9825,0.9874);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(202333,202972,0.9921,0.9921,0.9949,0.9949,0.9888,0.9994,0.9857,0.992,0.9812,0.9861);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(202973,203002,0.9916,0.9916,0.9944,0.9944,0.9857,0.9962,0.986,0.9924,0.9816,0.9865);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(203003,203852,0.993,0.993,0.9958,0.9958,0.9831,0.9936,0.9705,0.9769,0.9728,0.9778);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(203853,204099,0.9906,0.9906,0.9935,0.9935,0.9886,0.9991,0.9851,0.9914,0.9777,0.9827);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(204100,204562,0.991,0.991,0.9939,0.9939,0.9894,0.9999,0.9868,0.9931,0.9814,0.9864);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(204563,205085,0.991,0.991,0.9938,0.9938,0.9901,1.0006,0.9871,0.9934,0.9789,0.9839);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(205086,205310,0.9911,0.9911,0.9939,0.9939,0.9886,0.9991,0.9868,0.9932,0.9782,0.9832);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(205311,205617,0.9909,0.9909,0.9938,0.9938,0.989,0.9995,0.9852,0.9916,0.9766,0.9816);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(205618,205825,0.9914,0.9914,0.9942,0.9942,0.9877,0.9982,0.9866,0.993,0.9792,0.9842);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(205826,206207,0.9921,0.9921,0.9949,0.9949,0.9908,1.0012,0.9882,0.9945,0.9816,0.9866);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(206208,206389,0.9918,0.9918,0.9946,0.9946,0.989,0.9995,0.9842,0.9906,0.9794,0.9844);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(206390,206483,0.9917,0.9917,0.9945,0.9945,0.989,0.9995,0.9894,0.9957,0.9824,0.9873);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(206484,206597,0.9917,0.9917,0.9945,0.9945,0.9892,0.9997,0.9886,0.9949,0.9791,0.9841);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(206598,206896,0.9913,0.9913,0.9941,0.9941,0.9893,0.9998,0.9846,0.991,0.9789,0.9839);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(206897,207220,0.9923,0.9923,0.9952,0.9952,0.9892,0.9997,0.9872,0.9935,0.9813,0.9863);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(207221,207315,0.9921,0.9921,0.9949,0.9949,0.9901,1.0006,0.987,0.9934,0.9798,0.9848);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(207316,207489,0.9922,0.9922,0.995,0.995,0.9896,1.0001,0.9881,0.9944,0.9791,0.9841);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(207490,207919,0.9922,0.9922,0.9951,0.9951,0.9901,1.0006,0.9866,0.993,0.9769,0.9819);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(207920,208351,0.9919,0.9919,0.9947,0.9947,0.9895,1,0.9867,0.993,0.9813,0.9863);
  photcicnoeleveto->AddEnCorrPerRun2012HCP(208352,208686,0.9924,0.9924,0.9952,0.9952,0.9906,1.0011,0.9871,0.9934,0.9809,0.9859);

  //-----------------------------------------------------------------------------------------------------------------
  
  //photcicnoeleveto->SetDoMCR9Scaling(kTRUE);
  //photcicnoeleveto->SetMCR9Scale(1.0035, 1.0035);
  photcicnoeleveto->SetDoShowerShapeScaling(kFALSE);
  photcicnoeleveto->SetShowerShapeType("2012ShowerShape");
  //photcicnoeleveto->SetDoMCErrScaling(kTRUE);
  //photcicnoeleveto->SetMCErrScale(1.07, 1.045);    
  //photcicnoeleveto->SetMCErrScale(1, 1);    
  photcicnoeleveto->SetApplyEleVeto(kFALSE);
  photcicnoeleveto->SetInvertElectronVeto(kTRUE);
  photcicnoeleveto->SetJetsName(jetCorr->GetOutputName());  
  //photcicnoeleveto->SetRescaledBeamspotWidth(5.0);
  photcicnoeleveto->SetIsData(isData);
  photcicnoeleveto->SetApplyLeptonTag(kTRUE);
  photcicnoeleveto->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  photcicnoeleveto->SetLeptonTagMuonsName("HggLeptonTagMuons");  
  photcicnoeleveto->Set2012HCP(kTRUE);  
  photcicnoeleveto->SetIdMVAType("2013FinalIdMVA_8TeV");
 
  PhotonPairSelector         *photpresel = new PhotonPairSelector("PhotonPairSelectorPresel");
  photpresel->SetOutputName("GoodPhotonsPresel");
  photpresel->SetPhotonSelType("MITPFSelectionNoEcal");
  photpresel->SetVertexSelType("CiCMVA2012Selection");
  photpresel->SetUseSingleLegConversions(kFALSE);
  photpresel->SetIdMVAType("2013FinalIdMVA_8TeV");
  photpresel->DoMCSmear(kTRUE);
  //photpresel->DoDataEneCorr(kFALSE);
  photpresel->DoDataEneCorr(kTRUE);
  photpresel->SetPhotonsFromBranch(kFALSE);
  photpresel->SetInputPhotonsName(photreg->GetOutputName());
  //------------------------------------------2012 HCP--------------------------------------------------------------
  photpresel->SetMCSmearFactors2012HCP(0.0075,0.0075,0.0086,0.0086,0.0122,0.0188,0.0163,0.0198,0.0186,0.0192);
  photpresel->SetMCSmearFactors2012HCPMVA(0.0075,0.0075,0.0086,0.0086,0.0122,0.0188,0.0163,0.0198,0.0186,0.0192);
  photpresel->UseSpecialSmearForDPMVA(true);

  photpresel->AddEnCorrPerRun2012HCP(190645,190781,0.9894,0.9894,0.9922,0.9922,0.9876,0.9982,0.9856,0.9919,0.9812,0.9862);
  photpresel->AddEnCorrPerRun2012HCP(190782,191042,0.9961,0.9961,0.9989,0.9989,0.9909,1.0014,0.9854,0.9917,0.976,0.981);
  photpresel->AddEnCorrPerRun2012HCP(191043,191720,0.9902,0.9902,0.9931,0.9931,0.9858,0.9963,0.988,0.9944,0.9795,0.9845);
  photpresel->AddEnCorrPerRun2012HCP(191721,193833,0.9894,0.9894,0.9922,0.9922,0.9876,0.9982,0.9853,0.9917,0.9791,0.9841);
  photpresel->AddEnCorrPerRun2012HCP(193834,194116,0.99,0.99,0.9929,0.9929,0.9865,0.997,0.9842,0.9906,0.9793,0.9843);
  photpresel->AddEnCorrPerRun2012HCP(194117,194427,0.9907,0.9907,0.9935,0.9935,0.987,0.9975,0.9855,0.9918,0.9805,0.9855);
  photpresel->AddEnCorrPerRun2012HCP(194428,194618,0.9901,0.9901,0.9929,0.9929,0.9867,0.9973,0.9841,0.9905,0.9786,0.9836);
  photpresel->AddEnCorrPerRun2012HCP(194619,194789,0.9904,0.9904,0.9932,0.9932,0.9874,0.9979,0.9876,0.994,0.9788,0.9838);
  photpresel->AddEnCorrPerRun2012HCP(194790,195111,0.991,0.991,0.9938,0.9938,0.9887,0.9992,0.99,0.9963,0.9803,0.9853);
  photpresel->AddEnCorrPerRun2012HCP(195112,195377,0.9912,0.9912,0.994,0.994,0.9871,0.9976,0.9885,0.9948,0.9817,0.9866);
  photpresel->AddEnCorrPerRun2012HCP(195378,195398,0.9903,0.9903,0.9931,0.9931,0.9863,0.9968,0.9867,0.993,0.9801,0.9851);
  photpresel->AddEnCorrPerRun2012HCP(195399,195657,0.9908,0.9908,0.9936,0.9936,0.9888,0.9993,0.9866,0.9929,0.9806,0.9856);
  photpresel->AddEnCorrPerRun2012HCP(195658,195918,0.9914,0.9914,0.9942,0.9942,0.9878,0.9983,0.9898,0.9962,0.9797,0.9846);
  photpresel->AddEnCorrPerRun2012HCP(195919,196198,0.9908,0.9908,0.9936,0.9936,0.9874,0.998,0.9878,0.9942,0.9796,0.9846);
  photpresel->AddEnCorrPerRun2012HCP(196199,196356,0.9915,0.9915,0.9943,0.9943,0.9877,0.9983,0.988,0.9944,0.9781,0.9831);
  photpresel->AddEnCorrPerRun2012HCP(196357,198115,0.991,0.991,0.9938,0.9938,0.9871,0.9977,0.9859,0.9923,0.9791,0.9841);
  photpresel->AddEnCorrPerRun2012HCP(198116,198940,0.9906,0.9906,0.9934,0.9934,0.9864,0.997,0.9877,0.994,0.9818,0.9868);
  photpresel->AddEnCorrPerRun2012HCP(198941,199317,0.9907,0.9907,0.9936,0.9936,0.9865,0.997,0.9856,0.992,0.9805,0.9855);
  photpresel->AddEnCorrPerRun2012HCP(199318,199428,0.9904,0.9904,0.9933,0.9933,0.9868,0.9973,0.9857,0.9921,0.9808,0.9858);
  photpresel->AddEnCorrPerRun2012HCP(199429,199697,0.9907,0.9907,0.9935,0.9935,0.987,0.9976,0.9853,0.9916,0.9803,0.9853);
  photpresel->AddEnCorrPerRun2012HCP(199698,199832,0.991,0.991,0.9938,0.9938,0.9879,0.9985,0.9845,0.9909,0.9803,0.9853);
  photpresel->AddEnCorrPerRun2012HCP(199833,199960,0.9912,0.9912,0.994,0.994,0.9877,0.9983,0.9875,0.9938,0.9823,0.9873);
  photpresel->AddEnCorrPerRun2012HCP(199961,200151,0.9913,0.9913,0.9942,0.9942,0.9865,0.997,0.986,0.9924,0.9817,0.9867);
  photpresel->AddEnCorrPerRun2012HCP(200152,200490,0.9912,0.9912,0.994,0.994,0.9879,0.9984,0.9878,0.9942,0.9809,0.9858);
  photpresel->AddEnCorrPerRun2012HCP(200491,200991,0.9919,0.9919,0.9947,0.9947,0.9866,0.9971,0.9861,0.9924,0.9814,0.9864);
  photpresel->AddEnCorrPerRun2012HCP(200992,201201,0.9909,0.9909,0.9937,0.9937,0.9864,0.9969,0.985,0.9914,0.9833,0.9883);
  photpresel->AddEnCorrPerRun2012HCP(201202,201624,0.9915,0.9915,0.9943,0.9943,0.9886,0.9992,0.9857,0.992,0.9827,0.9876);
  photpresel->AddEnCorrPerRun2012HCP(201625,201707,0.9917,0.9917,0.9945,0.9945,0.9878,0.9983,0.9856,0.9919,0.982,0.987);
  photpresel->AddEnCorrPerRun2012HCP(201708,202059,0.9916,0.9916,0.9944,0.9944,0.9869,0.9974,0.9866,0.9929,0.9833,0.9883);
  photpresel->AddEnCorrPerRun2012HCP(202060,202204,0.9919,0.9919,0.9947,0.9947,0.9871,0.9977,0.985,0.9914,0.9827,0.9877);
  photpresel->AddEnCorrPerRun2012HCP(202205,202332,0.9923,0.9923,0.9951,0.9951,0.9883,0.9989,0.9887,0.995,0.9825,0.9874);
  photpresel->AddEnCorrPerRun2012HCP(202333,202972,0.9921,0.9921,0.9949,0.9949,0.9888,0.9994,0.9857,0.992,0.9812,0.9861);
  photpresel->AddEnCorrPerRun2012HCP(202973,203002,0.9916,0.9916,0.9944,0.9944,0.9857,0.9962,0.986,0.9924,0.9816,0.9865);
  photpresel->AddEnCorrPerRun2012HCP(203003,203852,0.993,0.993,0.9958,0.9958,0.9831,0.9936,0.9705,0.9769,0.9728,0.9778);
  photpresel->AddEnCorrPerRun2012HCP(203853,204099,0.9906,0.9906,0.9935,0.9935,0.9886,0.9991,0.9851,0.9914,0.9777,0.9827);
  photpresel->AddEnCorrPerRun2012HCP(204100,204562,0.991,0.991,0.9939,0.9939,0.9894,0.9999,0.9868,0.9931,0.9814,0.9864);
  photpresel->AddEnCorrPerRun2012HCP(204563,205085,0.991,0.991,0.9938,0.9938,0.9901,1.0006,0.9871,0.9934,0.9789,0.9839);
  photpresel->AddEnCorrPerRun2012HCP(205086,205310,0.9911,0.9911,0.9939,0.9939,0.9886,0.9991,0.9868,0.9932,0.9782,0.9832);
  photpresel->AddEnCorrPerRun2012HCP(205311,205617,0.9909,0.9909,0.9938,0.9938,0.989,0.9995,0.9852,0.9916,0.9766,0.9816);
  photpresel->AddEnCorrPerRun2012HCP(205618,205825,0.9914,0.9914,0.9942,0.9942,0.9877,0.9982,0.9866,0.993,0.9792,0.9842);
  photpresel->AddEnCorrPerRun2012HCP(205826,206207,0.9921,0.9921,0.9949,0.9949,0.9908,1.0012,0.9882,0.9945,0.9816,0.9866);
  photpresel->AddEnCorrPerRun2012HCP(206208,206389,0.9918,0.9918,0.9946,0.9946,0.989,0.9995,0.9842,0.9906,0.9794,0.9844);
  photpresel->AddEnCorrPerRun2012HCP(206390,206483,0.9917,0.9917,0.9945,0.9945,0.989,0.9995,0.9894,0.9957,0.9824,0.9873);
  photpresel->AddEnCorrPerRun2012HCP(206484,206597,0.9917,0.9917,0.9945,0.9945,0.9892,0.9997,0.9886,0.9949,0.9791,0.9841);
  photpresel->AddEnCorrPerRun2012HCP(206598,206896,0.9913,0.9913,0.9941,0.9941,0.9893,0.9998,0.9846,0.991,0.9789,0.9839);
  photpresel->AddEnCorrPerRun2012HCP(206897,207220,0.9923,0.9923,0.9952,0.9952,0.9892,0.9997,0.9872,0.9935,0.9813,0.9863);
  photpresel->AddEnCorrPerRun2012HCP(207221,207315,0.9921,0.9921,0.9949,0.9949,0.9901,1.0006,0.987,0.9934,0.9798,0.9848);
  photpresel->AddEnCorrPerRun2012HCP(207316,207489,0.9922,0.9922,0.995,0.995,0.9896,1.0001,0.9881,0.9944,0.9791,0.9841);
  photpresel->AddEnCorrPerRun2012HCP(207490,207919,0.9922,0.9922,0.9951,0.9951,0.9901,1.0006,0.9866,0.993,0.9769,0.9819);
  photpresel->AddEnCorrPerRun2012HCP(207920,208351,0.9919,0.9919,0.9947,0.9947,0.9895,1,0.9867,0.993,0.9813,0.9863);
  photpresel->AddEnCorrPerRun2012HCP(208352,208686,0.9924,0.9924,0.9952,0.9952,0.9906,1.0011,0.9871,0.9934,0.9809,0.9859);

  
  //-----------------------------------------------------------------------------------------------------------------
  //photpresel->SetDoMCR9Scaling(kTRUE);
  photpresel->SetDoShowerShapeScaling(kFALSE);
  photpresel->SetShowerShapeType("2012ShowerShape");
  //photpresel->SetDoMCErrScaling(kTRUE);
  //photpresel->SetMCErrScale(1.07, 1.045); 
  //photpresel->SetMCErrScale(1, 1); 
  photpresel->SetJetsName(jetCorr->GetOutputName()); 
  //photpresel->SetRescaledBeamspotWidth(5.0);
  photpresel->SetIsData(isData);
  photpresel->SetApplyLeptonTag(kTRUE);
  photpresel->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  photpresel->SetLeptonTagMuonsName("HggLeptonTagMuons");  
  photpresel->Set2012HCP(kTRUE); 
  //

 
  PhotonPairSelector         *photpreselinverteleveto = new PhotonPairSelector("PhotonPairSelectorPreselInvertEleVeto");
  photpreselinverteleveto->SetOutputName("GoodPhotonsPreselInvertEleVeto");
  photpreselinverteleveto->SetOutputVtxName("OutVtxPreselInvertEleVeto");    
  photpreselinverteleveto->SetPhotonSelType("MITPFSelectionNoEcal");
  photpreselinverteleveto->SetIdMVAType("2013FinalIdMVA_8TeV");
  //photpreselinverteleveto->SetVertexSelType("CiCMVA2012Selection");//ming change
  photpreselinverteleveto->SetUseSingleLegConversions(kFALSE);
  photpreselinverteleveto->DoMCSmear(kTRUE);
  photpreselinverteleveto->DoDataEneCorr(kFALSE);
  photpreselinverteleveto->SetPhotonsFromBranch(kFALSE);
  photpreselinverteleveto->SetInputPhotonsName(photreg->GetOutputName());
  //------------------------------------------2012 HCP--------------------------------------------------------------
  photpreselinverteleveto->SetMCSmearFactors2012HCP(0.0075,0.0075,0.0086,0.0086,0.0122,0.0188,0.0163,0.0198,0.0186,0.0192);
  photpreselinverteleveto->SetMCSmearFactors2012HCPMVA(0.0075,0.0075,0.0086,0.0086,0.0122,0.0188,0.0163,0.0198,0.0186,0.0192);
  photpreselinverteleveto->UseSpecialSmearForDPMVA(true);

  photpreselinverteleveto->AddEnCorrPerRun2012HCP(190645,190781,0.9894,0.9894,0.9922,0.9922,0.9876,0.9982,0.9856,0.9919,0.9812,0.9862);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(190782,191042,0.9961,0.9961,0.9989,0.9989,0.9909,1.0014,0.9854,0.9917,0.976,0.981);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(191043,191720,0.9902,0.9902,0.9931,0.9931,0.9858,0.9963,0.988,0.9944,0.9795,0.9845);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(191721,193833,0.9894,0.9894,0.9922,0.9922,0.9876,0.9982,0.9853,0.9917,0.9791,0.9841);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(193834,194116,0.99,0.99,0.9929,0.9929,0.9865,0.997,0.9842,0.9906,0.9793,0.9843);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(194117,194427,0.9907,0.9907,0.9935,0.9935,0.987,0.9975,0.9855,0.9918,0.9805,0.9855);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(194428,194618,0.9901,0.9901,0.9929,0.9929,0.9867,0.9973,0.9841,0.9905,0.9786,0.9836);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(194619,194789,0.9904,0.9904,0.9932,0.9932,0.9874,0.9979,0.9876,0.994,0.9788,0.9838);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(194790,195111,0.991,0.991,0.9938,0.9938,0.9887,0.9992,0.99,0.9963,0.9803,0.9853);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(195112,195377,0.9912,0.9912,0.994,0.994,0.9871,0.9976,0.9885,0.9948,0.9817,0.9866);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(195378,195398,0.9903,0.9903,0.9931,0.9931,0.9863,0.9968,0.9867,0.993,0.9801,0.9851);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(195399,195657,0.9908,0.9908,0.9936,0.9936,0.9888,0.9993,0.9866,0.9929,0.9806,0.9856);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(195658,195918,0.9914,0.9914,0.9942,0.9942,0.9878,0.9983,0.9898,0.9962,0.9797,0.9846);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(195919,196198,0.9908,0.9908,0.9936,0.9936,0.9874,0.998,0.9878,0.9942,0.9796,0.9846);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(196199,196356,0.9915,0.9915,0.9943,0.9943,0.9877,0.9983,0.988,0.9944,0.9781,0.9831);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(196357,198115,0.991,0.991,0.9938,0.9938,0.9871,0.9977,0.9859,0.9923,0.9791,0.9841);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(198116,198940,0.9906,0.9906,0.9934,0.9934,0.9864,0.997,0.9877,0.994,0.9818,0.9868);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(198941,199317,0.9907,0.9907,0.9936,0.9936,0.9865,0.997,0.9856,0.992,0.9805,0.9855);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(199318,199428,0.9904,0.9904,0.9933,0.9933,0.9868,0.9973,0.9857,0.9921,0.9808,0.9858);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(199429,199697,0.9907,0.9907,0.9935,0.9935,0.987,0.9976,0.9853,0.9916,0.9803,0.9853);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(199698,199832,0.991,0.991,0.9938,0.9938,0.9879,0.9985,0.9845,0.9909,0.9803,0.9853);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(199833,199960,0.9912,0.9912,0.994,0.994,0.9877,0.9983,0.9875,0.9938,0.9823,0.9873);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(199961,200151,0.9913,0.9913,0.9942,0.9942,0.9865,0.997,0.986,0.9924,0.9817,0.9867);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(200152,200490,0.9912,0.9912,0.994,0.994,0.9879,0.9984,0.9878,0.9942,0.9809,0.9858);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(200491,200991,0.9919,0.9919,0.9947,0.9947,0.9866,0.9971,0.9861,0.9924,0.9814,0.9864);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(200992,201201,0.9909,0.9909,0.9937,0.9937,0.9864,0.9969,0.985,0.9914,0.9833,0.9883);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(201202,201624,0.9915,0.9915,0.9943,0.9943,0.9886,0.9992,0.9857,0.992,0.9827,0.9876);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(201625,201707,0.9917,0.9917,0.9945,0.9945,0.9878,0.9983,0.9856,0.9919,0.982,0.987);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(201708,202059,0.9916,0.9916,0.9944,0.9944,0.9869,0.9974,0.9866,0.9929,0.9833,0.9883);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(202060,202204,0.9919,0.9919,0.9947,0.9947,0.9871,0.9977,0.985,0.9914,0.9827,0.9877);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(202205,202332,0.9923,0.9923,0.9951,0.9951,0.9883,0.9989,0.9887,0.995,0.9825,0.9874);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(202333,202972,0.9921,0.9921,0.9949,0.9949,0.9888,0.9994,0.9857,0.992,0.9812,0.9861);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(202973,203002,0.9916,0.9916,0.9944,0.9944,0.9857,0.9962,0.986,0.9924,0.9816,0.9865);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(203003,203852,0.993,0.993,0.9958,0.9958,0.9831,0.9936,0.9705,0.9769,0.9728,0.9778);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(203853,204099,0.9906,0.9906,0.9935,0.9935,0.9886,0.9991,0.9851,0.9914,0.9777,0.9827);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(204100,204562,0.991,0.991,0.9939,0.9939,0.9894,0.9999,0.9868,0.9931,0.9814,0.9864);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(204563,205085,0.991,0.991,0.9938,0.9938,0.9901,1.0006,0.9871,0.9934,0.9789,0.9839);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(205086,205310,0.9911,0.9911,0.9939,0.9939,0.9886,0.9991,0.9868,0.9932,0.9782,0.9832);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(205311,205617,0.9909,0.9909,0.9938,0.9938,0.989,0.9995,0.9852,0.9916,0.9766,0.9816);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(205618,205825,0.9914,0.9914,0.9942,0.9942,0.9877,0.9982,0.9866,0.993,0.9792,0.9842);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(205826,206207,0.9921,0.9921,0.9949,0.9949,0.9908,1.0012,0.9882,0.9945,0.9816,0.9866);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(206208,206389,0.9918,0.9918,0.9946,0.9946,0.989,0.9995,0.9842,0.9906,0.9794,0.9844);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(206390,206483,0.9917,0.9917,0.9945,0.9945,0.989,0.9995,0.9894,0.9957,0.9824,0.9873);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(206484,206597,0.9917,0.9917,0.9945,0.9945,0.9892,0.9997,0.9886,0.9949,0.9791,0.9841);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(206598,206896,0.9913,0.9913,0.9941,0.9941,0.9893,0.9998,0.9846,0.991,0.9789,0.9839);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(206897,207220,0.9923,0.9923,0.9952,0.9952,0.9892,0.9997,0.9872,0.9935,0.9813,0.9863);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(207221,207315,0.9921,0.9921,0.9949,0.9949,0.9901,1.0006,0.987,0.9934,0.9798,0.9848);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(207316,207489,0.9922,0.9922,0.995,0.995,0.9896,1.0001,0.9881,0.9944,0.9791,0.9841);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(207490,207919,0.9922,0.9922,0.9951,0.9951,0.9901,1.0006,0.9866,0.993,0.9769,0.9819);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(207920,208351,0.9919,0.9919,0.9947,0.9947,0.9895,1,0.9867,0.993,0.9813,0.9863);
  photpreselinverteleveto->AddEnCorrPerRun2012HCP(208352,208686,0.9924,0.9924,0.9952,0.9952,0.9906,1.0011,0.9871,0.9934,0.9809,0.9859);
  
  //-----------------------------------------------------------------------------------------------------------------
  photpreselinverteleveto->SetShowerShapeType("2012ShowerShape");
  photpreselinverteleveto->SetDoShowerShapeScaling(kFALSE);  
  //photpreselinverteleveto->SetDoMCErrScaling(kTRUE);
  //photpreselinverteleveto->SetMCErrScale(1.07, 1.045);  
  //photpreselinverteleveto->SetMCErrScale(1, 1);  
  photpreselinverteleveto->SetApplyEleVeto(kFALSE);
  photpreselinverteleveto->SetInvertElectronVeto(kTRUE);
  photpreselinverteleveto->SetJetsName(jetCorr->GetOutputName());  
  //photpreselinverteleveto->SetRescaledBeamspotWidth(5.0);   
  photpreselinverteleveto->SetIsData(isData);
  photpreselinverteleveto->SetApplyLeptonTag(kTRUE);
  photpreselinverteleveto->SetLeptonTagElectronsName("HggLeptonTagElectrons");
  photpreselinverteleveto->SetLeptonTagMuonsName("HggLeptonTagMuons");    
  photpreselinverteleveto->Set2012HCP(kTRUE);  
  photpreselinverteleveto->SetLeadingPtMin(20.);
  photpreselinverteleveto->SetTrailingPtMin(20.);

   PhotonPairSelector         *photpreselnosmear = new PhotonPairSelector("PhotonPairSelectorPreselNoSmear");
   photpreselnosmear->SetOutputName("GoodPhotonsPreselNoSmear");
   photpreselnosmear->SetPhotonSelType("MITPFSelectionNoEcal");
   photpreselnosmear->SetVertexSelType("CiCMVA2012Selection");
   photpreselnosmear->SetUseSingleLegConversions(kFALSE);
   photpreselnosmear->SetIdMVAType("2013FinalIdMVA_8TeV");
   photpreselnosmear->SetShowerShapeType("2012ShowerShape");
   photpreselnosmear->SetDoShowerShapeScaling(kFALSE);
   photpreselnosmear->SetPhotonsFromBranch(kFALSE);
   photpreselnosmear->SetInputPhotonsName(photreg->GetOutputName());
   photpreselnosmear->SetMCSmearFactors2012HCP(0.0075,0.0075,0.0086,0.0086,0.0122,0.0188,0.0163,0.0198,0.0186,0.0192);
   photpreselnosmear->SetMCSmearFactors2012HCPMVA(0.0075,0.0075,0.0086,0.0086,0.0122,0.0188,0.0163,0.0198,0.0186,0.0192);
   photpresel->UseSpecialSmearForDPMVA(true);
   photpreselnosmear->SetJetsName(jetCorr->GetOutputName());
   photpreselnosmear->SetOutputVtxName("OutVtxNoSmear");
   photpreselnosmear->SetLeadingPtMin(30.);
   photpreselnosmear->SetTrailingPtMin(22.);
   //photpreselnosmear->SetRescaledBeamspotWidth(5.0);
   photpreselnosmear->SetIsData(isData);
   photpreselnosmear->SetApplyLeptonTag(kTRUE);
   photpreselnosmear->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   photpreselnosmear->SetLeptonTagMuonsName("HggLeptonTagMuons");
   photpreselnosmear->Set2012HCP(kTRUE);
   photpreselnosmear->DoMCSmear(kFALSE);
   photpreselnosmear->DoMCEneSmear(kFALSE);
   photpreselnosmear->DoEneErrSmear(kTRUE);
   photpreselnosmear->DoDataEneCorr(kFALSE);

   PhotonPairSelector         *photcicnosmear = new PhotonPairSelector("PhotonPairSelectorCiCNoSmear");
   photcicnosmear->SetOutputName("GoodPhotonsCICNoSmear");
   photcicnosmear->SetOutputVtxName("OutVtxCiCNoSmear");        
   photcicnosmear->SetPhotonSelType("CiCPFSelection");
   photcicnosmear->SetVertexSelType("CiCMVA2012Selection");
   photcicnosmear->SetUseSingleLegConversions(kFALSE);
   photcicnosmear->DoMCSmear(kFALSE);
   photcicnosmear->DoDataEneCorr(kFALSE);
   photcicnosmear->SetPhotonsFromBranch(kFALSE);
   photcicnosmear->SetInputPhotonsName(photreg->GetOutputName());
   //  photcicnosmear->SetDoMCR9Scaling(kTRUE);
   //  photcicnosmear->SetMCR9Scale(1.0035, 1.0035);
   photcicnosmear->SetDoShowerShapeScaling(kFALSE); 
   photcicnosmear->SetShowerShapeType("2012ShowerShape");
   //photcicnosmear->SetDoMCErrScaling(kTRUE);
   //photcicnosmear->SetMCErrScale(1.07, 1.045); 
   //photcicnosmear->SetMCErrScale(1, 1); //ming:scale(sigE/E)
   photcicnosmear->SetJetsName(jetCorr->GetOutputName());    
   //photcicnosmear->SetRescaledBeamspotWidth(5.0);
   photcicnosmear->SetIsData(isData);
   photcicnosmear->SetApplyLeptonTag(kTRUE);
   photcicnosmear->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   photcicnosmear->SetLeptonTagMuonsName("HggLeptonTagMuons");  
   photcicnosmear->Set2012HCP(kTRUE);
   
   PhotonPairSelector         *photpreselinvertelevetonosmear = new PhotonPairSelector("PhotonPairSelectorPreselInvertEleVetoNoSmear");
   photpreselinvertelevetonosmear->SetOutputName("GoodPhotonsPreselInvertEleVetoNoSmear");
   photpreselinvertelevetonosmear->SetPhotonSelType("MITPFSelectionNoEcal");
   //photpreselinvertelevetonosmear->SetVertexSelType("CiCMVA2012Selection");//ming change
   photpreselinvertelevetonosmear->SetUseSingleLegConversions(kFALSE);  
   photpreselinvertelevetonosmear->SetIdMVAType("2013FinalIdMVA_8TeV");
   photpreselinvertelevetonosmear->SetShowerShapeType("2012ShowerShape");
   photpreselinvertelevetonosmear->SetDoShowerShapeScaling(kFALSE);  
   photpreselinvertelevetonosmear->SetPhotonsFromBranch(kFALSE);
   photpreselinvertelevetonosmear->SetInputPhotonsName(photreg->GetOutputName());
   photpreselinvertelevetonosmear->SetJetsName(jetCorr->GetOutputName());  
   photpreselinvertelevetonosmear->SetOutputVtxName("OutVtxPreselInvertEleVetoNoSmear");  
   //photpreselinvertelevetonosmear->SetLeadingPtMin(30.);
   photpreselinvertelevetonosmear->SetLeadingPtMin(20.);//ming
   //photpreselinvertelevetonosmear->SetTrailingPtMin(22.);
   photpreselinvertelevetonosmear->SetTrailingPtMin(20.);//ming
   //photpreselinvertelevetonosmear->SetRescaledBeamspotWidth(5.0);      
   photpreselinvertelevetonosmear->SetIsData(isData);  
   photpreselinvertelevetonosmear->SetApplyLeptonTag(kTRUE);
   photpreselinvertelevetonosmear->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   photpreselinvertelevetonosmear->SetLeptonTagMuonsName("HggLeptonTagMuons");
   photpreselinvertelevetonosmear->Set2012HCP(kTRUE); 
   photpreselinvertelevetonosmear->DoMCSmear(kFALSE);
   photpreselinvertelevetonosmear->DoDataEneCorr(kFALSE);
   photpreselinvertelevetonosmear->SetApplyEleVeto(kFALSE);
   photpreselinvertelevetonosmear->SetInvertElectronVeto(kTRUE);
   
   PhotonTreeWriter *phottreecic = new PhotonTreeWriter("PhotonTreeWriterCiC");
   phottreecic->SetPhotonsFromBranch(kFALSE);
   phottreecic->SetInputPhotonsName(photcic->GetOutputName());
   phottreecic->SetEnableJets(kTRUE);
   phottreecic->SetApplyJetId(kTRUE);
   phottreecic->SetPFJetsFromBranch(kFALSE);
   phottreecic->SetPFJetName(jetCorr->GetOutputName());
   phottreecic->SetExcludeDoublePrompt(excludedoubleprompt);
   phottreecic->SetIsData(isData);
   if (is25) phottreecic->SetEnablePFPhotons(kFALSE);
   phottreecic->SetBeamspotWidth(5.0);
   phottreecic->SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
   phottreecic->SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
   phottreecic->SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
   phottreecic->SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
   phottreecic->SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
   phottreecic->SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml"))));
   //phottreecic->SetApplyLeptonTag(kTRUE);
   //phottreecic->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   //phottreecic->SetLeptonTagMuonsName("HggLeptonTagMuons");
   //phottreecic->SetApplyVBFTag(kTRUE);
   phottreecic->SetApplyPFMetCorr(kTRUE);
   phottreecic->SetApplyVHLepTag(kTRUE);
   phottreecic->SetLeptonTagSoftElectronsName("HggLeptonTagSoftElectrons");
   phottreecic->SetLeptonTagSoftMuonsName("HggLeptonTagSoftMuons");
   phottreecic->SetApplyVHHadTag(kTRUE);
   phottreecic->SetApplyVBFTag(kTRUE);
   phottreecic->SetApplyTTHTag(kTRUE);
   
   PhotonTreeWriter *phottreecicnoeleveto = new PhotonTreeWriter("PhotonTreeWriterCiCInvertEleVeto");
   phottreecicnoeleveto->SetPhotonsFromBranch(kFALSE);
   phottreecicnoeleveto->SetInputPhotonsName(photcicnoeleveto->GetOutputName());
   phottreecicnoeleveto->SetEnableJets(kTRUE);
   phottreecicnoeleveto->SetApplyJetId(kTRUE);  
   phottreecicnoeleveto->SetPFJetsFromBranch(kFALSE);
   phottreecicnoeleveto->SetPFJetName(jetCorr->GetOutputName());
   phottreecicnoeleveto->SetApplyElectronVeto(kFALSE);
   phottreecicnoeleveto->SetExcludeDoublePrompt(excludedoubleprompt);
   phottreecicnoeleveto->SetIsData(isData);  
   if (is25) phottreecicnoeleveto->SetEnablePFPhotons(kFALSE);
   phottreecicnoeleveto->SetBeamspotWidth(5.0);
   phottreecicnoeleveto->SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
   phottreecicnoeleveto->SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
   phottreecicnoeleveto->SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
   phottreecicnoeleveto->SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
   phottreecicnoeleveto->SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
   phottreecicnoeleveto->SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml"))));
   //phottreecicnoeleveto->SetApplyLeptonTag(kTRUE);
   //phottreecicnoeleveto->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   //phottreecicnoeleveto->SetLeptonTagMuonsName("HggLeptonTagMuons");
   //phottreecicnoeleveto->SetApplyVBFTag(kTRUE);
   phottreecicnoeleveto->SetApplyPFMetCorr(kTRUE);
   phottreecicnoeleveto->SetApplyVHLepTag(kTRUE);
   phottreecicnoeleveto->SetLeptonTagSoftElectronsName("HggLeptonTagSoftElectrons");
   phottreecicnoeleveto->SetLeptonTagSoftMuonsName("HggLeptonTagSoftMuons");
   phottreecicnoeleveto->SetApplyVHHadTag(kTRUE);
   phottreecicnoeleveto->SetApplyVBFTag(kTRUE);
   phottreecicnoeleveto->SetApplyTTHTag(kTRUE);

   PhotonTreeWriter *phottreepresel = new PhotonTreeWriter("PhotonTreeWriterPresel");
   phottreepresel->SetPhotonsFromBranch(kFALSE);
   phottreepresel->SetInputPhotonsName(photpresel->GetOutputName());
   phottreepresel->SetEnableJets(kTRUE);
   phottreepresel->SetApplyJetId(kTRUE);    
   phottreepresel->SetPFJetsFromBranch(kFALSE);
   phottreepresel->SetPFJetName(jetCorr->GetOutputName());  
   phottreepresel->SetExcludeDoublePrompt(excludedoubleprompt);  
   phottreepresel->SetIsData(isData);  
   if (is25) phottreepresel->SetEnablePFPhotons(kFALSE);
   phottreepresel->SetBeamspotWidth(5.0);
   phottreepresel->SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
   phottreepresel->SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
   phottreepresel->SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
   phottreepresel->SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
   phottreepresel->SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
   phottreepresel->SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml"))));
   phottreepresel->SetDoSynching(kTRUE);
   //
   //phottreepresel->SetApplyLeptonTag(kTRUE);
   //phottreepresel->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   //phottreepresel->SetLeptonTagMuonsName("HggLeptonTagMuons");
   //phottreepresel->SetApplyVBFTag(kTRUE);
   phottreepresel->SetApplyPFMetCorr(kTRUE);
   phottreepresel->SetApplyVHLepTag(kTRUE);
   phottreepresel->SetLeptonTagSoftElectronsName("HggLeptonTagSoftElectrons");
   phottreepresel->SetLeptonTagSoftMuonsName("HggLeptonTagSoftMuons");
   phottreepresel->SetApplyVHHadTag(kTRUE);
   phottreepresel->SetApplyVBFTag(kTRUE);
   phottreepresel->SetApplyTTHTag(kTRUE);
 
   PhotonTreeWriter *phottreepreselinverteleveto = new PhotonTreeWriter("PhotonTreeWriterPreselInvertEleVeto");
   phottreepreselinverteleveto->SetPhotonsFromBranch(kFALSE);
   phottreepreselinverteleveto->SetInputPhotonsName(photpreselinverteleveto->GetOutputName());
   phottreepreselinverteleveto->SetEnableJets(kTRUE);
   phottreepreselinverteleveto->SetApplyJetId(kTRUE);    
   phottreepreselinverteleveto->SetPFJetsFromBranch(kFALSE);
   phottreepreselinverteleveto->SetPFJetName(jetCorr->GetOutputName());  
   phottreepreselinverteleveto->SetApplyElectronVeto(kFALSE);  
   phottreepreselinverteleveto->SetExcludeDoublePrompt(excludedoubleprompt);    
   phottreepreselinverteleveto->SetIsData(isData); 
   if (is25) phottreepreselinverteleveto->SetEnablePFPhotons(kFALSE);  
   phottreepreselinverteleveto->SetBeamspotWidth(5.0);
  phottreepreselinverteleveto->SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
  phottreepreselinverteleveto->SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
  phottreepreselinverteleveto->SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
  phottreepreselinverteleveto->SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
  phottreepreselinverteleveto->SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
  phottreepreselinverteleveto->SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml"))));
   // phottreepreselinverteleveto->SetApplyLeptonTag(kTRUE);
   // phottreepreselinverteleveto->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   // phottreepreselinverteleveto->SetLeptonTagMuonsName("HggLeptonTagMuons");
   // phottreepreselinverteleveto->SetApplyVBFTag(kTRUE);
    phottreepreselinverteleveto->SetApplyPFMetCorr(kTRUE);
    phottreepreselinverteleveto->SetApplyVHLepTag(kTRUE);
    phottreepreselinverteleveto->SetLeptonTagSoftElectronsName("HggLeptonTagSoftElectrons");
    phottreepreselinverteleveto->SetLeptonTagSoftMuonsName("HggLeptonTagSoftMuons");
    phottreepreselinverteleveto->SetApplyVHHadTag(kTRUE);
    phottreepreselinverteleveto->SetApplyVBFTag(kTRUE);
    phottreepreselinverteleveto->SetApplyTTHTag(kTRUE);

   PhotonTreeWriter *phottreepreselnosmear = new PhotonTreeWriter("PhotonTreeWriterPreselNoSmear");
   phottreepreselnosmear->SetPhotonsFromBranch(kFALSE);
   phottreepreselnosmear->SetInputPhotonsName(photpreselnosmear->GetOutputName());
   phottreepreselnosmear->SetEnableJets(kTRUE);
   phottreepreselnosmear->SetApplyJetId(kTRUE);  
   phottreepreselnosmear->SetPFJetsFromBranch(kFALSE);
   phottreepreselnosmear->SetPFJetName(jetCorr->GetOutputName());  
   phottreepreselnosmear->SetExcludeDoublePrompt(excludedoubleprompt);  
   phottreepreselnosmear->SetIsData(isData);    
   if (is25) phottreepreselnosmear->SetEnablePFPhotons(kFALSE); 
   phottreepreselnosmear->SetBeamspotWidth(5.0); 
   phottreepreselnosmear->SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
   phottreepreselnosmear->SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
   phottreepreselnosmear->SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
   phottreepreselnosmear->SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
   phottreepreselnosmear->SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
   phottreepreselnosmear->SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml"))));
   phottreepreselnosmear->SetDoSynching(kTRUE);
   // phottreepreselnosmear->SetApplyLeptonTag(kTRUE);
   // phottreepreselnosmear->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   // phottreepreselnosmear->SetLeptonTagMuonsName("HggLeptonTagMuons");
   // phottreepreselnosmear->SetApplyVBFTag(kTRUE);
    phottreepreselnosmear->SetApplyPFMetCorr(kTRUE);
    phottreepreselnosmear->SetApplyVHLepTag(kTRUE);
    phottreepreselnosmear->SetLeptonTagSoftElectronsName("HggLeptonTagSoftElectrons");
    phottreepreselnosmear->SetLeptonTagSoftMuonsName("HggLeptonTagSoftMuons");
    phottreepreselnosmear->SetApplyVHHadTag(kTRUE);
    phottreepreselnosmear->SetApplyVBFTag(kTRUE);
    phottreepreselnosmear->SetApplyTTHTag(kTRUE);

   PhotonTreeWriter *phottreecicnosmear = new PhotonTreeWriter("PhotonTreeWriterCiCNoSmear");
   phottreecicnosmear->SetPhotonsFromBranch(kFALSE);
   phottreecicnosmear->SetInputPhotonsName(photcicnosmear->GetOutputName());
   phottreecicnosmear->SetEnableJets(kTRUE);
   phottreecicnosmear->SetApplyJetId(kTRUE);
   phottreecicnosmear->SetPFJetsFromBranch(kFALSE);
   phottreecicnosmear->SetPFJetName(jetCorr->GetOutputName());
   phottreecicnosmear->SetExcludeDoublePrompt(excludedoubleprompt);
   phottreecicnosmear->SetIsData(isData);
   if (is25) phottreecicnosmear->SetEnablePFPhotons(kFALSE);
   phottreecicnosmear->SetBeamspotWidth(5.0);
   phottreecicnosmear->SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
   phottreecicnosmear->SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
   phottreecicnosmear->SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
   phottreecicnosmear->SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
   phottreecicnosmear->SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
   phottreecicnosmear->SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml"))));
   //  phottreecicnosmear->SetApplyLeptonTag(kTRUE);
   //  phottreecicnosmear->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   //  phottreecicnosmear->SetLeptonTagMuonsName("HggLeptonTagMuons");
   //  phottreecicnosmear->SetApplyVBFTag(kTRUE);
    phottreecicnosmear->SetApplyPFMetCorr(kTRUE);
    phottreecicnosmear->SetApplyVHLepTag(kTRUE);
    phottreecicnosmear->SetLeptonTagSoftElectronsName("HggLeptonTagSoftElectrons");
    phottreecicnosmear->SetLeptonTagSoftMuonsName("HggLeptonTagSoftMuons");
    phottreecicnosmear->SetApplyVHHadTag(kTRUE);
    phottreecicnosmear->SetApplyVBFTag(kTRUE);
    phottreecicnosmear->SetApplyTTHTag(kTRUE);
   
   PhotonTreeWriter *phottreepreselinvertelevetonosmear = new PhotonTreeWriter("PhotonTreeWriterPreselInvertEleVetoNoSmear");
   phottreepreselinvertelevetonosmear->SetPhotonsFromBranch(kFALSE);
   phottreepreselinvertelevetonosmear->SetInputPhotonsName(photpreselinvertelevetonosmear->GetOutputName());
   phottreepreselinvertelevetonosmear->SetEnableJets(kTRUE);
   phottreepreselinvertelevetonosmear->SetApplyJetId(kTRUE);  
   phottreepreselinvertelevetonosmear->SetPFJetsFromBranch(kFALSE);
   phottreepreselinvertelevetonosmear->SetPFJetName(jetCorr->GetOutputName());  
   phottreepreselinvertelevetonosmear->SetExcludeDoublePrompt(excludedoubleprompt);  
   phottreepreselinvertelevetonosmear->SetIsData(isData);    
   if (is25) phottreepreselinvertelevetonosmear->SetEnablePFPhotons(kFALSE); 
   phottreepreselinvertelevetonosmear->SetBeamspotWidth(5.0); 
   phottreepreselinvertelevetonosmear->SetElectronMVAWeightsSubdet0Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml"))));
   phottreepreselinvertelevetonosmear->SetElectronMVAWeightsSubdet1Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml"))));  
   phottreepreselinvertelevetonosmear->SetElectronMVAWeightsSubdet2Pt10To20(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml"))));  
   phottreepreselinvertelevetonosmear->SetElectronMVAWeightsSubdet0Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml")))); 
   phottreepreselinvertelevetonosmear->SetElectronMVAWeightsSubdet1Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml")))); 
   phottreepreselinvertelevetonosmear->SetElectronMVAWeightsSubdet2Pt20ToInf(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml"))));
   phottreepreselinvertelevetonosmear->SetApplyElectronVeto(kFALSE); 
   //   phottreepreselinvertelevetonosmear->SetApplyLeptonTag(kTRUE);
   //   phottreepreselinvertelevetonosmear->SetLeptonTagElectronsName("HggLeptonTagElectrons");
   //   phottreepreselinvertelevetonosmear->SetLeptonTagMuonsName("HggLeptonTagMuons");
   //   phottreepreselinvertelevetonosmear->SetApplyVBFTag(kTRUE);
    phottreepreselinvertelevetonosmear->SetApplyPFMetCorr(kTRUE);
    phottreepreselinvertelevetonosmear->SetApplyVHLepTag(kTRUE);
    phottreepreselinvertelevetonosmear->SetLeptonTagSoftElectronsName("HggLeptonTagSoftElectrons");
    phottreepreselinvertelevetonosmear->SetLeptonTagSoftMuonsName("HggLeptonTagSoftMuons");
    phottreepreselinvertelevetonosmear->SetApplyVHHadTag(kTRUE);
    phottreepreselinvertelevetonosmear->SetApplyVBFTag(kTRUE);
    phottreepreselinvertelevetonosmear->SetApplyTTHTag(kTRUE);
   
   PhotonIDMod         *photidpresel = new PhotonIDMod("PhotonIDModPresel");
   photidpresel->SetPtMin(25.0);
   photidpresel->SetOutputName("GoodPhotonsPreselid");
   photidpresel->SetIDType("MITPFSelection");
   photidpresel->SetApplyElectronVeto(kTRUE);
   photidpresel->SetIsData(isData);
   
  PhotonIDMod         *photidpreselinvert = new PhotonIDMod("PhotonIDModPreselInvert");
  photidpreselinvert->SetPtMin(25.0);
  photidpreselinvert->SetOutputName("GoodPhotonsPreselidInvert");
  photidpreselinvert->SetIDType("MITPFSelection");
  photidpreselinvert->SetApplyElectronVeto(kFALSE);
  photidpreselinvert->SetInvertElectronVeto(kTRUE);
  photidpreselinvert->SetIsData(isData);  
  
  PhotonTreeWriter *phottreesingle = new PhotonTreeWriter("PhotonTreeWriterSingle");
  phottreesingle->SetWriteDiphotonTree(kFALSE);
  phottreesingle->SetPhotonsFromBranch(kFALSE);
  phottreesingle->SetInputPhotonsName(photidpresel->GetOutputName());  
  phottreesingle->SetEnableJets(kTRUE);
  phottreesingle->SetPFJetsFromBranch(kFALSE);
  phottreesingle->SetPFJetName(jetCorr->GetOutputName());      
  phottreesingle->SetBeamspotWidth(5.0);
  phottreesingle->SetIsData(isData);
  if (is25) phottreesingle->SetEnablePFPhotons(kFALSE);  
  
  PhotonTreeWriter *phottreesingleinvert = new PhotonTreeWriter("PhotonTreeWriterSingleInvert");
  phottreesingleinvert->SetWriteDiphotonTree(kTRUE);
  phottreesingleinvert->SetPhotonsFromBranch(kFALSE);
  phottreesingleinvert->SetInputPhotonsName(photidpreselinvert->GetOutputName());  
  phottreesingleinvert->SetEnableJets(kTRUE);
  phottreesingleinvert->SetPFJetsFromBranch(kFALSE);
  phottreesingleinvert->SetPFJetName(pubJetOpen->GetOutputName());
  phottreesingleinvert->SetApplyElectronVeto(kFALSE);
  phottreesingleinvert->SetBeamspotWidth(5.0);
  phottreesingleinvert->SetIsData(isData);  
  if (is25) phottreesingleinvert->SetEnablePFPhotons(kFALSE);    
  

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel      ->Add(mcselmod);
    
  if (TString(dataset).Contains("-h")) {    
    mcselmod        ->Add(sysMod);
  }
  
  // high level trigger is always first
  //mcselmod         ->Add(hltModES);
  //mcselmod         ->Add(hltModE);
  
  if (!TString(dataset).Contains("meridiani")) {      
    mcselmod         ->Add(hltModP);
    hltModP         ->Add(goodPVFilterMod);
  }
  else {
    mcselmod->Add(goodPVFilterMod);
  }
  //hltModE         ->Add(goodPVFilterModE);
  //hltModES        ->Add(goodPVFilterModES);
  
  //goodPVFilterMod ->Add(muonId);
  goodPVFilterMod->Add(photreg);
  photreg->Add(pubJet);
  pubJet->Add(jetCorr);
  
  // simple object id modules
  //goodPVFilterModE -> Add(pubJetOpen);
  //pubJetOpen       -> Add(photidpreselinvert);
  //goodPVFilterModES ->Add(elecIdS);
  //jetCorr          ->Add(SepPUMod); 
  //SepPUMod         ->Add(eleIdMod);
  //eleIdMod         ->Add(muonIdMod);

  jetCorr          ->Add(SepPUMod); 
  SepPUMod         ->Add(eleIdMod);
  eleIdMod         ->Add(softEleIdMod);
  softEleIdMod     ->Add(softMuonIdMod);
  softMuonIdMod    ->Add(muonIdMod);
   
  muonIdMod          ->Add(photcic);
  muonIdMod          ->Add(photcicnoeleveto);  
  muonIdMod          ->Add(photpresel);  
  muonIdMod          ->Add(photpreselinverteleveto);  
  muonIdMod          ->Add(photpreselnosmear);  
  muonIdMod          ->Add(photpreselinvertelevetonosmear);  
  muonIdMod          ->Add(photcicnosmear);  

  photcic         ->Add(phottreecic);
  photcicnoeleveto       ->Add(phottreecicnoeleveto);
  photpresel    ->Add(phottreepresel);
  photpreselinverteleveto    ->Add(phottreepreselinverteleveto);
  photpreselnosmear    ->Add(phottreepreselnosmear);
  photpreselinvertelevetonosmear    ->Add(phottreepreselinvertelevetonosmear);
  photcicnosmear    ->Add(phottreecicnosmear);

  //jetCorr          ->Add(photidpresel);
  //photidpresel       ->Add(phottreesingle);
  
  //photidpreselinvert       ->Add(phottreesingleinvert);  
  //elecIdS->Add(phottreeES);
  

  //TFile::SetOpenTimeout(0);
  //TFile::SetCacheFileDir("./rootfilecache",kTRUE,kTRUE);
  TFile::SetReadaheadSize(128*1024*1024);
  
  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kTRUE);
  ana->SetSuperModule(runLumiSel);
  ana->SetPrintScale(100);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);

  if (TString(dataset).Contains("meridiani")) {    
    ana->SetUseHLT(kFALSE);
  }

 
  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  Catalog *c = new Catalog(catalogDir);
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  TString bookstr = book;
  Bool_t caching = kTRUE;
  //Bool_t caching = kFALSE;
  //if (TString(dataset).Contains("s11-h")) bookstr.ReplaceAll("local","t2mit");
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(bookstr,dataset,fileset,caching);
  else 
    d = c->FindDataset(bookstr,skimdataset.Data(),fileset,caching);

  ana->AddDataset(d);
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12a-pho-j13-v1/385B77DE-58D0-E111-B925-001E67396928.root");//79737729
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12b-dph-j13-v1/FE308E0A-23D2-E111-8B2D-00266CFAE7D0.root");//871378986
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/54E239DC-2725-E211-A52C-003048D373F6.root");//528937923
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/EC80B121-7E30-E211-BFFF-5404A640A63D.root");//261543921
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/66CE9443-472C-E211-B08F-001D09F2437B.root");//483561562
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/F03B81B2-0B29-E211-B28F-BCAEC518FF67.root");
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/A6EC91DD-CB39-E211-8897-BCAEC53296F8.root");
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12c-dph-pr-v2/64F7EAB0-CBEA-E111-B8AD-003048F117B6.root");
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12c-dph-pr-v2/0615B16F-ECDA-E111-A266-003048F117B6.root"); 
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12a-pho-j13-v1/385B77DE-58D0-E111-B925-001E67396928.root"); 
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12c-dph-pr-v2/068B257D-81D0-E111-A835-5404A63886B4.root");//369441614 
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/5443DBCB-562C-E211-8042-00237DDBE0E2.root");//333643114
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/B2698ACF-843A-E211-8D9E-001D09F244DE.root");//89022540
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12b-dph-j13-v1/1E050D78-6DD4-E111-81D2-00266CFAE8D0.root");//223740859
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12b-dph-j13-v1/2CD64E47-79D3-E111-90C5-00A0D1EE8E60.root");//541603559
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12b-dph-j13-v1/147EF338-B2D3-E111-B1A6-00266CF256CC.root");//8983064 
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12d-dph-pr-v1/8843BD9F-DC3E-E211-9E2F-003048678110.root");//876316897 
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/r12c-dph-pr-v2/CAA12077-F9F6-E111-8DF2-001D09F28D54.root");
 
  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  ana->SetOutputName(rootFile.Data());
  //ana->SetCacheSize(64*1024*1024);
  ana->SetCacheSize(0);
  
  //------------------------------------------------------------------------------------------------
  // Say what we are doing
  //------------------------------------------------------------------------------------------------
  printf("\n==== PARAMETER SUMMARY FOR THIS JOB ====\n");
  printf("\n JSON file: %s\n  and overlap cut: %f (%s)\n",jsonFile.Data(),overlapCut,overlap);
  printf("\n Rely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n",book,dataset,skim,fileset);
  printf("\n Root output: %s\n\n",rootFile.Data());  
  printf("\n========================================\n");

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());

  return;
}
