//--------------------------------------------------------------------------------------------------
// $Id: LeptonPairPhotonTreeWriter.h,v 1.0 2012/06/23 21:25:01 auhess ksingh
//
// LeptonPairPhotonTreeWriter
//
// Authors: A. Hess & K. Singh
//--------------------------------------------------------------------------------------------------
#include "MitPhysics/Utils/interface/ZGTools.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitPhysics/Utils/interface/MVATools.h"
#include "MitPhysics/Mods/interface/LeptonPairPhotonTreeWriter.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/StableParticle.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/PFMetCorrectionTools.h"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/Muon.h"
#include "TDataMember.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <TNtuple.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <fstream>


using namespace mithep;
//mithep::ElectronIDMVA *eleIDMVA; // The electron MVA 
//MVATools fTool; // The electron MVA tools
ClassImp(mithep::LeptonPairPhotonTreeWriter)
ClassImp(mithep::LeptonPairPhotonEvent)

//
// correction header is passed in as a define so as not to make 
// others checkout Si's UserCode. Set the env var using this syntax : 
//     export ELECTRON_CORRECTIONS_HEADER='\"UserCode/sixie/HiggsAna/Utils/LeptonScaleCorrections.hh\"'
// proper quoting is important ...
#ifdef ELECTRON_CORRECTIONS_HEADER
#include ELECTRON_CORRECTIONS_HEADER
#endif


#ifdef PHOSPHOR_CORRECTIONS_SRC
//ClassImp(zgamma::PhosphorCorrectionFunctor)
#include PHOSPHOR_CORRECTIONS_SRC
#endif

//--------------------------------------------------------------------------------------------------
  LeptonPairPhotonTreeWriter::LeptonPairPhotonTreeWriter(const char *name, const char *title) : 
    BaseMod                 (name,title),
    // define all the Branches to load
    fPhotonBranchName       (Names::gkPhotonBrn),
    fGoodElectronName       (Names::gkElectronBrn),  
    fGoodMuonName	    (Names::gkMuonBrn),

    // KH for sync
    //    fPVName                 (Names::gkPVBeamSpotBrn), 
    fPVName                 (Names::gkPVBrn), 

    fPFCandName		    (Names::gkPFCandidatesBrn),
    fTrackName              (Names::gkTrackBrn),
    fPileUpDenName          (Names::gkPileupEnergyDensityBrn),
    fPileUpName             (Names::gkPileupInfoBrn),
    fMCParticleName         (Names::gkMCPartBrn),
    fConversionName         (Names::gkMvfConversionBrn),  
    fBeamSpotName           (Names::gkBeamSpotBrn),
    fIsData                 (false),
    YEAR		    (2011),
//    fDoBlinding             (false),
    fPhotonsFromBranch      (kTRUE),  
    fPVFromBranch           (kTRUE),
    fGoodElectronsFromBranch(kTRUE),
    fGoodMuonsFromBranch    (kTRUE),
    fPFCandsFromBranch      (kTRUE),
    fTracksFromBranch       (kTRUE),
    fPileUpDenFromBranch    (kTRUE),
    fApplyElectronVeto      (kTRUE),

    // ----------------------------------------
    // collections....
    fPhotons                (NULL),
    fGoodElectrons          (NULL),
    fPV                     (NULL),
    fGoodMuons		    (NULL),
    fPFCands		    (NULL),
    fTracks		    (NULL),
    fPileUpDen              (NULL),
    fPileUp                 (0),
    fMCParticles            (0),
    fConversions	    (0),
    fBeamSpot               (0),
    fTupleName              ("h2LepPhotonTree"),

    verbose                 (false),
    _do_ElectronChannel     (true),
    _do_MuonChannel         (true),
    phosphorDataFile        ("PHOSPHOR_NUMBERS_EXPFIT.txt")
  
    //Photon MVA Variables
//    fVariableType_2011             (10),
//    fEndcapWeights_2011            (gSystem->Getenv("CMSSW_BASE")+
//				    TString("/src/MitPhysics/data/TMVAClassificationPhotonID_")+
//				    TString("Endcap_PassPreSel_Variable_10_BDTnCuts2000_BDT.")+
//				    TString("weights.xml")),
//    fBarrelWeights_2011            (gSystem->Getenv("CMSSW_BASE")+
//				    TString("/src/MitPhysics/data/TMVAClassificationPhotonID_")+
//				    TString("Barrel_PassPreSel_Variable_10_BDTnCuts2000_BDT.")+
//				    TString("weights.xml"))
    //2012 Photon MVA Variables not currently used
   // fVariableType_2012_globe       (1201),
   // fEndcapWeights_2012_globe      (gSystem->Getenv("CMSSW_BASE")+
//				    TString("/src/MitPhysics/data/")+
//				    TString("TMVA_EEpf_BDT_globe.")+
//				    TString("weights.xml")),
   // fBarrelWeights_2012_globe      (gSystem->Getenv("CMSSW_BASE")+
//				    TString("/src/MitPhysics/data/")+
//				    TString("TMVA_EBpf_BDT_globe.")+
//				    TString("weights.xml"))

{
  // Constructor
  rmcor = new rochcor();
#ifdef PHOSPHOR_CORRECTIONS_HEADER
  phosphor = new PhosphorCorrectionFunctor(phosphorDataFile.Data(), true);
#endif
}

LeptonPairPhotonTreeWriter::~LeptonPairPhotonTreeWriter()
{
  // Destructor
}

//--------------------------------------------------------------------------------------------------
void LeptonPairPhotonTreeWriter::fillEle1Variables(const mithep::Electron * ele1)  
//--------------------------------------------------------------------------------------------------
{
  fLeptonPairPhotonEvent->ele1charge = ele1->Charge();
  fLeptonPairPhotonEvent->ele1energy = ele1->E();
  fLeptonPairPhotonEvent->ele1px = ele1->Px();
  fLeptonPairPhotonEvent->ele1py = ele1->Py();
  fLeptonPairPhotonEvent->ele1pz = ele1->Pz();
  fLeptonPairPhotonEvent->ele1pt = ele1->Pt();
  fLeptonPairPhotonEvent->ele1eta = ele1->Eta();
  fLeptonPairPhotonEvent->ele1mass = ele1->Mass();
  fLeptonPairPhotonEvent->ele1phi = ele1->Phi();
   								    

}

//--------------------------------------------------------------------------------------------------
void LeptonPairPhotonTreeWriter::fillEle2Variables(const mithep::Electron * ele2)  
//--------------------------------------------------------------------------------------------------
{
  fLeptonPairPhotonEvent->ele2charge = ele2->Charge();
  fLeptonPairPhotonEvent->ele2energy = ele2->E();
  fLeptonPairPhotonEvent->ele2px = ele2->Px();
  fLeptonPairPhotonEvent->ele2py = ele2->Py();
  fLeptonPairPhotonEvent->ele2pz = ele2->Pz();
  fLeptonPairPhotonEvent->ele2pt = ele2->Pt();
  fLeptonPairPhotonEvent->ele2eta = ele2->Eta();
  fLeptonPairPhotonEvent->ele2mass = ele2->Mass();
  fLeptonPairPhotonEvent->ele2phi = ele2->Phi();
}

//--------------------------------------------------------------------------------------------------
void LeptonPairPhotonTreeWriter::resetTreeVariables()  
//--------------------------------------------------------------------------------------------------
{
  fLeptonPairPhotonEvent->electronZmass = -99.;
  fLeptonPairPhotonEvent->mllg = -99.;
  fLeptonPairPhotonEvent->mllgCorr = -99.;
  fLeptonPairPhotonEvent->ele1MVA = -99.;
  fLeptonPairPhotonEvent->ele2MVA = -99.;
  
  fLeptonPairPhotonEvent->ele1charge = -99.;
  fLeptonPairPhotonEvent->ele1energy = -99.;
  fLeptonPairPhotonEvent->ele1px = -99.;
  fLeptonPairPhotonEvent->ele1py = -99.;
  fLeptonPairPhotonEvent->ele1pz = -99.;
  fLeptonPairPhotonEvent->ele1pt = -99.;
  fLeptonPairPhotonEvent->ele1eta = -99.;
  fLeptonPairPhotonEvent->ele1mass = -99.;
  fLeptonPairPhotonEvent->ele1phi = -99.; 
  fLeptonPairPhotonEvent->ele1RegressionEnergyV0 = -99.;
  fLeptonPairPhotonEvent->ele1RegressionEnergyV1 = -99.;
  fLeptonPairPhotonEvent->ele1RegressionEnergyV2 = -99.;
  fLeptonPairPhotonEvent->ele1RegressionEnergyErrorV0 = -99.;
  fLeptonPairPhotonEvent->ele1RegressionEnergyErrorV1 = -99.;
  fLeptonPairPhotonEvent->ele1RegressionEnergyErrorV2 = -99.;
  fLeptonPairPhotonEvent->ele1energyCorr = -99.;
  fLeptonPairPhotonEvent->ele1pxCorr = -99.;
  fLeptonPairPhotonEvent->ele1pyCorr = -99.;
  fLeptonPairPhotonEvent->ele1pzCorr = -99.;
  fLeptonPairPhotonEvent->ele1ptCorr = -99.;

  fLeptonPairPhotonEvent->ele2charge = -99.;
  fLeptonPairPhotonEvent->ele2energy = -99.;
  fLeptonPairPhotonEvent->ele2px = -99.;
  fLeptonPairPhotonEvent->ele2py = -99.;
  fLeptonPairPhotonEvent->ele2pz = -99.;
  fLeptonPairPhotonEvent->ele2pt = -99.;
  fLeptonPairPhotonEvent->ele2eta = -99.;
  fLeptonPairPhotonEvent->ele2mass = -99.;
  fLeptonPairPhotonEvent->ele2phi = -99.;
  fLeptonPairPhotonEvent->ele2RegressionEnergyV0 = -99.;
  fLeptonPairPhotonEvent->ele2RegressionEnergyV1 = -99.;
  fLeptonPairPhotonEvent->ele2RegressionEnergyV2 = -99.;
  fLeptonPairPhotonEvent->ele2RegressionEnergyErrorV0 = -99.;
  fLeptonPairPhotonEvent->ele2RegressionEnergyErrorV1 = -99.;
  fLeptonPairPhotonEvent->ele2RegressionEnergyErrorV2 = -99.;
  fLeptonPairPhotonEvent->ele2energyCorr = -99.;
  fLeptonPairPhotonEvent->ele2pxCorr = -99.;
  fLeptonPairPhotonEvent->ele2pyCorr = -99.;
  fLeptonPairPhotonEvent->ele2pzCorr = -99.;
  fLeptonPairPhotonEvent->ele2ptCorr = -99.;

  fLeptonPairPhotonEvent->ele1dEtaIn = -99.;
  fLeptonPairPhotonEvent->ele1dPhiIn = -99.;
  fLeptonPairPhotonEvent->ele1sigmaIEtaIEta = -99.;
  fLeptonPairPhotonEvent->ele1HadOverEm = -99.;
  fLeptonPairPhotonEvent->ele1D0 = -99.;
  fLeptonPairPhotonEvent->ele1DZ = -99.;
  fLeptonPairPhotonEvent->ele1OneOverEMinusOneOverP = -99.;
  fLeptonPairPhotonEvent->ele1PFIsoOverPt = -99.;
  fLeptonPairPhotonEvent->ele1Conversion = kFALSE;
  fLeptonPairPhotonEvent->ele1missinghits = -99.;

  fLeptonPairPhotonEvent->ele2dEtaIn = -99.;
  fLeptonPairPhotonEvent->ele2dPhiIn = -99.;
  fLeptonPairPhotonEvent->ele2sigmaIEtaIEta = -99.;
  fLeptonPairPhotonEvent->ele2HadOverEm = -99.;
  fLeptonPairPhotonEvent->ele2D0 = -99.;
  fLeptonPairPhotonEvent->ele2DZ = -99.;
  fLeptonPairPhotonEvent->ele2OneOverEMinusOneOverP = -99.;
  fLeptonPairPhotonEvent->ele2PFIsoOverPt = -99.;
  fLeptonPairPhotonEvent->ele2Conversion = kFALSE;
  fLeptonPairPhotonEvent->ele2missinghits = -99.;

  fLeptonPairPhotonEvent->muonZgVeto = -99.;
  fLeptonPairPhotonEvent->m1E = -99.;
  fLeptonPairPhotonEvent->m1Pt = -99.;
  fLeptonPairPhotonEvent->m1Mass = -99.;
  fLeptonPairPhotonEvent->m1Px = -99.; 
  fLeptonPairPhotonEvent->m1Py = -99.;
  fLeptonPairPhotonEvent->m1Pz = -99.;
  fLeptonPairPhotonEvent->m1Eta = -99.;
  fLeptonPairPhotonEvent->m1Phi = -99.;
  fLeptonPairPhotonEvent->m1Charge = -99.;
  fLeptonPairPhotonEvent->m1PtErr = -99.;
  fLeptonPairPhotonEvent->m1ECorr = -99.;
  fLeptonPairPhotonEvent->m1PtCorr = -99.;
  fLeptonPairPhotonEvent->m1PxCorr = -99.; 
  fLeptonPairPhotonEvent->m1PyCorr = -99.;
  fLeptonPairPhotonEvent->m1PzCorr = -99.;

  fLeptonPairPhotonEvent->m2E = -99.;
  fLeptonPairPhotonEvent->m2Pt = -99.;
  fLeptonPairPhotonEvent->m2Mass = -99.;
  fLeptonPairPhotonEvent->m2Px = -99.;
  fLeptonPairPhotonEvent->m2Py = -99.;
  fLeptonPairPhotonEvent->m2Pz = -99.;
  fLeptonPairPhotonEvent->m2Eta = -99.;
  fLeptonPairPhotonEvent->m2Phi = -99.;
  fLeptonPairPhotonEvent->m2Charge = -99.;
  fLeptonPairPhotonEvent->m2PtErr = -99.;
  fLeptonPairPhotonEvent->m2ECorr = -99.;
  fLeptonPairPhotonEvent->m2PtCorr = -99.;
  fLeptonPairPhotonEvent->m2PxCorr = -99.; 
  fLeptonPairPhotonEvent->m2PyCorr = -99.;
  fLeptonPairPhotonEvent->m2PzCorr = -99.;

  fLeptonPairPhotonEvent->photonidmva = -99.;
  fLeptonPairPhotonEvent->photonenergy = -99.;
  fLeptonPairPhotonEvent->photonpx = -99.;
  fLeptonPairPhotonEvent->photonpy = -99.;
  fLeptonPairPhotonEvent->photonpz = -99.;
  fLeptonPairPhotonEvent->photonpt = -99.;
  fLeptonPairPhotonEvent->photoneta = -99.;
  fLeptonPairPhotonEvent->photonmass = -99.;
  fLeptonPairPhotonEvent->ele2phi = -99.;
  fLeptonPairPhotonEvent->photonr9 = -99.;
  fLeptonPairPhotonEvent->NPu = -99.;
  fLeptonPairPhotonEvent->NPuMinus = -99.;
  fLeptonPairPhotonEvent->NPuPlus = -99.;
  fLeptonPairPhotonEvent->photonmatchmc = -1;
  fLeptonPairPhotonEvent->photonenergyerror = -99.;

  fLeptonPairPhotonEvent->chargediso_ele1 = -99.;
  fLeptonPairPhotonEvent->gammaiso_ele1 = -99.;
  fLeptonPairPhotonEvent->neutraliso_ele1 = -99.; 
  fLeptonPairPhotonEvent->rho_ele1 = -99.;
  fLeptonPairPhotonEvent->effectivearea_ele1 = -99.;
  fLeptonPairPhotonEvent->chargediso_ele2 = -99.;
  fLeptonPairPhotonEvent->gammaiso_ele2 = -99.;
  fLeptonPairPhotonEvent->neutraliso_ele2 = -99.;
  fLeptonPairPhotonEvent->rho_ele2 = -99.;
  fLeptonPairPhotonEvent->effectivearea_ele2 = -99.;
  fLeptonPairPhotonEvent->muonZgVeto = kFALSE;
  fLeptonPairPhotonEvent->muonZmass = -99.;
  fLeptonPairPhotonEvent->costheta_lm_electrons = -99.;
  fLeptonPairPhotonEvent->costheta_lp_electrons = -99.;
  fLeptonPairPhotonEvent->phi_electrons = -99.;
  fLeptonPairPhotonEvent->cosTheta_electrons = -99.;
  fLeptonPairPhotonEvent->cosThetaG_electrons = -99.;
  fLeptonPairPhotonEvent->costheta_lm_muons = -99.;
  fLeptonPairPhotonEvent->costheta_lp_muons = -99.;
  fLeptonPairPhotonEvent->phi_muons = -99.;
  fLeptonPairPhotonEvent->cosTheta_muons = -99.;
  fLeptonPairPhotonEvent->cosThetaG_muons = -99.;
}

//--------------------------------------------------------------------------------------------------
void LeptonPairPhotonTreeWriter::regressEle1(const mithep::Electron               *ele1,
					     const mithep::PileupEnergyDensityCol *fPileUpDen,
					     const mithep::VertexCol              *fPV   )  
//--------------------------------------------------------------------------------------------------
{
  double ele1spp = ((!isnan(float(ele1->SCluster()->Seed()->CoviPhiiPhi()))) ? sqrt(ele1->SCluster()->Seed()->CoviPhiiPhi()) : 0.0);
  double ele1sep;
  if (ele1->CoviEtaiEta()*ele1spp > 0) {
    ele1sep = ele1->SCluster()->Seed()->CoviEtaiPhi()/(ele1->CoviEtaiEta()*ele1spp);
  } else if (ele1->SCluster()->Seed()->CoviEtaiPhi()) {
    ele1sep = 1.0; 
  } else {
    ele1sep = -1.0; 
  }
  
  fLeptonPairPhotonEvent->ele1RegressionEnergyV0 = 
    eleRegressionEvaluator_V0->regressionValueNoTrkVar(
						       ele1->SCluster()->RawEnergy(),
						       ele1->SCluster()->Eta(),
						       ele1->SCluster()->Phi(),
						       ele1->SCluster()->R9(),
						       ele1->SCluster()->EtaWidth(),
						       ele1->SCluster()->PhiWidth(),
						       ele1->NumberOfClusters(),
						       ele1->HadronicOverEm(),
						       fPileUpDen->At(0)->RhoKt6PFJets(),       
						       fPV->GetEntries(),   
						       ele1->SCluster()->Seed()->Eta(),
						       ele1->SCluster()->Seed()->Phi(),
						       ele1->SCluster()->Seed()->Energy(),
						       ele1->SCluster()->Seed()->E3x3(),
						       ele1->SCluster()->Seed()->E5x5(),
						       ele1->CoviEtaiEta(),
						       ele1spp,
						       ele1sep,
						       ele1->SCluster()->Seed()->EMax(),
						       ele1->SCluster()->Seed()->E2nd(),
						       ele1->SCluster()->Seed()->ETop(),
						       ele1->SCluster()->Seed()->EBottom(),
						       ele1->SCluster()->Seed()->ELeft(),
						       ele1->SCluster()->Seed()->ERight(),
						       ele1->SCluster()->Seed()->E2x5Max(),
						       ele1->SCluster()->Seed()->E2x5Top(),
						       ele1->SCluster()->Seed()->E2x5Bottom(),
						       ele1->SCluster()->Seed()->E2x5Left(),
						       ele1->SCluster()->Seed()->E2x5Right(),
						       ele1->SCluster()->Seed()->IEta(),
						       ele1->SCluster()->Seed()->IPhi(),
						       ele1->SCluster()->Seed()->EtaCry(),
						       ele1->SCluster()->Seed()->PhiCry(),
						       ele1->SCluster()->PreshowerEnergy() / ele1->SCluster()->RawEnergy(), 
						       false );
  
  fLeptonPairPhotonEvent->ele1RegressionEnergyV1 = 
    eleRegressionEvaluator_V1->regressionValueWithTrkVarV1(
							   ele1->SCluster()->RawEnergy(),
							   ele1->SCluster()->Eta(),
							   ele1->SCluster()->Phi(),
							   ele1->SCluster()->R9(),
							   ele1->SCluster()->EtaWidth(),
							   ele1->SCluster()->PhiWidth(),
							   ele1->NumberOfClusters(),
							   ele1->HadronicOverEm(),
							   fPileUpDen->At(0)->RhoKt6PFJets(),       
							   fPV->GetEntries(),   
							   ele1->SCluster()->Seed()->Eta(),
							   ele1->SCluster()->Seed()->Phi(),
							   ele1->SCluster()->Seed()->Energy(),
							   ele1->SCluster()->Seed()->E3x3(),
							   ele1->SCluster()->Seed()->E5x5(),
							   ele1->CoviEtaiEta(),
							   ele1spp,
							   ele1sep,
							   ele1->SCluster()->Seed()->EMax(),
							   ele1->SCluster()->Seed()->E2nd(),
							   ele1->SCluster()->Seed()->ETop(),
							   ele1->SCluster()->Seed()->EBottom(),
							   ele1->SCluster()->Seed()->ELeft(),
							   ele1->SCluster()->Seed()->ERight(),
							   ele1->SCluster()->Seed()->E2x5Max(),
							   ele1->SCluster()->Seed()->E2x5Top(),
							   ele1->SCluster()->Seed()->E2x5Bottom(),
							   ele1->SCluster()->Seed()->E2x5Left(),
							   ele1->SCluster()->Seed()->E2x5Right(),
							   ele1->SCluster()->Seed()->IEta(),
							   ele1->SCluster()->Seed()->IPhi(),
							   ele1->SCluster()->Seed()->EtaCry(),
							   ele1->SCluster()->Seed()->PhiCry(),
							   ele1->SCluster()->PreshowerEnergy() / ele1->SCluster()->RawEnergy(), 
							   ele1->PIn(),
							   ele1->FBrem(),
							   ele1->Charge(),
							   ele1->ESuperClusterOverP(),
							   fmin(ele1->TrackMomentumError(),500.0),
							   false );
  
  
  std::vector<double> inputvarsEle1;
  inputvarsEle1.push_back(ele1->SCluster()->RawEnergy());
  inputvarsEle1.push_back(ele1->SCluster()->Eta());
  inputvarsEle1.push_back(ele1->SCluster()->Phi());
  inputvarsEle1.push_back(ele1->SCluster()->R9());
  inputvarsEle1.push_back(ele1->SCluster()->EtaWidth());
  inputvarsEle1.push_back(ele1->SCluster()->PhiWidth());
  inputvarsEle1.push_back(ele1->NumberOfClusters());
  inputvarsEle1.push_back(ele1->HadronicOverEm());
  inputvarsEle1.push_back(fPileUpDen->At(0)->RhoKt6PFJets());       
  inputvarsEle1.push_back(fPV->GetEntries());   
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->Eta());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->Phi());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->Energy());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->E3x3());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->E5x5());
  inputvarsEle1.push_back(ele1->CoviEtaiEta());
  inputvarsEle1.push_back(ele1spp);
  inputvarsEle1.push_back(ele1sep);
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->EMax());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->E2nd());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->ETop());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->EBottom());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->ELeft());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->ERight());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->E2x5Max());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->E2x5Top());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->E2x5Bottom());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->E2x5Left());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->E2x5Right());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->IEta());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->IPhi());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->EtaCry());
  inputvarsEle1.push_back(ele1->SCluster()->Seed()->PhiCry());
  inputvarsEle1.push_back(ele1->SCluster()->PreshowerEnergy() / ele1->SCluster()->RawEnergy()); 
  inputvarsEle1.push_back(ele1->PIn());
  inputvarsEle1.push_back(ele1->FBrem());
  inputvarsEle1.push_back(ele1->Charge());
  inputvarsEle1.push_back(ele1->ESuperClusterOverP());
  inputvarsEle1.push_back(fmin(ele1->TrackMomentumError(),500.0));
  inputvarsEle1.push_back(fmin(fabs(ele1->DeltaEtaSuperClusterTrackAtVtx()), 0.6));
  inputvarsEle1.push_back(ele1->DeltaPhiSuperClusterTrackAtVtx());
  inputvarsEle1.push_back(fmin(ele1->DeltaEtaSeedClusterTrackAtCalo(), 0.2));
  inputvarsEle1.push_back(ele1->DeltaPhiSeedClusterTrackAtCalo());
  inputvarsEle1.push_back(fmin(ele1->BestTrk()->Chi2() / ele1->BestTrk()->Ndof(),200));
  inputvarsEle1.push_back(ele1->CTFTrkNLayersWithMeasurement());
  inputvarsEle1.push_back(fmin(ele1->EEleClusterOverPout(),20.0));
  fLeptonPairPhotonEvent->ele1RegressionEnergyV2 = 
    eleRegressionEvaluator_V2->regressionValueWithTrkVarV2(
							   inputvarsEle1,      
							   false );
	  

  fLeptonPairPhotonEvent->ele1RegressionEnergyErrorV0 = 
    eleRegressionEvaluator_V0->regressionUncertaintyNoTrkVar(
							     ele1->SCluster()->RawEnergy(),
							     ele1->SCluster()->Eta(),
							     ele1->SCluster()->Phi(),
							     ele1->SCluster()->R9(),
							     ele1->SCluster()->EtaWidth(),
							     ele1->SCluster()->PhiWidth(),
							     ele1->NumberOfClusters(),
							     ele1->HadronicOverEm(),
							     fPileUpDen->At(0)->RhoKt6PFJets(),       
							     fPV->GetEntries(),   
							     ele1->SCluster()->Seed()->Eta(),
							     ele1->SCluster()->Seed()->Phi(),
							     ele1->SCluster()->Seed()->Energy(),
							     ele1->SCluster()->Seed()->E3x3(),
							     ele1->SCluster()->Seed()->E5x5(),
							     ele1->CoviEtaiEta(),
							     ele1spp,
							     ele1sep,
							     ele1->SCluster()->Seed()->EMax(),
							     ele1->SCluster()->Seed()->E2nd(),
							     ele1->SCluster()->Seed()->ETop(),
							     ele1->SCluster()->Seed()->EBottom(),
							     ele1->SCluster()->Seed()->ELeft(),
							     ele1->SCluster()->Seed()->ERight(),
							     ele1->SCluster()->Seed()->E2x5Max(),
							     ele1->SCluster()->Seed()->E2x5Top(),
							     ele1->SCluster()->Seed()->E2x5Bottom(),
							     ele1->SCluster()->Seed()->E2x5Left(),
							     ele1->SCluster()->Seed()->E2x5Right(),
							     ele1->SCluster()->Seed()->IEta(),
							     ele1->SCluster()->Seed()->IPhi(),
							     ele1->SCluster()->Seed()->EtaCry(),
							     ele1->SCluster()->Seed()->PhiCry(),
							     ele1->SCluster()->PreshowerEnergy() / ele1->SCluster()->RawEnergy(), 
							     false );

  fLeptonPairPhotonEvent->ele1RegressionEnergyErrorV1 = 
    eleRegressionEvaluator_V1->regressionUncertaintyWithTrkVarV1(
								 ele1->SCluster()->RawEnergy(),
								 ele1->SCluster()->Eta(),
								 ele1->SCluster()->Phi(),
								 ele1->SCluster()->R9(),
								 ele1->SCluster()->EtaWidth(),
								 ele1->SCluster()->PhiWidth(),
								 ele1->NumberOfClusters(),
								 ele1->HadronicOverEm(),
								 fPileUpDen->At(0)->RhoKt6PFJets(),       
								 fPV->GetEntries(),   
								 ele1->SCluster()->Seed()->Eta(),
								 ele1->SCluster()->Seed()->Phi(),
								 ele1->SCluster()->Seed()->Energy(),
								 ele1->SCluster()->Seed()->E3x3(),
								 ele1->SCluster()->Seed()->E5x5(),
								 ele1->CoviEtaiEta(),
								 ele1spp,
								 ele1sep,
								 ele1->SCluster()->Seed()->EMax(),
								 ele1->SCluster()->Seed()->E2nd(),
								 ele1->SCluster()->Seed()->ETop(),
								 ele1->SCluster()->Seed()->EBottom(),
								 ele1->SCluster()->Seed()->ELeft(),
								 ele1->SCluster()->Seed()->ERight(),
								 ele1->SCluster()->Seed()->E2x5Max(),
								 ele1->SCluster()->Seed()->E2x5Top(),
								 ele1->SCluster()->Seed()->E2x5Bottom(),
								 ele1->SCluster()->Seed()->E2x5Left(),
								 ele1->SCluster()->Seed()->E2x5Right(),
								 ele1->SCluster()->Seed()->IEta(),
								 ele1->SCluster()->Seed()->IPhi(),
								 ele1->SCluster()->Seed()->EtaCry(),
								 ele1->SCluster()->Seed()->PhiCry(),
								 ele1->SCluster()->PreshowerEnergy()/ele1->SCluster()->RawEnergy(), 
								 ele1->PIn(),
								 ele1->FBrem(),
								 ele1->Charge(),
								 ele1->ESuperClusterOverP(),
								 fmin(ele1->TrackMomentumError(),500.0),
								 false );
	  
  fLeptonPairPhotonEvent->ele1RegressionEnergyErrorV2 
    = eleRegressionEvaluator_V2->regressionUncertaintyWithTrkVarV2(
								   inputvarsEle1,
								   false );
}


//--------------------------------------------------------------------------------------------------
void LeptonPairPhotonTreeWriter::regressEle2(const mithep::Electron               *ele2,
					     const mithep::PileupEnergyDensityCol *fPileUpDen,
					     const mithep::VertexCol              *fPV   )  
//--------------------------------------------------------------------------------------------------
{
  double ele2spp = ((!isnan(float(ele2->SCluster()->Seed()->CoviPhiiPhi()))) ? sqrt(ele2->SCluster()->Seed()->CoviPhiiPhi()) : 0.0);
  double ele2sep;
  if (ele2->CoviEtaiEta()*ele2spp > 0) {
    ele2sep = ele2->SCluster()->Seed()->CoviEtaiPhi()/(ele2->CoviEtaiEta()*ele2spp);
  } else if (ele2->SCluster()->Seed()->CoviEtaiPhi()) {
    ele2sep = 1.0; 
  } else {
    ele2sep = -1.0; 
  }

  fLeptonPairPhotonEvent->ele2RegressionEnergyV0 = 
    eleRegressionEvaluator_V0->regressionValueNoTrkVar(
						       ele2->SCluster()->RawEnergy(),
						       ele2->SCluster()->Eta(),
						       ele2->SCluster()->Phi(),
						       ele2->SCluster()->R9(),
						       ele2->SCluster()->EtaWidth(),
						       ele2->SCluster()->PhiWidth(),
						       ele2->NumberOfClusters(),
						       ele2->HadronicOverEm(),
						       fPileUpDen->At(0)->RhoKt6PFJets(),       
						       fPV->GetEntries(),   
						       ele2->SCluster()->Seed()->Eta(),
						       ele2->SCluster()->Seed()->Phi(),
						       ele2->SCluster()->Seed()->Energy(),
						       ele2->SCluster()->Seed()->E3x3(),
						       ele2->SCluster()->Seed()->E5x5(),
						       ele2->CoviEtaiEta(),
						       ele2spp,
						       ele2sep,
						       ele2->SCluster()->Seed()->EMax(),
						       ele2->SCluster()->Seed()->E2nd(),
						       ele2->SCluster()->Seed()->ETop(),
						       ele2->SCluster()->Seed()->EBottom(),
						       ele2->SCluster()->Seed()->ELeft(),
						       ele2->SCluster()->Seed()->ERight(),
						       ele2->SCluster()->Seed()->E2x5Max(),
						       ele2->SCluster()->Seed()->E2x5Top(),
						       ele2->SCluster()->Seed()->E2x5Bottom(),
						       ele2->SCluster()->Seed()->E2x5Left(),
						       ele2->SCluster()->Seed()->E2x5Right(),
						       ele2->SCluster()->Seed()->IEta(),
						       ele2->SCluster()->Seed()->IPhi(),
						       ele2->SCluster()->Seed()->EtaCry(),
						       ele2->SCluster()->Seed()->PhiCry(),
						       ele2->SCluster()->PreshowerEnergy() / ele2->SCluster()->RawEnergy(), 
						       false );

  fLeptonPairPhotonEvent->ele2RegressionEnergyV1 = 
    eleRegressionEvaluator_V1->regressionValueWithTrkVarV1(
							   ele2->SCluster()->RawEnergy(),
							   ele2->SCluster()->Eta(),
							   ele2->SCluster()->Phi(),
							   ele2->SCluster()->R9(),
							   ele2->SCluster()->EtaWidth(),
							   ele2->SCluster()->PhiWidth(),
							   ele2->NumberOfClusters(),
							   ele2->HadronicOverEm(),
							   fPileUpDen->At(0)->RhoKt6PFJets(),       
							   fPV->GetEntries(),   
							   ele2->SCluster()->Seed()->Eta(),
							   ele2->SCluster()->Seed()->Phi(),
							   ele2->SCluster()->Seed()->Energy(),
							   ele2->SCluster()->Seed()->E3x3(),
							   ele2->SCluster()->Seed()->E5x5(),
							   ele2->CoviEtaiEta(),
							   ele2spp,
							   ele2sep,
							   ele2->SCluster()->Seed()->EMax(),
							   ele2->SCluster()->Seed()->E2nd(),
							   ele2->SCluster()->Seed()->ETop(),
							   ele2->SCluster()->Seed()->EBottom(),
							   ele2->SCluster()->Seed()->ELeft(),
							   ele2->SCluster()->Seed()->ERight(),
							   ele2->SCluster()->Seed()->E2x5Max(),
							   ele2->SCluster()->Seed()->E2x5Top(),
							   ele2->SCluster()->Seed()->E2x5Bottom(),
							   ele2->SCluster()->Seed()->E2x5Left(),
							   ele2->SCluster()->Seed()->E2x5Right(),
							   ele2->SCluster()->Seed()->IEta(),
							   ele2->SCluster()->Seed()->IPhi(),
							   ele2->SCluster()->Seed()->EtaCry(),
							   ele2->SCluster()->Seed()->PhiCry(),
							   ele2->SCluster()->PreshowerEnergy() / ele2->SCluster()->RawEnergy(), 
							   ele2->PIn(),
							   ele2->FBrem(),
							   ele2->Charge(),
							   ele2->ESuperClusterOverP(),
							   fmin(ele2->TrackMomentumError(),500.0),
							   false );

  std::vector<double> inputvarsEle2;
  inputvarsEle2.push_back(ele2->SCluster()->RawEnergy());
  inputvarsEle2.push_back(ele2->SCluster()->Eta());
  inputvarsEle2.push_back(ele2->SCluster()->Phi());
  inputvarsEle2.push_back(ele2->SCluster()->R9());
  inputvarsEle2.push_back(ele2->SCluster()->EtaWidth());
  inputvarsEle2.push_back(ele2->SCluster()->PhiWidth());
  inputvarsEle2.push_back(ele2->NumberOfClusters());
  inputvarsEle2.push_back(ele2->HadronicOverEm());
  inputvarsEle2.push_back(fPileUpDen->At(0)->RhoKt6PFJets());       
  inputvarsEle2.push_back(fPV->GetEntries());   
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->Eta());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->Phi());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->Energy());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->E3x3());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->E5x5());
  inputvarsEle2.push_back(ele2->CoviEtaiEta());
  inputvarsEle2.push_back(ele2spp);
  inputvarsEle2.push_back(ele2sep);
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->EMax());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->E2nd());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->ETop());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->EBottom());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->ELeft());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->ERight());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->E2x5Max());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->E2x5Top());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->E2x5Bottom());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->E2x5Left());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->E2x5Right());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->IEta());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->IPhi());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->EtaCry());
  inputvarsEle2.push_back(ele2->SCluster()->Seed()->PhiCry());
  inputvarsEle2.push_back(ele2->SCluster()->PreshowerEnergy() / ele2->SCluster()->RawEnergy()); 
  inputvarsEle2.push_back(ele2->PIn());
  inputvarsEle2.push_back(ele2->FBrem());
  inputvarsEle2.push_back(ele2->Charge());
  inputvarsEle2.push_back(ele2->ESuperClusterOverP());
  inputvarsEle2.push_back(fmin(ele2->TrackMomentumError(),500.0));
  inputvarsEle2.push_back(fmin(fabs(ele2->DeltaEtaSuperClusterTrackAtVtx()), 0.6));
  inputvarsEle2.push_back(ele2->DeltaPhiSuperClusterTrackAtVtx());
  inputvarsEle2.push_back(fmin(ele2->DeltaEtaSeedClusterTrackAtCalo(), 0.2));
  inputvarsEle2.push_back(ele2->DeltaPhiSeedClusterTrackAtCalo());
  inputvarsEle2.push_back(fmin(ele2->BestTrk()->Chi2() / ele2->BestTrk()->Ndof(),200));
  inputvarsEle2.push_back(ele2->CTFTrkNLayersWithMeasurement());
  inputvarsEle2.push_back(fmin(ele2->EEleClusterOverPout(),20.0));
  fLeptonPairPhotonEvent->ele2RegressionEnergyV2 = 
    eleRegressionEvaluator_V2->regressionValueWithTrkVarV2(
							   inputvarsEle2,         
							   false );
  

  fLeptonPairPhotonEvent->ele2RegressionEnergyErrorV0 = 
    eleRegressionEvaluator_V0->regressionUncertaintyNoTrkVar(
							     ele2->SCluster()->RawEnergy(),
							     ele2->SCluster()->Eta(),
							     ele2->SCluster()->Phi(),
							     ele2->SCluster()->R9(),
							     ele2->SCluster()->EtaWidth(),
							     ele2->SCluster()->PhiWidth(),
							     ele2->NumberOfClusters(),
							     ele2->HadronicOverEm(),
							     fPileUpDen->At(0)->RhoKt6PFJets(),       
							     fPV->GetEntries(),   
							     ele2->SCluster()->Seed()->Eta(),
							     ele2->SCluster()->Seed()->Phi(),
							     ele2->SCluster()->Seed()->Energy(),
							     ele2->SCluster()->Seed()->E3x3(),
							     ele2->SCluster()->Seed()->E5x5(),
							     ele2->CoviEtaiEta(),
							     ele2spp,
							     ele2sep,
							     ele2->SCluster()->Seed()->EMax(),
							     ele2->SCluster()->Seed()->E2nd(),
							     ele2->SCluster()->Seed()->ETop(),
							     ele2->SCluster()->Seed()->EBottom(),
							     ele2->SCluster()->Seed()->ELeft(),
							     ele2->SCluster()->Seed()->ERight(),
							     ele2->SCluster()->Seed()->E2x5Max(),
							     ele2->SCluster()->Seed()->E2x5Top(),
							     ele2->SCluster()->Seed()->E2x5Bottom(),
							     ele2->SCluster()->Seed()->E2x5Left(),
							     ele2->SCluster()->Seed()->E2x5Right(),
							     ele2->SCluster()->Seed()->IEta(),
							     ele2->SCluster()->Seed()->IPhi(),
							     ele2->SCluster()->Seed()->EtaCry(),
							     ele2->SCluster()->Seed()->PhiCry(),
							     ele2->SCluster()->PreshowerEnergy() / ele2->SCluster()->RawEnergy(), 
							     false );

  fLeptonPairPhotonEvent->ele2RegressionEnergyErrorV1 = 
    eleRegressionEvaluator_V1->regressionUncertaintyWithTrkVarV1(
								 ele2->SCluster()->RawEnergy(),
								 ele2->SCluster()->Eta(),
								 ele2->SCluster()->Phi(),
								 ele2->SCluster()->R9(),
								 ele2->SCluster()->EtaWidth(),
								 ele2->SCluster()->PhiWidth(),
								 ele2->NumberOfClusters(),
								 ele2->HadronicOverEm(),
								 fPileUpDen->At(0)->RhoKt6PFJets(),       
								 fPV->GetEntries(),   
								 ele2->SCluster()->Seed()->Eta(),
								 ele2->SCluster()->Seed()->Phi(),
								 ele2->SCluster()->Seed()->Energy(),
								 ele2->SCluster()->Seed()->E3x3(),
								 ele2->SCluster()->Seed()->E5x5(),
								 ele2->CoviEtaiEta(),
								 ele2spp,
								 ele2sep,
								 ele2->SCluster()->Seed()->EMax(),
								 ele2->SCluster()->Seed()->E2nd(),
								 ele2->SCluster()->Seed()->ETop(),
								 ele2->SCluster()->Seed()->EBottom(),
								 ele2->SCluster()->Seed()->ELeft(),
								 ele2->SCluster()->Seed()->ERight(),
								 ele2->SCluster()->Seed()->E2x5Max(),
								 ele2->SCluster()->Seed()->E2x5Top(),
								 ele2->SCluster()->Seed()->E2x5Bottom(),
								 ele2->SCluster()->Seed()->E2x5Left(),
								 ele2->SCluster()->Seed()->E2x5Right(),
								 ele2->SCluster()->Seed()->IEta(),
								 ele2->SCluster()->Seed()->IPhi(),
								 ele2->SCluster()->Seed()->EtaCry(),
								 ele2->SCluster()->Seed()->PhiCry(),
								 ele2->SCluster()->PreshowerEnergy() / ele2->SCluster()->RawEnergy(), 
								 ele2->PIn(),
								 ele2->FBrem(),
								 ele2->Charge(),
								 ele2->ESuperClusterOverP(),
								 fmin(ele2->TrackMomentumError(),500.0),
								 false );

  fLeptonPairPhotonEvent->ele2RegressionEnergyErrorV2 = 
    eleRegressionEvaluator_V2->regressionUncertaintyWithTrkVarV2(
								 inputvarsEle2,
								 false );
  
}

//--------------------------------------------------------------------------------------------------
void LeptonPairPhotonTreeWriter::Process()
//--------------------------------------------------------------------------------------------------
{
  Bool_t electronZgfilled = kFALSE;  // needed if we pass both channels ... 
  bool store_event_ele = false, store_event_mu = false;

  // Process entries of the tree. 
  LoadEventObject(fPhotonBranchName,   fPhotons);
  LoadEventObject(fGoodElectronName,   fGoodElectrons);
  LoadEventObject(fPVName,             fPV);
  LoadEventObject(fPFCandName,         fPFCands);
  LoadEventObject(fTrackName,          fTracks);
  LoadEventObject(fPileUpDenName,      fPileUpDen);
  LoadEventObject(fGoodMuonName,       fGoodMuons);
  LoadEventObject(fConversionName,     fConversions);
  LoadEventObject(fBeamSpotName,       fBeamSpot);  
  if (!fIsData){
    LoadBranch(fPileUpName);
    LoadBranch(fMCParticleName);
  }


  //Initialize all tree leaf entries to -99.;
  resetTreeVariables();
  fLeptonPairPhotonEvent->run = GetEventHeader()->RunNum();
  fLeptonPairPhotonEvent->lumi = GetEventHeader()->LumiSec();
  fLeptonPairPhotonEvent->event = GetEventHeader()->EvtNum();



  // ************************************************************
  // KH check for good PV
  // ************************************************************
  const mithep::Vertex *bestPV = 0;
  const UInt_t   fMinNTracksFit = 0;
  const Double_t fMinNdof       = 4;
  const Double_t fMaxAbsZ       = 24;
  const Double_t fMaxRho        = 2;

  bool pv_found=false;
  for( unsigned i=0; i<fPV->GetEntries(); i++ ) {
    const mithep::Vertex * pv = fPV->At(i);
    if(!pv->IsValid())                                continue;
    if(pv->Ndof()	          < fMinNdof)	      continue;
    if(fabs(pv->Z()) > fMaxAbsZ)	              continue;
    if(pv->Position().Rho()   > fMaxRho)	      continue;
    
    if( !pv_found ) { 
      bestPV = pv;
      pv_found = true;
      //      cout << "\t ^^^ this PV selected ... ^^^ " << endl;
    }
  }
  if(!pv_found)
    return;
  else if( verbose ) { 
    cout << "pass VTX" << endl;
    cout <<"-----------> best PV <----------------" << endl;
    cout << "X: " << bestPV->X() << "\t"
	 << "Y: " << bestPV->Y() << "\t"
	 << "Z: " << bestPV->Z() << endl;
  }	
  // ************************************************************  

  Int_t _numPU      = -99.;        // some sensible default values....
  Int_t _numPUminus = -99.;        // some sensible default values....
  Int_t _numPUplus  = -99.; 
  
  if (!fIsData){  
    for (UInt_t i=0; i<fPileUp->GetEntries(); ++i) {
      const PileupInfo *puinfo = fPileUp->At(i);
      if      (puinfo->GetBunchCrossing() ==  0) _numPU      = puinfo->GetPU_NumInteractions();
      else if (puinfo->GetBunchCrossing() == -1) _numPUminus = puinfo->GetPU_NumInteractions();
      else if (puinfo->GetBunchCrossing() ==  1) _numPUplus  = puinfo->GetPU_NumInteractions();
    }
  }
  
  fLeptonPairPhotonEvent->NPu      = (float) _numPU;
  fLeptonPairPhotonEvent->NPuMinus = (float) _numPUminus;
  fLeptonPairPhotonEvent->NPuPlus  = (float) _numPUplus;

  // *******************************************************************************************
  //
  // Here we go ...
  //
  // *******************************************************************************************
  if ((fGoodElectrons->GetEntries() > 1) || (fGoodMuons->GetEntries() > 1)){
    
    if( verbose )    cout << "goodEle: " << fGoodElectrons->GetEntries() << "\t"
			  << "goodMu: " << fGoodMuons->GetEntries() << "\t"
			  << "gamma: " << fPhotons->GetEntries() << endl;
    
    //
    // first deal w/ PFnoPU
    //
    vector<bool> pfNoPileUpflag;
    for(UInt_t i = 0; i < fPFCands->GetEntries(); i++) {
      const PFCandidate *pf = fPFCands->At(i);
      assert(pf);
      if(pf->PFType() == PFCandidate::eHadron) {
        if(pf->HasTrackerTrk() && fPV->At(0)->HasTrack(pf->TrackerTrk()) && fPV->At(0)->TrackWeight(pf->TrackerTrk()) > 0) {
          pfNoPileUpflag.push_back(1);
        }
        else {
          Bool_t vertexFound = kFALSE;
          const Vertex *closestVtx = 0;
          Double_t dzmin = 10000;

          for(UInt_t j = 0; j < fPV->GetEntries(); j++) {
            const Vertex *vtx = fPV->At(j);
            assert(vtx);

            if(pf->HasTrackerTrk() && vtx->HasTrack(pf->TrackerTrk()) && vtx->TrackWeight(pf->TrackerTrk()) > 0) {
              vertexFound = kTRUE;
              closestVtx = vtx;
              break;
            }
            Double_t dz = fabs(pf->SourceVertex().Z() - vtx->Z());
            if(dz < dzmin) {
              closestVtx = vtx;
              dzmin = dz;
            }
          }
          if(vertexFound || closestVtx != fPV->At(0)) {
            pfNoPileUpflag.push_back(0);
          }
          else {
            pfNoPileUpflag.push_back(1);
          }
        }
      }else {
        pfNoPileUpflag.push_back(1);
      }
    }
    

    // --------------------------------------------------------------------
    // Photons :: just make a list of those passing ID/Iso here
    //            additional selection after dileptons ...
    // --------------------------------------------------------------------
    vector<mithep::Photon const *>  selected_photons;
    if (fPhotons->GetEntries() >= 1 ){
      for (UInt_t i = 0; i < fPhotons->GetEntries(); ++i){
	const Photon *tmppho = fPhotons->At(i);	
	if( verbose ) 
	  cout << "gamma :: pt: " << tmppho->Pt() << "\teta: " << tmppho->Eta() << "(" << tmppho->SCluster()->Eta() << ")";
	if (tmppho->Pt() > 10 && 
	    (fabs(tmppho->SCluster()->Eta()) < 1.4442 || fabs(tmppho->SCluster()->Eta()) > 1.566) && fabs(tmppho->SCluster()->Eta()) < 2.5)
	  {
	    if( verbose ) cout << "\tpresel";
	    if (ZGTools::photonCutBasedMedium2012ID(tmppho, bestPV, fPFCands, fPileUpDen, fGoodElectrons, fConversions,fPV) ){
	      if( verbose )cout << "\tID";
	      if (ZGTools::photonCutBasedMedium2012Isolation(tmppho, fPV, fPFCands, fPileUpDen, pfNoPileUpflag) ){
		if( verbose ) cout << "\tIso";
		selected_photons.push_back(tmppho);
	      }
	    }
	  }
	if( verbose ) cout << endl;
      }
    }
    // sort selected photons with decreasing pt
    sort( selected_photons.begin(), selected_photons.end(), PhotonPtComparison() );
    const Photon *pho; // this will be the selected Zg photon, assigned below    

    // --------------------------------------------------------------------
    //
    // Dielectrons
    if( _do_ElectronChannel ) { 
    //
    // --------------------------------------------------------------------
      const Electron *ele1=0, *ele2=0; // these will be the selected Z eles, assigned below
      if (fGoodElectrons->GetEntries() > 1){
	vector<bool> electronpass, electronpassID;
	for (UInt_t j = 0; j < fGoodElectrons->GetEntries();++j){
	  const Electron *tmpele = fGoodElectrons->At(j);
	  if( verbose ) cout << "ele :: pt: " << tmpele->Pt() << "\teta: " << tmpele->Eta();
	  if (tmpele->Pt() > 7 && fabs(tmpele->SCluster()->Eta()) < 2.5  ){ //
	    if( verbose )cout << "\tpreSel";
	    bool ele_is_clean=true;
	    for (UInt_t k = 0; k < fGoodMuons->GetEntries(); ++k){
	      if (fabs(fGoodMuons->At(k)->Eta()) < 2.4 && fGoodMuons->At(k)->IsGlobalMuon() &&
		  mithep::MathUtils::DeltaR(tmpele->Phi(),tmpele->Eta(), 
					    fGoodMuons->At(k)->Phi(), fGoodMuons->At(k)->Eta()) < 0.05) {
		ele_is_clean = false;
		break;
	      }
	    }
	    if( ele_is_clean ) { if( verbose )cout << "\tclean"; }
	    if( ele_is_clean && tmpele->Pt() > 10 && ZGTools::electronCutBasedIDLoose(tmpele, bestPV,fConversions,YEAR)){
	      if( verbose )cout << "\tID";
	      if( (YEAR == 2011 && 
		   ZGTools::electronPFIso04(tmpele,fPV->At(0),fPFCands,fPileUpDen, mithep::ElectronTools::kEleEAData2011, pfNoPileUpflag,YEAR)) || 
		  (YEAR == 2012 && 
		   ZGTools::electronPFIso04(tmpele,fPV->At(0),fPFCands,fPileUpDen, mithep::ElectronTools::kEleEAData2012, pfNoPileUpflag,YEAR)) )
		{
		  if( verbose ) cout << "\tIso";
		  electronpass.push_back(1);
		}
	      else  electronpass.push_back(0);
	      electronpassID.push_back(1);
	    } // ID & clean
	    else { electronpass.push_back(0); electronpassID.push_back(0); }
	  } // presel
	  else {electronpass.push_back(0); electronpassID.push_back(0);}
	  if( verbose )cout << endl;
	} // loop over electons
      

	int nIDPairEL=0, nIDIsoPairEL=0;	
	for (UInt_t j = 0; j < fGoodElectrons->GetEntries(); ++j){
	  for (UInt_t k = 0; k < j; ++k){    
	    if (j == k) continue;
	    // KH : no charge here for resync, moved to Z selection
	    //	  if (ele1->Charge() != ele2->Charge() && electronpassID[j] && electronpassID[k]){ 
	    if (electronpassID[j] && electronpassID[k]){ 
	      nIDPairEL++;
	    
	      if (electronpass[j] && electronpass[k] ) { 
		nIDIsoPairEL++;
	      
	      }
	    }
	  }
	}
	if( nIDPairEL >= 1 &&  verbose ){cout << "PASSID :: " << GetEventHeader()->RunNum() 
					      << "\t" << GetEventHeader()->EvtNum() << endl;}
	if( nIDIsoPairEL >= 1 && verbose ){cout << "PASSIDISO :: " << GetEventHeader()->RunNum() 
						<< "\t" << GetEventHeader()->EvtNum() << endl;}
      

	// 
	// now make the ee and eeg systems ...
	// 
	if( nIDPairEL >= 1 && nIDIsoPairEL >= 1 ) {
	
	  // find the e-pair that guves the best Z mass ...
	  bool found_Z=false;
	  Float_t zdifference = 999;
	  UInt_t electron1,electron2;
	  for( unsigned j=0; j<fGoodElectrons->GetEntries(); j++ ) {
	    if(!electronpass[j]) continue;
	    ele1 = fGoodElectrons->At(j);
	    for( unsigned  k=j+1; k<fGoodElectrons->GetEntries(); k++ ) {
	      if(!electronpass[k]) continue;
	      ele2 = fGoodElectrons->At(k);

	      float tmp_mll = (ele1->Mom() + ele2->Mom()).M();	    
	      if( verbose) cout << "j,k:"  << j<< ","<<k << "\ttmp_mll: " << tmp_mll << endl;
	      if( tmp_mll < 50. 
		  || (ele1->Mom().Pt() < 10 || ele2->Mom().Pt() < 10. )  
		  || (ele1->Mom().Pt() < 20 && ele2->Mom().Pt() < 20. )  
		  || (ele1->Charge() == ele2->Charge())  ) { // KH , charge here for resync
		continue;
	      }

	      if( fabs(91.19 - tmp_mll) < zdifference){
		found_Z = true;
		zdifference = fabs(91.19 - tmp_mll);
		electron1 = k;
		electron2 = j;
	      }
	    }
	  }
	
	  if( found_Z ) {  

	    ele1 = fGoodElectrons->At(electron1);
	    ele2 = fGoodElectrons->At(electron2); 
	    float best_mll = (ele1->Mom() + ele2->Mom()).M();
  
	    if( verbose ) cout << "GOODZ :: run " << GetEventHeader()->RunNum() 
			       << "\t" << GetEventHeader()->EvtNum() 
			       << "\t" << best_mll
			       << endl;	  
 
	    for (UInt_t i = 0; i < selected_photons.size(); ++i) {
	      const Photon *tmppho = selected_photons[i];	
	      float mllg = (ele1->Mom() + ele2->Mom() + tmppho->Mom()).M();
	      if( tmppho->Mom().Pt() > 15 && (tmppho->Mom().Pt()/mllg) > 15./110. ) 
		{ 
		  if( verbose ) cout << "PASSPHOTON :: run " << GetEventHeader()->RunNum() << "\t" << GetEventHeader()->EvtNum() << endl;
		
		  if( mithep::MathUtils::DeltaR(ele1->Phi(),ele1->Eta(), tmppho->Phi(), tmppho->Eta()) > 0.4 && 
		      mithep::MathUtils::DeltaR(ele2->Phi(),ele2->Eta(), tmppho->Phi(), tmppho->Eta()) > 0.4 ) {
		    if( verbose )cout << "PASS_DR_LEP_PHO :: run " << GetEventHeader()->RunNum() 
				      << "\t" << GetEventHeader()->EvtNum() 
				      << "\tmllg: " << mllg << endl;
		    if( mllg > 115 ) { 
		      if( verbose ) cout << "MLLG 115 :: run " << GetEventHeader()->RunNum() << "\t" << GetEventHeader()->EvtNum() << endl;
		      if( mllg < 180 ) { 
			if( verbose ) cout << "MLLG 180 :: run " << GetEventHeader()->RunNum() << "\t" << GetEventHeader()->EvtNum() << endl;
			store_event_ele = true;		  
			pho = tmppho;
			break;
		      } // m(llg)<180
		    } // m(llg)>115
		  } // dR(l,g)
		} // pho kinematics
	    } // loop over phos
	  
	  } // goodZ
	} // good dielectron pair
      } // >1 electrons

      //Compute Photon MVA Value and photon quantities
      //	Float_t rho = -99.;
      //	if (fPileUpDen->GetEntries() > 0) rho = (Double_t) fPileUpDen->At(0)->RhoRandomLowEta();
      //      			fLeptonPairPhotonEvent->photonidmva = fTool.GetMVAbdtValue_2011(pho,fPV->At(0),fTracks,fPV,rho,fGoodElectrons,kTRUE);
      if( store_event_ele ) { 
	if( verbose ) cout << "storing electron event ..." << endl;

	electronZgfilled = kTRUE;
	fLeptonPairPhotonEvent->mllg = (ele1->Mom() + ele2->Mom() + pho->Mom()).M();

	fLeptonPairPhotonEvent->photonenergy = pho->E();
	fLeptonPairPhotonEvent->photonpx = pho->Px();
	fLeptonPairPhotonEvent->photonpy = pho->Py();
	fLeptonPairPhotonEvent->photonpz = pho->Pz();
	fLeptonPairPhotonEvent->photonpt = pho->Pt();
	fLeptonPairPhotonEvent->photoneta = pho->Eta();
	fLeptonPairPhotonEvent->photonmass = pho->Mass();
	fLeptonPairPhotonEvent->photonphi = pho->Phi();
	fLeptonPairPhotonEvent->photonr9 = pho->R9();	
	fLeptonPairPhotonEvent->photonenergyerror = pho->EnergyErr();
	// check if pho matched
	const MCParticle *phgen = NULL;
	if (!fIsData){
	  phgen = PhotonTools::MatchMC(pho,fMCParticles,!fApplyElectronVeto);
	  // 	cout << "It went into Photon Tools phgen loop." << endl;
	}
	if (phgen != NULL){
	  fLeptonPairPhotonEvent->photonmatchmc = 1;
	  //	cout << "photonmatchmc written to be kTRUE" << endl;	
	}
	else fLeptonPairPhotonEvent->photonmatchmc = 0;	


	//
	// photon phosphor corrections not applied right now ...
	//
	TLorentzVector gcorr;
	gcorr.SetPtEtaPhiM(fLeptonPairPhotonEvent->photonpt,
			   fLeptonPairPhotonEvent->photoneta,
			   fLeptonPairPhotonEvent->photonphi,
			   fLeptonPairPhotonEvent->photonmass);
	if( phgen != NULL ) {  // this also means it's MC ...

          float eCorr = 0;
#ifdef PHOSPHOR_CORRECTIONS_HEADER
	  eCorr = phosphor->GetCorrEnergy(pho->R9(),YEAR,pho->Pt(),pho->SCluster()->Eta(),phgen->E());
#endif
	  /*
	  float xMC = (pho->E()/phgen->E())-1;
	  std::pair<float,float> mcScaleRes 
	    = ZGTools::getPhosphorScaleRes(YEAR, true, pho->SCluster()->Eta(), pho->Pt(), pho->R9());
	  std::pair<float,float> dataScaleRes 
	    = ZGTools::getPhosphorScaleRes(YEAR, false, pho->SCluster()->Eta(), pho->Pt(), pho->R9());
	  float rData_over_rMC = dataScaleRes.second/mcScaleRes.second;
	  float sMC = mcScaleRes.first;
	  float xCorrSmear = rData_over_rMC*(xMC-0.01*sMC);
	  float eCorr = (1+xCorrSmear)*phgen->E();
	  */
	  float sf = eCorr/pho->E();
	  gcorr.SetE(eCorr);
	  gcorr.SetPx(sf*pho->Px());
	  gcorr.SetPy(sf*pho->Py());
	  gcorr.SetPz(sf*pho->Pz());
	}
	if( fIsData ) { 
	  float eCorr = 0;
#ifdef PHOSPHOR_CORRECTIONS_HEADER
	  phosphor->GetCorrEnergy(pho->R9(),YEAR,pho->Pt(),pho->SCluster()->Eta());
#endif
	  /*
	  std::pair<float,float> dataScaleRes 
	    = ZGTools::getPhosphorScaleRes(YEAR, false, pho->SCluster()->Eta(), pho->Pt(), pho->R9());
	  float eCorr = pho->E()/(1.+dataScaleRes.first);
	  */
	  float sf = eCorr/pho->E();
	  gcorr.SetE(eCorr);
	  gcorr.SetPx(sf*pho->Px());
	  gcorr.SetPy(sf*pho->Py());
	  gcorr.SetPz(sf*pho->Pz());
	}
	fLeptonPairPhotonEvent->photonenergyCorr = gcorr.E();
	fLeptonPairPhotonEvent->photonpxCorr = gcorr.Px();
	fLeptonPairPhotonEvent->photonpyCorr = gcorr.Py();
	fLeptonPairPhotonEvent->photonpzCorr = gcorr.Pz();
	fLeptonPairPhotonEvent->photonptCorr = gcorr.Pt();	
	//


	fillEle1Variables(ele1);
	fillEle2Variables(ele2);
	regressEle1(ele1,fPileUpDen,fPV);
	regressEle2(ele2,fPileUpDen,fPV);

     	//
	// electron scale/res corrections
	//
#if xELECTRON_CORRECTIONS_HEADER != x 
	char tmpbuf[256];
	if( YEAR == 2011 ) sprintf( tmpbuf, "2011");
	else sprintf(tmpbuf, "HCP2012");
	string erastr( tmpbuf);

	float ele1CorrFactor = correctedElectronEnergy( ele1->E(),
							ele1->SCluster()->Eta(),
							ele1->SCluster()->R9(),
							GetEventHeader()->RunNum(), 
							0, 
							erastr,
							!(fIsData),
							&rand) / ele1->E();
	fLeptonPairPhotonEvent->ele1energyCorr = ele1CorrFactor*fLeptonPairPhotonEvent->ele1energy;
	fLeptonPairPhotonEvent->ele1pxCorr     = ele1CorrFactor*fLeptonPairPhotonEvent->ele1px;
	fLeptonPairPhotonEvent->ele1pyCorr     = ele1CorrFactor*fLeptonPairPhotonEvent->ele1py;
	fLeptonPairPhotonEvent->ele1pzCorr     = ele1CorrFactor*fLeptonPairPhotonEvent->ele1pz;
	fLeptonPairPhotonEvent->ele1ptCorr     = ele1CorrFactor*fLeptonPairPhotonEvent->ele1pt;

	float ele2CorrFactor = correctedElectronEnergy( ele2->E(),
							ele2->SCluster()->Eta(),
							ele2->SCluster()->R9(),
							GetEventHeader()->RunNum(), 
							0, 
							erastr,
							!(fIsData),
							&rand) / ele2->E();
	fLeptonPairPhotonEvent->ele2energyCorr = ele2CorrFactor*fLeptonPairPhotonEvent->ele2energy;
	fLeptonPairPhotonEvent->ele2pxCorr     = ele2CorrFactor*fLeptonPairPhotonEvent->ele2px;
	fLeptonPairPhotonEvent->ele2pyCorr     = ele2CorrFactor*fLeptonPairPhotonEvent->ele2py;
	fLeptonPairPhotonEvent->ele2pzCorr     = ele2CorrFactor*fLeptonPairPhotonEvent->ele2pz;
	fLeptonPairPhotonEvent->ele2ptCorr     = ele2CorrFactor*fLeptonPairPhotonEvent->ele2pt;
	//
	//
	//



	TLorentzVector e1corr, e2corr;
	e1corr.SetPxPyPzE(fLeptonPairPhotonEvent->ele1pxCorr,
			  fLeptonPairPhotonEvent->ele1pyCorr,
			  fLeptonPairPhotonEvent->ele1pzCorr,
			  fLeptonPairPhotonEvent->ele1energyCorr);
	e2corr.SetPxPyPzE(fLeptonPairPhotonEvent->ele2pxCorr,
			  fLeptonPairPhotonEvent->ele2pyCorr,
			  fLeptonPairPhotonEvent->ele2pzCorr,
			  fLeptonPairPhotonEvent->ele2energyCorr);

	fLeptonPairPhotonEvent->mllgCorr = (e1corr + e2corr + gcorr).M();
#endif // electron corrections

	ZGLabVectors l;
	ZGAngles b;
      
	l.vecg.SetPxPyPzE(pho->Px(),pho->Py(),pho->Pz(),pho->E());
	l.veclp.SetPxPyPzE(ele1->Px(),ele1->Py(),ele1->Pz(),ele1->E());
	l.veclm.SetPxPyPzE(ele2->Px(),ele2->Py(),ele2->Pz(),ele2->E());
	l.vecz = (l.veclp+l.veclm);
	l.veczg = (l.vecg+l.veclp+l.veclm);
      
	b = ZGTools::getZGAngles(l,kFALSE);
      
	fLeptonPairPhotonEvent->costheta_lm_electrons = b.costheta_lm;
	fLeptonPairPhotonEvent->costheta_lp_electrons = b.costheta_lp;
	fLeptonPairPhotonEvent->phi_electrons = b.phi;
	fLeptonPairPhotonEvent->cosTheta_electrons = b.cosTheta;

	fLeptonPairPhotonEvent->electronZmass = (ele1->Mom() + ele2->Mom()).M();
      
	fLeptonPairPhotonEvent->ele1dEtaIn = TMath::Abs(ele1->DeltaEtaSuperClusterTrackAtVtx());
	fLeptonPairPhotonEvent->ele1dPhiIn = TMath::Abs(ele1->DeltaPhiSuperClusterTrackAtVtx());
	fLeptonPairPhotonEvent->ele1sigmaIEtaIEta = ele1->CoviEtaiEta();
	fLeptonPairPhotonEvent->ele1HadOverEm = ele1->HadronicOverEm();
	fLeptonPairPhotonEvent->ele1D0 = fabs(ele1->BestTrk()->D0Corrected(*bestPV));
	fLeptonPairPhotonEvent->ele1DZ = fabs(ele1->BestTrk()->DzCorrected(*bestPV));
	fLeptonPairPhotonEvent->ele1OneOverEMinusOneOverP = fabs((1 - ele1->ESuperClusterOverP())/(ele1->EcalEnergy()));
	fLeptonPairPhotonEvent->ele1Conversion = mithep::ElectronTools::PassConversionFilter(ele1,fConversions,bestPV, 0, 1e-6, 2.0, kTRUE, kFALSE, 7);
	fLeptonPairPhotonEvent->ele1missinghits = ele1->CorrectedNExpectedHitsInner();

	fLeptonPairPhotonEvent->ele2dEtaIn = TMath::Abs(ele2->DeltaEtaSuperClusterTrackAtVtx());
	fLeptonPairPhotonEvent->ele2dPhiIn = TMath::Abs(ele2->DeltaPhiSuperClusterTrackAtVtx());
	fLeptonPairPhotonEvent->ele2sigmaIEtaIEta = ele2->CoviEtaiEta();
	fLeptonPairPhotonEvent->ele2HadOverEm = ele2->HadronicOverEm();
	fLeptonPairPhotonEvent->ele2D0 = fabs(ele2->BestTrk()->D0Corrected(*bestPV));
	fLeptonPairPhotonEvent->ele2DZ = fabs(ele2->BestTrk()->DzCorrected(*bestPV));
	fLeptonPairPhotonEvent->ele2OneOverEMinusOneOverP = fabs((1 - ele2->ESuperClusterOverP())/(ele2->EcalEnergy()));
	fLeptonPairPhotonEvent->ele2Conversion = mithep::ElectronTools::PassConversionFilter(ele2,fConversions,bestPV, 0, 1e-6, 2.0, kTRUE, kFALSE,7);
	fLeptonPairPhotonEvent->ele2missinghits = ele2->CorrectedNExpectedHitsInner();

      } // store event
    } // doElectronChannel

    // --------------------------------------------------------------------
    //
    // Dimuons
    if( _do_MuonChannel ) { 
    //
    // --------------------------------------------------------------------
      const Muon *muona=0, *muonb=0; // these will be the selected Z muons, assigned below
      if (fGoodMuons->GetEntries() > 1) {
	vector<bool> muonpass, muonpassID;
	for (UInt_t j = 0; j < fGoodMuons->GetEntries();++j){
	  const Muon *tmpmu  = fGoodMuons->At(j);
	  if( verbose ) cout << "mu :: pt: " << tmpmu->Pt() << "\teta: " << tmpmu->Eta();
	  if (tmpmu->Pt() > 10 && fabs(tmpmu->Eta()) < 2.4 && tmpmu->IsGlobalMuon() == true){
	    if( verbose ) cout << "\tpresel";
	    if (ZGTools::muonIDPOGTightSelection(YEAR,tmpmu,bestPV,fPFCands)){
	      if( verbose )cout << "\tID";
	      if ((YEAR == 2011 && ZGTools::muonPFIso04(tmpmu,fPV->At(0),fPFCands,fPileUpDen, mithep::MuonTools::kMuEAData2011, pfNoPileUpflag,YEAR)) || 
		  (YEAR == 2012 && ZGTools::muonPFIso04(tmpmu,fPV->At(0),fPFCands,fPileUpDen, mithep::MuonTools::kMuEAData2012, pfNoPileUpflag,YEAR))){
		muonpass.push_back(1);
		if( verbose ) cout << "\tIso";
	      } else muonpass.push_back(0);
	      muonpassID.push_back(1);
	    }
	    else { 
	      muonpass.push_back(0);
	      muonpassID.push_back(0);
	    }
	  }
	  else { 
	    muonpass.push_back(0);
	    muonpassID.push_back(0);
	  }
	  if( verbose ) cout << endl;
	} // loop over muons
      
     
	
	int nIDPairMU=0, nIDIsoPairMU=0;
	for (UInt_t j = 0; j < fGoodMuons->GetEntries(); ++j){
	  for (UInt_t k = 0; k < j; ++k){
	    if (j == k) continue;
	    muona = fGoodMuons->At(j);
	    muonb = fGoodMuons->At(k);
	    if (muonpassID[j] && muonpassID[k]){ // KH, no charge here in resync muona->Charge() != muonb->Charge() && 
	      nIDPairMU++;
	      if (muonpass[j] && muonpass[k]){ // KH, no charge here in resync muona->Charge() != muonb->Charge() && 
		nIDIsoPairMU++;
	      }
	    }
	  }
	}
      
	if( nIDPairMU >= 1 &&  verbose ){ cout << "PASSID :: " << GetEventHeader()->RunNum() << "\t" << GetEventHeader()->EvtNum() << endl; }
	if( nIDIsoPairMU >= 1 &&  verbose ){ cout << "PASSIDISO :: " << GetEventHeader()->RunNum() << "\t" << GetEventHeader()->EvtNum() << endl; }


	// 
	// now make the mm and mmg systems ...
	//       
	if( nIDPairMU >= 1 && nIDIsoPairMU >= 1 ) {


	  // find the mu-pair that guves the best Z mass ...
	  bool found_Z=false;
	  Float_t zdifference = 999;
	  UInt_t mindex1,mindex2;
	  for( unsigned j=0; j<fGoodMuons->GetEntries(); j++ ) {
	    if(!muonpass[j]) continue;
	    muona = fGoodMuons->At(j);
	    for( unsigned  k=j+1; k<fGoodMuons->GetEntries(); k++ ) {
	      if(!muonpass[k]) continue;
	      muonb = fGoodMuons->At(k);

	      float tmp_mll = (muona->Mom() + muonb->Mom()).M();	    
	      if( verbose ) cout << "j,k:"  << j<< ","<<k << "\ttmp_mll: " << tmp_mll << endl;
	      if( tmp_mll < 50. 
		  || (muona->Mom().Pt() < 10 || muonb->Mom().Pt() < 10. )  
		  || (muona->Mom().Pt() < 20 && muonb->Mom().Pt() < 20. )  
		  || (muona->Charge() == muonb->Charge())  ) { // KH , charge here for resync
		continue;
	      }

	      if( fabs(91.19 - tmp_mll) < zdifference){
		found_Z = true;
		zdifference = fabs(91.19 - tmp_mll);
		mindex1 = k;
		mindex2 = j;
	      }
	    }
	  }
	
	  if( found_Z ) {  

	    muona = fGoodMuons->At(mindex1);
	    muonb = fGoodMuons->At(mindex2);
	    float best_mll = (muona->Mom() + muonb->Mom()).M();
	    if( verbose ) cout << "GOODZ :: run " << GetEventHeader()->RunNum() 
			       << "\t" << GetEventHeader()->EvtNum() 
			       << "\t" << best_mll
			       << endl;	  
	    
	    for (UInt_t i = 0; i < selected_photons.size(); ++i) {
	      const Photon *tmppho = selected_photons[i];	
	      float mllg = (muona->Mom() + muonb->Mom() + tmppho->Mom()).M();
	      if( tmppho->Mom().Pt() > 15 && (tmppho->Mom().Pt()/mllg) > 15./110. ) {  // KH , add scaled pT cut
		if( verbose ) cout << "PASSPHOTON :: run " << GetEventHeader()->RunNum() << "\t" << GetEventHeader()->EvtNum() 
				   << "\t" << tmppho->Mom().Pt() << "\t" << tmppho->Mom().Eta() << endl;
		if( mithep::MathUtils::DeltaR(muona->Phi(),muona->Eta(), tmppho->Phi(), tmppho->Eta()) > 0.4 && 
		    mithep::MathUtils::DeltaR(muonb->Phi(),muonb->Eta(), tmppho->Phi(), tmppho->Eta()) > 0.4 ) {
		  if( verbose )cout << "PASS_DR_LEP_PHO :: run " << GetEventHeader()->RunNum() 
				    << "\t" << GetEventHeader()->EvtNum() 
				    << "\tmllg: " << mllg << endl;
		  if( mllg > 115 ) { 
		    if( verbose ) cout << "MLLG 115 :: run " << GetEventHeader()->RunNum() << "\t" << GetEventHeader()->EvtNum() << endl;
		    if( mllg < 180 ) { 
		      if( verbose ) cout << "MLLG 180 :: run " << GetEventHeader()->RunNum() << "\t" << GetEventHeader()->EvtNum() << endl;
		      store_event_mu = true;		  
		      pho = tmppho;
		      break;
		    } // m(llg)<180
		  } // m(llg)>115
		} // dR(l,g)
	      } // pho kinematics
	    } // loop over phos
	  
	  } // goodZ
	} // good dimuon pair
      } // >1 muons



      if( store_event_mu ) { 

	fLeptonPairPhotonEvent->mllg = (pho->Mom()+muona->Mom()+muonb->Mom()).M();

	fLeptonPairPhotonEvent->photonenergy = pho->E();
	fLeptonPairPhotonEvent->photonpx = pho->Px();
	fLeptonPairPhotonEvent->photonpy = pho->Py();
	fLeptonPairPhotonEvent->photonpz = pho->Pz();
	fLeptonPairPhotonEvent->photonpt = pho->Pt();
	fLeptonPairPhotonEvent->photoneta = pho->Eta();
	fLeptonPairPhotonEvent->photonmass = pho->Mass();
	fLeptonPairPhotonEvent->photonphi = pho->Phi();
	fLeptonPairPhotonEvent->photonr9 = pho->R9();	
	fLeptonPairPhotonEvent->photonenergyerror = pho->EnergyErr();
	// check if pho matched
	const MCParticle *phgen = NULL;
	if (!fIsData){
	  phgen = PhotonTools::MatchMC(pho,fMCParticles,!fApplyElectronVeto);
	  // 	cout << "It went into Photon Tools phgen loop." << endl;
	}
	if (phgen != NULL){
	  fLeptonPairPhotonEvent->photonmatchmc = 1;
	  //	cout << "photonmatchmc written to be kTRUE" << endl;	
	}
	else fLeptonPairPhotonEvent->photonmatchmc = 0;	
      

	//
	// photon phosphor corrections not applied right now ...
	//
	TLorentzVector gcorr;
	gcorr.SetPtEtaPhiM(fLeptonPairPhotonEvent->photonpt,
			   fLeptonPairPhotonEvent->photoneta,
			   fLeptonPairPhotonEvent->photonphi,
			   fLeptonPairPhotonEvent->photonmass);
	if( phgen != NULL ) {  // also means it's MC ...
	  float xMC = (pho->E()/phgen->E())-1;
	  std::pair<float,float> mcScaleRes 
	    = ZGTools::getPhosphorScaleRes(YEAR, true, pho->SCluster()->Eta(), pho->Pt(), pho->R9());
	  std::pair<float,float> dataScaleRes 
	    = ZGTools::getPhosphorScaleRes(YEAR, false, pho->SCluster()->Eta(), pho->Pt(), pho->R9());
	  float rData_over_rMC = dataScaleRes.second/mcScaleRes.second;
	  float sMC = mcScaleRes.first;
	  float xCorrSmear = rData_over_rMC*(xMC-0.01*sMC);
	  float eCorr = (1+xCorrSmear)*phgen->E();
	  assert(eCorr>0);
	  float sf = eCorr/pho->E();
	  gcorr.SetE(eCorr);
	  gcorr.SetPx(sf*pho->Px());
	  gcorr.SetPy(sf*pho->Py());
	  gcorr.SetPz(sf*pho->Pz());
	}
	if( fIsData ) { 
	  std::pair<float,float> dataScaleRes 
	    = ZGTools::getPhosphorScaleRes(YEAR, false, pho->SCluster()->Eta(), pho->Pt(), pho->R9());
	  float eCorr = pho->E()/(1.+dataScaleRes.first);
	  float sf = eCorr/pho->E();
	  gcorr.SetE(eCorr);
	  gcorr.SetPx(sf*pho->Px());
	  gcorr.SetPy(sf*pho->Py());
	  gcorr.SetPz(sf*pho->Pz());
	}
	fLeptonPairPhotonEvent->photonenergyCorr = gcorr.E();
	fLeptonPairPhotonEvent->photonpxCorr = gcorr.Px();
	fLeptonPairPhotonEvent->photonpyCorr = gcorr.Py();
	fLeptonPairPhotonEvent->photonpzCorr = gcorr.Pz();
	fLeptonPairPhotonEvent->photonptCorr = gcorr.Pt();	
	//

	if (electronZgfilled == kTRUE) fLeptonPairPhotonEvent->muonZgVeto = kTRUE;

	ZGLabVectors l;
	ZGAngles b;
	l.vecg.SetPxPyPzE(pho->Px(),pho->Py(),pho->Pz(),pho->E());
	l.veclp.SetPxPyPzE(muona->Px(),muona->Py(),muona->Pz(),muona->E());
	l.veclm.SetPxPyPzE(muonb->Px(),muonb->Py(),muonb->Pz(),muonb->E());
	l.vecz = (l.veclp+l.veclm);
	l.veczg = (l.vecg+l.veclp+l.veclm);
      
	b = ZGTools::getZGAngles(l,kFALSE);
      
	fLeptonPairPhotonEvent->costheta_lm_muons = b.costheta_lm;
	fLeptonPairPhotonEvent->costheta_lp_muons = b.costheta_lp;
	fLeptonPairPhotonEvent->phi_muons = b.phi;
	fLeptonPairPhotonEvent->cosTheta_muons = b.cosTheta;

	fLeptonPairPhotonEvent->muonZmass = (muona->Mom() + muonb->Mom()).M();
	fLeptonPairPhotonEvent->m1E = muona->E();
	fLeptonPairPhotonEvent->m1Pt = muona->Pt();
	fLeptonPairPhotonEvent->m1Mass = muona->Mass();
	fLeptonPairPhotonEvent->m1Px = muona->Px();
	fLeptonPairPhotonEvent->m1Py = muona->Py();
	fLeptonPairPhotonEvent->m1Pz = muona->Pz();
	fLeptonPairPhotonEvent->m1Eta = muona->Eta();
	fLeptonPairPhotonEvent->m1Phi = muona->Phi();
	fLeptonPairPhotonEvent->m1Charge = muona->Charge();
	if (muona->TrackerTrk()) fLeptonPairPhotonEvent->m1PtErr = muona->TrackerTrk()->PtErr();
	else fLeptonPairPhotonEvent->m1PtErr = muona->BestTrk()->PtErr();
	fLeptonPairPhotonEvent->m2E = muonb->E();
	fLeptonPairPhotonEvent->m2Pt = muonb->Pt();
	fLeptonPairPhotonEvent->m2Mass = muonb->Mass();
	fLeptonPairPhotonEvent->m2Px = muonb->Px();
	fLeptonPairPhotonEvent->m2Py = muonb->Py();
	fLeptonPairPhotonEvent->m2Pz = muonb->Pz();
	fLeptonPairPhotonEvent->m2Eta = muonb->Eta();
	fLeptonPairPhotonEvent->m2Phi = muonb->Phi();
	fLeptonPairPhotonEvent->m2Charge = muonb->Charge();
	if (muonb->TrackerTrk()) fLeptonPairPhotonEvent->m2PtErr = muonb->TrackerTrk()->PtErr();
	else fLeptonPairPhotonEvent->m2PtErr = muonb->BestTrk()->PtErr();

	//
	// scale/res corrections
	//
	TLorentzVector m1corr, m2corr;
	m1corr.SetPtEtaPhiM(muona->Pt(),muona->Eta(),muona->Phi(),muona->Mass());
	m2corr.SetPtEtaPhiM(muonb->Pt(),muonb->Eta(),muonb->Phi(),muonb->Mass());
	if(fIsData) { 
	  rmcor->momcor_data(m1corr, muona->Charge(), 0, 0); 
	  rmcor->momcor_data(m2corr, muonb->Charge(), 0, 0); 
	} else { 
	  rmcor->momcor_mc(m1corr, muona->Charge(), 0, 0); 
	  rmcor->momcor_mc(m2corr, muonb->Charge(), 0, 0); 
	}
	fLeptonPairPhotonEvent->m1ECorr  = m1corr.E();
	fLeptonPairPhotonEvent->m1PtCorr = m1corr.Pt();
	fLeptonPairPhotonEvent->m1PxCorr = m1corr.Px();
	fLeptonPairPhotonEvent->m1PyCorr = m1corr.Py();
	fLeptonPairPhotonEvent->m1PzCorr = m1corr.Pz();
	fLeptonPairPhotonEvent->m2ECorr  = m2corr.E();
	fLeptonPairPhotonEvent->m2PtCorr = m2corr.Pt();
	fLeptonPairPhotonEvent->m2PxCorr = m2corr.Px();
	fLeptonPairPhotonEvent->m2PyCorr = m2corr.Py();
	fLeptonPairPhotonEvent->m2PzCorr = m2corr.Pz();


	fLeptonPairPhotonEvent->mllgCorr = (m1corr+m2corr+gcorr).M();
	//
	//
	//

      } // store event   
    } // doMuonChannel

    if( store_event_ele || store_event_mu ) ZgllTuple->Fill();
  }
}


//--------------------------------------------------------------------------------------------------
void LeptonPairPhotonTreeWriter::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the branch collections.

  ReqEventObject(fPhotonBranchName,	fPhotons,      	fPhotonsFromBranch);
  ReqEventObject(fGoodElectronName,	fGoodElectrons,	fGoodElectronsFromBranch);  
  ReqEventObject(fPVName,             	fPV, 		fPVFromBranch);
  ReqEventObject(fPFCandName,         	fPFCands, 	fPFCandsFromBranch);
  ReqEventObject(fTrackName,    	fTracks, 	fTracksFromBranch);
  ReqEventObject(fPileUpDenName,      	fPileUpDen, 	fPileUpDenFromBranch);
  ReqEventObject(fGoodMuonName,    	fGoodMuons,    	fGoodMuonsFromBranch); 
  if (!fIsData){
    ReqBranch(fPileUpName,         fPileUp);
    ReqBranch(fMCParticleName,     fMCParticles);
  }
  ReqEventObject(fConversionName,  fConversions,  true);
  ReqEventObject(fBeamSpotName,    fBeamSpot,     true);

  fLeptonPairPhotonEvent = new LeptonPairPhotonEvent;//Declares all tree leaves  
  ZgllTuple = new TTree(fTupleName.Data(),fTupleName.Data());//Declares tree name  
  ZgllTuple->SetAutoSave(300e9);
  

  ZgllTuple->Branch("run",&fLeptonPairPhotonEvent->run,"run/i");
  ZgllTuple->Branch("lumi",&fLeptonPairPhotonEvent->lumi,"lumi/i");
  ZgllTuple->Branch("event",&fLeptonPairPhotonEvent->event,"event/i");

  ZgllTuple->Branch("electronZmass",&fLeptonPairPhotonEvent->electronZmass,"electronZmass/F");
  ZgllTuple->Branch("mllg",&fLeptonPairPhotonEvent->mllg,"mllg/F");
  ZgllTuple->Branch("mllgCorr",&fLeptonPairPhotonEvent->mllgCorr,"mllgCorr/F");
  ZgllTuple->Branch("ele1MVA",&fLeptonPairPhotonEvent->ele1MVA,"ele1MVA/F");
  ZgllTuple->Branch("ele2MVA",&fLeptonPairPhotonEvent->ele2MVA,"ele2MVA/F");
  
  ZgllTuple->Branch("ele1charge",&fLeptonPairPhotonEvent->ele1charge,"ele1charge/F");
  ZgllTuple->Branch("ele1energy",&fLeptonPairPhotonEvent->ele1energy,"ele1energy/F");
  ZgllTuple->Branch("ele1px",&fLeptonPairPhotonEvent->ele1px,"ele1px/F");
  ZgllTuple->Branch("ele1py",&fLeptonPairPhotonEvent->ele1py,"ele1py/F");
  ZgllTuple->Branch("ele1pz",&fLeptonPairPhotonEvent->ele1pz,"ele1pz/F");
  ZgllTuple->Branch("ele1pt",&fLeptonPairPhotonEvent->ele1pt,"ele1pt/F");
  ZgllTuple->Branch("ele1eta",&fLeptonPairPhotonEvent->ele1eta,"ele1eta/F");
  ZgllTuple->Branch("ele1mass",&fLeptonPairPhotonEvent->ele1mass,"ele1mass/F");
  ZgllTuple->Branch("ele1phi",&fLeptonPairPhotonEvent->ele1phi,"ele1phi/F");
  ZgllTuple->Branch("ele1RegressionEnergyV0",&fLeptonPairPhotonEvent->ele1RegressionEnergyV0,"ele1RegressionEnergyV0/F");
  ZgllTuple->Branch("ele1RegressionEnergyV1",&fLeptonPairPhotonEvent->ele1RegressionEnergyV1,"ele1RegressionEnergyV1/F");
  ZgllTuple->Branch("ele1RegressionEnergyV2",&fLeptonPairPhotonEvent->ele1RegressionEnergyV2,"ele1RegressionEnergyV2/F");
  ZgllTuple->Branch("ele1RegressionEnergyErrorV0",&fLeptonPairPhotonEvent->ele1RegressionEnergyErrorV0,"ele1RegressionEnergyErrorV0/F");
  ZgllTuple->Branch("ele1RegressionEnergyErrorV1",&fLeptonPairPhotonEvent->ele1RegressionEnergyErrorV1,"ele1RegressionEnergyErrorV1/F");
  ZgllTuple->Branch("ele1RegressionEnergyErrorV2",&fLeptonPairPhotonEvent->ele1RegressionEnergyErrorV2,"ele1RegressionEnergyErrorV2/F");

  // scale/res corrected 
  ZgllTuple->Branch("ele1energyCorr",&fLeptonPairPhotonEvent->ele1energyCorr,"ele1energyCorr/F");
  ZgllTuple->Branch("ele1pxCorr",&fLeptonPairPhotonEvent->ele1pxCorr,"ele1pxCorr/F");
  ZgllTuple->Branch("ele1pyCorr",&fLeptonPairPhotonEvent->ele1pyCorr,"ele1pyCorr/F");
  ZgllTuple->Branch("ele1pzCorr",&fLeptonPairPhotonEvent->ele1pzCorr,"ele1pzCorr/F");
  ZgllTuple->Branch("ele1ptCorr",&fLeptonPairPhotonEvent->ele1ptCorr,"ele1ptCorr/F");


  ZgllTuple->Branch("ele2charge",&fLeptonPairPhotonEvent->ele2charge,"ele2charge/F");
  ZgllTuple->Branch("ele2energy",&fLeptonPairPhotonEvent->ele2energy,"ele2energy/F");
  ZgllTuple->Branch("ele2px",&fLeptonPairPhotonEvent->ele2px,"ele2px/F");
  ZgllTuple->Branch("ele2py",&fLeptonPairPhotonEvent->ele2py,"ele2py/F");
  ZgllTuple->Branch("ele2pz",&fLeptonPairPhotonEvent->ele2pz,"ele2pz/F");
  ZgllTuple->Branch("ele2pt",&fLeptonPairPhotonEvent->ele2pt,"ele2pt/F");
  ZgllTuple->Branch("ele2eta",&fLeptonPairPhotonEvent->ele2eta,"ele2eta/F");
  ZgllTuple->Branch("ele2mass",&fLeptonPairPhotonEvent->ele2mass,"ele2mass/F");
  ZgllTuple->Branch("ele2phi",&fLeptonPairPhotonEvent->ele2phi,"ele2phi/F");
  ZgllTuple->Branch("ele2RegressionEnergyV0",&fLeptonPairPhotonEvent->ele2RegressionEnergyV0,"ele2RegressionEnergyV0/F");
  ZgllTuple->Branch("ele2RegressionEnergyV1",&fLeptonPairPhotonEvent->ele2RegressionEnergyV1,"ele2RegressionEnergyV1/F");
  ZgllTuple->Branch("ele2RegressionEnergyV2",&fLeptonPairPhotonEvent->ele2RegressionEnergyV2,"ele2RegressionEnergyV2/F");
  ZgllTuple->Branch("ele2RegressionEnergyErrorV0",&fLeptonPairPhotonEvent->ele2RegressionEnergyErrorV0,"ele2RegressionEnergyErrorV0/F");
  ZgllTuple->Branch("ele2RegressionEnergyErrorV1",&fLeptonPairPhotonEvent->ele2RegressionEnergyErrorV1,"ele2RegressionEnergyErrorV1/F");
  ZgllTuple->Branch("ele2RegressionEnergyErrorV2",&fLeptonPairPhotonEvent->ele2RegressionEnergyErrorV2,"ele2RegressionEnergyErrorV2/F");

  // scale/res corrected 
  ZgllTuple->Branch("ele2energyCorr",&fLeptonPairPhotonEvent->ele2energyCorr,"ele2energyCorr/F");
  ZgllTuple->Branch("ele2pxCorr",&fLeptonPairPhotonEvent->ele2pxCorr,"ele2pxCorr/F");
  ZgllTuple->Branch("ele2pyCorr",&fLeptonPairPhotonEvent->ele2pyCorr,"ele2pyCorr/F");
  ZgllTuple->Branch("ele2pzCorr",&fLeptonPairPhotonEvent->ele2pzCorr,"ele2pzCorr/F");
  ZgllTuple->Branch("ele2ptCorr",&fLeptonPairPhotonEvent->ele2ptCorr,"ele2ptCorr/F");

  ZgllTuple->Branch("ele1dEtaIn",&fLeptonPairPhotonEvent->ele1dEtaIn,"ele1dEtaIn/F");
  ZgllTuple->Branch("ele1dPhiIn",&fLeptonPairPhotonEvent->ele1dPhiIn,"ele1dPhiIn/F");
  ZgllTuple->Branch("ele1sigmaIEtaIEta",&fLeptonPairPhotonEvent->ele1sigmaIEtaIEta, "ele1sigmaIEtaIEta/F");
  ZgllTuple->Branch("ele1HadOverEm",&fLeptonPairPhotonEvent->ele1HadOverEm,"ele1HadOverEm/F");
  ZgllTuple->Branch("ele1D0",&fLeptonPairPhotonEvent->ele1D0,"ele1D0/F");
  ZgllTuple->Branch("ele1DZ",&fLeptonPairPhotonEvent->ele1DZ,"ele1DZ/F");
  ZgllTuple->Branch("ele1OneOverEMinusOneOverP",&fLeptonPairPhotonEvent->ele1OneOverEMinusOneOverP,"ele1OneOverEMinusOneOverP/F");
  ZgllTuple->Branch("ele1PFIsoOverPt",&fLeptonPairPhotonEvent->ele1PFIsoOverPt,"ele1PFIsoOverPt/F");
  ZgllTuple->Branch("ele1Conversion",&fLeptonPairPhotonEvent->ele1Conversion,"ele1Conversion/O");
  ZgllTuple->Branch("ele1missinghits",&fLeptonPairPhotonEvent->ele1missinghits,"ele1missinghits/F");

  ZgllTuple->Branch("ele2dEtaIn",&fLeptonPairPhotonEvent->ele2dEtaIn,"ele2dEtaIn/F");
  ZgllTuple->Branch("ele2dPhiIn",&fLeptonPairPhotonEvent->ele2dPhiIn,"ele2dPhiIn/F");
  ZgllTuple->Branch("ele2sigmaIEtaIEta",&fLeptonPairPhotonEvent->ele2sigmaIEtaIEta, "ele2sigmaIEtaIEta/F");
  ZgllTuple->Branch("ele2HadOverEm",&fLeptonPairPhotonEvent->ele2HadOverEm,"ele2HadOverEm/F");
  ZgllTuple->Branch("ele2D0",&fLeptonPairPhotonEvent->ele2D0,"ele2D0/F");
  ZgllTuple->Branch("ele2DZ",&fLeptonPairPhotonEvent->ele2DZ,"ele2DZ/F");
  ZgllTuple->Branch("ele2OneOverEMinusOneOverP",&fLeptonPairPhotonEvent->ele2OneOverEMinusOneOverP,"ele2OneOverEMinusOneOverP/F");
  ZgllTuple->Branch("ele2PFIsoOverPt",&fLeptonPairPhotonEvent->ele2PFIsoOverPt,"ele2PFIsoOverPt/F");
  ZgllTuple->Branch("ele2Conversion",&fLeptonPairPhotonEvent->ele2Conversion,"ele2Conversion/O");
  ZgllTuple->Branch("ele2missinghits",&fLeptonPairPhotonEvent->ele2missinghits,"ele2missinghits/F");

  ZgllTuple->Branch("chargediso_ele1",&fLeptonPairPhotonEvent->chargediso_ele1,"chargediso_ele1/F");
  ZgllTuple->Branch("gammaiso_ele1",&fLeptonPairPhotonEvent->gammaiso_ele1,"gammaiso_ele1/F");
  ZgllTuple->Branch("neutraliso_ele1",&fLeptonPairPhotonEvent->neutraliso_ele1,"neutraliso_ele1/F");
  ZgllTuple->Branch("rho_ele1",&fLeptonPairPhotonEvent->rho_ele1,"rho_ele1/F");
  ZgllTuple->Branch("effectivearea_ele1",&fLeptonPairPhotonEvent->effectivearea_ele1,"effectivearea_ele1/F");
  ZgllTuple->Branch("chargediso_ele2",&fLeptonPairPhotonEvent->chargediso_ele2,"chargediso_ele2/F");
  ZgllTuple->Branch("gammaiso_ele2",&fLeptonPairPhotonEvent->gammaiso_ele2,"gammaiso_ele2/F");
  ZgllTuple->Branch("neutraliso_ele2",&fLeptonPairPhotonEvent->neutraliso_ele2,"neutraliso_ele2/F");
  ZgllTuple->Branch("rho_ele2",&fLeptonPairPhotonEvent->rho_ele2,"rho_ele2/F");
  ZgllTuple->Branch("effectivearea_ele2",&fLeptonPairPhotonEvent->effectivearea_ele2,"effectivearea_ele2/F");

  ZgllTuple->Branch("costheta_lm_electrons",&fLeptonPairPhotonEvent->costheta_lm_electrons,"costheta_lm_electrons/F");
  ZgllTuple->Branch("costheta_lp_electrons",&fLeptonPairPhotonEvent->costheta_lp_electrons,"costheta_lp_electrons/F");
  ZgllTuple->Branch("phi_electrons",&fLeptonPairPhotonEvent->phi_electrons,"phi_electrons/F");
  ZgllTuple->Branch("cosTheta_electrons",&fLeptonPairPhotonEvent->cosTheta_electrons,"cosTheta_electrons/F");
  ZgllTuple->Branch("cosThetaG_electrons",&fLeptonPairPhotonEvent->cosThetaG_electrons,"cosThetaG_electrons");
  ZgllTuple->Branch("costheta_lm_muons",&fLeptonPairPhotonEvent->costheta_lm_muons,"costheta_lm_muons/F");
  ZgllTuple->Branch("costheta_lp_muons",&fLeptonPairPhotonEvent->costheta_lp_muons,"costheta_lp_muons/F");
  ZgllTuple->Branch("phi_muons",&fLeptonPairPhotonEvent->phi_muons,"phi_muons/F");
  ZgllTuple->Branch("cosTheta_muons",&fLeptonPairPhotonEvent->cosTheta_muons,"cosTheta_muons/F");
  ZgllTuple->Branch("cosThetaG_muons",&fLeptonPairPhotonEvent->cosThetaG_muons,"cosThetaG_muons/F");

  ZgllTuple->Branch("muonZgVeto",&fLeptonPairPhotonEvent->muonZgVeto,"muonZgVeto/O");
  ZgllTuple->Branch("muonZmass",&fLeptonPairPhotonEvent->muonZmass,"muonZmass/F");
  ZgllTuple->Branch("m1E",&fLeptonPairPhotonEvent->m1E,"m1E/F");
  ZgllTuple->Branch("m1Pt",&fLeptonPairPhotonEvent->m1Pt,"m1Pt/F");
  ZgllTuple->Branch("m1Mass",&fLeptonPairPhotonEvent->m1Mass,"m1Mass/F");
  ZgllTuple->Branch("m1Px",&fLeptonPairPhotonEvent->m1Px,"m1Px/F");
  ZgllTuple->Branch("m1Py",&fLeptonPairPhotonEvent->m1Py,"m1Py/F");
  ZgllTuple->Branch("m1Pz",&fLeptonPairPhotonEvent->m1Pz,"m1Pz/F");
  ZgllTuple->Branch("m1Eta",&fLeptonPairPhotonEvent->m1Eta,"m1Eta/F");
  ZgllTuple->Branch("m1Phi",&fLeptonPairPhotonEvent->m1Phi,"m1Phi/F");
  ZgllTuple->Branch("m1Charge",&fLeptonPairPhotonEvent->m1Charge,"m1Charge/F");
  ZgllTuple->Branch("m1PtErr",&fLeptonPairPhotonEvent->m1PtErr,"m1PtErr/F");
  ZgllTuple->Branch("m2E",&fLeptonPairPhotonEvent->m2E,"m2E/F");
  ZgllTuple->Branch("m2Pt",&fLeptonPairPhotonEvent->m2Pt,"m2Pt/F");
  ZgllTuple->Branch("m2Mass",&fLeptonPairPhotonEvent->m2Mass,"m2Mass/F");
  ZgllTuple->Branch("m2Px",&fLeptonPairPhotonEvent->m2Px,"m2Px/F");
  ZgllTuple->Branch("m2Py",&fLeptonPairPhotonEvent->m2Py,"m2Py/F");
  ZgllTuple->Branch("m2Pz",&fLeptonPairPhotonEvent->m2Pz,"m2Pz/F");
  ZgllTuple->Branch("m2Eta",&fLeptonPairPhotonEvent->m2Eta,"m2Eta/F");
  ZgllTuple->Branch("m2Phi",&fLeptonPairPhotonEvent->m2Phi,"m2Phi/F");
  ZgllTuple->Branch("m2Charge",&fLeptonPairPhotonEvent->m2Charge,"m2Charge/F");
  ZgllTuple->Branch("m2PtErr",&fLeptonPairPhotonEvent->m2PtErr,"m2PtErr/F");

  ZgllTuple->Branch("m1ECorr",&fLeptonPairPhotonEvent->m1ECorr,"m1ECorr/F");
  ZgllTuple->Branch("m1PtCorr",&fLeptonPairPhotonEvent->m1PtCorr,"m1PtCorr/F");
  ZgllTuple->Branch("m1PxCorr",&fLeptonPairPhotonEvent->m1PxCorr,"m1PxCorr/F");
  ZgllTuple->Branch("m1PyCorr",&fLeptonPairPhotonEvent->m1PyCorr,"m1PyCorr/F");
  ZgllTuple->Branch("m1PzCorr",&fLeptonPairPhotonEvent->m1PzCorr,"m1PzCorr/F");

  ZgllTuple->Branch("m2ECorr",&fLeptonPairPhotonEvent->m2ECorr,"m2ECorr/F");
  ZgllTuple->Branch("m2PtCorr",&fLeptonPairPhotonEvent->m2PtCorr,"m2PtCorr/F");
  ZgllTuple->Branch("m2PxCorr",&fLeptonPairPhotonEvent->m2PxCorr,"m2PxCorr/F");
  ZgllTuple->Branch("m2PyCorr",&fLeptonPairPhotonEvent->m2PyCorr,"m2PyCorr/F");
  ZgllTuple->Branch("m2PzCorr",&fLeptonPairPhotonEvent->m2PzCorr,"m2PzCorr/F");

  ZgllTuple->Branch("photonidmva",&fLeptonPairPhotonEvent->photonidmva,"photonidmva/F");
  ZgllTuple->Branch("photonr9",&fLeptonPairPhotonEvent->photonr9,"photonr9/F");
  ZgllTuple->Branch("photonenergy",&fLeptonPairPhotonEvent->photonenergy,"photonenergy/F");
  ZgllTuple->Branch("photonpx",&fLeptonPairPhotonEvent->photonpx,"photonpx/F");
  ZgllTuple->Branch("photonpy",&fLeptonPairPhotonEvent->photonpy,"photonpy/F");
  ZgllTuple->Branch("photonpz",&fLeptonPairPhotonEvent->photonpz,"photonpz/F");
  ZgllTuple->Branch("photonpt",&fLeptonPairPhotonEvent->photonpt,"photonpt/F");
  ZgllTuple->Branch("photoneta",&fLeptonPairPhotonEvent->photoneta,"photoneta/F");
  ZgllTuple->Branch("photonmass",&fLeptonPairPhotonEvent->photonmass,"photonmass/F");
  ZgllTuple->Branch("photonphi",&fLeptonPairPhotonEvent->photonphi,"photonphi/F");
  ZgllTuple->Branch("photonenergyerror",&fLeptonPairPhotonEvent->photonenergyerror,"photonenergyerror/F");
 
  ZgllTuple->Branch("photonenergyCorr",&fLeptonPairPhotonEvent->photonenergyCorr,"photonenergyCorr/F");
  ZgllTuple->Branch("photonpxCorr",&fLeptonPairPhotonEvent->photonpxCorr,"photonpxCorr/F");
  ZgllTuple->Branch("photonpyCorr",&fLeptonPairPhotonEvent->photonpyCorr,"photonpyCorr/F");
  ZgllTuple->Branch("photonpzCorr",&fLeptonPairPhotonEvent->photonpzCorr,"photonpzCorr/F");
  ZgllTuple->Branch("photonptCorr",&fLeptonPairPhotonEvent->photonptCorr,"photonptCorr/F");

  ZgllTuple->Branch("NPu",&fLeptonPairPhotonEvent->NPu,"NPu/F");
  ZgllTuple->Branch("NPuPlus",&fLeptonPairPhotonEvent->NPuPlus,"NPuPlus/F");
  ZgllTuple->Branch("NPuMinus",&fLeptonPairPhotonEvent->NPuMinus,"NPuMinus/F");
 
  ZgllTuple->Branch("photonmatchmc",&fLeptonPairPhotonEvent->photonmatchmc,"photonmatchmc/F");
  //Initialize Electron MVA
  //eleIDMVA = new mithep::ElectronIDMVA();

  // vector<string> weightFiles;
  // weightFiles.push_back("/home/ksingh/cms/cmssw/026/CMSSW_5_2_3/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml");
  // weightFiles.push_back("/home/ksingh/cms/cmssw/026/CMSSW_5_2_3/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml");
  // weightFiles.push_back("/home/ksingh/cms/cmssw/026/CMSSW_5_2_3/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml");
  // weightFiles.push_back("/home/ksingh/cms/cmssw/026/CMSSW_5_2_3/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml");
  // weightFiles.push_back("/home/ksingh/cms/cmssw/026/CMSSW_5_2_3/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml");
  // weightFiles.push_back("/home/ksingh/cms/cmssw/026/CMSSW_5_2_3/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml");
  // eleIDMVA->Initialize( "ElectronIDMVA",
  //                      mithep::ElectronIDMVA::kIDEGamma2012NonTrigV0,
  //                      kTRUE, weightFiles);
  //fTool.InitializeMVA(fVariableType_2011,fEndcapWeights_2011,fBarrelWeights_2011);
  

  //Set up Electron Regression Evalatuator
  eleRegressionEvaluator_V0 = new ElectronEnergyRegression();
  eleRegressionEvaluator_V1 = new ElectronEnergyRegression();
  eleRegressionEvaluator_V2 = new ElectronEnergyRegression();

  eleRegressionEvaluator_V0->initialize ("/afs/cern.ch/user/s/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronRegressionWeights/weightFile_V00.root",
                                         ElectronEnergyRegression::kNoTrkVar);
  eleRegressionEvaluator_V1->initialize   ("/afs/cern.ch/user/s/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronRegressionWeights/weightFile_V01.root",
                                           ElectronEnergyRegression::kWithTrkVar);
  eleRegressionEvaluator_V2->initialize   ("/afs/cern.ch/user/s/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronRegressionWeights/weightFile_V02.root",
                                           ElectronEnergyRegression::kWithTrkVar);
  assert(eleRegressionEvaluator_V0->isInitialized());
  assert(eleRegressionEvaluator_V1->isInitialized());
  assert(eleRegressionEvaluator_V2->isInitialized());


  //Add Output Tree
  AddOutput(ZgllTuple);

}


