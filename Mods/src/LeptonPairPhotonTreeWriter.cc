//--------------------------------------------------------------------------------------------------
// $Id: LeptonPairPhotonTreeWriter.h,v 1.0 2012/06/23 21:25:01 auhess ksingh
//
// LeptonPairPhotonTreeWriter
//
// Authors: A. Hess & K. Singh
//--------------------------------------------------------------------------------------------------
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
#include <TNtuple.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

using namespace mithep;
mithep::ElectronIDMVA *eleIDMVA; // The electron MVA 
MVATools fTool; // The electron MVA tools

ClassImp(mithep::LeptonPairPhotonTreeWriter)
  ClassImp(mithep::LeptonPairPhotonEvent)
//--------------------------------------------------------------------------------------------------
  LeptonPairPhotonTreeWriter::LeptonPairPhotonTreeWriter(const char *name, const char *title) : 
    BaseMod                 (name,title),
    // define all the Branches to load
    fPhotonBranchName       (Names::gkPhotonBrn),
    fGoodElectronName       (Names::gkElectronBrn),  
    fGoodMuonName	    (Names::gkMuonBrn),
    fPVName                 (Names::gkPVBeamSpotBrn), 
    fVertexName             (ModNames::gkGoodVertexesName),
    fPFCandName		    (Names::gkPFCandidatesBrn),
    fTrackName              (Names::gkTrackBrn),
    fPileUpDenName          (Names::gkPileupEnergyDensityBrn),
    fPileUpName             (Names::gkPileupInfoBrn),
  
    fIsData                 (false),
    fElectronMuonFlag	    (1), // 1 if for electrons, 2 if for muons

    fDoBlinding             (false),

    fPhotonsFromBranch      (kTRUE),  
    fPVFromBranch           (kTRUE),
    fGoodElectronsFromBranch(kTRUE),
    fGoodMuonsFromBranch    (kTRUE),
    fPFCandsFromBranch      (kTRUE),
    fTracksFromBranch       (kTRUE),
    fPileUpDenFromBranch    (kTRUE),

    // ----------------------------------------
    // collections....
    fPhotons                (NULL),
    fGoodElectrons          (NULL),
    fPV                     (NULL),
    fGoodMuons		    (NULL),
    fVertices               (NULL),
    fPFCands		    (NULL),
    fTracks		    (NULL),
    fPileUpDen              (NULL),
    fPileUp                 (0),
  
    fTupleName              ("h2LepPhotonTree"),
  
    //Photon MVA Variables
    fVariableType_2011             (10),
    fEndcapWeights_2011            (gSystem->Getenv("CMSSW_BASE")+
				    TString("/src/MitPhysics/data/TMVAClassificationPhotonID_")+
				    TString("Endcap_PassPreSel_Variable_10_BDTnCuts2000_BDT.")+
				    TString("weights.xml")),
    fBarrelWeights_2011            (gSystem->Getenv("CMSSW_BASE")+
				    TString("/src/MitPhysics/data/TMVAClassificationPhotonID_")+
				    TString("Barrel_PassPreSel_Variable_10_BDTnCuts2000_BDT.")+
				    TString("weights.xml")),
    //2012 Photon MVA Variables not currently used
    fVariableType_2012_globe       (1201),
    fEndcapWeights_2012_globe      (gSystem->Getenv("CMSSW_BASE")+
				    TString("/src/MitPhysics/data/")+
				    TString("TMVA_EEpf_BDT_globe.")+
				    TString("weights.xml")),
    fBarrelWeights_2012_globe      (gSystem->Getenv("CMSSW_BASE")+
				    TString("/src/MitPhysics/data/")+
				    TString("TMVA_EBpf_BDT_globe.")+
				    TString("weights.xml"))

{
  // Constructor
}

LeptonPairPhotonTreeWriter::~LeptonPairPhotonTreeWriter()
{
  // Destructor
}

//--------------------------------------------------------------------------------------------------
void LeptonPairPhotonTreeWriter::Process()
{
  // ------------------------------------------------------------  
  // Process entries of the tree. 
  LoadEventObject(fPhotonBranchName,   fPhotons);
  LoadEventObject(fGoodElectronName,   fGoodElectrons);
  LoadEventObject(fPVName,             fPV);
  LoadEventObject(fPFCandName,         fPFCands);
  LoadEventObject(fTrackName,          fTracks);
  LoadEventObject(fPileUpDenName,      fPileUpDen);
  LoadEventObject(fGoodMuonName,       fGoodMuons);

  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);

  if (!fIsData){
    LoadBranch(fPileUpName);
  }

  //Initialize all tree leaf entries to -99
  fLeptonPairPhotonEvent->Meeg = -99.;
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
  
  fLeptonPairPhotonEvent->ele2charge = -99.;
  fLeptonPairPhotonEvent->ele2energy = -99.;
  fLeptonPairPhotonEvent->ele2px = -99.;
  fLeptonPairPhotonEvent->ele2py = -99.;
  fLeptonPairPhotonEvent->ele2pz = -99.;
  fLeptonPairPhotonEvent->ele2pt = -99.;
  fLeptonPairPhotonEvent->ele2eta = -99.;
  fLeptonPairPhotonEvent->ele2mass = -99.;
  fLeptonPairPhotonEvent->ele2phi = -99.;

  fLeptonPairPhotonEvent->Mmmg = -99;
  fLeptonPairPhotonEvent->m1E = -99;
  fLeptonPairPhotonEvent->m1Pt = -99;
  fLeptonPairPhotonEvent->m1Mass = -99;
  fLeptonPairPhotonEvent->m1Px = -99; 
  fLeptonPairPhotonEvent->m1Py = -99;
  fLeptonPairPhotonEvent->m1Pz = -99;
  fLeptonPairPhotonEvent->m1Eta = -99;
  fLeptonPairPhotonEvent->m1Phi = -99;
  fLeptonPairPhotonEvent->m1Charge = -99;
  fLeptonPairPhotonEvent->m2E = -99;
  fLeptonPairPhotonEvent->m2Pt = -99;
  fLeptonPairPhotonEvent->m2Mass = -99;
  fLeptonPairPhotonEvent->m2Px = -99;
  fLeptonPairPhotonEvent->m2Py = -99;
  fLeptonPairPhotonEvent->m2Pz = -99;
  fLeptonPairPhotonEvent->m2Eta = -99;
  fLeptonPairPhotonEvent->m2Phi = -99;
  fLeptonPairPhotonEvent->m2Charge = -99;

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

  Int_t _numPU      = -1.;        // some sensible default values....
  Int_t _numPUminus = -1.;        // some sensible default values....
  Int_t _numPUplus  = -1.; 
    
  if (!fIsData){  
    for (UInt_t i=0; i<fPileUp->GetEntries(); ++i) {
      const PileupInfo *puinfo = fPileUp->At(i);
      if      (puinfo->GetBunchCrossing() ==  0) _numPU      = puinfo->GetPU_NumInteractions();
      else if (puinfo->GetBunchCrossing() == -1) _numPUminus = puinfo->GetPU_NumInteractions();
      else if (puinfo->GetBunchCrossing() ==  1) _numPUplus  = puinfo->GetPU_NumInteractions();
    }
  }
  
  fLeptonPairPhotonEvent->NPu      = _numPU;
  fLeptonPairPhotonEvent->NPuMinus = _numPUminus;
  fLeptonPairPhotonEvent->NPuPlus  = _numPUplus;
  
  //This for loop computes dielectron quantities
  if (fElectronMuonFlag == 1){
    if (fGoodElectrons->GetEntries() > 1 && fPhotons->GetEntries() > 0){
      
      const Photon *pho;
      pho = 0;
      pho = fPhotons->At(0);
      
      //Compute Photon MVA Value and photon quantities
      Float_t rho = -99.;
      if (fPileUpDen->GetEntries() > 0)
	rho = (Double_t) fPileUpDen->At(0)->RhoRandomLowEta();
      
      fLeptonPairPhotonEvent->photonidmva = fTool.GetMVAbdtValue_2011(pho,fPV->At(0),fTracks,fPV,rho,fGoodElectrons,kTRUE);
      fLeptonPairPhotonEvent->photonenergy = pho->E();
      fLeptonPairPhotonEvent->photonpx = pho->Px();
      fLeptonPairPhotonEvent->photonpy = pho->Py();
      fLeptonPairPhotonEvent->photonpz = pho->Pz();
      fLeptonPairPhotonEvent->photonpt = pho->Pt();
      fLeptonPairPhotonEvent->photoneta = pho->Eta();
      fLeptonPairPhotonEvent->photonmass = pho->Mass();
      fLeptonPairPhotonEvent->photonphi = pho->Phi();
      fLeptonPairPhotonEvent->photonr9 = pho->R9();
  
      //Compute Dielectron Quantities (Currently takes the two highest Pt electrons
      // and the highest Pt photon
      const Electron *ele1;
      const Electron *ele2;
      ele1 = 0;
      ele2 = 0;
  
      ele1 = fGoodElectrons->At(0);
      ele2 = fGoodElectrons->At(1);
      if(ele1->Charge() != ele2->Charge()){
  
	fLeptonPairPhotonEvent->ele1charge = ele1->Charge();
	fLeptonPairPhotonEvent->ele1energy = ele1->E();
	fLeptonPairPhotonEvent->ele1px = ele1->Px();
	fLeptonPairPhotonEvent->ele1py = ele1->Py();
	fLeptonPairPhotonEvent->ele1pz = ele1->Pz();
	fLeptonPairPhotonEvent->ele1pt = ele1->Pt();
	fLeptonPairPhotonEvent->ele1eta = ele1->Eta();
	fLeptonPairPhotonEvent->ele1mass = ele1->Mass();
	fLeptonPairPhotonEvent->ele1phi = ele1->Phi();

	fLeptonPairPhotonEvent->ele2charge = ele2->Charge();
	fLeptonPairPhotonEvent->ele2energy = ele2->E();
	fLeptonPairPhotonEvent->ele2px = ele2->Px();
	fLeptonPairPhotonEvent->ele2py = ele2->Py();
	fLeptonPairPhotonEvent->ele2pz = ele2->Pz();
	fLeptonPairPhotonEvent->ele2pt = ele2->Pt();
	fLeptonPairPhotonEvent->ele2eta = ele2->Eta();
	fLeptonPairPhotonEvent->ele2mass = ele2->Mass();
	fLeptonPairPhotonEvent->ele2phi = ele2->Phi();

	fLeptonPairPhotonEvent->Meeg = (ele1->Mom() + ele2->Mom() + pho->Mom()).M();

	//Compute Electron MVA Values
	double _fbrem = max(double(ele1->FBrem()),-1.0);
	double _kftrk_chisq = (ele1->HasTrackerTrk()) ? ele1->TrackerTrk()->Chi2() / ele1->TrackerTrk()->Ndof() : 0;
	double _kftrk_nhits = (ele1->HasTrackerTrk()) ? (Float_t)ele1->TrackerTrk()->NHits() : -1; 
	double _gsftrk_chisq = (Double_t)min(double(ele1->BestTrk()->Chi2() / ele1->BestTrk()->Ndof()),200.0);
	double seedE1x5OverE = ele1->SCluster()->Seed()->E1x5() / ele1->SCluster()->Seed()->Energy();
	double seedE5x5OverE = ele1->SCluster()->Seed()->E5x5() / ele1->SCluster()->Seed()->Energy();
	double _e1x5e5x5 = min(max(1 - double(seedE1x5OverE/seedE5x5OverE),-1.0),2.0);
	double _r9 = min(double(ele1->SCluster()->R9()),5.0);
	double _e_o_p = min(double(ele1->ESuperClusterOverP()), 20.0);
	double _eseed_o_pout = min(double(ele1->ESeedClusterOverPout()),20.0);
	double _IoEmIoP =  (Double_t)(1 - ele1->ESuperClusterOverP())/(ele1->SCluster()->Et()*TMath::CosH(ele1->SCluster()->Eta()));
	double _epreoraw = (Float_t)ele1->SCluster()->PreshowerEnergy() / ele1->SCluster()->RawEnergy();
  
	fLeptonPairPhotonEvent->ele1MVA = eleIDMVA->MVAValue_IDNonTrig(ele1->Pt(),
								       ele1->SCluster()->Eta(),
								       _fbrem,
								       _kftrk_chisq,
								       _kftrk_nhits,
								       _gsftrk_chisq,
								       fabs(ele1->DeltaEtaSuperClusterTrackAtVtx()),
								       fabs(ele1->DeltaPhiSuperClusterTrackAtVtx()),
								       fabs(ele1->DeltaEtaSeedClusterTrackAtCalo()),
								       ele1->CoviEtaiEta(),
								       sqrt(ele1->SCluster()->Seed()->CoviPhiiPhi()),
								       ele1->SCluster()->EtaWidth(), 
								       ele1->SCluster()->PhiWidth(),
								       _e1x5e5x5,
								       _r9,
								       ele1->HadronicOverEm(),
								       _e_o_p,
								       _IoEmIoP,
								       _eseed_o_pout,
								       _epreoraw,
								       kFALSE );


	double _fbrem_2 = max(double(ele2->FBrem()),-1.0);
	double _kftrk_chisq_2 = (ele2->HasTrackerTrk()) ? ele2->TrackerTrk()->Chi2() / ele2->TrackerTrk()->Ndof() : 0;
	double _kftrk_nhits_2 = (ele2->HasTrackerTrk()) ? (Float_t)ele2->TrackerTrk()->NHits() : -1; 
	double _gsftrk_chisq_2 = (Double_t)min(double(ele2->BestTrk()->Chi2() / ele2->BestTrk()->Ndof()),200.0);
	double seedE1x5OverE_2 = ele2->SCluster()->Seed()->E1x5() / ele2->SCluster()->Seed()->Energy();
	double seedE5x5OverE_2 = ele2->SCluster()->Seed()->E5x5() / ele2->SCluster()->Seed()->Energy();
	double _e1x5e5x5_2 = min(max(1 - double(seedE1x5OverE_2/seedE5x5OverE_2),-1.0),2.0);
	double _r9_2 = min(double(ele2->SCluster()->R9()),5.0);
	double _e_o_p_2 = min(double(ele2->ESuperClusterOverP()), 20.0);
	double _eseed_o_pout_2 = min(double(ele2->ESeedClusterOverPout()),20.0);
	double _IoEmIoP_2 =  (Double_t)(1 - ele2->ESuperClusterOverP())/(ele2->SCluster()->Et()*TMath::CosH(ele2->SCluster()->Eta()));
	double _epreoraw_2 = (Float_t)ele2->SCluster()->PreshowerEnergy() / ele2->SCluster()->RawEnergy();
  
  
	fLeptonPairPhotonEvent->ele2MVA = eleIDMVA->MVAValue_IDNonTrig(ele2->Pt(),
								       ele2->SCluster()->Eta(),
								       _fbrem_2,
								       _kftrk_chisq_2,
								       _kftrk_nhits_2,
								       _gsftrk_chisq_2,
								       fabs(ele2->DeltaEtaSuperClusterTrackAtVtx()),
								       fabs(ele2->DeltaPhiSuperClusterTrackAtVtx()),
								       fabs(ele2->DeltaEtaSeedClusterTrackAtCalo()),
								       ele2->CoviEtaiEta(),
								       sqrt(ele2->SCluster()->Seed()->CoviPhiiPhi()),
								       ele2->SCluster()->EtaWidth(), 
								       ele2->SCluster()->PhiWidth(),
								       _e1x5e5x5_2,
								       _r9_2,
								       ele2->HadronicOverEm(),
								       _e_o_p_2,
								       _IoEmIoP_2,
								       _eseed_o_pout_2,
								       _epreoraw_2,
								       kFALSE );

      }
    }//Set Blinding Range
    if ( fDoBlinding && !fIsData &&  (fLeptonPairPhotonEvent->Meeg < 110 || fLeptonPairPhotonEvent->Meeg > 140) ){
      ZgllTuple->Fill();
    }
    
  }

  //Calculations for muons
  if (fElectronMuonFlag == 2){
    if (fGoodMuons->GetEntries()>1 && fPhotons->GetEntries() > 0){
      //Takes first (highest Pt) photon
      const Photon *pho;
      pho = 0;
      pho = fPhotons->At(0);
      //Photon MVA again
      Float_t rho = -99.;
      if (fPileUpDen->GetEntries() > 0)
	rho = (Double_t) fPileUpDen->At(0)->RhoRandomLowEta();

      fLeptonPairPhotonEvent->photonidmva = fTool.GetMVAbdtValue_2011(pho,fPV->At(0),fTracks,fPV,rho,fGoodElectrons,kTRUE);
      fLeptonPairPhotonEvent->photonenergy = pho->E();
      fLeptonPairPhotonEvent->photonpx = pho->Px();
      fLeptonPairPhotonEvent->photonpy = pho->Py();
      fLeptonPairPhotonEvent->photonpz = pho->Pz();
      fLeptonPairPhotonEvent->photonpt = pho->Pt();
      fLeptonPairPhotonEvent->photoneta = pho->Eta();
      fLeptonPairPhotonEvent->photonmass = pho->Mass();
      fLeptonPairPhotonEvent->photonphi = pho->Phi();
      fLeptonPairPhotonEvent->photonr9 = pho->R9();

      //Take 2 highest Pt muons 
      const Muon *muona = NULL;
      const Muon *muonb = NULL;
      muona = fGoodMuons->At(0);
      muonb = fGoodMuons->At(1);
 
      if(muona->Charge() != muonb->Charge()){  
	fLeptonPairPhotonEvent->Mmmg = (pho->Mom()+muona->Mom()+muonb->Mom()).M();
 
	fLeptonPairPhotonEvent->m1E = muona->E();
	fLeptonPairPhotonEvent->m1Pt = muona->Pt();
	fLeptonPairPhotonEvent->m1Mass = muona->Mass();
	fLeptonPairPhotonEvent->m1Px = muona->Px();
	fLeptonPairPhotonEvent->m1Py = muona->Py();
	fLeptonPairPhotonEvent->m1Pz = muona->Pz();
	fLeptonPairPhotonEvent->m1Eta = muona->Eta();
	fLeptonPairPhotonEvent->m1Phi = muona->Phi();
	fLeptonPairPhotonEvent->m1Charge = muona->Charge();
	fLeptonPairPhotonEvent->m2E = muonb->E();
	fLeptonPairPhotonEvent->m2Pt = muonb->Pt();
	fLeptonPairPhotonEvent->m2Mass = muonb->Mass();
	fLeptonPairPhotonEvent->m2Px = muonb->Px();
	fLeptonPairPhotonEvent->m2Py = muonb->Py();
	fLeptonPairPhotonEvent->m2Pz = muonb->Pz();
	fLeptonPairPhotonEvent->m2Eta = muonb->Eta();
	fLeptonPairPhotonEvent->m2Phi = muonb->Phi();
	fLeptonPairPhotonEvent->m2Charge = muonb->Charge();

      }
    }
    if ( fDoBlinding && !fIsData && (fLeptonPairPhotonEvent->Mmmg < 110 || fLeptonPairPhotonEvent->Mmmg > 140) ){
      ZgllTuple->Fill();
    }
  }//ElectronsMuons
  
}



//--------------------------------------------------------------------------------------------------
void LeptonPairPhotonTreeWriter::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the branch collections.

  ReqEventObject(fPhotonBranchName,fPhotons,      fPhotonsFromBranch);
  ReqEventObject(fGoodElectronName,fGoodElectrons,fGoodElectronsFromBranch);  
  ReqEventObject(fPVName,             fPV, fPVFromBranch);
  ReqEventObject(fPFCandName,         fPFCands, fPFCandsFromBranch);
  ReqEventObject(fTrackName,    fTracks, fTracksFromBranch);
  ReqEventObject(fPileUpDenName,      fPileUpDen, fPileUpDenFromBranch);
  ReqEventObject(fGoodMuonName,    fGoodMuons,    fGoodMuonsFromBranch); 
  if (!fIsData){
    ReqBranch(fPileUpName,         fPileUp);
  }
 
  fLeptonPairPhotonEvent = new LeptonPairPhotonEvent;//Declares all tree leaves  
  ZgllTuple = new TTree(fTupleName.Data(),fTupleName.Data());//Declares tree name  
  ZgllTuple->SetAutoSave(300e9);
  
  ZgllTuple->Branch("Meeg",&fLeptonPairPhotonEvent->Meeg,"F");
  ZgllTuple->Branch("ele1MVA",&fLeptonPairPhotonEvent->ele1MVA,"F");
  ZgllTuple->Branch("ele2MVA",&fLeptonPairPhotonEvent->ele2MVA,"F");
  
  ZgllTuple->Branch("ele1charge",&fLeptonPairPhotonEvent->ele1charge,"F");
  ZgllTuple->Branch("ele1energy",&fLeptonPairPhotonEvent->ele1energy,"F");
  ZgllTuple->Branch("ele1px",&fLeptonPairPhotonEvent->ele1px,"F");
  ZgllTuple->Branch("ele1py",&fLeptonPairPhotonEvent->ele1py,"F");
  ZgllTuple->Branch("ele1pz",&fLeptonPairPhotonEvent->ele1pz,"F");
  ZgllTuple->Branch("ele1pt",&fLeptonPairPhotonEvent->ele1pt,"F");
  ZgllTuple->Branch("ele1eta",&fLeptonPairPhotonEvent->ele1eta,"F");
  ZgllTuple->Branch("ele1mass",&fLeptonPairPhotonEvent->ele1mass,"F");
  ZgllTuple->Branch("ele1phi",&fLeptonPairPhotonEvent->ele1phi,"F");

  ZgllTuple->Branch("ele2charge",&fLeptonPairPhotonEvent->ele2charge,"F");
  ZgllTuple->Branch("ele2energy",&fLeptonPairPhotonEvent->ele2energy,"F");
  ZgllTuple->Branch("ele2px",&fLeptonPairPhotonEvent->ele2px,"F");
  ZgllTuple->Branch("ele2py",&fLeptonPairPhotonEvent->ele2py,"F");
  ZgllTuple->Branch("ele2pz",&fLeptonPairPhotonEvent->ele2pz,"F");
  ZgllTuple->Branch("ele2pt",&fLeptonPairPhotonEvent->ele2pt,"F");
  ZgllTuple->Branch("ele2eta",&fLeptonPairPhotonEvent->ele2eta,"F");
  ZgllTuple->Branch("ele2mass",&fLeptonPairPhotonEvent->ele2mass,"F");
  ZgllTuple->Branch("ele2phi",&fLeptonPairPhotonEvent->ele2phi,"F");

  ZgllTuple->Branch("Mmmg",&fLeptonPairPhotonEvent->Mmmg,"F");
  ZgllTuple->Branch("m1E",&fLeptonPairPhotonEvent->m1E,"F");
  ZgllTuple->Branch("m1Pt",&fLeptonPairPhotonEvent->m1Pt,"F");
  ZgllTuple->Branch("m1Mass",&fLeptonPairPhotonEvent->m1Mass,"F");
  ZgllTuple->Branch("m1Px",&fLeptonPairPhotonEvent->m1Px,"F");
  ZgllTuple->Branch("m1Py",&fLeptonPairPhotonEvent->m1Py,"F");
  ZgllTuple->Branch("m1Pz",&fLeptonPairPhotonEvent->m1Pz,"F");
  ZgllTuple->Branch("m1Eta",&fLeptonPairPhotonEvent->m1Eta,"F");
  ZgllTuple->Branch("m1Phi",&fLeptonPairPhotonEvent->m1Phi,"F");
  ZgllTuple->Branch("m1Charge",&fLeptonPairPhotonEvent->m1Charge,"F");
  ZgllTuple->Branch("m2E",&fLeptonPairPhotonEvent->m2E,"F");
  ZgllTuple->Branch("m2Pt",&fLeptonPairPhotonEvent->m2Pt,"F");
  ZgllTuple->Branch("m2Mass",&fLeptonPairPhotonEvent->m2Mass,"F");
  ZgllTuple->Branch("m2Px",&fLeptonPairPhotonEvent->m2Px,"F");
  ZgllTuple->Branch("m2Py",&fLeptonPairPhotonEvent->m2Py,"F");
  ZgllTuple->Branch("m2Pz",&fLeptonPairPhotonEvent->m2Pz,"F");
  ZgllTuple->Branch("m2Eta",&fLeptonPairPhotonEvent->m2Eta,"F");
  ZgllTuple->Branch("m2Phi",&fLeptonPairPhotonEvent->m2Phi,"F");
  ZgllTuple->Branch("m2Charge",&fLeptonPairPhotonEvent->m2Charge,"F");

  ZgllTuple->Branch("photonidmva",&fLeptonPairPhotonEvent->photonidmva,"F");
  ZgllTuple->Branch("photonr9",&fLeptonPairPhotonEvent->photonr9,"F");
  ZgllTuple->Branch("photonenergy",&fLeptonPairPhotonEvent->photonenergy,"F");
  ZgllTuple->Branch("photonpx",&fLeptonPairPhotonEvent->photonpx,"F");
  ZgllTuple->Branch("photonpy",&fLeptonPairPhotonEvent->photonpy,"F");
  ZgllTuple->Branch("photonpz",&fLeptonPairPhotonEvent->photonpz,"F");
  ZgllTuple->Branch("photonpt",&fLeptonPairPhotonEvent->photonpt,"F");
  ZgllTuple->Branch("photoneta",&fLeptonPairPhotonEvent->photoneta,"F");
  ZgllTuple->Branch("photonmass",&fLeptonPairPhotonEvent->photonmass,"F");
  ZgllTuple->Branch("photonphi",&fLeptonPairPhotonEvent->photonphi,"F");
 
  ZgllTuple->Branch("NPu",&fLeptonPairPhotonEvent->NPu,"I");
  ZgllTuple->Branch("NPuPlus",&fLeptonPairPhotonEvent->NPuPlus,"I");
  ZgllTuple->Branch("NPuMinus",&fLeptonPairPhotonEvent->NPuMinus,"I");
 
  //Initialize Electron MVA
  eleIDMVA = new mithep::ElectronIDMVA();

  vector<string> weightFiles;
  weightFiles.push_back("/home/ksingh/cms/cmssw/026/CMSSW_5_2_3/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml");
  weightFiles.push_back("/home/ksingh/cms/cmssw/026/CMSSW_5_2_3/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml");
  weightFiles.push_back("/home/ksingh/cms/cmssw/026/CMSSW_5_2_3/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml");
  weightFiles.push_back("/home/ksingh/cms/cmssw/026/CMSSW_5_2_3/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml");
  weightFiles.push_back("/home/ksingh/cms/cmssw/026/CMSSW_5_2_3/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml");
  weightFiles.push_back("/home/ksingh/cms/cmssw/026/CMSSW_5_2_3/src/MitPhysics/data/ElectronMVAWeights/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml");
  eleIDMVA->Initialize( "ElectronIDMVA",
                        mithep::ElectronIDMVA::kIDEGamma2012NonTrigV0,
                        kTRUE, weightFiles);
  fTool.InitializeMVA(fVariableType_2011,fEndcapWeights_2011,fBarrelWeights_2011);
  
  //Add Output Tree
  AddOutput(ZgllTuple);
}
