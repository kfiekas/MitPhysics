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
//--------------------------------------------------------------------------------------------------
  LeptonPairPhotonTreeWriter::LeptonPairPhotonTreeWriter(const char *name, const char *title) : 
    BaseMod                 (name,title),
    // define all the Branches to load
    fPhotonBranchName       (Names::gkPhotonBrn),
    fGoodElectronName       (Names::gkElectronBrn),  
    fGoodMuonName	    (Names::gkMuonBrn),
    fPVName                 (Names::gkPVBeamSpotBrn), 
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
    fTupleName              ("h2LepPhotonTree")
  
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
  LoadEventObject(fConversionName,     fConversions);
  LoadEventObject(fBeamSpotName,       fBeamSpot);  
  if (!fIsData){
    LoadBranch(fPileUpName);
    LoadBranch(fMCParticleName);
  }
  //Initialize all tree leaf entries to -99.;
  fLeptonPairPhotonEvent->electronZmass = -99.;
  fLeptonPairPhotonEvent->mllg = -99.;
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
  fLeptonPairPhotonEvent->m2E = -99.;
  fLeptonPairPhotonEvent->m2Pt = -99.;
  fLeptonPairPhotonEvent->m2Mass = -99.;
  fLeptonPairPhotonEvent->m2Px = -99.;
  fLeptonPairPhotonEvent->m2Py = -99.;
  fLeptonPairPhotonEvent->m2Pz = -99.;
  fLeptonPairPhotonEvent->m2Eta = -99.;
  fLeptonPairPhotonEvent->m2Phi = -99.;
  fLeptonPairPhotonEvent->m2Charge = -99.;

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

//This for loop computes dielectron quantities
if ((fGoodElectrons->GetEntries() > 1) || (fGoodMuons->GetEntries() > 1)){
	Bool_t electronZgfilled; 
	electronZgfilled = kFALSE;
	bool photonpass = false;
	UInt_t photonindex = 0;
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


	if (fPhotons->GetEntries() >= 1 ){
		for (UInt_t i = 0; i < fPhotons->GetEntries(); ++i){
			const Photon *pho;
			pho = fPhotons->At(i);	
			if (pho->Pt() > 10 && (fabs(pho->SCluster()->Eta()) < 1.4442 || (fabs(pho->SCluster()->Eta()) > 1.566 && fabs(pho->SCluster()->Eta()) < 2.5))){
				if (ZGTools::photonCutBasedLoose2012ID(pho, fPV->At(0), fPFCands, fPileUpDen, fGoodElectrons, fConversions,fPV) ){
					if (ZGTools::photonCutBasedLoose2012Isolation(pho, fPV, fPFCands, fPileUpDen,pfNoPileUpflag) ){
						photonpass = true;
						photonindex = i;
						break;
					}
				}
			}
		}
		
    		if (photonpass){
			const Photon *pho;
			pho = fPhotons->At(photonindex);
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

      			//Compute Photon MVA Value and photon quantities
      		//	Float_t rho = -99.;
      		//	if (fPileUpDen->GetEntries() > 0) rho = (Double_t) fPileUpDen->At(0)->RhoRandomLowEta();
//      			fLeptonPairPhotonEvent->photonidmva = fTool.GetMVAbdtValue_2011(pho,fPV->At(0),fTracks,fPV,rho,fGoodElectrons,kTRUE);
			fLeptonPairPhotonEvent->photonenergy = pho->E();
      			fLeptonPairPhotonEvent->photonpx = pho->Px();
      			fLeptonPairPhotonEvent->photonpy = pho->Py();
      			fLeptonPairPhotonEvent->photonpz = pho->Pz();
      			fLeptonPairPhotonEvent->photonpt = pho->Pt();
      			fLeptonPairPhotonEvent->photoneta = pho->Eta();
      			fLeptonPairPhotonEvent->photonmass = pho->Mass();
      			fLeptonPairPhotonEvent->photonphi = pho->Phi();
      			fLeptonPairPhotonEvent->photonr9 = pho->R9();	
    		}
   	}
   //Compute Dielectron Quantities
   if (fGoodElectrons->GetEntries() > 1){
	vector<bool> electronpass;
	for (UInt_t j = 0; j < fGoodElectrons->GetEntries();++j){
		const Electron *test;
		test = fGoodElectrons->At(j);
		if (test->Pt() > 7 && fabs(test->Eta()) < 2.5){
			if (ZGTools::electronCutBasedIDLoose(test, fPV->At(0),fConversions,YEAR)){
				if ((YEAR == 2011 && ZGTools::electronPFIso04(test,fPV->At(0),fPFCands,fPileUpDen, mithep::ElectronTools::kEleEAData2011, pfNoPileUpflag,YEAR)) || (YEAR == 2012 && ZGTools::electronPFIso04(test,fPV->At(0),fPFCands,fPileUpDen, mithep::ElectronTools::kEleEAData2012, pfNoPileUpflag,YEAR))){
					electronpass.push_back(1);
				}
				else electronpass.push_back(0);
			}
			else electronpass.push_back(0);
		}
		else electronpass.push_back(0);

	}
	const Electron *ele1;
	const Electron *ele2;
      	ele1 = 0;
      	ele2 = 0;
      	Float_t zdifference = 999;
      	UInt_t electron1 = 0;
      	UInt_t electron2 = 1;
	
	for (UInt_t j = 0; j < fGoodElectrons->GetEntries(); ++j){
		for (UInt_t k = 0; k < j; ++k){    
		if (j == k) continue;
		ele1 = fGoodElectrons->At(j);
		ele2 = fGoodElectrons->At(k);
			if (ele1->Charge() != ele2->Charge() && electronpass[j] && electronpass[k]){ 
				if (fabs(91.19 - (ele1->Mom() + ele2->Mom()).M()) < zdifference){
					zdifference = fabs(91.19 - (ele1->Mom() + ele2->Mom()).M());
					electron1 = k;
					electron2 = j;
				}
			}	
		}
	}
 	
      	ele1 = fGoodElectrons->At(electron1);
      	ele2 = fGoodElectrons->At(electron2); 
	if (ele1->Charge() != ele2->Charge() && electronpass[electron1] && electronpass[electron2]){
		Bool_t cleaning = kTRUE;
		for (UInt_t i = 0; i < fGoodMuons->GetEntries(); ++i){
			if (fabs(fGoodMuons->At(i)->Eta()) < 2.4 && fGoodMuons->At(i)->IsGlobalMuon() &&
                        	mithep::MathUtils::DeltaR(ele1->Phi(),ele1->Eta(), fGoodMuons->At(i)->Phi(), fGoodMuons->At(i)->Eta()) < 0.05) cleaning = kFALSE;
	
			if (fabs(fGoodMuons->At(i)->Eta()) < 2.4 && fGoodMuons->At(i)->IsGlobalMuon() &&
                        	mithep::MathUtils::DeltaR(ele2->Phi(),ele2->Eta(), fGoodMuons->At(i)->Phi(), fGoodMuons->At(i)->Eta()) < 0.05) cleaning = kFALSE;
		}
		if (cleaning){	


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

	if (photonpass) {
		const Photon *pho;
		pho = fPhotons->At(photonindex);
		fLeptonPairPhotonEvent->mllg = (ele1->Mom() + ele2->Mom() + pho->Mom()).M();
		electronZgfilled = kTRUE;
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
	
	}
	fLeptonPairPhotonEvent->electronZmass = (ele1->Mom() + ele2->Mom()).M();

        fLeptonPairPhotonEvent->ele1dEtaIn = TMath::Abs(ele1->DeltaEtaSuperClusterTrackAtVtx());
        fLeptonPairPhotonEvent->ele1dPhiIn = TMath::Abs(ele1->DeltaPhiSuperClusterTrackAtVtx());
        fLeptonPairPhotonEvent->ele1sigmaIEtaIEta = ele1->CoviEtaiEta();
        fLeptonPairPhotonEvent->ele1HadOverEm = ele1->HadronicOverEm();
        fLeptonPairPhotonEvent->ele1D0 = fabs(ele1->BestTrk()->D0Corrected(*fPV->At(0)));
        fLeptonPairPhotonEvent->ele1DZ = fabs(ele1->BestTrk()->DzCorrected(*fPV->At(0)));
        fLeptonPairPhotonEvent->ele1OneOverEMinusOneOverP = fabs((1 - ele1->ESuperClusterOverP())/(ele1->EcalEnergy()));
	const BaseVertex *bsp = dynamic_cast<const BaseVertex*>(fBeamSpot->At(0));
	fLeptonPairPhotonEvent->ele1Conversion = mithep::ElectronTools::PassConversionFilter(ele1,fConversions,bsp, 0, 1e-6, 2.0, kTRUE, kFALSE, 7);
	fLeptonPairPhotonEvent->ele1missinghits = ele1->CorrectedNExpectedHitsInner();

        fLeptonPairPhotonEvent->ele2dEtaIn = TMath::Abs(ele2->DeltaEtaSuperClusterTrackAtVtx());
        fLeptonPairPhotonEvent->ele2dPhiIn = TMath::Abs(ele2->DeltaPhiSuperClusterTrackAtVtx());
        fLeptonPairPhotonEvent->ele2sigmaIEtaIEta = ele2->CoviEtaiEta();
        fLeptonPairPhotonEvent->ele2HadOverEm = ele2->HadronicOverEm();
        fLeptonPairPhotonEvent->ele2D0 = fabs(ele2->BestTrk()->D0Corrected(*fPV->At(0)));
        fLeptonPairPhotonEvent->ele2DZ = fabs(ele2->BestTrk()->DzCorrected(*fPV->At(0)));
        fLeptonPairPhotonEvent->ele2OneOverEMinusOneOverP = fabs((1 - ele2->ESuperClusterOverP())/(ele2->EcalEnergy()));
	fLeptonPairPhotonEvent->ele2Conversion = mithep::ElectronTools::PassConversionFilter(ele2,fConversions,bsp, 0, 1e-6, 2.0, kTRUE, kFALSE,7);
	fLeptonPairPhotonEvent->ele2missinghits = ele2->CorrectedNExpectedHitsInner();
	
	//Compute Electron MVA Values
	//double _fbrem = max(double(ele1->FBrem()),-1.0);
	//double _kftrk_chisq = (ele1->HasTrackerTrk()) ? ele1->TrackerTrk()->Chi2() / ele1->TrackerTrk()->Ndof() : 0;
	//double _kftrk_nhits = (ele1->HasTrackerTrk()) ? (Float_t)ele1->TrackerTrk()->NHits() : -1; 
	//double _gsftrk_chisq = (Double_t)min(double(ele1->BestTrk()->Chi2() / ele1->BestTrk()->Ndof()),200.0);
	//double seedE1x5OverE = ele1->SCluster()->Seed()->E1x5() / ele1->SCluster()->Seed()->Energy();
	//double seedE5x5OverE = ele1->SCluster()->Seed()->E5x5() / ele1->SCluster()->Seed()->Energy();
	//double _e1x5e5x5 = min(max(1 - double(seedE1x5OverE/seedE5x5OverE),-1.0),2.0);
	//double _r9 = min(double(ele1->SCluster()->R9()),5.0);
	//double _e_o_p = min(double(ele1->ESuperClusterOverP()), 20.0);
	//double _eseed_o_pout = min(double(ele1->ESeedClusterOverPout()),20.0);
	//double _IoEmIoP =  (Double_t)(1 - ele1->ESuperClusterOverP())/(ele1->SCluster()->Et()*TMath::CosH(ele1->SCluster()->Eta()));
	//double _epreoraw = (Float_t)ele1->SCluster()->PreshowerEnergy() / ele1->SCluster()->RawEnergy();
  
	//fLeptonPairPhotonEvent->ele1MVA = eleIDMVA->MVAValue_IDNonTrig(ele1->Pt(),
	//							       ele1->SCluster()->Eta(),
	//								       _fbrem,
	//								       _kftrk_chisq,
	//							       _kftrk_nhits,
	//							       _gsftrk_chisq,
	//							       fabs(ele1->DeltaEtaSuperClusterTrackAtVtx()),
	//							       fabs(ele1->DeltaPhiSuperClusterTrackAtVtx()),
	//							       fabs(ele1->DeltaEtaSeedClusterTrackAtCalo()),
	//							       ele1->CoviEtaiEta(),
	//							       sqrt(ele1->SCluster()->Seed()->CoviPhiiPhi()),
	//							       ele1->SCluster()->EtaWidth(), 
	//							       ele1->SCluster()->PhiWidth(),
	//							       _e1x5e5x5,
	//							       _r9,
	//							       ele1->HadronicOverEm(),
	//							       _e_o_p,
	//							       _IoEmIoP,
	//							       _eseed_o_pout,
	//							       _epreoraw,
	//							       kFALSE );


	//double _fbrem_2 = max(double(ele2->FBrem()),-1.0);
	//double _kftrk_chisq_2 = (ele2->HasTrackerTrk()) ? ele2->TrackerTrk()->Chi2() / ele2->TrackerTrk()->Ndof() : 0;
	//double _kftrk_nhits_2 = (ele2->HasTrackerTrk()) ? (Float_t)ele2->TrackerTrk()->NHits() : -1; 
	//double _gsftrk_chisq_2 = (Double_t)min(double(ele2->BestTrk()->Chi2() / ele2->BestTrk()->Ndof()),200.0);
	//double seedE1x5OverE_2 = ele2->SCluster()->Seed()->E1x5() / ele2->SCluster()->Seed()->Energy();
	//double seedE5x5OverE_2 = ele2->SCluster()->Seed()->E5x5() / ele2->SCluster()->Seed()->Energy();
	//double _e1x5e5x5_2 = min(max(1 - double(seedE1x5OverE_2/seedE5x5OverE_2),-1.0),2.0);
	//double _r9_2 = min(double(ele2->SCluster()->R9()),5.0);
	//double _e_o_p_2 = min(double(ele2->ESuperClusterOverP()), 20.0);
	//double _eseed_o_pout_2 = min(double(ele2->ESeedClusterOverPout()),20.0);
	//double _IoEmIoP_2 =  (Double_t)(1 - ele2->ESuperClusterOverP())/(ele2->SCluster()->Et()*TMath::CosH(ele2->SCluster()->Eta()));
	//double _epreoraw_2 = (Float_t)ele2->SCluster()->PreshowerEnergy() / ele2->SCluster()->RawEnergy();
  
  
	//fLeptonPairPhotonEvent->ele2MVA = eleIDMVA->MVAValue_IDNonTrig(ele2->Pt(),
	//							       ele2->SCluster()->Eta(),
	//							       _fbrem_2,
	//							       _kftrk_chisq_2,
	//							       _kftrk_nhits_2,
	//							       _gsftrk_chisq_2,
	//							       fabs(ele2->DeltaEtaSuperClusterTrackAtVtx()),
	//							       fabs(ele2->DeltaPhiSuperClusterTrackAtVtx()),
	//							       fabs(ele2->DeltaEtaSeedClusterTrackAtCalo()),
	//							       ele2->CoviEtaiEta(),
	//							       sqrt(ele2->SCluster()->Seed()->CoviPhiiPhi()),
	//							       ele2->SCluster()->EtaWidth(), 
	//							       ele2->SCluster()->PhiWidth(),
	//							       _e1x5e5x5_2,
	//							       _r9_2,
	//							       ele2->HadronicOverEm(),
	//							       _e_o_p_2,
	//							       _IoEmIoP_2,
	//							       _eseed_o_pout_2,
	//							       _epreoraw_2,
	//							       kFALSE );
   			
   		}  
	}
   }
   
   if (fGoodMuons->GetEntries() > 1){
      	vector<bool> muonpass;
        for (UInt_t j = 0; j < fGoodMuons->GetEntries();++j){
                const Muon *test;
		test = fGoodMuons->At(j);
		if (test->Pt() > 10 && fabs(test->Eta()) < 2.4 && test->IsGlobalMuon() == true){
                        if (ZGTools::muonIDPOGTightSelection(YEAR,test, fPV->At(0),fPFCands)){
			        if ((YEAR == 2011 && ZGTools::muonPFIso04(test,fPV->At(0),fPFCands,fPileUpDen, mithep::MuonTools::kMuEAData2011, pfNoPileUpflag,YEAR)) || (YEAR == 2012 && ZGTools::muonPFIso04(test,fPV->At(0),fPFCands,fPileUpDen, mithep::MuonTools::kMuEAData2012, pfNoPileUpflag,YEAR))){
				         muonpass.push_back(1);
                                }
                                else muonpass.push_back(0);
                        }
                        else muonpass.push_back(0);
                }
                else muonpass.push_back(0);

        }
        const Muon *muona;
        const Muon *muonb;
        muona = 0;
        muonb = 0;
        Float_t zdifference = 999;
        UInt_t muon1 = 0;
        UInt_t muon2 = 1;

        for (UInt_t j = 0; j < fGoodMuons->GetEntries(); ++j){
                for (UInt_t k = 0; k < j; ++k){
                if (j == k) continue;
                muona = fGoodMuons->At(j);
                muonb = fGoodMuons->At(k);
                        if (muona->Charge() != muonb->Charge() && muonpass[j] && muonpass[k]){
                                if (fabs(91.19 - (muona->Mom() + muonb->Mom()).M()) < zdifference){
                                        zdifference = fabs(91.19 - (muona->Mom() + muonb->Mom()).M());
                                        muon1 = k;
                                        muon2 = j;
                                }
                        }
                }
        }

        muona = fGoodMuons->At(muon1);
        muonb = fGoodMuons->At(muon2);
        if (muona->Charge() != muonb->Charge() && muonpass[muon1] && muonpass[muon2]){
	
        	if (photonpass){
			const Photon *pho;
			pho = fPhotons->At(photonindex);
			fLeptonPairPhotonEvent->mllg = (pho->Mom()+muona->Mom()+muonb->Mom()).M();
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
		}
	

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

	ZgllTuple->Fill();
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
  
  ZgllTuple->Branch("electronZmass",&fLeptonPairPhotonEvent->electronZmass,"electronZmass/F");
  ZgllTuple->Branch("mllg",&fLeptonPairPhotonEvent->mllg,"mllg/F");
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

  ZgllTuple->Branch("ele2charge",&fLeptonPairPhotonEvent->ele2charge,"ele2charge/F");
  ZgllTuple->Branch("ele2energy",&fLeptonPairPhotonEvent->ele2energy,"ele2energy/F");
  ZgllTuple->Branch("ele2px",&fLeptonPairPhotonEvent->ele2px,"ele2px/F");
  ZgllTuple->Branch("ele2py",&fLeptonPairPhotonEvent->ele2py,"ele2py/F");
  ZgllTuple->Branch("ele2pz",&fLeptonPairPhotonEvent->ele2pz,"ele2pz/F");
  ZgllTuple->Branch("ele2pt",&fLeptonPairPhotonEvent->ele2pt,"ele2pt/F");
  ZgllTuple->Branch("ele2eta",&fLeptonPairPhotonEvent->ele2eta,"ele2eta/F");
  ZgllTuple->Branch("ele2mass",&fLeptonPairPhotonEvent->ele2mass,"ele2mass/F");
  ZgllTuple->Branch("ele2phi",&fLeptonPairPhotonEvent->ele2phi,"ele2phi/F");

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
  ZgllTuple->Branch("m2E",&fLeptonPairPhotonEvent->m2E,"m2E/F");
  ZgllTuple->Branch("m2Pt",&fLeptonPairPhotonEvent->m2Pt,"m2Pt/F");
  ZgllTuple->Branch("m2Mass",&fLeptonPairPhotonEvent->m2Mass,"m2Mass/F");
  ZgllTuple->Branch("m2Px",&fLeptonPairPhotonEvent->m2Px,"m2Px/F");
  ZgllTuple->Branch("m2Py",&fLeptonPairPhotonEvent->m2Py,"m2Py/F");
  ZgllTuple->Branch("m2Pz",&fLeptonPairPhotonEvent->m2Pz,"m2Pz/F");
  ZgllTuple->Branch("m2Eta",&fLeptonPairPhotonEvent->m2Eta,"m2Eta/F");
  ZgllTuple->Branch("m2Phi",&fLeptonPairPhotonEvent->m2Phi,"m2Phi/F");
  ZgllTuple->Branch("m2Charge",&fLeptonPairPhotonEvent->m2Charge,"m2Charge/F");

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
  
  //Add Output Tree
  AddOutput(ZgllTuple);
}

