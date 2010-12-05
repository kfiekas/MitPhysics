 // $Id $

#include "MitPhysics/SelMods/interface/HwwExampleAnalysisMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "TFile.h"
#include "TTree.h"

using namespace mithep;
ClassImp(mithep::HwwExampleAnalysisMod)

//--------------------------------------------------------------------------------------------------
HwwExampleAnalysisMod::HwwExampleAnalysisMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fMuonBranchName(Names::gkMuonBrn),
  fMetName("NoDefaultNameSet"),
  fCleanJetsName("NoDefaultNameSet"),
  fCleanJetsNoPtCutName("NoDefaultNameSet"),
  fCaloJetName0("AKt5Jets"),
  fVertexName(ModNames::gkGoodVertexesName),
  fMuons(0),
  fMet(0),
  fVertices(0),
  fCaloJet0(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void HwwExampleAnalysisMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void HwwExampleAnalysisMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  // Load Branches
  ReqBranch(fMuonBranchName,  fMuons);
  ReqBranch(fCaloJetName0,    fCaloJet0);

  //Create your histograms here

  //*************************************************************************************************
  // Selection Histograms
  //*************************************************************************************************
  AddTH1(fHWWSelection,"hHWWSelection", ";Cut Number;Number of Events", 11, -1.5, 9.5);
  AddTH1(fHWWToEESelection,"hHWWToEESelection", ";Cut Number;Number of Events", 11, -1.5, 9.5);
  AddTH1(fHWWToMuMuSelection,"hHWWToMuMuSelection", ";Cut Number;Number of Events", 11, -1.5, 9.5);
  AddTH1(fHWWToEMuSelection,"hHWWToEMuSelection", ";Cut Number;Number of Events", 11, -1.5, 9.5);

  //***********************************************************************************************
  // Histograms after preselection
  //***********************************************************************************************
  AddTH1(fLeptonEta          ,"hLeptonEta",";LeptonEta;Number of Events",100,-5.,5.0);
  AddTH1(fLeptonPtMax        ,"hLeptonPtMax",";Lepton P_t Max;Number of Events",150,0.,150.);
  AddTH1(fLeptonPtMin        ,"hLeptonPtMin",";Lepton P_t Min;Number of Events",150,0.,150.);
  AddTH1(fMetPtHist          ,"hMetPtHist",";Met;Number of Events",150,0.,300.);  
  AddTH1(fMetPhiHist         ,"hMetPhiHist",";#phi;Number of Events",28,-3.5,3.5);
  AddTH1(fUncorrMetPtHist    ,"hUncorrMetPtHist",";Met;Number of Events",150,0.,300.);  
  AddTH1(fUncorrMetPhiHist   ,"hUncorrMetPhiHist",";#phi;Number of Events",28,-3.5,3.5);
  AddTH1(fDeltaPhiLeptons    ,"hDeltaPhiLeptons",";#Delta#phi_{ll};Number of Events",90,0,180);
  AddTH1(fDeltaEtaLeptons    ,"hDeltaEtaLeptons",";#Delta#eta_{ll};Number of Events",100,-5.,5.0);
  AddTH1(fDileptonMass       ,"hDileptonMass",";Mass_{ll};Number of Events",150,0.,300.);
 
  //***********************************************************************************************
  // N-1 Histograms
  //***********************************************************************************************
  //All events
  AddTH1(fLeptonPtMax_NMinusOne            ,"hLeptonPtMax_NMinusOne",
                                            ";Lepton P_t Max;Number of Events",150,0.,150.);
  AddTH1(fLeptonPtMin_NMinusOne            ,"hLeptonPtMin_NMinusOne",
                                            ";Lepton P_t Min;Number of Events",150,0.,150.);
  AddTH1(fMetPtHist_NMinusOne              ,"hMetPtHist_NMinusOne",
                                            ";Met;Number of Events",150,0.,300.);  
  AddTH1(fMetPhiHist_NMinusOne             ,"hMetPhiHist_NMinusOne",
                                            ";#phi;Number of Events",28,-3.5,3.5);
  AddTH1(fMETdeltaPhilEtHist_NMinusOne     ,"hMETdeltaPhilEtHist_NMinusOne",
                                            ";METdeltaPhilEtHist;Number of Events",150,0.,300.);
  AddTH1(fNCentralJets_NMinusOne           ,"hNCentralJets_NMinusOne",
                                            ";Number of Central Jets;Number of Events",6,-0.5,5.5);
  AddTH1(fNSoftMuonsHist_NMinusOne        ,"hNSoftMuonsHist_NMinusOne",
                                            ";Number of Dirty Muons; Number of Events",6,-0.5,5.5);
  AddTH1(fDeltaPhiLeptons_NMinusOne        ,"hDeltaPhiLeptons_NMinusOne",
                                            ";#Delta#phi_{ll};Number of Events",90,0,180);
  AddTH1(fDeltaEtaLeptons_NMinusOne        ,"hDeltaEtaLeptons_NMinusOne",
                                            ";#Delta#eta_{ll};Number of Events",100,-5.0,5.0);
  AddTH1(fDileptonMass_NMinusOne           ,"hDileptonMass_NMinusOne",
                                            ";Mass_{ll};Number of Events",150,0.,300.);
  AddTH1(fMinDeltaPhiLeptonMet_NMinusOne   ,"hMinDeltaPhiLeptonMet_NMinusOne", 
                                            ";Min #Delta#phi_{l,Met};Number of Events",90,0.,180);

  //***********************************************************************************************
  // After all cuts Histograms
  //***********************************************************************************************
  AddTH1(fMinDeltaPhiLeptonMet_afterCuts    ,"hMinDeltaPhiLeptonMet_afterCuts", 
                                             ";Min #Delta#phi_{l,Met};Number of Events",90,0.,180);
  AddTH1(fMtLepton1_afterCuts               ,"hMtLepton1_afterCuts",
                                             ";M_t (Lepton1,Met);Number of Events",100,0.,200.);
  AddTH1(fMtLepton2_afterCuts               ,"hMtLepton2_afterCuts",
                                             ";M_t (Lepton2,Met);Number of Events",100,0.,200.);
  AddTH1(fMtHiggs_afterCuts                 ,"hMtHiggs_afterCuts",
                                             ";M_t (l1+l2+Met);Number of Events",150,0.,300.);
  AddTH1(fLeptonPtPlusMet_afterCuts         ,"hLeptonPtPlusMet_afterCuts",
                                             ";LeptonPtPlusMet;Number of Events",150,0., 300.);

}

//--------------------------------------------------------------------------------------------------
void HwwExampleAnalysisMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  LoadBranch(fMuonBranchName);
  LoadBranch(fCaloJetName0);
 
  //Obtain all the good objects from the event cleaning module
  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
  ObjArray<Muon> *CleanMuons = dynamic_cast<ObjArray<Muon>* >(FindObjThisEvt(ModNames::gkCleanMuonsName));
  ObjArray<Electron> *CleanElectrons = dynamic_cast<ObjArray<Electron>* >(FindObjThisEvt(ModNames::gkCleanElectronsName));
  ParticleOArr *CleanLeptons = dynamic_cast<mithep::ParticleOArr*>
    (FindObjThisEvt(ModNames::gkMergedLeptonsName));
  ObjArray<Jet> *CleanJets = dynamic_cast<ObjArray<Jet>* >
    (FindObjThisEvt(fCleanJetsName.Data()));
  ObjArray<Jet> *CleanJetsNoPtCut = dynamic_cast<ObjArray<Jet>* >
    (FindObjThisEvt(fCleanJetsNoPtCutName.Data()));
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  MetCol *met = dynamic_cast<ObjArray<Met>* >(FindObjThisEvt(fMetName));
  const Met *caloMet = 0;
  if (met) {
    caloMet = met->At(0);
  } else {
    cout << "Error: Met Collection " << fMetName << " could not be loaded.\n";
    return;
  }

  //***********************************************************************************************
  //Kinematic PreSelection
  //***********************************************************************************************
  // At least two leptons in the event
  if (CleanLeptons->GetEntries() < 2) return;
  // Pt1 > 20 && Pt2 > 10
  if(CleanLeptons->At(0)->Pt() <= 20 || CleanLeptons->At(1)->Pt() <= 10) return;
  // opposite charge leptons
  if(CleanLeptons->At(0)->Charge() * CleanLeptons->At(1)->Charge() > 0) return;
    
  CompositeParticle *dilepton = new CompositeParticle();
  dilepton->AddDaughter(CleanLeptons->At(0));
  dilepton->AddDaughter(CleanLeptons->At(1));
   
  //***********************************************************************************************
  //Get Dirty Muons: Non-isolated Muons (exclude the clean muons)
  //***********************************************************************************************
  ObjArray<Muon> *SoftMuons = new ObjArray<Muon>;
  for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i);
    if(!MuonTools::PassSoftMuonCut(mu, fVertices)) continue;
    
    bool isCleanMuon = kFALSE;
    for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
      if(fMuons->At(i) == CleanMuons->At(j) &&
  	 CleanMuons->At(j)->Pt() > 10) isCleanMuon = kTRUE;
    }
    if(isCleanMuon == kFALSE) SoftMuons->Add(mu);
  }

  //***********************************************************************************************
  //|Z_vert-Z_l| maximum
  //***********************************************************************************************
  double zDiffMax = 0.0;
  if(fVertices->GetEntries() > 0) {
    for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
      double pDz = 0.0;
      for(uint i0 = 0; i0 < fVertices->GetEntries(); i0++) {
        if(fVertices->At(i0)->NTracks() > 0){
	  pDz = TMath::Abs(CleanMuons->At(j)->BestTrk()->DzCorrected(*fVertices->At(i0)));
          break;
        }
      }
      if(pDz > zDiffMax) zDiffMax = pDz;
    }
    for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {   
      double pDz = 0.0;
      for(uint i0 = 0; i0 < fVertices->GetEntries(); i0++) {
        if(fVertices->At(i0)->NTracks() > 0){
	  pDz = TMath::Abs(CleanElectrons->At(j)->GsfTrk()->DzCorrected(*fVertices->At(i0)));
          break;
        }
      }
      if(pDz > zDiffMax) zDiffMax = pDz;
    }
  }

  //***********************************************************************************************
  //Define Event Variables
  //***********************************************************************************************
  //delta phi between the 2 leptons in degrees
  double deltaPhiLeptons = MathUtils::DeltaPhi(CleanLeptons->At(0)->Phi(), 
                                               CleanLeptons->At(1)->Phi())* 180.0 / TMath::Pi();

  double deltaEtaLeptons = CleanLeptons->At(0)->Eta() - CleanLeptons->At(1)->Eta();

  double deltaPhiDileptonMet = MathUtils::DeltaPhi(caloMet->Phi(), 
                                                   dilepton->Phi())*180.0 / TMath::Pi();

  double mtHiggs = TMath::Sqrt(2.0*dilepton->Pt() * caloMet->Pt()*
			       (1.0 - cos(deltaPhiDileptonMet * TMath::Pi() / 180.0)));

  //angle between MET and closest lepton
  double deltaPhiMetLepton[2] = {MathUtils::DeltaPhi(caloMet->Phi(), CleanLeptons->At(0)->Phi()),
                                 MathUtils::DeltaPhi(caloMet->Phi(), CleanLeptons->At(1)->Phi())};
  
  double mTW[2] = {TMath::Sqrt(2.0*CleanLeptons->At(0)->Pt()*caloMet->Pt()*
                               (1.0 - cos(deltaPhiMetLepton[0]))),
		   TMath::Sqrt(2.0*CleanLeptons->At(1)->Pt()*caloMet->Pt()*
                               (1.0 - cos(deltaPhiMetLepton[1])))};

  double minDeltaPhiMetLepton = (deltaPhiMetLepton[0] < deltaPhiMetLepton[1])?
    deltaPhiMetLepton[0]:deltaPhiMetLepton[1];

  double METdeltaPhilEt = caloMet->Pt();
  if(minDeltaPhiMetLepton < TMath::Pi()/2.)
      METdeltaPhilEt = METdeltaPhilEt * sin(minDeltaPhiMetLepton);

  //count the number of central Jets for vetoing and b-tagging
  vector<Jet*> sortedJetsAll;
  vector<Jet*> sortedJets;
  vector<Jet*> sortedJetsLowPt;
  for(UInt_t i=0; i<CleanJetsNoPtCut->GetEntries(); i++){
    Jet* jet_a = new Jet(CleanJetsNoPtCut->At(i)->Px(),
   			 CleanJetsNoPtCut->At(i)->Py(),
   			 CleanJetsNoPtCut->At(i)->Pz(),
   			 CleanJetsNoPtCut->At(i)->E() );

    int nCloseStdJet = -1;
    double deltaRMin = 999.;
    for(UInt_t nj=0; nj<fCaloJet0->GetEntries(); nj++){
      const CaloJet *jet = fCaloJet0->At(nj);
      Double_t deltaR = MathUtils::DeltaR(jet_a->Mom(),jet->Mom());
      if(deltaR < deltaRMin) {
   	nCloseStdJet = nj;
   	deltaRMin = deltaR;
      }
    }
    if(nCloseStdJet >= 0 && deltaRMin < 0.5){
      jet_a->SetMatchedMCFlavor(fCaloJet0->At(nCloseStdJet)->MatchedMCFlavor());
      jet_a->SetCombinedSecondaryVertexBJetTagsDisc(fCaloJet0->At(nCloseStdJet)->CombinedSecondaryVertexBJetTagsDisc());
      jet_a->SetCombinedSecondaryVertexMVABJetTagsDisc(fCaloJet0->At(nCloseStdJet)->CombinedSecondaryVertexMVABJetTagsDisc());
      jet_a->SetJetProbabilityBJetTagsDisc(fCaloJet0->At(nCloseStdJet)->JetProbabilityBJetTagsDisc());
      jet_a->SetJetBProbabilityBJetTagsDisc(fCaloJet0->At(nCloseStdJet)->JetBProbabilityBJetTagsDisc());
      jet_a->SetTrackCountingHighEffBJetTagsDisc(fCaloJet0->At(nCloseStdJet)->TrackCountingHighEffBJetTagsDisc());
      jet_a->SetTrackCountingHighPurBJetTagsDisc(fCaloJet0->At(nCloseStdJet)->TrackCountingHighPurBJetTagsDisc());
      jet_a->SetSimpleSecondaryVertexBJetTagsDisc(fCaloJet0->At(nCloseStdJet)->SimpleSecondaryVertexBJetTagsDisc());
      jet_a->SetSimpleSecondaryVertexHighEffBJetTagsDisc(fCaloJet0->At(nCloseStdJet)->SimpleSecondaryVertexHighEffBJetTagsDisc());
      jet_a->SetSimpleSecondaryVertexHighPurBJetTagsDisc(fCaloJet0->At(nCloseStdJet)->SimpleSecondaryVertexHighPurBJetTagsDisc());
    }
    else {
      jet_a->SetMatchedMCFlavor(CleanJetsNoPtCut->At(i)->MatchedMCFlavor());
      jet_a->SetCombinedSecondaryVertexBJetTagsDisc(CleanJetsNoPtCut->At(i)->CombinedSecondaryVertexBJetTagsDisc());
      jet_a->SetCombinedSecondaryVertexMVABJetTagsDisc(CleanJetsNoPtCut->At(i)->CombinedSecondaryVertexMVABJetTagsDisc());
      jet_a->SetJetProbabilityBJetTagsDisc(CleanJetsNoPtCut->At(i)->JetProbabilityBJetTagsDisc());
      jet_a->SetJetBProbabilityBJetTagsDisc(CleanJetsNoPtCut->At(i)->JetBProbabilityBJetTagsDisc());
      jet_a->SetTrackCountingHighEffBJetTagsDisc(CleanJetsNoPtCut->At(i)->TrackCountingHighEffBJetTagsDisc());
      jet_a->SetTrackCountingHighPurBJetTagsDisc(CleanJetsNoPtCut->At(i)->TrackCountingHighPurBJetTagsDisc());
      jet_a->SetSimpleSecondaryVertexBJetTagsDisc(CleanJetsNoPtCut->At(i)->SimpleSecondaryVertexBJetTagsDisc());
      jet_a->SetSimpleSecondaryVertexHighEffBJetTagsDisc(CleanJetsNoPtCut->At(i)->SimpleSecondaryVertexHighEffBJetTagsDisc());
      jet_a->SetSimpleSecondaryVertexHighPurBJetTagsDisc(CleanJetsNoPtCut->At(i)->SimpleSecondaryVertexHighPurBJetTagsDisc());
    }
    sortedJetsAll.push_back(jet_a);
  }

  for(UInt_t i=0; i<CleanJets->GetEntries(); i++){
    if(TMath::Abs(CleanJets->At(i)->Eta()) < 5.0 &&
       CleanJets->At(i)->Pt() > 25.0){
      Jet* jet_b = new Jet(CleanJets->At(i)->Px(),
     			   CleanJets->At(i)->Py(),
   			   CleanJets->At(i)->Pz(),
   			   CleanJets->At(i)->E() );
      sortedJets.push_back(jet_b);
    }
  }

  for(UInt_t i=0; i<sortedJetsAll.size(); i++){
    bool overlap = kFALSE;
    for(UInt_t j=0; j<sortedJets.size(); j++){
      if(sortedJetsAll[i]->Pt() == sortedJets[j]->Pt() ||
        (sortedJetsAll[i]->CombinedSecondaryVertexBJetTagsDisc() == sortedJets[j]->CombinedSecondaryVertexBJetTagsDisc() &&
	 sortedJetsAll[i]->JetBProbabilityBJetTagsDisc()	 == sortedJets[j]->JetBProbabilityBJetTagsDisc() &&
	 sortedJetsAll[i]->TrackCountingHighPurBJetTagsDisc()	 == sortedJets[j]->TrackCountingHighPurBJetTagsDisc())
        ) {
        sortedJets[j]->SetMatchedMCFlavor(sortedJetsAll[i]->MatchedMCFlavor());
        sortedJets[j]->SetCombinedSecondaryVertexBJetTagsDisc(sortedJetsAll[i]->CombinedSecondaryVertexBJetTagsDisc());
        sortedJets[j]->SetCombinedSecondaryVertexMVABJetTagsDisc(sortedJetsAll[i]->CombinedSecondaryVertexMVABJetTagsDisc());
        sortedJets[j]->SetJetProbabilityBJetTagsDisc(sortedJetsAll[i]->JetProbabilityBJetTagsDisc());
        sortedJets[j]->SetJetBProbabilityBJetTagsDisc(sortedJetsAll[i]->JetBProbabilityBJetTagsDisc());
        sortedJets[j]->SetTrackCountingHighEffBJetTagsDisc(sortedJetsAll[i]->TrackCountingHighEffBJetTagsDisc());
        sortedJets[j]->SetTrackCountingHighPurBJetTagsDisc(sortedJetsAll[i]->TrackCountingHighPurBJetTagsDisc());
        sortedJets[j]->SetSimpleSecondaryVertexBJetTagsDisc(sortedJetsAll[i]->SimpleSecondaryVertexBJetTagsDisc());
        sortedJets[j]->SetSimpleSecondaryVertexHighEffBJetTagsDisc(sortedJetsAll[i]->SimpleSecondaryVertexHighEffBJetTagsDisc());
        sortedJets[j]->SetSimpleSecondaryVertexHighPurBJetTagsDisc(sortedJetsAll[i]->SimpleSecondaryVertexHighPurBJetTagsDisc());        
  	overlap = kTRUE;
        break;
      }
    }
    if(overlap == kFALSE){
      sortedJetsLowPt.push_back(sortedJetsAll[i]);
    }
  }
  double maxBtag = -99999.;
  double imaxBtag = -1;
  for(UInt_t i=0; i<sortedJetsLowPt.size(); i++){
    if(sortedJetsLowPt[i]->TrackCountingHighEffBJetTagsDisc() > maxBtag){
      maxBtag  = sortedJetsLowPt[i]->TrackCountingHighEffBJetTagsDisc();
      imaxBtag = i;
    }
  }

  //Lepton Type
  int finalstateType = -1;
  if (CleanLeptons->At(0)->ObjType() == kMuon && CleanLeptons->At(1)->ObjType() == kMuon ){ // mumu
    finalstateType = 10;
  } else if(CleanLeptons->At(0)->ObjType() == kElectron && CleanLeptons->At(1)->ObjType() == kElectron ){ // ee
    finalstateType = 11;
  } else if((CleanLeptons->At(0)->ObjType() == kElectron && CleanLeptons->At(1)->ObjType() == kMuon) ||
            (CleanLeptons->At(1)->ObjType() == kElectron && CleanLeptons->At(0)->ObjType() == kMuon)) {
    finalstateType = 12;
  } else {
    cerr << "Error: finalstate lepton type not supported\n";
  }
                        
  //*********************************************************************************************
  //Define Cuts
  //*********************************************************************************************
  const int nCuts = 10;
  bool passCut[nCuts] = {false, false, false, false, false,
                         false, false, false, false, false};
  
  if(CleanLeptons->At(0)->Pt() >  20.0 &&
     CleanLeptons->At(1)->Pt() >= 20.0) passCut[0] = true;
  
  if(zDiffMax < 1.0)                    passCut[1] = true;
  
  if(caloMet->Pt()    > 20.0)           passCut[2] = true;
  
  if(dilepton->Mass() > 12.0)           passCut[3] = true;
  
  if(sortedJets.size() < 1)             passCut[6] = true;

  if(SoftMuons->GetEntries() == 0)      passCut[7] = true;

  if(CleanLeptons->GetEntries() == 2)   passCut[8] = true;

  if(maxBtag < 2.1)                     passCut[9] = true;

  if (finalstateType == 10 || finalstateType == 11){ // mumu/ee
    if(fabs(dilepton->Mass()-91.1876)   > 15.0)   passCut[4] = true;
    if(METdeltaPhilEt > 35) passCut[5] = true;
  }
  else if(finalstateType == 12) { // emu
    passCut[4] = true;
    if(METdeltaPhilEt > 20) passCut[5] = true;
  }
  
  //*********************************************************************************************
  //Make Selection Histograms. Number of events passing each level of cut
  //*********************************************************************************************  
  bool passAllCuts = true;
  for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
    
  //Cut Selection Histograms
  fHWWSelection->Fill(-1,NNLOWeight->GetVal());
  if (finalstateType == 10 )
    fHWWToMuMuSelection->Fill(-1,NNLOWeight->GetVal());
  else if(finalstateType == 11 )
    fHWWToEESelection->Fill(-1,NNLOWeight->GetVal());
  else if(finalstateType == 12 )
    fHWWToEMuSelection->Fill(-1,NNLOWeight->GetVal());

  for (int k=0;k<nCuts;k++) {
    bool pass = true;
    bool passPreviousCut = true;
    for (int p=0;p<=k;p++) {
      pass = (pass && passCut[p]);
      if (p<k)
        passPreviousCut = (passPreviousCut&& passCut[p]);
    }
    
    if (pass) {
      fHWWSelection->Fill(k,NNLOWeight->GetVal());
      if (finalstateType == 10 )
        fHWWToMuMuSelection->Fill(k,NNLOWeight->GetVal());
      else if(finalstateType == 11)
        fHWWToEESelection->Fill(k,NNLOWeight->GetVal());
      else if(finalstateType == 12)
        fHWWToEMuSelection->Fill(k,NNLOWeight->GetVal());
    }
  }
  
  //*****************************************************************************************
  //Make Preselection Histograms  
  //*****************************************************************************************
  fLeptonEta->Fill(CleanLeptons->At(0)->Eta(),NNLOWeight->GetVal()); 
  fLeptonEta->Fill(CleanLeptons->At(1)->Eta(),NNLOWeight->GetVal());
  fLeptonPtMax->Fill(CleanLeptons->At(0)->Pt(),NNLOWeight->GetVal());
  fLeptonPtMin->Fill(CleanLeptons->At(1)->Pt(),NNLOWeight->GetVal());
  fMetPtHist->Fill(caloMet->Pt(),NNLOWeight->GetVal());                             
  fMetPhiHist->Fill(caloMet->Phi(),NNLOWeight->GetVal());                            
  fDeltaPhiLeptons->Fill(deltaPhiLeptons,NNLOWeight->GetVal());
  fDeltaEtaLeptons->Fill(deltaEtaLeptons,NNLOWeight->GetVal());
  fDileptonMass->Fill(dilepton->Mass(),NNLOWeight->GetVal());    

  //*********************************************************************************************
  // N-1 Histograms
  //*********************************************************************************************
  bool pass;;
  
  //N Jet Veto  
  pass = true;
  for (int k=0;k<nCuts;k++) {
    if (k != 6) {
      pass = (pass && passCut[k]);      
    }
  }
  if (pass) {
    fNCentralJets_NMinusOne->Fill(sortedJets.size(),NNLOWeight->GetVal());
  }     
  
  // Final Met Cut
  pass = true;
  for (int k=0;k<nCuts;k++) {
    if (k != 5) {
      pass = (pass && passCut[k]);      
    }
  }
  if (pass) {
    fMetPtHist_NMinusOne->Fill(caloMet->Pt(),NNLOWeight->GetVal());  
  }

  // dilepton mass
  pass = true;
  for (int k=0;k<nCuts;k++) {
    if (k != 3 && k !=  4)
      pass = (pass && passCut[k]);    
  }
  if (pass) {
    fDileptonMass_NMinusOne->Fill(dilepton->Mass(),NNLOWeight->GetVal());
  }
  
  // Lepton Pt Max, Lepton Pt Min, DeltaPhiLeptons
  pass = true;
  for (int k=0;k<nCuts;k++) {
    pass = (pass && passCut[k]);      
  }
  if (pass) {
    fLeptonPtMax_NMinusOne->Fill(CleanLeptons->At(0)->Pt(),NNLOWeight->GetVal());
    fLeptonPtMin_NMinusOne->Fill(CleanLeptons->At(1)->Pt(),NNLOWeight->GetVal());
    fDeltaPhiLeptons_NMinusOne->Fill(deltaPhiLeptons,NNLOWeight->GetVal()); 
  }
  
  // NSoftMuons
  pass = true;
  for (int k=0;k<nCuts;k++) {
    if (k != 7)
      pass = (pass && passCut[k]);    
  }
  if (pass) {
    fNSoftMuonsHist_NMinusOne->Fill(SoftMuons->GetEntries(),NNLOWeight->GetVal());
  }

  //*********************************************************************************************
  //Plots after all Cuts
  //*********************************************************************************************
  if (passAllCuts) {
    fMinDeltaPhiLeptonMet_afterCuts->Fill(minDeltaPhiMetLepton,NNLOWeight->GetVal());
    fMtLepton1_afterCuts->Fill(mTW[0],NNLOWeight->GetVal());
    fMtLepton2_afterCuts->Fill(mTW[1],NNLOWeight->GetVal());
    fMtHiggs_afterCuts->Fill(mtHiggs,NNLOWeight->GetVal());
    fLeptonPtPlusMet_afterCuts->Fill(CleanLeptons->At(0)->Pt()+CleanLeptons->At(1)->Pt()+caloMet->Pt(),NNLOWeight->GetVal());
  }
  
  delete dilepton;
  delete SoftMuons;
  for(UInt_t i=0; i<sortedJets.size();      i++) delete sortedJets[i];
  for(UInt_t i=0; i<sortedJetsAll.size();   i++) delete sortedJetsAll[i];
  return;
}

//--------------------------------------------------------------------------------------------------
void HwwExampleAnalysisMod::SlaveTerminate()
{
  
  // Run finishing code on the computer (slave) that did the analysis. For this
  // module, we dont do anything here.

} 
//--------------------------------------------------------------------------------------------------
void HwwExampleAnalysisMod::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.

}
