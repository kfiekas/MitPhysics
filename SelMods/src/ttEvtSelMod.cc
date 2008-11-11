 // $Id: ttEvtSelMod.cc,v 1.2 2008/10/23 12:21:37 ceballos Exp $

#include "MitPhysics/SelMods/interface/ttEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/MuonTools.h"

using namespace mithep;
ClassImp(mithep::ttEvtSelMod)

//--------------------------------------------------------------------------------------------------
ttEvtSelMod::ttEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(false),
  fMetName(Names::gkCaloMetBrn),
  fMuonName(Names::gkMuonBrn),
  fTrackName(Names::gkTrackBrn),
  fVertexName(string("PrimaryVertexesBeamSpot").c_str()),
  fCleanJetsName(Names::gkCleanJetsName),
  fMCLeptonsName(Names::gkMCLeptonsName),
  fMCAllLeptonsName(Names::gkMCAllLeptonsName),
  fMCQuarksName(Names::gkMCQuarksName),
  fMet(0),
  fMuons(0),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void ttEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void ttEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  if (fNEventsProcessed % 1000000 == 0 || fPrintDebug) {
    time_t systime;
    systime = time(NULL);
    cerr << endl << "ttEvtSelMod : Process Event " << fNEventsProcessed << "  Time: " << ctime(&systime) << endl;
  }

  //Get Generator Level information for matching
  //ObjArray<MCParticle> *GenLeptons    = dynamic_cast<ObjArray<MCParticle>* > (FindObjThisEvt(fMCLeptonsName.Data()));
  ObjArray<MCParticle> *GenAllLeptons = dynamic_cast<ObjArray<MCParticle>* > (FindObjThisEvt(fMCAllLeptonsName.Data()));
  ObjArray<MCParticle> *GenQuarks     = dynamic_cast<ObjArray<MCParticle>* > (FindObjThisEvt(fMCQuarksName.Data()));

  //Obtain all the good objects from the event cleaning module
  ObjArray<Electron> *CleanElectrons = dynamic_cast<ObjArray<Electron>* >
    (FindObjThisEvt(Names::gkCleanElectronsName));
  ObjArray<Muon> *CleanMuons = dynamic_cast<ObjArray<Muon>* >
    (FindObjThisEvt(Names::gkCleanMuonsName));
  ObjArray<Jet> *CleanJets = dynamic_cast<ObjArray<Jet>* >
    (FindObjThisEvt(fCleanJetsName.Data()));

  LoadBranch(fMetName);
  Met *caloMet = fMet->At(0);

  vector<ChargedParticle*> leptons;
  vector<string> leptonType; 
  double zAverage = 0.0;
  MuonTools myMuonTools;

  LoadBranch(fMuonName);
  ObjArray<Muon> *DirtyMuons = new ObjArray<Muon>;
  for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    Muon *mu = fMuons->At(i);
    if(!mu->GlobalTrk()) continue;
    if(mu->Pt() < 5.0)   continue;

    bool isCleanMuon = false;
    for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
      if(fMuons->At(i) == CleanMuons->At(j) &&
         CleanMuons->At(j)->Pt() > 15) isCleanMuon = true;
    }
    if(isCleanMuon == false) DirtyMuons->Add(mu);
  }

  // Make lepton vector from muons and electrons
  for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
    if(CleanMuons->At(j)->Pt() <= 15) continue;
    zAverage = zAverage + CleanMuons->At(j)->BestTrk()->Z0();
    leptons.push_back(CleanMuons->At(j));
    leptonType.push_back("mu");
  }
  for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {   
    if(CleanElectrons->At(j)->Pt() <= 15) continue;
    zAverage = zAverage + CleanElectrons->At(j)->BestTrk()->Z0();
    leptons.push_back(CleanElectrons->At(j));
    leptonType.push_back("e");
  }

  // Computing Z average (our primary vertex)
  if(leptons.size() > 0) zAverage = zAverage / leptons.size();

  // Sort the Leptons by Pt   
  for(UInt_t i=0; i<leptons.size(); i++){
    for(UInt_t j=i+1; j<leptons.size(); j++){
      if(leptons[i]->Pt() < leptons[j]->Pt()) {
	//swap i and j
	ChargedParticle* templepton = leptons[i];
	leptons[i] = leptons[j];
	leptons[j] = templepton;
	string templeptonType = leptonType[i];
	leptonType[i] = leptonType[j];
	leptonType[j] = templeptonType;	 
      }
    }
  }
  if (fPrintDebug) {
    cout << "Check Lepton Sort\n";
    for(UInt_t i=0; i<leptons.size(); i++)
      cout << leptons[i]->Pt() << endl;
  }

  // Minimun Pt, Nleptons>=2 requirements
  if (leptons.size() >= 2 &&
      leptons[0]->Pt() > 20 && leptons[1]->Pt() > 15){

    CompositeParticle *dilepton = new CompositeParticle();
    dilepton->AddDaughter(leptons[0]);
    dilepton->AddDaughter(leptons[1]);

    // Sort and count the number of central Jets for vetoing
    vector<Jet*> sortedJets;
    int nCentralJets = 0;
    for(UInt_t i=0; i<CleanJets->GetEntries(); i++){
      if(fabs(CleanJets->At(i)->Eta()) < 2.5){
        nCentralJets++;
        Jet* jet_f = new Jet(CleanJets->At(i)->Px()*CleanJets->At(i)->L2RelativeCorrectionScale()*CleanJets->At(i)->L3AbsoluteCorrectionScale(),
    			     CleanJets->At(i)->Py()*CleanJets->At(i)->L2RelativeCorrectionScale()*CleanJets->At(i)->L3AbsoluteCorrectionScale(),
    			     CleanJets->At(i)->Pz()*CleanJets->At(i)->L2RelativeCorrectionScale()*CleanJets->At(i)->L3AbsoluteCorrectionScale(),
        		     CleanJets->At(i)->E() *CleanJets->At(i)->L2RelativeCorrectionScale()*CleanJets->At(i)->L3AbsoluteCorrectionScale());
        sortedJets.push_back(jet_f);
      }
    }
    for(UInt_t i=0; i<sortedJets.size(); i++){
      for(UInt_t j=i+1; j<sortedJets.size(); j++){
        if(sortedJets[i]->Pt() < sortedJets[j]->Pt()) {
          //swap i and j
    	  Jet* tempjet = sortedJets[i];
    	  sortedJets[i] = sortedJets[j];
    	  sortedJets[j] = tempjet;    
        }
      }
    }

    // Delta phi between the 2 leptons in degrees
    double deltaPhiLeptons = MathUtils::DeltaPhi(leptons[0]->Phi(), 
    						 leptons[1]->Phi())* 180./TMath::Pi();

    double deltaEtaLeptons = abs(leptons[0]->Eta() - leptons[1]->Eta());

    double deltaPhiDileptonMet = MathUtils::DeltaPhi(caloMet->Phi(), 
    						     dilepton->Phi())* 180./TMath::Pi();

    double mtHiggs = TMath::Sqrt(2.0*dilepton->Pt() * caloMet->Pt()*
        			 (1.0 - cos(deltaPhiDileptonMet * TMath::Pi()/180.0)));

    // Angle between MET and closest lepton
    double deltaPhiMetLepton[2] = {MathUtils::DeltaPhi(caloMet->Phi(), leptons[0]->Phi()),
    				   MathUtils::DeltaPhi(caloMet->Phi(), leptons[1]->Phi())};

    double minDeltaPhiMetLepton = (deltaPhiMetLepton[0] < deltaPhiMetLepton[1])?
      deltaPhiMetLepton[0]:deltaPhiMetLepton[1];
    minDeltaPhiMetLepton = minDeltaPhiMetLepton * 180./TMath::Pi();

    double METdeltaPhilEt = caloMet->Pt();
    if(minDeltaPhiMetLepton < 90.0) 
      METdeltaPhilEt = METdeltaPhilEt * sin(minDeltaPhiMetLepton * TMath::Pi()/180.);

    double mTW[2] = {TMath::Sqrt(2.0*leptons[0]->Pt()*caloMet->Pt()*
    				 (1.0 - cos(deltaPhiMetLepton[0]))),
        	     TMath::Sqrt(2.0*leptons[1]->Pt()*caloMet->Pt()*
    				 (1.0 - cos(deltaPhiMetLepton[1])))};

    int pairType = -1;
    if (leptonType[0] == "mu" && leptonType[1] == "mu" )
      pairType = 0;
    else if(leptonType[0] == "e" && leptonType[1] == "e")
      pairType = 1;
    else if((leptonType[0] == "e" && leptonType[1] == "mu") || 
            (leptonType[0] == "mu" && leptonType[1] == "e"))
      pairType = 2;
    else {
      cout << "Hey, this is not possible, leptonTypes: "
    	   << leptonType[0] << " - " 
           << leptonType[1] << endl;
    }

    // No more than 2 isolated good leptons
    if (leptons.size() == 2){
      if(dilepton->Charge() != 0) hDttPresel[ 0+100*pairType]->Fill(dilepton->Mass());
      else			  hDttPresel[ 1+100*pairType]->Fill(dilepton->Mass());
      LoadBranch(fVertexName);

     // Njets == 0, Preselection level
     if(nCentralJets == 0 &&
       dilepton->Mass() > 12 && dilepton->Charge() == 0 &&
       caloMet->Pt() > 40){
       LoadBranch(fTrackName);
       int nCleanTracks[3] = {0, 0, 0};    
       for(UInt_t i=0; i<fTracks->GetEntries(); ++i){
    	 bool trackIsLepton = false;
    	 if(MathUtils::DeltaR(fTracks->At(i)->Eta(), fTracks->At(i)->Phi(),
    			      leptons[0]->Eta(), leptons[0]->Phi()) < 0.03 ||
    	    MathUtils::DeltaR(fTracks->At(i)->Eta(),fTracks->At(i)->Phi(),
    			       leptons[1]->Eta(), leptons[1]->Phi()) < 0.03)
    	   trackIsLepton = true;

    	 if(!trackIsLepton && fTracks->At(i)->Pt() > 2.0){
    	   hDttJetVeto[0]->Fill(TMath::Min((double)fTracks->At(i)->NHits(),29.499));
    	   hDttJetVeto[1]->Fill(TMath::Min(fabs(zAverage-fTracks->At(i)->Z0()),1.999));
    	   hDttJetVeto[2]->Fill(TMath::Min(fTracks->At(i)->Pt(),49.999));
    	 }
    	 if(!trackIsLepton && fTracks->At(i)->Pt() > 3.0 && fTracks->At(i)->NHits() >= 8 &&
    	    fabs(zAverage-fTracks->At(i)->Z0()) < 0.5){
    	   nCleanTracks[0]++;
    	 }
    	 if(!trackIsLepton && fTracks->At(i)->Pt() > 3.5 && fTracks->At(i)->NHits() >= 8 &&
    	    fabs(zAverage-fTracks->At(i)->Z0()) < 0.5){
    	   nCleanTracks[1]++;
	   double d0_real = DecayXY(fTracks->At(i), fVertices);
	   hDttVert[0]->Fill(TMath::Min(fabs(d0_real),0.999));
	   hDttVert[1]->Fill(TMath::Min(fabs(d0_real)/fTracks->At(i)->D0Err(),19.999));
	   hDttSelD0Phi->Fill(fTracks->At(i)->Phi0() * 180./TMath::Pi(), d0_real);
    	 }
    	  if(!trackIsLepton && fTracks->At(i)->Pt() > 4.0 && fTracks->At(i)->NHits() >= 8 &&
    	    fabs(zAverage-fTracks->At(i)->Z0()) < 0.5){
    	   nCleanTracks[2]++;
    	 }
       }
       hDttJetVeto[3]->Fill(TMath::Min((double)nCleanTracks[0],19.4999));
       hDttJetVeto[4]->Fill(TMath::Min((double)nCleanTracks[1],19.4999));
       hDttJetVeto[5]->Fill(TMath::Min((double)nCleanTracks[2],19.4999));
       if(DirtyMuons->GetEntries() == 0){
         hDttJetVeto[6]->Fill(TMath::Min((double)nCleanTracks[0],19.4999));
         hDttJetVeto[7]->Fill(TMath::Min((double)nCleanTracks[1],19.4999));
         hDttJetVeto[8]->Fill(TMath::Min((double)nCleanTracks[2],19.4999));
       }
     } // Njets == 0

      // Z veto
      if((pairType == 2 || fabs(dilepton->Mass()-91.1876) > 20) &&
         dilepton->Mass() > 12 && dilepton->Charge() == 0){
    	// MET requirements
    	hDttPresel[ 4+100*pairType]->Fill(caloMet->Pt());
    	if((pairType == 2 && caloMet->Pt() > 40) || caloMet->Pt() > 50){
          hDttPresel[17+100*pairType]->Fill(180.-deltaPhiDileptonMet);
    	  // deltaPhiDileptonMet requirement
    	  if(180.-deltaPhiDileptonMet > 15.0){
            hDttPresel[14+100*pairType]->Fill(TMath::Min(mTW[0],mTW[1]));
            // mTW min requirement
    	    if(TMath::Min(mTW[0],mTW[1]) > 5.0){
              hDttPresel[ 2+100*pairType]->Fill(leptons[0]->Pt());
              hDttPresel[ 3+100*pairType]->Fill(leptons[1]->Pt());
              hDttPresel[ 5+100*pairType]->Fill((double)nCentralJets);
              hDttPresel[ 6+100*pairType]->Fill(deltaEtaLeptons);
              hDttPresel[ 7+100*pairType]->Fill(deltaPhiLeptons);
              hDttPresel[ 8+100*pairType]->Fill(mtHiggs);
              hDttPresel[ 9+100*pairType]->Fill(minDeltaPhiMetLepton);
              hDttPresel[10+100*pairType]->Fill(METdeltaPhilEt);
              hDttPresel[11+100*pairType]->Fill(caloMet->MetSig());
              hDttPresel[12+100*pairType]->Fill(caloMet->SumEt());
              hDttPresel[13+100*pairType]->Fill(TMath::Max(mTW[0],mTW[1]));
              hDttPresel[15+100*pairType]->Fill(leptons[0]->Pt()+leptons[1]->Pt()+caloMet->Pt());
              hDttPresel[16+100*pairType]->Fill(TMath::Min(caloMet->Pt()/dilepton->Pt(),3.99));
              hDttPreselAlphaEP->Fill(180.-deltaPhiDileptonMet,TMath::Min(caloMet->Pt()/dilepton->Pt(),3.99));
              for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
        	hDttPresel[18+100*pairType]->Fill(TMath::Min(fabs(CleanMuons->At(j)->GlobalTrk()->D0()),0.3999));
              }
              for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {	
        	hDttPresel[19+100*pairType]->Fill(TMath::Min(fabs(CleanElectrons->At(j)->GsfTrk()->D0()),0.3999));
              }

              // Study for genuine tt->WbWb->2l2nubb
              if(nCentralJets >= 2){
        	hDttPresel[20+100*pairType]->Fill(TMath::Max(fabs(sortedJets[0]->Eta()),
    							     fabs(sortedJets[1]->Eta())));
        	hDttPresel[21+100*pairType]->Fill(TMath::Min(fabs(sortedJets[0]->Eta()),
    							     fabs(sortedJets[1]->Eta())));
        	hDttPresel[22+100*pairType]->Fill(MathUtils::DeltaPhi(sortedJets[0]->Phi(),sortedJets[1]->Phi())
    						  * 180.0 / TMath::Pi());
        	hDttPresel[23+100*pairType]->Fill(fabs(sortedJets[0]->Eta()-sortedJets[1]->Eta()));
        	hDttPresel[24+100*pairType]->Fill(MathUtils::DeltaR(sortedJets[0]->Phi(),sortedJets[0]->Eta(),
    								    sortedJets[1]->Phi(),sortedJets[1]->Eta()));
        	hDttPresel[25+100*pairType]->Fill(TMath::Max(fabs(sortedJets[0]->Pt()),
    							     fabs(sortedJets[1]->Pt())));
        	hDttPresel[26+100*pairType]->Fill(TMath::Min(fabs(sortedJets[0]->Pt()),
    							     fabs(sortedJets[1]->Pt())));
        	CompositeParticle *dijet = new CompositeParticle();
        	dijet->AddDaughter(sortedJets[0]);
        	dijet->AddDaughter(sortedJets[1]);
        	hDttPresel[27+100*pairType]->Fill(dijet->Mass());
        	hDttPresel[28+100*pairType]->Fill(dijet->Pt());

        	// Study energy scale corrections
        	if(GenQuarks->GetEntries() >= 2){
        	  double deltarAux[2];
        	  double deltaPt;
        	  deltarAux[0] = MathUtils::DeltaR(sortedJets[0]->Phi(),sortedJets[0]->Eta(),
    						   GenQuarks->At(0)->Phi(),GenQuarks->At(0)->Eta());
        	  deltarAux[1] = MathUtils::DeltaR(sortedJets[0]->Phi(),sortedJets[0]->Eta(),
    						   GenQuarks->At(1)->Phi(),GenQuarks->At(1)->Eta());
        	  hDttPresel[29+100*pairType]->Fill(TMath::Min(deltarAux[0],deltarAux[1]));
        	  deltaPt = -999;
        	  if	 (deltarAux[0] < deltarAux[1] && deltarAux[0] < 1.0){
        	    deltaPt = (sortedJets[0]->Pt()-GenQuarks->At(0)->Pt())/sortedJets[0]->Pt();
        	  }
        	  else if(deltarAux[1] < deltarAux[0] && deltarAux[1] < 1.0){
        	    deltaPt = (sortedJets[0]->Pt()-GenQuarks->At(1)->Pt())/sortedJets[0]->Pt();
        	  }
        	  if(deltaPt != -999) hDttPresel[30+100*pairType]->Fill(deltaPt);
        	  deltarAux[0] = MathUtils::DeltaR(sortedJets[1]->Phi(),sortedJets[1]->Eta(),
    						   GenQuarks->At(0)->Phi(),GenQuarks->At(0)->Eta());
        	  deltarAux[1] = MathUtils::DeltaR(sortedJets[1]->Phi(),sortedJets[1]->Eta(),
    						   GenQuarks->At(1)->Phi(),GenQuarks->At(1)->Eta());
        	  hDttPresel[31+100*pairType]->Fill(TMath::Min(deltarAux[0],deltarAux[1]));
        	  deltaPt = -999;
        	  if	 (deltarAux[0] < deltarAux[1] && deltarAux[0] < 1.0){
        	    deltaPt = (sortedJets[1]->Pt()-GenQuarks->At(0)->Pt())/sortedJets[1]->Pt();
        	  }
        	  else if(deltarAux[1] < deltarAux[0] && deltarAux[1] < 1.0){
        	    deltaPt = (sortedJets[1]->Pt()-GenQuarks->At(1)->Pt())/sortedJets[1]->Pt();
        	  }
        	  if(deltaPt != -999) hDttPresel[32+100*pairType]->Fill(deltaPt);
        	}
        	delete dijet;
              } // Njets >=2

              // Study for additional muons
              hDttPresel[33+100*pairType]->Fill((double)DirtyMuons->GetEntries());
              for (UInt_t i=0; i<DirtyMuons->GetEntries(); ++i) {
        	Muon *mu = DirtyMuons->At(i);
        	hDttPresel[34+100*pairType]->Fill(TMath::Min(mu->Pt(),49.999));
        	hDttPresel[35+100*pairType]->Fill(TMath::Min(fabs(mu->GlobalTrk()->D0()),0.3999));
    		bool isGenLepton = false;
    		for (UInt_t j=0; j<GenAllLeptons->GetEntries(); ++j) {
    		  MCParticle *gen = GenAllLeptons->At(j);
    		  if(MathUtils::DeltaR(gen->Eta(), gen->Phi(), mu->Eta(), mu->Phi()) < 0.1)
    		    isGenLepton = true;
    		}
		if(isGenLepton == true){
        	  hDttPresel[36+100*pairType]->Fill(TMath::Min(mu->Pt(),49.999));
		  if(myMuonTools.isGood(mu, MuonTools::TMOneStationLoose))
		    hDttPresel[37+100*pairType]->Fill(TMath::Min(mu->Pt(),49.999));
		  if(myMuonTools.isGood(mu, MuonTools::TMOneStationTight))
		    hDttPresel[38+100*pairType]->Fill(TMath::Min(mu->Pt(),49.999));
		  if(myMuonTools.isGood(mu, MuonTools::TM2DCompatibilityLoose))
		    hDttPresel[39+100*pairType]->Fill(TMath::Min(mu->Pt(),49.999));
		  if(myMuonTools.isGood(mu, MuonTools::TM2DCompatibilityTight))
		    hDttPresel[40+100*pairType]->Fill(TMath::Min(mu->Pt(),49.999));
		  if(myMuonTools.isGood(mu, MuonTools::TMOneStationLoose) &&
		     myMuonTools.isGood(mu, MuonTools::TM2DCompatibilityLoose))
		    hDttPresel[41+100*pairType]->Fill(TMath::Min(mu->Pt(),49.999));
		  if(myMuonTools.isGood(mu, MuonTools::TMOneStationTight) &&
		     myMuonTools.isGood(mu, MuonTools::TM2DCompatibilityTight))
		    hDttPresel[42+100*pairType]->Fill(TMath::Min(mu->Pt(),49.999));
	        } else {
        	  hDttPresel[43+100*pairType]->Fill(TMath::Min(mu->Pt(),49.999));
		  if(myMuonTools.isGood(mu, MuonTools::TMOneStationLoose))
		    hDttPresel[44+100*pairType]->Fill(TMath::Min(mu->Pt(),49.999));
		  if(myMuonTools.isGood(mu, MuonTools::TMOneStationTight))
		    hDttPresel[45+100*pairType]->Fill(TMath::Min(mu->Pt(),49.999));
		  if(myMuonTools.isGood(mu, MuonTools::TM2DCompatibilityLoose))
		    hDttPresel[46+100*pairType]->Fill(TMath::Min(mu->Pt(),49.999));
		  if(myMuonTools.isGood(mu, MuonTools::TM2DCompatibilityTight))
		    hDttPresel[47+100*pairType]->Fill(TMath::Min(mu->Pt(),49.999));
		  if(myMuonTools.isGood(mu, MuonTools::TMOneStationLoose) &&
		     myMuonTools.isGood(mu, MuonTools::TM2DCompatibilityLoose))
		    hDttPresel[48+100*pairType]->Fill(TMath::Min(mu->Pt(),49.999));
		  if(myMuonTools.isGood(mu, MuonTools::TMOneStationTight) &&
		     myMuonTools.isGood(mu, MuonTools::TM2DCompatibilityTight))
		    hDttPresel[49+100*pairType]->Fill(TMath::Min(mu->Pt(),49.999));
		}
	        double d0_real = DecayXY(DirtyMuons->At(i)->GlobalTrk(), fVertices);
	        hDttVert[2]->Fill(TMath::Min(fabs(d0_real),0.999));
	        hDttVert[3]->Fill(TMath::Min(fabs(d0_real)/DirtyMuons->At(i)->GlobalTrk()->D0Err(),19.999));
              }
              if(DirtyMuons->GetEntries() > 0){
        	hDttPresel[50+100*pairType]->Fill((double)nCentralJets);
              }

              // Study dijet mass
              if(nCentralJets >= 2){
        	double dijetMass = 999.;
        	int indDiJet[2] = {-1, -1};
        	for(UInt_t i=0; i<sortedJets.size(); i++){
        	  for(UInt_t j=i+1; j<sortedJets.size(); j++){
        	    CompositeParticle *dijetTemp = new CompositeParticle();
        	    dijetTemp->AddDaughter(sortedJets[i]);
        	    dijetTemp->AddDaughter(sortedJets[j]);
        	    if(fabs(dijetTemp->Mass()-91.1876) < fabs(dijetMass-91.1876)){
        	      dijetMass = dijetTemp->Mass();
        	      indDiJet[0] = i;
        	      indDiJet[1] = j;
        	    }
        	    delete dijetTemp;
        	  }
        	}
        	hDttDiJetSel[ 0+100*pairType]->Fill(dijetMass);
        	if(indDiJet[0] == 0 && indDiJet[1] == 1)
        	  hDttDiJetSel[ 1+100*pairType]->Fill(dijetMass);
        	else
        	  hDttPresel[ 2+100*pairType]->Fill(dijetMass);
        	hDttDiJetSel[ 3+100*pairType]->Fill((double)(indDiJet[0]+10*indDiJet[1]));
              } // Study dijet mass

    	      // Njets == 0
              if(nCentralJets == 0){
    		LoadBranch(fTrackName);
    		// Get all Tracks not corresponding to leptons
        	// Applying additional cuts
    		//ObjArray<Track> *CleanTracks = new ObjArray<Track>;
    		int nCleanTracks = 0;	 
    		for(UInt_t i=0; i<fTracks->GetEntries(); ++i){
    		  bool trackIsLepton = false;
    		  if(MathUtils::DeltaR(fTracks->At(i)->Eta(), fTracks->At(i)->Phi(),
        			       leptons[0]->Eta(), leptons[0]->Phi()) < 0.03 ||
    		     MathUtils::DeltaR(fTracks->At(i)->Eta(),fTracks->At(i)->Phi(),
        				leptons[1]->Eta(), leptons[1]->Phi()) < 0.03)
        	    trackIsLepton = true;

    		  if(!trackIsLepton && fTracks->At(i)->Pt() > 2.0){
        	    hDttNoJetsSel[ 0+100*pairType]->Fill(TMath::Min((double)fTracks->At(i)->NHits(),29.499));
        	    hDttNoJetsSel[ 1+100*pairType]->Fill(TMath::Min(fabs(zAverage-fTracks->At(i)->Z0()),1.999));
        	    hDttNoJetsSel[ 2+100*pairType]->Fill(TMath::Min(fTracks->At(i)->Pt(),49.999));
    		  }
    		  if(!trackIsLepton && fTracks->At(i)->Pt() > 3.0 && fTracks->At(i)->NHits() >= 8 &&
        	     fabs(zAverage-fTracks->At(i)->Z0()) < 0.5){
    		    //CleanTracks->Add(fTracks->At(i));
        	    nCleanTracks++;
    		  }
    		}
        	hDttNoJetsSel[ 3+100*pairType]->Fill(TMath::Min((double)nCleanTracks,19.4999));
        	if(DirtyMuons->GetEntries() == 0)
		     hDttNoJetsSel[ 4+100*pairType]->Fill(TMath::Min((double)nCleanTracks,19.4999));
    		//delete CleanTracks;
    	      } // Njets == 0

            } // mTW min requirement
    	  } // deltaPhiDileptonMet requirement
    	} // MET requirements
      } // Z veto && q1+q2==0 && mass_ll>12
    } // Nleptons==2
    delete dilepton;
    for(UInt_t i=0; i<sortedJets.size(); i++) delete sortedJets[i];
  } // Minimun Pt, Nleptons>=2 requirements
  delete DirtyMuons;
  leptons.clear();
}
//--------------------------------------------------------------------------------------------------
void ttEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fMetName,    fMet);
  ReqBranch(fMuonName,   fMuons);
  ReqBranch(fTrackName,  fTracks);
  ReqBranch(fVertexName, fVertices);

  char sb[200];
  for(int j=0; j<3; j++){
    int ind = 100 * j;
    sprintf(sb,"hDttPresel_%d",ind+0);  hDttPresel[ind+0]  = new TH1D(sb,sb,150,0.0,300.); 
    sprintf(sb,"hDttPresel_%d",ind+1);  hDttPresel[ind+1]  = new TH1D(sb,sb,150,0.0,300.);
    sprintf(sb,"hDttPresel_%d",ind+2);  hDttPresel[ind+2]  = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDttPresel_%d",ind+3);  hDttPresel[ind+3]  = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDttPresel_%d",ind+4);  hDttPresel[ind+4]  = new TH1D(sb,sb,150,0.0,300.);
    sprintf(sb,"hDttPresel_%d",ind+5);  hDttPresel[ind+5]  = new TH1D(sb,sb,10,-0.5,9.5); 
    sprintf(sb,"hDttPresel_%d",ind+6);  hDttPresel[ind+6]  = new TH1D(sb,sb,50,0.0,5.);  
    sprintf(sb,"hDttPresel_%d",ind+7);  hDttPresel[ind+7]  = new TH1D(sb,sb,90,0.0,180.); 
    sprintf(sb,"hDttPresel_%d",ind+8);  hDttPresel[ind+8]  = new TH1D(sb,sb,150,0.0,300.);
    sprintf(sb,"hDttPresel_%d",ind+9);  hDttPresel[ind+9]  = new TH1D(sb,sb,90,0.0,180.); 
    sprintf(sb,"hDttPresel_%d",ind+10); hDttPresel[ind+10] = new TH1D(sb,sb,150,0.0,300.);
    sprintf(sb,"hDttPresel_%d",ind+11); hDttPresel[ind+11] = new TH1D(sb,sb,100,0.0,20.); 
    sprintf(sb,"hDttPresel_%d",ind+12); hDttPresel[ind+12] = new TH1D(sb,sb,200,0.0,800.);
    sprintf(sb,"hDttPresel_%d",ind+13); hDttPresel[ind+13] = new TH1D(sb,sb,300,0.0,300.);
    sprintf(sb,"hDttPresel_%d",ind+14); hDttPresel[ind+14] = new TH1D(sb,sb,300,0.0,300.);
    sprintf(sb,"hDttPresel_%d",ind+15); hDttPresel[ind+15] = new TH1D(sb,sb,250,0.0,500.);
    sprintf(sb,"hDttPresel_%d",ind+16); hDttPresel[ind+16] = new TH1D(sb,sb,200,0.0,4.0); 
    sprintf(sb,"hDttPresel_%d",ind+17); hDttPresel[ind+17] = new TH1D(sb,sb,90,0.0,180.); 
    sprintf(sb,"hDttPresel_%d",ind+18); hDttPresel[ind+18] = new TH1D(sb,sb,200,0.0,0.4); 
    sprintf(sb,"hDttPresel_%d",ind+19); hDttPresel[ind+19] = new TH1D(sb,sb,200,0.0,0.4); 
    sprintf(sb,"hDttPresel_%d",ind+20); hDttPresel[ind+20] = new TH1D(sb,sb,25,0.0,2.5);
    sprintf(sb,"hDttPresel_%d",ind+21); hDttPresel[ind+21] = new TH1D(sb,sb,25,0.0,2.5);
    sprintf(sb,"hDttPresel_%d",ind+22); hDttPresel[ind+22] = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDttPresel_%d",ind+23); hDttPresel[ind+23] = new TH1D(sb,sb,50,0.0, 5.0);
    sprintf(sb,"hDttPresel_%d",ind+24); hDttPresel[ind+24] = new TH1D(sb,sb,50,0.0,10.0);
    sprintf(sb,"hDttPresel_%d",ind+25); hDttPresel[ind+25] = new TH1D(sb,sb,100,0.0,200.0);
    sprintf(sb,"hDttPresel_%d",ind+26); hDttPresel[ind+26] = new TH1D(sb,sb,100,0.0,200.0);
    sprintf(sb,"hDttPresel_%d",ind+27); hDttPresel[ind+27] = new TH1D(sb,sb,100,0.0,400.0);
    sprintf(sb,"hDttPresel_%d",ind+28); hDttPresel[ind+28] = new TH1D(sb,sb,100,0.0,400.0);
    sprintf(sb,"hDttPresel_%d",ind+29); hDttPresel[ind+29] = new TH1D(sb,sb,200,0.0,10.0);
    sprintf(sb,"hDttPresel_%d",ind+30); hDttPresel[ind+30] = new TH1D(sb,sb,80,-2.0,2.0);
    sprintf(sb,"hDttPresel_%d",ind+31); hDttPresel[ind+31] = new TH1D(sb,sb,200,0.0,10.0);
    sprintf(sb,"hDttPresel_%d",ind+32); hDttPresel[ind+32] = new TH1D(sb,sb,80,-2.0,2.0);
    sprintf(sb,"hDttPresel_%d",ind+33); hDttPresel[ind+33] = new TH1D(sb,sb,10,-0.5,9.5); 
    sprintf(sb,"hDttPresel_%d",ind+34); hDttPresel[ind+34] = new TH1D(sb,sb,200,0.0,50.0);
    sprintf(sb,"hDttPresel_%d",ind+35); hDttPresel[ind+35] = new TH1D(sb,sb,200,0.0,0.4); 
    sprintf(sb,"hDttPresel_%d",ind+36); hDttPresel[ind+36] = new TH1D(sb,sb,200,0.0,50.0);
    sprintf(sb,"hDttPresel_%d",ind+37); hDttPresel[ind+37] = new TH1D(sb,sb,200,0.0,50.0);
    sprintf(sb,"hDttPresel_%d",ind+38); hDttPresel[ind+38] = new TH1D(sb,sb,200,0.0,50.0);
    sprintf(sb,"hDttPresel_%d",ind+39); hDttPresel[ind+39] = new TH1D(sb,sb,200,0.0,50.0);
    sprintf(sb,"hDttPresel_%d",ind+40); hDttPresel[ind+40] = new TH1D(sb,sb,200,0.0,50.0);
    sprintf(sb,"hDttPresel_%d",ind+41); hDttPresel[ind+41] = new TH1D(sb,sb,200,0.0,50.0);
    sprintf(sb,"hDttPresel_%d",ind+42); hDttPresel[ind+42] = new TH1D(sb,sb,200,0.0,50.0);
    sprintf(sb,"hDttPresel_%d",ind+43); hDttPresel[ind+43] = new TH1D(sb,sb,200,0.0,50.0);
    sprintf(sb,"hDttPresel_%d",ind+44); hDttPresel[ind+44] = new TH1D(sb,sb,200,0.0,50.0);
    sprintf(sb,"hDttPresel_%d",ind+45); hDttPresel[ind+45] = new TH1D(sb,sb,200,0.0,50.0);
    sprintf(sb,"hDttPresel_%d",ind+46); hDttPresel[ind+46] = new TH1D(sb,sb,200,0.0,50.0);
    sprintf(sb,"hDttPresel_%d",ind+47); hDttPresel[ind+47] = new TH1D(sb,sb,200,0.0,50.0);
    sprintf(sb,"hDttPresel_%d",ind+48); hDttPresel[ind+48] = new TH1D(sb,sb,200,0.0,50.0);
    sprintf(sb,"hDttPresel_%d",ind+49); hDttPresel[ind+49] = new TH1D(sb,sb,200,0.0,50.0);
    sprintf(sb,"hDttPresel_%d",ind+50); hDttPresel[ind+50] = new TH1D(sb,sb,10,-0.5,9.5); 
  }

  for(int i=0; i<51; i++){
    for(int j=0; j<3; j++){
      AddOutput(hDttPresel[i+j*100]);
    }
  }

  for(int j=0; j<3; j++){
    int ind = 100 * j;
    sprintf(sb,"hDttDiJetSel_%d",ind+0); hDttDiJetSel[ind+0] = new TH1D(sb,sb,100,0.0,400.0);
    sprintf(sb,"hDttDiJetSel_%d",ind+1); hDttDiJetSel[ind+1] = new TH1D(sb,sb,100,0.0,400.0);
    sprintf(sb,"hDttDiJetSel_%d",ind+2); hDttDiJetSel[ind+2] = new TH1D(sb,sb,100,0.0,400.0);
    sprintf(sb,"hDttDiJetSel_%d",ind+3); hDttDiJetSel[ind+3] = new TH1D(sb,sb,100,-0.5,99.5); 
  }

  for(int i=0; i<4; i++){
    for(int j=0; j<3; j++){
      AddOutput(hDttDiJetSel[i+j*100]);
    }
  }

  for(int j=0; j<3; j++){
    int ind = 100 * j;
    sprintf(sb,"hDttNoJetsSel_%d",ind+0); hDttNoJetsSel[ind+0] = new TH1D(sb,sb,30,-0.5,29.5); 
    sprintf(sb,"hDttNoJetsSel_%d",ind+1); hDttNoJetsSel[ind+1] = new TH1D(sb,sb,200,0.0,2.0);
    sprintf(sb,"hDttNoJetsSel_%d",ind+2); hDttNoJetsSel[ind+2] = new TH1D(sb,sb,100,0.0,50.0);
    sprintf(sb,"hDttNoJetsSel_%d",ind+3); hDttNoJetsSel[ind+3] = new TH1D(sb,sb,20,-0.5,19.5); 
    sprintf(sb,"hDttNoJetsSel_%d",ind+4); hDttNoJetsSel[ind+4] = new TH1D(sb,sb,20,-0.5,19.5); 
  }

  for(int i=0; i<5; i++){
    for(int j=0; j<3; j++){
      AddOutput(hDttNoJetsSel[i+j*100]);
    }
  }

  sprintf(sb,"hDttPreselAlphaEP");
       hDttPreselAlphaEP = new TH2D(sb,sb,36,0.0,180.0,40,0.0,4.0);  
  AddOutput(hDttPreselAlphaEP);

  sprintf(sb,"hDttJetVeto_%d",0); hDttJetVeto[0] = new TH1D(sb,sb,30,-0.5,29.5); 
  sprintf(sb,"hDttJetVeto_%d",1); hDttJetVeto[1] = new TH1D(sb,sb,200,0.0,2.0);
  sprintf(sb,"hDttJetVeto_%d",2); hDttJetVeto[2] = new TH1D(sb,sb,100,0.0,50.0);
  sprintf(sb,"hDttJetVeto_%d",3); hDttJetVeto[3] = new TH1D(sb,sb,20,-0.5,19.5); 
  sprintf(sb,"hDttJetVeto_%d",4); hDttJetVeto[4] = new TH1D(sb,sb,20,-0.5,19.5);
  sprintf(sb,"hDttJetVeto_%d",5); hDttJetVeto[5] = new TH1D(sb,sb,20,-0.5,19.5); 
  sprintf(sb,"hDttJetVeto_%d",6); hDttJetVeto[6] = new TH1D(sb,sb,20,-0.5,19.5); 
  sprintf(sb,"hDttJetVeto_%d",7); hDttJetVeto[7] = new TH1D(sb,sb,20,-0.5,19.5); 
  sprintf(sb,"hDttJetVeto_%d",8); hDttJetVeto[8] = new TH1D(sb,sb,20,-0.5,19.5); 
  for(int i=0; i<9; i++){
    AddOutput(hDttJetVeto[i]);
  }

  sprintf(sb,"hDttVert_%d",0);  hDttVert[0]  = new TH1D(sb,sb,1000,0.0,1.);
  sprintf(sb,"hDttVert_%d",1);  hDttVert[1]  = new TH1D(sb,sb,200,0.0,20.);
  sprintf(sb,"hDttVert_%d",2);  hDttVert[2]  = new TH1D(sb,sb,1000,0.0,1.);
  sprintf(sb,"hDttVert_%d",3);  hDttVert[3]  = new TH1D(sb,sb,200,0.0,20.);
  for(int i=0; i<4; i++){
    AddOutput(hDttVert[i]);
  }
  sprintf(sb,"hDttSelD0Phi");
       hDttSelD0Phi = new TH2D(sb,sb,90,-180.0,180.0,200,-0.2,0.2);  
  AddOutput(hDttSelD0Phi);
}

//--------------------------------------------------------------------------------------------------
void ttEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void ttEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}

//--------------------------------------------------------------------------------------------------
double ttEvtSelMod::DecayXY(const mithep::Track *lTrack, mithep::Vertex *iVertex) {
  if(lTrack == 0) return 999999;

  double lXM =  -sin(lTrack->Phi()) * (lTrack->D0());
  double lYM =  cos(lTrack->Phi()) * (lTrack->D0());
  double lDX = (lXM + iVertex->X()); 
  double lDY = (lYM + iVertex->Y());
  return (lTrack->Px()*lDY - lTrack->Py()*lDX) / lTrack->Pt();
}

//--------------------------------------------------------------------------------------------------
double ttEvtSelMod::DecayXY(const mithep::Track *lTrack, mithep::VertexCol *iVertices) {
  double lD0 = 10000;

  for(uint i0 = 0; i0 < iVertices->GetEntries(); i0++) {
    double pD0 = DecayXY(lTrack,iVertices->At(i0));
    if(fabs(pD0) < fabs(lD0)) lD0 = pD0;
  }

  if(lD0 == 10000) return -1;
  return lD0;
}
