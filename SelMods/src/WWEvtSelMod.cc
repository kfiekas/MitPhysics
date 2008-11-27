 // $Id: WWEvtSelMod.cc,v 1.1 2008/11/11 21:22:54 ceballos Exp $

#include "MitPhysics/SelMods/interface/WWEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

using namespace mithep;
ClassImp(mithep::WWEvtSelMod)

//--------------------------------------------------------------------------------------------------
WWEvtSelMod::WWEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(false),
  fMetName(Names::gkCaloMetBrn),
  fMuonName(Names::gkMuonBrn),
  fTrackName(Names::gkTrackBrn),
  fVertexName(string("PrimaryVertexesBeamSpot").c_str()),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fMCLeptonsName(ModNames::gkMCLeptonsName),
  fMCNeutrinosName(ModNames::gkMCNeutrinosName),
  fMet(0),
  fMuons(0),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void WWEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void WWEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  if (fNEventsProcessed % 1000000 == 0 || fPrintDebug) {
    time_t systime;
    systime = time(NULL);
    cerr << endl << "WWEvtSelMod : Process Event " << fNEventsProcessed << "  Time: " << ctime(&systime) << endl;
  }

  //Get Generator Level information for matching
  //ObjArray<MCParticle> *GenLeptons   = dynamic_cast<ObjArray<MCParticle>* > (FindObjThisEvt(fMCLeptonsName.Data()));
  //ObjArray<MCParticle> *GenNeutrinos = dynamic_cast<ObjArray<MCParticle>* > (FindObjThisEvt(fMCNeutrinosName.Data()));

  //Obtain all the good objects from the event cleaning module
  ObjArray<Electron> *CleanElectrons = dynamic_cast<ObjArray<Electron>* >
    (FindObjThisEvt(ModNames::gkCleanElectronsName));
  ObjArray<Muon> *CleanMuons = dynamic_cast<ObjArray<Muon>* >
    (FindObjThisEvt(ModNames::gkCleanMuonsName));
  ObjArray<Jet> *CleanJets = dynamic_cast<ObjArray<Jet>* >
    (FindObjThisEvt(fCleanJetsName.Data()));

  LoadBranch(fMetName);
  Met *caloMet = fMet->At(0);

  vector<ChargedParticle*> leptons;
  vector<string> leptonType; 
  double zAverage = 0.0;

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

  hDwwPresel[0]->Fill(TMath::Min((double)leptons.size(),9.499));
  if(leptons.size() >= 1)
    hDwwPresel[1]->Fill(TMath::Min(leptons[0]->Pt(),199.999));
  if(leptons.size() >= 2)
    hDwwPresel[2]->Fill(TMath::Min(leptons[1]->Pt(),199.999));
  // Minimun Pt, Nleptons>=2 requirements
  if (leptons.size() == 2 &&
      leptons[0]->Pt() > 20 && leptons[1]->Pt() > 15){

    CompositeParticle *dilepton = new CompositeParticle();
    dilepton->AddDaughter(leptons[0]);
    dilepton->AddDaughter(leptons[1]);

    // Sort and count the number of central Jets for vetoing
    int nCentralJets = 0;
    for(UInt_t i=0; i<CleanJets->GetEntries(); i++){
      if(fabs(CleanJets->At(i)->Eta()) < 2.5){
        nCentralJets++;
      }
    }

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

    hDwwPresel[3]->Fill(TMath::Min(dilepton->Mass(),199.999));
    hDwwPresel[4]->Fill((double)nCentralJets);
    hDwwPresel[5]->Fill((double)dilepton->Charge());
    hDwwPresel[6]->Fill(TMath::Min(caloMet->Pt(),199.999));

    // Njets == 0, Preselection level
    if((pairType == 2 || fabs(dilepton->Mass()-91.1876) > 20) &&
      nCentralJets == 0 && caloMet->Pt() > 40 &&
      dilepton->Mass() > 12 && dilepton->Charge() == 0){
      LoadBranch(fVertexName);
      LoadBranch(fTrackName);
      int nCleanTracks = 0;
      for(UInt_t i=0; i<fTracks->GetEntries(); ++i){
        bool trackIsLepton = false;
        if(MathUtils::DeltaR(fTracks->At(i)->Eta(), fTracks->At(i)->Phi(),
          		     leptons[0]->Eta(), leptons[0]->Phi()) < 0.03 ||
           MathUtils::DeltaR(fTracks->At(i)->Eta(),fTracks->At(i)->Phi(),
        		     leptons[1]->Eta(), leptons[1]->Phi()) < 0.03)
          trackIsLepton = true;

        if(!trackIsLepton && fTracks->At(i)->Pt() > 3.5 && fTracks->At(i)->NHits() >= 8 &&
           fabs(zAverage-fTracks->At(i)->Z0()) < 0.5){
          nCleanTracks++;
        }
	double d0_real = DecayXY(fTracks->At(i), fVertices);
	hDwwVert[0]->Fill(TMath::Min(fabs(d0_real),0.999));
	hDwwVert[1]->Fill(TMath::Min(fabs(d0_real)/fTracks->At(i)->D0Err(),19.999));
	hDwwSelD0Phi->Fill(fTracks->At(i)->Phi0() * 180./TMath::Pi(), d0_real);
      }
      hDwwPresel[7]->Fill(TMath::Min((double)nCleanTracks,19.4999));
      hDwwPresel[8]->Fill((double)DirtyMuons->GetEntries());

      for(UInt_t i=0; i<DirtyMuons->GetEntries(); i++){
	double d0_real = DecayXY(DirtyMuons->At(i)->GlobalTrk(), fVertices);
	hDwwVert[2]->Fill(TMath::Min(fabs(d0_real),0.999));
	hDwwVert[3]->Fill(TMath::Min(fabs(d0_real)/DirtyMuons->At(i)->GlobalTrk()->D0Err(),19.999));
      }
      // Ntracks <= 3 and NDirtyMuons == 0
      if(nCleanTracks <= 3 && DirtyMuons->GetEntries() == 0){
        // Delta phi between the 2 leptons in degrees
        double deltaPhiLeptons = MathUtils::DeltaPhi(leptons[0]->Phi(), 
    	  					     leptons[1]->Phi())* 180./TMath::Pi();

        double deltaPhiDileptonMet = MathUtils::DeltaPhi(caloMet->Phi(), 
      						         dilepton->Phi())* 180./TMath::Pi();

        double mtHiggs = TMath::Sqrt(2.0*dilepton->Pt() * caloMet->Pt()*
        			    (1.0 - cos(deltaPhiDileptonMet * TMath::Pi()/180.0)));

        // Angle between MET and closest lepton
        double deltaPhiMetLepton[2] = {MathUtils::DeltaPhi(caloMet->Phi(), leptons[0]->Phi()),
    				       MathUtils::DeltaPhi(caloMet->Phi(), leptons[1]->Phi())};

        double mTW[2] = {TMath::Sqrt(2.0*leptons[0]->Pt()*caloMet->Pt()*
    				    (1.0 - cos(deltaPhiMetLepton[0]))),
        	         TMath::Sqrt(2.0*leptons[1]->Pt()*caloMet->Pt()*
    				    (1.0 - cos(deltaPhiMetLepton[1])))};

        double minDeltaPhiMetLepton = (deltaPhiMetLepton[0] < deltaPhiMetLepton[1])?
           deltaPhiMetLepton[0]:deltaPhiMetLepton[1];

        double METdeltaPhilEt = caloMet->Pt();
        if(minDeltaPhiMetLepton < 90.0) 
          METdeltaPhilEt = METdeltaPhilEt * sin(minDeltaPhiMetLepton);

        int metBins = 50;
        for(int i=0; i<=metBins; i++){
	  double metAux[2] = {caloMet->Pt() * i/metBins*1.0,
	                      caloMet->Pt() * (metBins-i)/metBins*1.0};
          double mTW[2] = {TMath::Sqrt(2.0*leptons[0]->Pt()*metAux[0]*
    		  		      (1.0 - cos(deltaPhiMetLepton[0]))),
        	           TMath::Sqrt(2.0*leptons[1]->Pt()*metAux[1]*
    				      (1.0 - cos(deltaPhiMetLepton[1])))};
          hDwwSel[ 0+100*pairType]->Fill(TMath::Min(TMath::Max(mTW[0],mTW[1]),199.999));
          hDwwSel[ 1+100*pairType]->Fill(TMath::Min(TMath::Min(mTW[0],mTW[1]),199.999));	  
	}
	double minMass = 9000.;
	double nBin[2] = {-1, -1};
        for(int i=0; i<=metBins; i++){
          for(int j=0; j<=metBins; j++){
	    double metAuxPx[2] = {caloMet->Px() * i/metBins*1.0,
	                          caloMet->Px() * (metBins-i)/metBins*1.0};
	    double metAuxPy[2] = {caloMet->Py() * j/metBins*1.0,
	                          caloMet->Py() * (metBins-j)/metBins*1.0};
	    double phiAux[2] = {atan2(metAuxPy[0],metAuxPx[0]),
	                        atan2(metAuxPy[1],metAuxPx[1])};
	    double ptAux[2] = {sqrt(metAuxPx[0]*metAuxPx[0]+metAuxPy[0]*metAuxPy[0]),
	                       sqrt(metAuxPx[1]*metAuxPx[1]+metAuxPy[1]*metAuxPy[1])};
            double deltaPhiMetLepton[2] = {MathUtils::DeltaPhi(phiAux[0], leptons[0]->Phi()),
    				           MathUtils::DeltaPhi(phiAux[1], leptons[1]->Phi())};
            double mTW[2] = {TMath::Sqrt(2.0*leptons[0]->Pt()*ptAux[0]*
    	  	  		        (1.0 - cos(deltaPhiMetLepton[0]))),
        	             TMath::Sqrt(2.0*leptons[1]->Pt()*ptAux[1]*
    				        (1.0 - cos(deltaPhiMetLepton[1])))};
            hDwwSel[ 2+100*pairType]->Fill(TMath::Min(TMath::Max(mTW[0],mTW[1]),199.999));
            hDwwSel[ 3+100*pairType]->Fill(TMath::Min(TMath::Min(mTW[0],mTW[1]),199.999));
	    if(minMass > TMath::Max(mTW[0],mTW[1])){
	      minMass = TMath::Max(mTW[0],mTW[1]);
	      nBin[0] = i;
	      nBin[1] = j;
	    }
          }
	}
        hDwwSel[ 4+100*pairType]->Fill(TMath::Min(minMass,199.999));

	double metAuxPx[2] = {caloMet->Px() * nBin[0]/metBins*1.0,
			      caloMet->Px() * (metBins-nBin[0])/metBins*1.0};
	double metAuxPy[2] = {caloMet->Py() * nBin[1]/metBins*1.0,
			      caloMet->Py() * (metBins-nBin[1])/metBins*1.0};
	//double phiAux[2] = {atan2(metAuxPy[0],metAuxPx[0]),
	//		      atan2(metAuxPy[1],metAuxPx[1])};
	//double ptAux[2] = {sqrt(metAuxPx[0]*metAuxPx[0]+metAuxPy[0]*metAuxPy[0]),
	//		     sqrt(metAuxPx[1]*metAuxPx[1]+metAuxPy[1]*metAuxPy[1])};
        //double deltaPhiMetLeptonAux[2] = {MathUtils::DeltaPhi(phiAux[0], leptons[0]->Phi()),
    	//			            MathUtils::DeltaPhi(phiAux[1], leptons[1]->Phi())};
        double T[2] = {80.40*80.40/2.0+leptons[0]->Px()*metAuxPx[0]+leptons[0]->Py()*metAuxPy[0],
		       80.40*80.40/2.0+leptons[1]->Px()*metAuxPx[1]+leptons[1]->Py()*metAuxPy[1]};
	double S[2] = {metAuxPx[0]*metAuxPx[0]+metAuxPy[0]*metAuxPy[0],
		       metAuxPx[1]*metAuxPx[1]+metAuxPy[1]*metAuxPy[1]};
	double B[2] = {T[0]*T[0]*leptons[0]->Pz()*leptons[0]->Pz()-(leptons[0]->P()*leptons[0]->P()*S[0]-T[0]*T[0])*(leptons[0]->P()*leptons[0]->P()-leptons[0]->Pz()*leptons[0]->Pz()),
		       T[1]*T[1]*leptons[1]->Pz()*leptons[1]->Pz()-(leptons[1]->P()*leptons[1]->P()*S[1]-T[1]*T[1])*(leptons[1]->P()*leptons[1]->P()-leptons[1]->Pz()*leptons[1]->Pz())};
	for(int i=0; i<2; i++) {if(B[i] > 0) B[i] = sqrt(B[i]); else B[i] = 0.0;}
	double PzNeuP[2] = {(T[0]*leptons[0]->Pz()+B[0])/(leptons[0]->P()*leptons[0]->P()-leptons[0]->Pz()*leptons[0]->Pz()),
	                    (T[1]*leptons[1]->Pz()+B[1])/(leptons[1]->P()*leptons[1]->P()-leptons[1]->Pz()*leptons[1]->Pz())};
	double PzNeuN[2] = {(T[0]*leptons[0]->Pz()-B[0])/(leptons[0]->P()*leptons[0]->P()-leptons[0]->Pz()*leptons[0]->Pz()),
	                    (T[1]*leptons[1]->Pz()-B[1])/(leptons[1]->P()*leptons[1]->P()-leptons[1]->Pz()*leptons[1]->Pz())};
	//for(int i=0; i<2; i++) cout << B[i] << endl;
	double np[2] = {sqrt(metAuxPx[0]*metAuxPx[0]+metAuxPy[0]*metAuxPy[0]+PzNeuP[0]*PzNeuP[0]),
			sqrt(metAuxPx[1]*metAuxPx[1]+metAuxPy[1]*metAuxPy[1]+PzNeuP[1]*PzNeuP[1])};
	double nn[2] = {sqrt(metAuxPx[0]*metAuxPx[0]+metAuxPy[0]*metAuxPy[0]+PzNeuN[0]*PzNeuN[0]),
			sqrt(metAuxPx[1]*metAuxPx[1]+metAuxPy[1]*metAuxPy[1]+PzNeuN[1]*PzNeuN[1])};
	double massH[4];
        massH[0]= (np[0]+leptons[0]->P()+np[1]+leptons[1]->P())*(np[0]+leptons[0]->P()+np[1]+leptons[1]->P())
	       -(metAuxPx[0]+leptons[0]->Px()+metAuxPx[1]+leptons[1]->Px())*(metAuxPx[0]+leptons[0]->Px()+metAuxPx[1]+leptons[1]->Px())
	       -(metAuxPy[0]+leptons[0]->Py()+metAuxPy[1]+leptons[1]->Py())*(metAuxPy[0]+leptons[0]->Py()+metAuxPy[1]+leptons[1]->Py())
	       -(PzNeuP[0]  +leptons[0]->Pz()+PzNeuP[1]  +leptons[1]->Pz())*(PzNeuP[0]  +leptons[0]->Py()+PzNeuP[1]  +leptons[1]->Pz());
        massH[1]= (np[0]+leptons[0]->P()+nn[1]+leptons[1]->P())*(np[0]+leptons[0]->P()+nn[1]+leptons[1]->P())
	       -(metAuxPx[0]+leptons[0]->Px()+metAuxPx[1]+leptons[1]->Px())*(metAuxPx[0]+leptons[0]->Px()+metAuxPx[1]+leptons[1]->Px())
	       -(metAuxPy[0]+leptons[0]->Py()+metAuxPy[1]+leptons[1]->Py())*(metAuxPy[0]+leptons[0]->Py()+metAuxPy[1]+leptons[1]->Py())
	       -(PzNeuP[0]  +leptons[0]->Pz()+PzNeuN[1]  +leptons[1]->Pz())*(PzNeuP[0]  +leptons[0]->Py()+PzNeuN[1]  +leptons[1]->Pz());
        massH[2]= (nn[0]+leptons[0]->P()+np[1]+leptons[1]->P())*(nn[0]+leptons[0]->P()+np[1]+leptons[1]->P())
	       -(metAuxPx[0]+leptons[0]->Px()+metAuxPx[1]+leptons[1]->Px())*(metAuxPx[0]+leptons[0]->Px()+metAuxPx[1]+leptons[1]->Px())
	       -(metAuxPy[0]+leptons[0]->Py()+metAuxPy[1]+leptons[1]->Py())*(metAuxPy[0]+leptons[0]->Py()+metAuxPy[1]+leptons[1]->Py())
	       -(PzNeuN[0]  +leptons[0]->Pz()+PzNeuP[1]  +leptons[1]->Pz())*(PzNeuN[0]  +leptons[0]->Py()+PzNeuP[1]  +leptons[1]->Pz());
        massH[3]= (np[0]+leptons[0]->P()+nn[1]+leptons[1]->P())*(np[0]+leptons[0]->P()+nn[1]+leptons[1]->P())
	       -(metAuxPx[0]+leptons[0]->Px()+metAuxPx[1]+leptons[1]->Px())*(metAuxPx[0]+leptons[0]->Px()+metAuxPx[1]+leptons[1]->Px())
	       -(metAuxPy[0]+leptons[0]->Py()+metAuxPy[1]+leptons[1]->Py())*(metAuxPy[0]+leptons[0]->Py()+metAuxPy[1]+leptons[1]->Py())
	       -(PzNeuN[0]  +leptons[0]->Pz()+PzNeuN[1]  +leptons[1]->Pz())*(PzNeuN[0]  +leptons[0]->Py()+PzNeuN[1]  +leptons[1]->Pz());
	for(int i=0; i<4; i++) {if(massH[i] > 0) massH[i] = sqrt(massH[i]); else massH[i] = 0.0;}
        hDwwSel[ 5+100*pairType]->Fill(TMath::Min(massH[0],499.999));
        hDwwSel[ 6+100*pairType]->Fill(TMath::Min(massH[1],499.999));
        hDwwSel[ 7+100*pairType]->Fill(TMath::Min(massH[2],499.999));
        hDwwSel[ 8+100*pairType]->Fill(TMath::Min(massH[3],499.999));
        hDwwPresel[9]->Fill(TMath::Min(massH[0],499.999));

	/*
	for(int i=0; i<2; i++) {
	  double mas0 = (np[i]+leptons[i]->P())*(np[i]+leptons[i]->P())-(metAuxPx[i]+leptons[i]->Px())*(metAuxPx[i]+leptons[i]->Px())
	                                                         -(metAuxPy[i]+leptons[i]->Py())*(metAuxPy[i]+leptons[i]->Py())
	                                                         -(PzNeuP[i]+leptons[i]->Pz())  *(PzNeuP[i]+leptons[i]->Pz());
	  double mas1 = (nn[i]+leptons[i]->P())*(nn[i]+leptons[i]->P())-(metAuxPx[i]+leptons[i]->Px())*(metAuxPx[i]+leptons[i]->Px())
	                                                         -(metAuxPy[i]+leptons[i]->Py())*(metAuxPy[i]+leptons[i]->Py())
	                                                         -(PzNeuN[i]+leptons[i]->Pz())  *(PzNeuN[i]+leptons[i]->Pz());
	  cout << mas0 << " " << mas1 << endl;
	}
	*/
        hDwwSel[ 9+100*pairType]->Fill(deltaPhiLeptons);
        hDwwSel[10+100*pairType]->Fill(TMath::Min(mtHiggs,299.999));
        hDwwSel[11+100*pairType]->Fill(TMath::Min(TMath::Max(mTW[0],mTW[1]),199.999));
        hDwwSel[12+100*pairType]->Fill(TMath::Min(TMath::Min(mTW[0],mTW[1]),199.999));
        hDwwSel[13+100*pairType]->Fill(minDeltaPhiMetLepton * 180./TMath::Pi());
        hDwwSel[14+100*pairType]->Fill(TMath::Min(METdeltaPhilEt,199.999));
        hDwwSel[15+100*pairType]->Fill(TMath::Min(dilepton->Mass(),499.999));
        hDwwSel[16+100*pairType]->Fill(TMath::Min(caloMet->Pt(),199.999));
        hDwwSel[17+100*pairType]->Fill(TMath::Min(leptons[0]->Pt(),199.999));
        hDwwSel[18+100*pairType]->Fill(TMath::Min(leptons[1]->Pt(),199.999));
        hDwwSel[19+100*pairType]->Fill(180.-deltaPhiDileptonMet);
        hDwwSel[20+100*pairType]->Fill(TMath::Min(caloMet->Pt()/dilepton->Pt(),3.99));
        hDwwSelAlphaEP->Fill(180.-deltaPhiDileptonMet,TMath::Min(caloMet->Pt()/dilepton->Pt(),3.99));
      } // Ntracks <= 3 and NDirtyMuons == 0
    } // Njets == 0, Preselection level
    delete dilepton;
  } // Minimun Pt, Nleptons==2 requirements
  delete DirtyMuons;
  leptons.clear();
}
//--------------------------------------------------------------------------------------------------
void WWEvtSelMod::SlaveBegin()
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
    sprintf(sb,"hDwwSel_%d",ind+0);  hDwwSel[ind+0]  = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+1);  hDwwSel[ind+1]  = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+2);  hDwwSel[ind+2]  = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+3);  hDwwSel[ind+3]  = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+4);  hDwwSel[ind+4]  = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+5);  hDwwSel[ind+5]  = new TH1D(sb,sb,100,0.0,500.);
    sprintf(sb,"hDwwSel_%d",ind+6);  hDwwSel[ind+6]  = new TH1D(sb,sb,100,0.0,500.);
    sprintf(sb,"hDwwSel_%d",ind+7);  hDwwSel[ind+7]  = new TH1D(sb,sb,100,0.0,500.);
    sprintf(sb,"hDwwSel_%d",ind+8);  hDwwSel[ind+8]  = new TH1D(sb,sb,100,0.0,500.);
    sprintf(sb,"hDwwSel_%d",ind+9);  hDwwSel[ind+9]  = new TH1D(sb,sb,90,0.0,180.); 
    sprintf(sb,"hDwwSel_%d",ind+10); hDwwSel[ind+10] = new TH1D(sb,sb,150,0.0,300.);
    sprintf(sb,"hDwwSel_%d",ind+11); hDwwSel[ind+11] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+12); hDwwSel[ind+12] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+13); hDwwSel[ind+13] = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDwwSel_%d",ind+14); hDwwSel[ind+14] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+15); hDwwSel[ind+15] = new TH1D(sb,sb,200,0.0,400.);
    sprintf(sb,"hDwwSel_%d",ind+16); hDwwSel[ind+16] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+17); hDwwSel[ind+17] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+18); hDwwSel[ind+18] = new TH1D(sb,sb,100,0.0,200.);
    sprintf(sb,"hDwwSel_%d",ind+19); hDwwSel[ind+19] = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDwwSel_%d",ind+20); hDwwSel[ind+20] = new TH1D(sb,sb,80,0.0,4.0);
  }

  for(int i=0; i<21; i++){
    for(int j=0; j<3; j++){
      AddOutput(hDwwSel[i+j*100]);
    }
  }

  sprintf(sb,"hDwwSelAlphaEP");
       hDwwSelAlphaEP = new TH2D(sb,sb,36,0.0,180.0,40,0.0,4.0);  
  AddOutput(hDwwSelAlphaEP);

  sprintf(sb,"hDwwPresel_%d",0); hDwwPresel[0] = new TH1D(sb,sb,10,-0.5,9.5); 
  sprintf(sb,"hDwwPresel_%d",1); hDwwPresel[1] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwPresel_%d",2); hDwwPresel[2] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwPresel_%d",3); hDwwPresel[3] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwPresel_%d",4); hDwwPresel[4] = new TH1D(sb,sb,10,-0.5,9.5); 
  sprintf(sb,"hDwwPresel_%d",5); hDwwPresel[5] = new TH1D(sb,sb,9,-4.5,4.5); 
  sprintf(sb,"hDwwPresel_%d",6); hDwwPresel[6] = new TH1D(sb,sb,200,0.0,200.0);
  sprintf(sb,"hDwwPresel_%d",7); hDwwPresel[7] = new TH1D(sb,sb,20,-0.5,19.5); 
  sprintf(sb,"hDwwPresel_%d",8); hDwwPresel[8] = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDwwPresel_%d",9); hDwwPresel[9] = new TH1D(sb,sb,100,0.0,500.);
  for(int i=0; i<10; i++){
    AddOutput(hDwwPresel[i]);
  }

  sprintf(sb,"hDwwVert_%d",0);  hDwwVert[0]  = new TH1D(sb,sb,1000,0.0,1.);
  sprintf(sb,"hDwwVert_%d",1);  hDwwVert[1]  = new TH1D(sb,sb,200,0.0,20.);
  sprintf(sb,"hDwwVert_%d",2);  hDwwVert[2]  = new TH1D(sb,sb,1000,0.0,1.);
  sprintf(sb,"hDwwVert_%d",3);  hDwwVert[3]  = new TH1D(sb,sb,200,0.0,20.);
  for(int i=0; i<4; i++){
    AddOutput(hDwwVert[i]);
  }
  sprintf(sb,"hDwwSelD0Phi");
       hDwwSelD0Phi = new TH2D(sb,sb,90,-180.0,180.0,200,-0.2,0.2);  
  AddOutput(hDwwSelD0Phi);
}

//--------------------------------------------------------------------------------------------------
void WWEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void WWEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}

//--------------------------------------------------------------------------------------------------
double WWEvtSelMod::DecayXY(const mithep::Track *lTrack, mithep::Vertex *iVertex) {
  if(lTrack == 0) return 999999;

  double lXM =  -sin(lTrack->Phi()) * (lTrack->D0());
  double lYM =  cos(lTrack->Phi()) * (lTrack->D0());
  double lDX = (lXM + iVertex->X()); 
  double lDY = (lYM + iVertex->Y());
  return (lTrack->Px()*lDY - lTrack->Py()*lDX) / lTrack->Pt();
}

//--------------------------------------------------------------------------------------------------
double WWEvtSelMod::DecayXY(const mithep::Track *lTrack, mithep::VertexCol *iVertices) {
  double lD0 = 10000;

  for(uint i0 = 0; i0 < iVertices->GetEntries(); i0++) {
    double pD0 = DecayXY(lTrack,iVertices->At(i0));
    if(fabs(pD0) < fabs(lD0)) lD0 = pD0;
  }

  if(lD0 == 10000) return -1;
  return lD0;
}
