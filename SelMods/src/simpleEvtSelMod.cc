 // $Id: simpleEvtSelMod.cc,v 1.1 2008/11/11 21:22:54 ceballos Exp $

#include "MitPhysics/SelMods/interface/simpleEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/MuonTools.h"

using namespace mithep;
ClassImp(mithep::simpleEvtSelMod)

//--------------------------------------------------------------------------------------------------
simpleEvtSelMod::simpleEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(false),
  fMetName(Names::gkCaloMetBrn),
  fMuonName(Names::gkMuonBrn),
  fTrackName(Names::gkTrackBrn),
  fVertexName("PrimaryVertexesBeamSpot"),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fMet(0),
  fMuons(0),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void simpleEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void simpleEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  if (fNEventsProcessed % 1000 == 0 || fPrintDebug) {
    time_t systime;
    systime = time(NULL);
    cerr << endl << "simpleEvtSelMod : Process Event " << fNEventsProcessed << "  Time: " << ctime(&systime) << endl;
  }

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
  MuonTools myMuonTools;

  LoadBranch(fMuonName);
  // DirtyMuons == soft muons!
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

  // Minimun Pt, Nleptons == 2 requirements
  if (leptons.size() == 2 &&
      leptons[0]->Pt() > 20 && leptons[1]->Pt() > 15){

    CompositeParticle *dilepton = new CompositeParticle();
    dilepton->AddDaughter(leptons[0]);
    dilepton->AddDaughter(leptons[1]);

    // Count the number of central Jets for vetoing
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

    LoadBranch(fVertexName);

    // Njets == 0, Preselection level
    if(nCentralJets == 0 &&
      dilepton->Mass() > 12 && dilepton->Charge() == 0 &&
      caloMet->Pt() > 40){
      LoadBranch(fTrackName);

      for(UInt_t i=0; i<fTracks->GetEntries(); ++i){
        bool trackIsLepton = false;
        if(MathUtils::DeltaR(fTracks->At(i)->Eta(), fTracks->At(i)->Phi(),
        		     leptons[0]->Eta(), leptons[0]->Phi()) < 0.03 ||
           MathUtils::DeltaR(fTracks->At(i)->Eta(),fTracks->At(i)->Phi(),
        		      leptons[1]->Eta(), leptons[1]->Phi()) < 0.03)
          trackIsLepton = true;

        if(!trackIsLepton && fTracks->At(i)->Pt() > 1.0 &&
           fTracks->At(i)->NHits() >= 8 &&
           fabs(zAverage-fTracks->At(i)->Z0()) < 0.5){
          double d0_real = DecayXY(fTracks->At(i), fVertices);
          hDSelVert[0]->Fill(TMath::Min(fabs(d0_real),0.999));
          hDSelVert[1]->Fill(TMath::Min(fabs(d0_real)/fTracks->At(i)->D0Err(),19.999));
          hDSelD0Phi->Fill(fTracks->At(i)->Phi0() * 180./TMath::Pi(), d0_real);
        }
      }
      for (UInt_t i=0; i<DirtyMuons->GetEntries(); ++i) {
    	Muon *mu = DirtyMuons->At(i);
        double d0_real = DecayXY(mu->GlobalTrk(), fVertices);
        hDSelVert[2]->Fill(TMath::Min(fabs(d0_real),0.999));
        hDSelVert[3]->Fill(TMath::Min(fabs(d0_real)/
        		  DirtyMuons->At(i)->GlobalTrk()->D0Err(),19.999));
      }
    } // Njets == 0

    double deltaPhiDileptonMet = MathUtils::DeltaPhi(caloMet->Phi(), 
    						     dilepton->Phi())* 180./TMath::Pi();

    double deltaPhiMetLepton[2] = {MathUtils::DeltaPhi(caloMet->Phi(), leptons[0]->Phi()),
    				   MathUtils::DeltaPhi(caloMet->Phi(), leptons[1]->Phi())};

    double mTW[2] = {TMath::Sqrt(2.0*leptons[0]->Pt()*caloMet->Pt()*
    				 (1.0 - cos(deltaPhiMetLepton[0]))),
        	     TMath::Sqrt(2.0*leptons[1]->Pt()*caloMet->Pt()*
    				 (1.0 - cos(deltaPhiMetLepton[1])))};
    // All ttbar cuts
    if((pairType == 2 || fabs(dilepton->Mass()-91.1876) > 20) &&
        dilepton->Mass() > 12 && dilepton->Charge() == 0 &&
	(pairType == 2 && caloMet->Pt() > 40) || caloMet->Pt() > 50 &&
        180.-deltaPhiDileptonMet > 15.0 && TMath::Min(mTW[0],mTW[1]) > 5.0){
      // To be filled with anything
    } // All ttbar cuts
    delete dilepton;
  } // Minimun Pt, Nleptons>=2 requirements
  delete DirtyMuons;
  leptons.clear();
}
//--------------------------------------------------------------------------------------------------
void simpleEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fMetName,    fMet);
  ReqBranch(fMuonName,   fMuons);
  ReqBranch(fTrackName,  fTracks);
  ReqBranch(fVertexName, fVertices);

  char sb[200];
  sprintf(sb,"hDSelVert_%d",0);  hDSelVert[0]  = new TH1D(sb,sb,1000,0.0,1.);
  sprintf(sb,"hDSelVert_%d",1);  hDSelVert[1]  = new TH1D(sb,sb,200,0.0,20.);
  sprintf(sb,"hDSelVert_%d",2);  hDSelVert[2]  = new TH1D(sb,sb,1000,0.0,1.);
  sprintf(sb,"hDSelVert_%d",3);  hDSelVert[3]  = new TH1D(sb,sb,200,0.0,20.);
  for(int i=0; i<4; i++){
    AddOutput(hDSelVert[i]);
  }
  sprintf(sb,"hDSelD0Phi");
       hDSelD0Phi = new TH2D(sb,sb,90,-180.0,180.0,200,-0.2,0.2);  
  AddOutput(hDSelD0Phi);
}

//--------------------------------------------------------------------------------------------------
void simpleEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void simpleEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}

//--------------------------------------------------------------------------------------------------
double simpleEvtSelMod::DecayXY(const mithep::Track *lTrack, mithep::Vertex *iVertex) {
  if(lTrack == 0) return 999999;

  double lXM =  -sin(lTrack->Phi()) * (lTrack->D0());
  double lYM =  cos(lTrack->Phi()) * (lTrack->D0());
  double lDX = (lXM + iVertex->X()); 
  double lDY = (lYM + iVertex->Y());
  return (lTrack->Px()*lDY - lTrack->Py()*lDX) / lTrack->Pt();
}

//--------------------------------------------------------------------------------------------------
double simpleEvtSelMod::DecayXY(const mithep::Track *lTrack, mithep::VertexCol *iVertices) {
  double lD0 = 10000;

  for(uint i0 = 0; i0 < iVertices->GetEntries(); i0++) {
    double pD0 = DecayXY(lTrack,iVertices->At(i0));
    if(fabs(pD0) < fabs(lD0)) lD0 = pD0;
  }

  if(lD0 == 10000) return -1;
  return lD0;
}
