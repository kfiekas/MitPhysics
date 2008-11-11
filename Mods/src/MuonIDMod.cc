// $Id: MuonIDMod.cc,v 1.3 2008/11/05 14:06:09 ceballos Exp $

#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/MuonTools.h"

using namespace mithep;

ClassImp(mithep::MuonIDMod)

//--------------------------------------------------------------------------------------------------
  MuonIDMod::MuonIDMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(false),
  fMuonName(Names::gkMuonBrn),
  fCleanMuonsName(Names::gkCleanMuonsName),  
  fMuonIDType("Tight"),
  fMuonIsoType("TrackCalo"),  
  fMuons(0),
  fTrackIsolationCut(3.0),
  fCaloIsolationCut(3.0),
  fCombIsolationCut(-1.0),
  fTMOneStationLooseCut(false),
  fTMOneStationTightCut	(false),  
  fTM2DCompatibilityLooseCut(false),
  fTM2DCompatibilityTightCut(false),
  fMuonPtMin(10),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void MuonIDMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void MuonIDMod::Process()
{
  // Process entries of the tree. 

  fNEventsProcessed++;
 
  if (fNEventsProcessed % 1000000 == 0 || fPrintDebug) {
    time_t systime;
    systime = time(NULL);

    cerr << endl << "MuonIDMod : Process Event " << fNEventsProcessed << "  Time: " << ctime(&systime) << endl;  
  }  

  MuonTools myMuonTools;

  //Get Muons
  LoadBranch(fMuonName);
  ObjArray<Muon> *CleanMuons = new ObjArray<Muon>; 
  for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    Muon *mu = fMuons->At(i);
  
    Double_t MuonClass = -1;    
    if (mu->GlobalTrk())      
      MuonClass = 0;
    else if (mu->StandaloneTrk())      
      MuonClass = 1;
    else if (mu->TrackerTrk())
      MuonClass = 2;

    bool allCuts = false;

    // We always want global muons
    if(MuonClass == 0) allCuts = true;

    // Isolation requirements
    if(fCombIsolationCut < 0.0){
      if(mu->IsoR03SumPt() >= fTrackIsolationCut) allCuts = false;
      if(mu->IsoR03EmEt() + 
         mu->IsoR03HadEt() >= fCaloIsolationCut) allCuts = false;
    }
    else {
      if(1.0 * mu->IsoR03SumPt() + 
         1.0 * mu->IsoR03EmEt() + 
         1.0 * mu->IsoR03HadEt() >= fCombIsolationCut) allCuts = false;
      
    }

    // Muon chambers and calo compatibility requirements
    if(fTMOneStationLooseCut == true &&
       myMuonTools.isGood(mu, MuonTools::TMOneStationLoose) == false)
      allCuts = false;

    if(fTMOneStationTightCut == true &&
       myMuonTools.isGood(mu, MuonTools::TMOneStationTight) == false)
      allCuts = false;

    if(fTM2DCompatibilityLooseCut == true &&
       myMuonTools.isGood(mu, MuonTools::TM2DCompatibilityLoose) == false)
      allCuts = false;

    if(fTM2DCompatibilityTightCut == true &&
       myMuonTools.isGood(mu, MuonTools::TM2DCompatibilityTight) == false)
      allCuts = false;

    // Min Pt requirement
    if(mu->Pt() <= fMuonPtMin) allCuts = false;
        
    if(allCuts) {     
      CleanMuons->Add(mu);
    }
  }

  //Final Summary Debug Output   
  if ( fPrintDebug ) {
    cerr << "Event Dump: " << fNEventsProcessed << endl;  
    cerr << "Muons" << endl;
    for (UInt_t i = 0; i < CleanMuons->GetEntries(); i++) {
      cerr << i << " " << CleanMuons->At(i)->Pt() << " " << CleanMuons->At(i)->Eta() 
           << " " << CleanMuons->At(i)->Phi() << endl;    
    }  
  }   
  
  //Save Objects for Other Modules to use
  AddObjThisEvt(CleanMuons, fCleanMuonsName.Data());  
}


//--------------------------------------------------------------------------------------------------
void MuonIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fMuonName,              fMuons);
}

//--------------------------------------------------------------------------------------------------
void MuonIDMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis. For this
  // module, we dont do anything here.

}

//--------------------------------------------------------------------------------------------------
void MuonIDMod::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.
}
