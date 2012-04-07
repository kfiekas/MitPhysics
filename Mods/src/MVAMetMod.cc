// $Id: JetIDMod.cc,v 1.27 2012/04/07 09:36:32 pharris Exp $

#include "MitPhysics/Mods/interface/MVAMetMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Init/interface/ModNames.h"


using namespace mithep;

ClassImp(mithep::MVAMetMod)

//--------------------------------------------------------------------------------------------------
MVAMetMod::MVAMetMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fMVAMetName("MVAMet"),  
  fJetsName  (ModNames::gkGoodJetsName),
  fPFCandName(Names::gkPFCandidatesBrn),
  fVertexName(ModNames::gkGoodVertexesName),
  fPFMetName ("PFMet"),
  fJets(0),
  fCands(0),
  fVertices(0),
  fPFMet(0),
  fMuons(0)
{
  // Constructor.
  fMVAMet    = new MVAMet();
  fMVAMet->Initialize(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_lowpt.weights.xml"))),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_highpt.weights.xml"))),
                      TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/Utils/python/JetIdParams_cfi.py")),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmet.root"))),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetphi.root")))
                      );

 }

//-------------------------------------------------------------------------------------------------
void MVAMetMod::Process()
{
  // Process entries of the tree. 

  fJets     = GetObjThisEvt<PFJetCol>       (fJetsName); //corrected Jets
  fCands    = GetObjThisEvt<PFCandidateCol>(fPFCandName);
  fVertices = GetObjThisEvt<VertexOArr>    (fVertexName);
  fPFMet    = GetObjThisEvt<PFMetCol>      (fPFMetName);
  //fMuons    = GetObjThisEvt<MuonCol>       (fMuons); //Random lepton selection
  if (!fJets || !fCands || !fVertices || !fPFMet || !fMuons) {
    SendError(kAbortModule, "Process", 
              "Pointer to input jet collection %s is null.",
              fJetsName.Data());
    return;
  }
  //Random lepton selection
  Float_t lPt0 = 0; Float_t lEta0 = 0; Float_t lPhi0 = 0;
  Float_t lPt1 = 0; Float_t lEta1 = 0; Float_t lPhi1 = 0;

  //const Muon *lMuon0 = 0; if(fMuons->GetEntries() > 0) lMuon0 = fMuons->At(0);
  //const Muon *lMuon1 = 0; if(fMuons->GetEntries() > 1) lMuon1 = fMuons->At(1);
  //lPt0  = lMuon0->Pt();
  //lEta0 = lMuon0->Eta();
  //lPhi0 = lMuon0->Phi();

  //lPt1  = lMuon1->Pt();
  //lEta1 = lMuon1->Eta();
  //lPhi1 = lMuon1->Phi();


  MetOArr *MVAMet = new MetOArr; 
  MVAMet->SetName(fMVAMetName);
  Met lMVAMet = fMVAMet->GetMet(  false,
                                  lPt0,lPhi0,lEta0,
                                  lPt1,lPhi1,lEta1,
                                  fPFMet->At(0),
                                  fCands,fVertices->At(0),fVertices,
                                  fJets,
                                  int(fVertices->GetEntries()));

  MVAMet->Add(&lMVAMet);
  
  // sort according to pt
  MVAMet->Sort();
  
  // add to event for other modules to use
  AddObjThisEvt(MVAMet);  
}

