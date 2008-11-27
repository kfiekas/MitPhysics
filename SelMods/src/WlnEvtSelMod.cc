// $Id: WlnEvtSelMod.cc,v 1.1 2008/11/11 21:22:54 ceballos Exp $

#include "MitPhysics/SelMods/interface/WlnEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/MuonTools.h"

using namespace mithep;
ClassImp(mithep::WlnEvtSelMod)

//--------------------------------------------------------------------------------------------------
WlnEvtSelMod::WlnEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(false),
  fMetName(Names::gkCaloMetBrn),
  fTrackName(Names::gkTrackBrn),
  fCaloTowerName(Names::gkCaloTowerBrn),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fMuonName(Names::gkMuonBrn),
  fElectronName(Names::gkElectronBrn),
  fMCLeptonsName(ModNames::gkMCLeptonsName),
  fMet(0),
  fNEventsProcessed(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void WlnEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void WlnEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  fNEventsProcessed++;

  if (fNEventsProcessed % 1000000 == 0 || fPrintDebug) {
    time_t systime;
    systime = time(NULL);
    cerr << endl << "WlnEvtSelMod : Process Event " << fNEventsProcessed << "  Time: " << ctime(&systime) << endl;
  }

  //Get Generator Level information for matching
  ObjArray<MCParticle> *GenLeptons = 
    dynamic_cast<ObjArray<MCParticle>*> (FindObjThisEvt(fMCLeptonsName.Data()));

  const TriggerObjectCol *objs = GetHLTObjects(fObjsName);

  if (!objs){
    //printf("TriggerObjectCol not found\n");
    //return;
  }

  vector<ChargedParticle*> leptons;
  vector<string> leptonType; 
  double zAverage = 0.0;
  double isoAux;
  MuonTools myMuonTools;

  // Make lepton vector from muons and electrons
  LoadBranch(fMuonName);
  ObjArray<Muon> *DirtyMuons = new ObjArray<Muon>;
  for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    Muon *mu = fMuons->At(i);
    if(!mu->GlobalTrk()) continue;
    bool isGenLepton = false;
    for (UInt_t j=0; j<GenLeptons->GetEntries(); ++j) {
      MCParticle *gen = GenLeptons->At(j);
      if(MathUtils::DeltaR(gen->Eta(), gen->Phi(), mu->Eta(), mu->Phi()) < 0.1)
        isGenLepton = true;
    }
    if(mu->Pt() > 5){
      if(isGenLepton == true){
	hDWlnSel[50]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kPromptTight))
          hDWlnSel[51]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kTMOneStationLoose))
          hDWlnSel[52]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kTMOneStationTight))
          hDWlnSel[53]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kTMLastStationLoose))
          hDWlnSel[54]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kTMLastStationTight))
          hDWlnSel[55]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kTM2DCompatibilityLoose))
          hDWlnSel[56]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kTM2DCompatibilityTight))
          hDWlnSel[57]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kTMOneStationLoose) &&
	   myMuonTools.IsGood(mu, MuonTools::kTM2DCompatibilityLoose))
          hDWlnSel[58]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kTMOneStationTight) &&
	   myMuonTools.IsGood(mu, MuonTools::kTM2DCompatibilityTight))
          hDWlnSel[59]->Fill(TMath::Min(mu->Pt(),99.999));
      } else {
	hDWlnSel[60]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kPromptTight))
          hDWlnSel[61]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kTMOneStationLoose))
          hDWlnSel[62]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kTMOneStationTight))
          hDWlnSel[63]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kTMLastStationLoose))
          hDWlnSel[64]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kTMLastStationTight))
          hDWlnSel[65]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kTM2DCompatibilityLoose))
          hDWlnSel[66]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kTM2DCompatibilityTight))
          hDWlnSel[67]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kTMOneStationLoose) &&
	   myMuonTools.IsGood(mu, MuonTools::kTM2DCompatibilityLoose))
          hDWlnSel[68]->Fill(TMath::Min(mu->Pt(),99.999));
	if(myMuonTools.IsGood(mu, MuonTools::kTMOneStationTight) &&
	   myMuonTools.IsGood(mu, MuonTools::kTM2DCompatibilityTight))
          hDWlnSel[69]->Fill(TMath::Min(mu->Pt(),99.999));
      }
      if(mu->IsoR03EmEt() + mu->IsoR03HadEt() < 3.0 &&
         mu->IsoR03SumPt() < 3.0){
        if(isGenLepton == true){
	  hDWlnSel[70]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kPromptTight))
            hDWlnSel[71]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kTMOneStationLoose))
            hDWlnSel[72]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kTMOneStationTight))
            hDWlnSel[73]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kTMLastStationLoose))
            hDWlnSel[74]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kTMLastStationTight))
            hDWlnSel[75]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kTM2DCompatibilityLoose))
            hDWlnSel[76]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kTM2DCompatibilityTight))
            hDWlnSel[77]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kTMOneStationLoose) &&
	     myMuonTools.IsGood(mu, MuonTools::kTM2DCompatibilityLoose))
            hDWlnSel[78]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kTMOneStationTight) &&
	     myMuonTools.IsGood(mu, MuonTools::kTM2DCompatibilityTight))
            hDWlnSel[79]->Fill(TMath::Min(mu->Pt(),99.999));
        } else {
	  hDWlnSel[80]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kPromptTight))
            hDWlnSel[81]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kTMOneStationLoose))
            hDWlnSel[82]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kTMOneStationTight))
            hDWlnSel[83]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kTMLastStationLoose))
            hDWlnSel[84]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kTMLastStationTight))
            hDWlnSel[85]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kTM2DCompatibilityLoose))
            hDWlnSel[86]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kTM2DCompatibilityTight))
            hDWlnSel[87]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kTMOneStationLoose) &&
	     myMuonTools.IsGood(mu, MuonTools::kTM2DCompatibilityLoose))
            hDWlnSel[88]->Fill(TMath::Min(mu->Pt(),99.999));
	  if(myMuonTools.IsGood(mu, MuonTools::kTMOneStationTight) &&
	     myMuonTools.IsGood(mu, MuonTools::kTM2DCompatibilityTight))
            hDWlnSel[89]->Fill(TMath::Min(mu->Pt(),99.999));
	}
      }
    }
    if(mu->Pt() <= 10) continue;
    if(mu->IsoR03EmEt() + mu->IsoR03HadEt() < 3.0 &&
       mu->IsoR03SumPt() < 3.0){
      DirtyMuons->Add(mu);
      zAverage = zAverage + mu->BestTrk()->Z0();
      leptons.push_back(mu);
      leptonType.push_back("mu");
    }

    if(mu->IsoR03EmEt() + mu->IsoR03HadEt() < 3.0){
      isoAux = TMath::Min(mu->IsoR03SumPt(),19.999);
      hDWlnSel[0]->Fill(isoAux);
      hDWlnSel[1]->Fill(TMath::Min(isoAux/mu->Pt(),1.999));
    }
    if(mu->IsoR03SumPt() < 3.0){
      isoAux = TMath::Min(mu->IsoR03EmEt() + mu->IsoR03HadEt(),19.999);
      hDWlnSel[2]->Fill(isoAux);
      hDWlnSel[3]->Fill(TMath::Min(isoAux/mu->Pt(),1.999));
    }
    isoAux = TMath::Min(mu->IsoR03SumPt() + mu->IsoR03EmEt() + mu->IsoR03HadEt(),19.999);
    hDWlnSel[4]->Fill(isoAux);
    hDWlnSel[5]->Fill(TMath::Min(isoAux/mu->Pt(),1.999));
    isoAux = TMath::Min(2.0*mu->IsoR03SumPt() + mu->IsoR03EmEt() + mu->IsoR03HadEt(),19.999);
    hDWlnSel[6]->Fill(isoAux);
    hDWlnSel[7]->Fill(TMath::Min(isoAux/mu->Pt(),1.999));
    isoAux = TMath::Min(2.0*mu->IsoR03SumPt() + 1.5*mu->IsoR03EmEt() + mu->IsoR03HadEt(),19.999);
    hDWlnSel[8]->Fill(isoAux);
    hDWlnSel[9]->Fill(TMath::Min(isoAux/mu->Pt(),1.999));
    isoAux = TMath::Min(mu->IsoR03SumPt() + 1.5*mu->IsoR03EmEt() + mu->IsoR03HadEt(),19.999);
    hDWlnSel[10]->Fill(isoAux);
    hDWlnSel[11]->Fill(TMath::Min(isoAux/mu->Pt(),1.999));
    isoAux = TMath::Min(mu->IsoR03SumPt() + mu->IsoR03EmEt() + mu->IsoR03HadEt(),19.999);
    if(mu->Pt() < 15){
      hDWlnSel[12]->Fill(isoAux);
      hDWlnSel[13]->Fill(TMath::Min(isoAux/mu->Pt(),1.999));
    } else if(mu->Pt() < 20){
      hDWlnSel[14]->Fill(isoAux);
      hDWlnSel[15]->Fill(TMath::Min(isoAux/mu->Pt(),1.999));
    } else if(mu->Pt() < 25){
      hDWlnSel[16]->Fill(isoAux);
      hDWlnSel[17]->Fill(TMath::Min(isoAux/mu->Pt(),1.999));
    } else if(mu->Pt() < 30){
      hDWlnSel[18]->Fill(isoAux);
      hDWlnSel[19]->Fill(TMath::Min(isoAux/mu->Pt(),1.999));
    } else if(mu->Pt() < 40){
      hDWlnSel[20]->Fill(isoAux);
      hDWlnSel[21]->Fill(TMath::Min(isoAux/mu->Pt(),1.999));
    } else if(mu->Pt() < 60){
      hDWlnSel[22]->Fill(isoAux);
      hDWlnSel[23]->Fill(TMath::Min(isoAux/mu->Pt(),1.999));
    } else {
      hDWlnSel[24]->Fill(isoAux);
      hDWlnSel[25]->Fill(TMath::Min(isoAux/mu->Pt(),1.999));
    }
    hDWlnSel[26]->Fill(TMath::Min(mu->Pt(),99.99));
    hDWlnSel2D[0]->Fill(isoAux,mu->Pt());
    if(isoAux > 0)
      hDWlnSel2D[1]->Fill(isoAux,mu->Pt());

    if(isoAux < 7.0){
      int theHisto = 27;
      if(isGenLepton == false) theHisto = 28;
      hDWlnSel[theHisto]->Fill(0.0);
      if((isoAux < (mu->Pt()-10.0)*7.0/5.0) ||
         (isoAux < 7.0 && mu->Pt() > 15)) hDWlnSel[theHisto]->Fill(1.0);
      if((isoAux < (mu->Pt()-10.0)*6.5/5.0) ||
         (isoAux < 6.5 && mu->Pt() > 15)) hDWlnSel[theHisto]->Fill(2.0);
      if((isoAux < (mu->Pt()-10.0)*6.0/5.0) ||
         (isoAux < 6.0 && mu->Pt() > 15)) hDWlnSel[theHisto]->Fill(3.0);
      if((isoAux < (mu->Pt()-10.0)*5.5/5.0) ||
         (isoAux < 5.5 && mu->Pt() > 15)) hDWlnSel[theHisto]->Fill(4.0);
      if((isoAux < (mu->Pt()-10.0)*5.0/5.0) ||
         (isoAux < 5.0 && mu->Pt() > 15)) hDWlnSel[theHisto]->Fill(5.0);
      if((isoAux < (mu->Pt()-10.0)*4.5/5.0) ||
         (isoAux < 4.5 && mu->Pt() > 15)) hDWlnSel[theHisto]->Fill(6.0);
      if((isoAux < (mu->Pt()-10.0)*4.0/5.0) ||
         (isoAux < 4.0 && mu->Pt() > 15)) hDWlnSel[theHisto]->Fill(7.0);
      if((isoAux < (mu->Pt()-10.0)*3.5/5.0) ||
         (isoAux < 3.0 && mu->Pt() > 15)) hDWlnSel[theHisto]->Fill(8.0);
      if((isoAux < (mu->Pt()-10.0)*3.0/5.0) ||
         (isoAux < 3.0 && mu->Pt() > 15)) hDWlnSel[theHisto]->Fill(9.0);
      if((isoAux < (mu->Pt()-10.0)*2.5/5.0) ||
         (isoAux < 2.5 && mu->Pt() > 15)) hDWlnSel[theHisto]->Fill(10.0);

      if((isoAux < (mu->Pt()-10.0)*7.0/10.0) ||
         (isoAux < 7.0 && mu->Pt() > 20)) hDWlnSel[theHisto]->Fill(11.0);
      if((isoAux < (mu->Pt()-10.0)*6.5/10.0) ||
         (isoAux < 6.5 && mu->Pt() > 20)) hDWlnSel[theHisto]->Fill(12.0);
      if((isoAux < (mu->Pt()-10.0)*6.0/10.0) ||
         (isoAux < 6.0 && mu->Pt() > 20)) hDWlnSel[theHisto]->Fill(13.0);
      if((isoAux < (mu->Pt()-10.0)*5.5/10.0) ||
         (isoAux < 5.5 && mu->Pt() > 20)) hDWlnSel[theHisto]->Fill(14.0);
      if((isoAux < (mu->Pt()-10.0)*5.0/10.0) ||
         (isoAux < 5.0 && mu->Pt() > 20)) hDWlnSel[theHisto]->Fill(15.0);
      if((isoAux < (mu->Pt()-10.0)*4.5/10.0) ||
         (isoAux < 4.5 && mu->Pt() > 20)) hDWlnSel[theHisto]->Fill(16.0);
      if((isoAux < (mu->Pt()-10.0)*4.0/10.0) ||
         (isoAux < 4.0 && mu->Pt() > 20)) hDWlnSel[theHisto]->Fill(17.0);
      if((isoAux < (mu->Pt()-10.0)*3.5/10.0) ||
         (isoAux < 3.0 && mu->Pt() > 20)) hDWlnSel[theHisto]->Fill(18.0);
      if((isoAux < (mu->Pt()-10.0)*3.0/10.0) ||
         (isoAux < 3.0 && mu->Pt() > 20)) hDWlnSel[theHisto]->Fill(19.0);
      if((isoAux < (mu->Pt()-10.0)*2.5/10.0) ||
         (isoAux < 2.5 && mu->Pt() > 20)) hDWlnSel[theHisto]->Fill(20.0);

      if((isoAux < (mu->Pt()-10.0)*7.0/15.0) ||
         (isoAux < 7.0 && mu->Pt() > 25)) hDWlnSel[theHisto]->Fill(21.0);
      if((isoAux < (mu->Pt()-10.0)*6.5/15.0) ||
         (isoAux < 6.5 && mu->Pt() > 25)) hDWlnSel[theHisto]->Fill(22.0);
      if((isoAux < (mu->Pt()-10.0)*6.0/15.0) ||
         (isoAux < 6.0 && mu->Pt() > 25)) hDWlnSel[theHisto]->Fill(23.0);
      if((isoAux < (mu->Pt()-10.0)*5.5/15.0) ||
         (isoAux < 5.5 && mu->Pt() > 25)) hDWlnSel[theHisto]->Fill(24.0);
      if((isoAux < (mu->Pt()-10.0)*5.0/15.0) ||
         (isoAux < 5.0 && mu->Pt() > 25)) hDWlnSel[theHisto]->Fill(25.0);
      if((isoAux < (mu->Pt()-10.0)*4.5/15.0) ||
         (isoAux < 4.5 && mu->Pt() > 25)) hDWlnSel[theHisto]->Fill(26.0);
      if((isoAux < (mu->Pt()-10.0)*4.0/15.0) ||
         (isoAux < 4.0 && mu->Pt() > 25)) hDWlnSel[theHisto]->Fill(27.0);
      if((isoAux < (mu->Pt()-10.0)*3.5/15.0) ||
         (isoAux < 3.0 && mu->Pt() > 25)) hDWlnSel[theHisto]->Fill(28.0);
      if((isoAux < (mu->Pt()-10.0)*3.0/15.0) ||
         (isoAux < 3.0 && mu->Pt() > 25)) hDWlnSel[theHisto]->Fill(29.0);
      if((isoAux < (mu->Pt()-10.0)*2.5/15.0) ||
         (isoAux < 2.5 && mu->Pt() > 25)) hDWlnSel[theHisto]->Fill(30.0);

      if((isoAux < (mu->Pt()-10.0)*7.0/20.0) ||
         (isoAux < 7.0 && mu->Pt() > 30)) hDWlnSel[theHisto]->Fill(31.0);
      if((isoAux < (mu->Pt()-10.0)*6.5/20.0) ||
         (isoAux < 6.5 && mu->Pt() > 30)) hDWlnSel[theHisto]->Fill(32.0);
      if((isoAux < (mu->Pt()-10.0)*6.0/20.0) ||
         (isoAux < 6.0 && mu->Pt() > 30)) hDWlnSel[theHisto]->Fill(33.0);
      if((isoAux < (mu->Pt()-10.0)*5.5/20.0) ||
         (isoAux < 5.5 && mu->Pt() > 30)) hDWlnSel[theHisto]->Fill(34.0);
      if((isoAux < (mu->Pt()-10.0)*5.0/20.0) ||
         (isoAux < 5.0 && mu->Pt() > 30)) hDWlnSel[theHisto]->Fill(35.0);
      if((isoAux < (mu->Pt()-10.0)*4.5/20.0) ||
         (isoAux < 4.5 && mu->Pt() > 30)) hDWlnSel[theHisto]->Fill(36.0);
      if((isoAux < (mu->Pt()-10.0)*4.0/20.0) ||
         (isoAux < 4.0 && mu->Pt() > 30)) hDWlnSel[theHisto]->Fill(37.0);
      if((isoAux < (mu->Pt()-10.0)*3.5/20.0) ||
         (isoAux < 3.0 && mu->Pt() > 30)) hDWlnSel[theHisto]->Fill(38.0);
      if((isoAux < (mu->Pt()-10.0)*3.0/20.0) ||
         (isoAux < 3.0 && mu->Pt() > 30)) hDWlnSel[theHisto]->Fill(39.0);
      if((isoAux < (mu->Pt()-10.0)*2.5/20.0) ||
         (isoAux < 2.5 && mu->Pt() > 30)) hDWlnSel[theHisto]->Fill(40.0);
    } // iso vs. Pt study
  }

  LoadBranch(fElectronName);
  ObjArray<Electron> *DirtyElectrons = new ObjArray<Electron>;
  for (UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    Electron *el = fElectrons->At(i);
    if(el->Pt() <= 10) continue;
    bool isGenLepton = false;
    for (UInt_t j=0; j<GenLeptons->GetEntries(); ++j) {
      MCParticle *gen = GenLeptons->At(j);
      if(MathUtils::DeltaR(gen->Eta(), gen->Phi(), el->Eta(), el->Phi()) < 0.1)
        isGenLepton = true;
    }
    if(el->TrackIsolation() < 5.0 &&
       el->TrackIsolation() < 5.0){
      if(isGenLepton == true){
	hDWlnSel[150]->Fill(TMath::Min(el->Pt(),99.999));
	if(el->PassTightID())
          hDWlnSel[151]->Fill(TMath::Min(el->Pt(),99.999));
	if(el->PassLooseID())
          hDWlnSel[152]->Fill(TMath::Min(el->Pt(),99.999));
	hDWlnSel[153]->Fill(TMath::Max(TMath::Min(el->IDLikelihood(),0.99999),0.));
      } else {
	hDWlnSel[154]->Fill(TMath::Min(el->Pt(),99.999));
	if(el->PassTightID())
          hDWlnSel[155]->Fill(TMath::Min(el->Pt(),99.999));
	if(el->PassLooseID())
          hDWlnSel[156]->Fill(TMath::Min(el->Pt(),99.999));
	hDWlnSel[157]->Fill(TMath::Max(TMath::Min(el->IDLikelihood(),0.99999),0.));
      }
    }
    if(!el->PassTightID()) continue;
    if(el->TrackIsolation() < 5.0 &&
       el->TrackIsolation() < 5.0){
      DirtyElectrons->Add(el);
      zAverage = zAverage + el->BestTrk()->Z0();
      leptons.push_back(el);
      leptonType.push_back("e");
    }

    if(el->CaloIsolation() < 5.0){
      isoAux = TMath::Min(el->TrackIsolation(),19.999);
      hDWlnSel[100]->Fill(isoAux);
      hDWlnSel[101]->Fill(TMath::Min(isoAux/el->Pt(),1.999));
    }
    if(el->TrackIsolation() < 5.0){
      isoAux = TMath::Min(el->CaloIsolation(),19.999);
      hDWlnSel[102]->Fill(isoAux);
      hDWlnSel[103]->Fill(TMath::Min(isoAux/el->Pt(),1.999));
      isoAux = TMath::Min(el->CaloTowerIsolation(),19.999);
      hDWlnSel[104]->Fill(isoAux);
      hDWlnSel[105]->Fill(TMath::Min(isoAux/el->Pt(),1.999));
    }
    isoAux = TMath::Min(el->TrackIsolation() + el->CaloIsolation(),19.999);
    hDWlnSel[106]->Fill(isoAux);
    hDWlnSel[107]->Fill(TMath::Min(isoAux/el->Pt(),1.999));
    isoAux = TMath::Min(2.0*el->TrackIsolation() + el->CaloIsolation(),19.999);
    hDWlnSel[108]->Fill(isoAux);
    hDWlnSel[109]->Fill(TMath::Min(isoAux/el->Pt(),1.999));
    isoAux = TMath::Min(el->TrackIsolation() + el->CaloTowerIsolation(),19.999);
    hDWlnSel[110]->Fill(isoAux);
    hDWlnSel[111]->Fill(TMath::Min(isoAux/el->Pt(),1.999));
    isoAux = TMath::Min(el->TrackIsolation() + el->CaloIsolation(),19.999);
    if(el->Pt() < 15){
      hDWlnSel[112]->Fill(isoAux);
      hDWlnSel[113]->Fill(TMath::Min(isoAux/el->Pt(),1.999));
    } else if(el->Pt() < 20){
      hDWlnSel[114]->Fill(isoAux);
      hDWlnSel[115]->Fill(TMath::Min(isoAux/el->Pt(),1.999));
    } else if(el->Pt() < 25){
      hDWlnSel[116]->Fill(isoAux);
      hDWlnSel[117]->Fill(TMath::Min(isoAux/el->Pt(),1.999));
    } else if(el->Pt() < 30){
      hDWlnSel[118]->Fill(isoAux);
      hDWlnSel[119]->Fill(TMath::Min(isoAux/el->Pt(),1.999));
    } else if(el->Pt() < 40){
      hDWlnSel[120]->Fill(isoAux);
      hDWlnSel[121]->Fill(TMath::Min(isoAux/el->Pt(),1.999));
    } else if(el->Pt() < 60){
      hDWlnSel[122]->Fill(isoAux);
      hDWlnSel[123]->Fill(TMath::Min(isoAux/el->Pt(),1.999));
    } else {
      hDWlnSel[124]->Fill(isoAux);
      hDWlnSel[125]->Fill(TMath::Min(isoAux/el->Pt(),1.999));
    }
    hDWlnSel[126]->Fill(TMath::Min(el->Pt(),99.99));
    hDWlnSel2D[2]->Fill(isoAux,el->Pt());
    if(isoAux > 0)
      hDWlnSel2D[3]->Fill(isoAux,el->Pt());

    if(isoAux < 7.0){
      int theHisto = 127;
      if(isGenLepton == false) theHisto = 128;
      hDWlnSel[theHisto]->Fill(0.0);
      if((isoAux < (el->Pt()-10.0)*7.0/5.0) ||
         (isoAux < 7.0 && el->Pt() > 15)) hDWlnSel[theHisto]->Fill(1.0);
      if((isoAux < (el->Pt()-10.0)*6.5/5.0) ||
         (isoAux < 6.5 && el->Pt() > 15)) hDWlnSel[theHisto]->Fill(2.0);
      if((isoAux < (el->Pt()-10.0)*6.0/5.0) ||
         (isoAux < 6.0 && el->Pt() > 15)) hDWlnSel[theHisto]->Fill(3.0);
      if((isoAux < (el->Pt()-10.0)*5.5/5.0) ||
         (isoAux < 5.5 && el->Pt() > 15)) hDWlnSel[theHisto]->Fill(4.0);
      if((isoAux < (el->Pt()-10.0)*5.0/5.0) ||
         (isoAux < 5.0 && el->Pt() > 15)) hDWlnSel[theHisto]->Fill(5.0);
      if((isoAux < (el->Pt()-10.0)*4.5/5.0) ||
         (isoAux < 4.5 && el->Pt() > 15)) hDWlnSel[theHisto]->Fill(6.0);
      if((isoAux < (el->Pt()-10.0)*4.0/5.0) ||
         (isoAux < 4.0 && el->Pt() > 15)) hDWlnSel[theHisto]->Fill(7.0);
      if((isoAux < (el->Pt()-10.0)*3.5/5.0) ||
         (isoAux < 3.0 && el->Pt() > 15)) hDWlnSel[theHisto]->Fill(8.0);
      if((isoAux < (el->Pt()-10.0)*3.0/5.0) ||
         (isoAux < 3.0 && el->Pt() > 15)) hDWlnSel[theHisto]->Fill(9.0);
      if((isoAux < (el->Pt()-10.0)*2.5/5.0) ||
         (isoAux < 2.5 && el->Pt() > 15)) hDWlnSel[theHisto]->Fill(10.0);

      if((isoAux < (el->Pt()-10.0)*7.0/10.0) ||
         (isoAux < 7.0 && el->Pt() > 20)) hDWlnSel[theHisto]->Fill(11.0);
      if((isoAux < (el->Pt()-10.0)*6.5/10.0) ||
         (isoAux < 6.5 && el->Pt() > 20)) hDWlnSel[theHisto]->Fill(12.0);
      if((isoAux < (el->Pt()-10.0)*6.0/10.0) ||
         (isoAux < 6.0 && el->Pt() > 20)) hDWlnSel[theHisto]->Fill(13.0);
      if((isoAux < (el->Pt()-10.0)*5.5/10.0) ||
         (isoAux < 5.5 && el->Pt() > 20)) hDWlnSel[theHisto]->Fill(14.0);
      if((isoAux < (el->Pt()-10.0)*5.0/10.0) ||
         (isoAux < 5.0 && el->Pt() > 20)) hDWlnSel[theHisto]->Fill(15.0);
      if((isoAux < (el->Pt()-10.0)*4.5/10.0) ||
         (isoAux < 4.5 && el->Pt() > 20)) hDWlnSel[theHisto]->Fill(16.0);
      if((isoAux < (el->Pt()-10.0)*4.0/10.0) ||
         (isoAux < 4.0 && el->Pt() > 20)) hDWlnSel[theHisto]->Fill(17.0);
      if((isoAux < (el->Pt()-10.0)*3.5/10.0) ||
         (isoAux < 3.0 && el->Pt() > 20)) hDWlnSel[theHisto]->Fill(18.0);
      if((isoAux < (el->Pt()-10.0)*3.0/10.0) ||
         (isoAux < 3.0 && el->Pt() > 20)) hDWlnSel[theHisto]->Fill(19.0);
      if((isoAux < (el->Pt()-10.0)*2.5/10.0) ||
         (isoAux < 2.5 && el->Pt() > 20)) hDWlnSel[theHisto]->Fill(20.0);

      if((isoAux < (el->Pt()-10.0)*7.0/15.0) ||
         (isoAux < 7.0 && el->Pt() > 25)) hDWlnSel[theHisto]->Fill(21.0);
      if((isoAux < (el->Pt()-10.0)*6.5/15.0) ||
         (isoAux < 6.5 && el->Pt() > 25)) hDWlnSel[theHisto]->Fill(22.0);
      if((isoAux < (el->Pt()-10.0)*6.0/15.0) ||
         (isoAux < 6.0 && el->Pt() > 25)) hDWlnSel[theHisto]->Fill(23.0);
      if((isoAux < (el->Pt()-10.0)*5.5/15.0) ||
         (isoAux < 5.5 && el->Pt() > 25)) hDWlnSel[theHisto]->Fill(24.0);
      if((isoAux < (el->Pt()-10.0)*5.0/15.0) ||
         (isoAux < 5.0 && el->Pt() > 25)) hDWlnSel[theHisto]->Fill(25.0);
      if((isoAux < (el->Pt()-10.0)*4.5/15.0) ||
         (isoAux < 4.5 && el->Pt() > 25)) hDWlnSel[theHisto]->Fill(26.0);
      if((isoAux < (el->Pt()-10.0)*4.0/15.0) ||
         (isoAux < 4.0 && el->Pt() > 25)) hDWlnSel[theHisto]->Fill(27.0);
      if((isoAux < (el->Pt()-10.0)*3.5/15.0) ||
         (isoAux < 3.0 && el->Pt() > 25)) hDWlnSel[theHisto]->Fill(28.0);
      if((isoAux < (el->Pt()-10.0)*3.0/15.0) ||
         (isoAux < 3.0 && el->Pt() > 25)) hDWlnSel[theHisto]->Fill(29.0);
      if((isoAux < (el->Pt()-10.0)*2.5/15.0) ||
         (isoAux < 2.5 && el->Pt() > 25)) hDWlnSel[theHisto]->Fill(30.0);

      if((isoAux < (el->Pt()-10.0)*7.0/20.0) ||
         (isoAux < 7.0 && el->Pt() > 30)) hDWlnSel[theHisto]->Fill(31.0);
      if((isoAux < (el->Pt()-10.0)*6.5/20.0) ||
         (isoAux < 6.5 && el->Pt() > 30)) hDWlnSel[theHisto]->Fill(32.0);
      if((isoAux < (el->Pt()-10.0)*6.0/20.0) ||
         (isoAux < 6.0 && el->Pt() > 30)) hDWlnSel[theHisto]->Fill(33.0);
      if((isoAux < (el->Pt()-10.0)*5.5/20.0) ||
         (isoAux < 5.5 && el->Pt() > 30)) hDWlnSel[theHisto]->Fill(34.0);
      if((isoAux < (el->Pt()-10.0)*5.0/20.0) ||
         (isoAux < 5.0 && el->Pt() > 30)) hDWlnSel[theHisto]->Fill(35.0);
      if((isoAux < (el->Pt()-10.0)*4.5/20.0) ||
         (isoAux < 4.5 && el->Pt() > 30)) hDWlnSel[theHisto]->Fill(36.0);
      if((isoAux < (el->Pt()-10.0)*4.0/20.0) ||
         (isoAux < 4.0 && el->Pt() > 30)) hDWlnSel[theHisto]->Fill(37.0);
      if((isoAux < (el->Pt()-10.0)*3.5/20.0) ||
         (isoAux < 3.0 && el->Pt() > 30)) hDWlnSel[theHisto]->Fill(38.0);
      if((isoAux < (el->Pt()-10.0)*3.0/20.0) ||
         (isoAux < 3.0 && el->Pt() > 30)) hDWlnSel[theHisto]->Fill(39.0);
      if((isoAux < (el->Pt()-10.0)*2.5/20.0) ||
         (isoAux < 2.5 && el->Pt() > 30)) hDWlnSel[theHisto]->Fill(40.0);
    } // iso vs. Pt study
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

  if(leptons.size() >= 1){
    LoadBranch(fMetName);
    Met *caloMet = fMet->At(0);

    double deltaPhiMetLepton = MathUtils::DeltaPhi(caloMet->Phi(), leptons[0]->Phi());

    double mTW = TMath::Sqrt(2.0*leptons[0]->Pt()*caloMet->Pt()*
    		            (1.0 - cos(deltaPhiMetLepton)));

    hDWlnSel[200]->Fill(TMath::Min(caloMet->Pt(),99.999));
    hDWlnSel[201]->Fill(TMath::Min((double)DirtyMuons->GetEntries(),9.499));
    hDWlnSel[202]->Fill(TMath::Min((double)DirtyElectrons->GetEntries(),9.499));
    hDWlnSel[203]->Fill(deltaPhiMetLepton * 180./TMath::Pi());
    hDWlnSel[204]->Fill(TMath::Min(mTW, 199.999));
  }
  delete DirtyMuons;
  delete DirtyElectrons;
  leptons.clear();
  delete objs;

}
//--------------------------------------------------------------------------------------------------
void WlnEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fMetName,       fMet);
  ReqBranch(fTrackName,     fTracks);
  ReqBranch(fCaloTowerName, fCaloTowers);
  ReqBranch(fMuonName,      fMuons);
  ReqBranch(fElectronName,  fElectrons);

  char sb[200];
  // Muons
  sprintf(sb,"hDWlnSel_%d", 0);  hDWlnSel[ 0]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d", 1);  hDWlnSel[ 1]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d", 2);  hDWlnSel[ 2]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d", 3);  hDWlnSel[ 3]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d", 4);  hDWlnSel[ 4]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d", 5);  hDWlnSel[ 5]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d", 6);  hDWlnSel[ 6]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d", 7);  hDWlnSel[ 7]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d", 8);  hDWlnSel[ 8]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d", 9);  hDWlnSel[ 9]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",10);  hDWlnSel[10]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d",11);  hDWlnSel[11]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",12);  hDWlnSel[12]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d",13);  hDWlnSel[13]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",14);  hDWlnSel[14]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d",15);  hDWlnSel[15]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",16);  hDWlnSel[16]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d",17);  hDWlnSel[17]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",18);  hDWlnSel[18]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d",19);  hDWlnSel[19]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",20);  hDWlnSel[20]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d",21);  hDWlnSel[21]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",22);  hDWlnSel[22]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d",23);  hDWlnSel[23]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",24);  hDWlnSel[24]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d",25);  hDWlnSel[25]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",26);  hDWlnSel[26]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",27);  hDWlnSel[27]  = new TH1D(sb,sb,41,-0.5,40.5);
  sprintf(sb,"hDWlnSel_%d",28);  hDWlnSel[28]  = new TH1D(sb,sb,41,-0.5,40.5);

  for(int i=0; i<29; i++){
    AddOutput(hDWlnSel[i]);
  }

  sprintf(sb,"hDWlnSel_%d",50);  hDWlnSel[50]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",51);  hDWlnSel[51]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",52);  hDWlnSel[52]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",53);  hDWlnSel[53]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",54);  hDWlnSel[54]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",55);  hDWlnSel[55]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",56);  hDWlnSel[56]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",57);  hDWlnSel[57]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",58);  hDWlnSel[58]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",59);  hDWlnSel[59]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",60);  hDWlnSel[60]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",61);  hDWlnSel[61]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",62);  hDWlnSel[62]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",63);  hDWlnSel[63]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",64);  hDWlnSel[64]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",65);  hDWlnSel[65]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",66);  hDWlnSel[66]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",67);  hDWlnSel[67]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",68);  hDWlnSel[68]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",69);  hDWlnSel[69]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",70);  hDWlnSel[70]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",71);  hDWlnSel[71]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",72);  hDWlnSel[72]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",73);  hDWlnSel[73]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",74);  hDWlnSel[74]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",75);  hDWlnSel[75]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",76);  hDWlnSel[76]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",77);  hDWlnSel[77]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",78);  hDWlnSel[78]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",79);  hDWlnSel[79]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",80);  hDWlnSel[80]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",81);  hDWlnSel[81]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",82);  hDWlnSel[82]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",83);  hDWlnSel[83]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",84);  hDWlnSel[84]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",85);  hDWlnSel[85]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",86);  hDWlnSel[86]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",87);  hDWlnSel[87]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",88);  hDWlnSel[88]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",89);  hDWlnSel[89]  = new TH1D(sb,sb,200,0.,100.0);

  for(int i=50; i<90; i++){
    AddOutput(hDWlnSel[i]);
  }

  // Electrons
  sprintf(sb,"hDWlnSel_%d", 0+100);  hDWlnSel[ 0+100]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d", 1+100);  hDWlnSel[ 1+100]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d", 2+100);  hDWlnSel[ 2+100]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d", 3+100);  hDWlnSel[ 3+100]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d", 4+100);  hDWlnSel[ 4+100]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d", 5+100);  hDWlnSel[ 5+100]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d", 6+100);  hDWlnSel[ 6+100]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d", 7+100);  hDWlnSel[ 7+100]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d", 8+100);  hDWlnSel[ 8+100]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d", 9+100);  hDWlnSel[ 9+100]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",10+100);  hDWlnSel[10+100]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d",11+100);  hDWlnSel[11+100]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",12+100);  hDWlnSel[12+100]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d",13+100);  hDWlnSel[13+100]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",14+100);  hDWlnSel[14+100]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d",15+100);  hDWlnSel[15+100]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",16+100);  hDWlnSel[16+100]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d",17+100);  hDWlnSel[17+100]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",18+100);  hDWlnSel[18+100]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d",19+100);  hDWlnSel[19+100]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",20+100);  hDWlnSel[20+100]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d",21+100);  hDWlnSel[21+100]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",22+100);  hDWlnSel[22+100]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d",23+100);  hDWlnSel[23+100]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",24+100);  hDWlnSel[24+100]  = new TH1D(sb,sb,200,0.,20.);
  sprintf(sb,"hDWlnSel_%d",25+100);  hDWlnSel[25+100]  = new TH1D(sb,sb,200,0.,2.0);
  sprintf(sb,"hDWlnSel_%d",26+100);  hDWlnSel[26+100]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",27+100);  hDWlnSel[27+100]  = new TH1D(sb,sb,41,-0.5,40.5);
  sprintf(sb,"hDWlnSel_%d",28+100);  hDWlnSel[28+100]  = new TH1D(sb,sb,41,-0.5,40.5);

  for(int i=0; i<29; i++){
    AddOutput(hDWlnSel[i+100]);
  }

  sprintf(sb,"hDWlnSel_%d",150);  hDWlnSel[150]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",151);  hDWlnSel[151]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",152);  hDWlnSel[152]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",153);  hDWlnSel[153]  = new TH1D(sb,sb,200,0.,1.0);
  sprintf(sb,"hDWlnSel_%d",154);  hDWlnSel[154]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",155);  hDWlnSel[155]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",156);  hDWlnSel[156]  = new TH1D(sb,sb,200,0.,100.0);
  sprintf(sb,"hDWlnSel_%d",157);  hDWlnSel[157]  = new TH1D(sb,sb,200,0.,1.0);

  for(int i=150; i<158; i++){
    AddOutput(hDWlnSel[i]);
  }
  // Event
  sprintf(sb,"hDWlnSel_%d",200);  hDWlnSel[200]  = new TH1D(sb,sb,100,0.,100.);
  sprintf(sb,"hDWlnSel_%d",201);  hDWlnSel[201]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDWlnSel_%d",202);  hDWlnSel[202]  = new TH1D(sb,sb,10,-0.5,9.5);
  sprintf(sb,"hDWlnSel_%d",203);  hDWlnSel[203]  = new TH1D(sb,sb,90,0.0,180.);
  sprintf(sb,"hDWlnSel_%d",204);  hDWlnSel[204]  = new TH1D(sb,sb,100,0.0,200.);

  for(int i=0; i<5; i++){
    AddOutput(hDWlnSel[i+200]);
  }

  // 2D
  sprintf(sb,"hDWlnSel2D_%d",0);   hDWlnSel2D[0] = new TH2D(sb,sb,50,0.0,10.0,50,0,100);
  sprintf(sb,"hDWlnSel2D_%d",1);   hDWlnSel2D[1] = new TH2D(sb,sb,50,0.0,10.0,50,0,100); 
  sprintf(sb,"hDWlnSel2D_%d",2);   hDWlnSel2D[2] = new TH2D(sb,sb,50,0.0,10.0,50,0,100); 
  sprintf(sb,"hDWlnSel2D_%d",3);   hDWlnSel2D[3] = new TH2D(sb,sb,50,0.0,10.0,50,0,100); 

  for(int i=0; i<4; i++){
    AddOutput(hDWlnSel2D[i]);
  }
}

//--------------------------------------------------------------------------------------------------
void WlnEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void WlnEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}

//--------------------------------------------------------------------------------------------------
void WlnEvtSelMod::getTrackIsolation(double theEta, double thePhi,
                                     vector<ChargedParticle*> leptons,
                                     mithep::TrackCol *iTracks,
				     double &sumPt, int &nTracks) {
  sumPt   = 0.0; nTracks = 0;
  if(leptons.size() < 2) return;

  for(unsigned int i = 0; i < iTracks->GetEntries(); i++) {
    mithep::Track* pTrack = iTracks->At(i);
    if(MathUtils::DeltaR(pTrack->Eta(), pTrack->Phi(),
    		         leptons[0]->Eta(), leptons[0]->Phi()) < 0.05 ||
       MathUtils::DeltaR(pTrack->Eta(), pTrack->Phi(),
    		         leptons[1]->Eta(), leptons[1]->Phi()) < 0.05 ||
       MathUtils::DeltaR(pTrack->Eta(), pTrack->Phi(),
    		         theEta, thePhi) > 0.3) continue;

    if(pTrack->Pt() <= 0.0) continue;
    sumPt = sumPt + pTrack->Pt(); nTracks++;
  }
}

//--------------------------------------------------------------------------------------------------
void WlnEvtSelMod::getCaloIsolation(double theEta, double thePhi,
                                    vector<ChargedParticle*> leptons,
                                    mithep::CaloTowerCol *iTowers,
				    double lIso[4]) {

  for(int i=0; i<4; i++) lIso[i] = 0.0;
  if(leptons.size() < 2) return;

  // input (to be included somewhere else)
  double fEth    = 0.2;
  double fHth    = 0.5;
  double fHoh    = 0.5;
    
  double fNoiseEB = 0.04;
  double fNoiseEE = 0.15;
  double fNoiseHB = 0.2;
  double fNoiseHO = 0.2;

  double fDRMax  = 0.3;
  double fVetoEm = 0.07;
  double fVetoHa = 0.1;
  double fVetoHo = 0.1;

  for (unsigned int i0 = 0; i0 < iTowers->GetEntries(); i0++){
    mithep::CaloTower *pTower = iTowers->At(i0);

    if(MathUtils::DeltaR(pTower->Eta(), pTower->Phi(),
    		         leptons[0]->Eta(), leptons[0]->Phi()) < 0.1 ||
       MathUtils::DeltaR(pTower->Eta(), pTower->Phi(),
    		         leptons[1]->Eta(), leptons[1]->Phi()) < 0.1) continue;

    double pDR = MathUtils::DeltaR(pTower->Eta(), pTower->Phi(), theEta, thePhi);

    double pEmEt = pTower->EmEt();
    double pEmE  = pTower->EmEnergy();
    bool   pEcal = pEmEt > fEth && pEmE > 3 * (fabs(theEta) > 1.479 ? fNoiseEE : fNoiseEB);
    
    double pHaEt = pTower->HadEt();
    double pHaE  = pTower->HadEnergy();
    bool   pHcal = pHaEt > fHth && pHaE > 3 * fNoiseHB;
    
    double pHoEt = pTower->OuterEt();
    double pHoE  = pTower->OuterEnergy();
    bool   pHoal = pHoEt > fHoh && pHoE > 3 * fNoiseHO;
    
    if (pDR > fDRMax) continue;
    
    if ((!pEcal) && (!pHcal) && (!pHoal)) continue;

    bool pEmVeto = false; bool pHaVeto = false; bool pHoVeto = false;
    if (pDR < fVetoEm ) pEmVeto = true;
    if (pDR < fVetoHa ) pHaVeto = true;
    if (pDR < fVetoHo ) pHoVeto = true;
    
    if (pEcal && pEmVeto == false) {lIso[1] += pEmEt; lIso[0] += pEmEt;}
    if (pHcal && pHaVeto == false) {lIso[2] += pHaEt; lIso[0] += pHaEt;} 
    if (pHoal && pHoVeto == false) {lIso[3] += pHoEt; lIso[0] += pHoEt;} 
  }
}
