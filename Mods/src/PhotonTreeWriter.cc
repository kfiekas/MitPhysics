#include "MitPhysics/Mods/interface/PhotonTreeWriter.h"
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
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "TDataMember.h"
#include "TFile.h"
#include <TNtuple.h>
#include <TRandom3.h>
#include <TSystem.h>

using namespace mithep;

ClassImp(mithep::PhotonTreeWriter)
templateClassImp(mithep::PhotonTreeWriterPhoton)
ClassImp(mithep::PhotonTreeWriterDiphotonEvent)

//--------------------------------------------------------------------------------------------------
PhotonTreeWriter::PhotonTreeWriter(const char *name, const char *title) : 
  // Base Module...
  BaseMod                 (name,title),
  // define all the Branches to load
  fPhotonBranchName       (Names::gkPhotonBrn),
  fPFPhotonName           ("PFPhotons"),
  fElectronName           (Names::gkElectronBrn),
  fGoodElectronName       (Names::gkElectronBrn),  
  fConversionName         (Names::gkMvfConversionBrn),  
  fTrackBranchName        (Names::gkTrackBrn),
  fPileUpDenName          (Names::gkPileupEnergyDensityBrn),
  fPVName                 (Names::gkPVBeamSpotBrn),
  fBeamspotName           (Names::gkBeamSpotBrn),
  fPFCandName             (Names::gkPFCandidatesBrn),
  fMCParticleName         (Names::gkMCPartBrn),
  fPileUpName             (Names::gkPileupInfoBrn),  
  fSuperClusterName       ("PFSuperClusters"),
  fPFMetName              ("PFMet"),
  fPFJetName              (Names::gkPFJetBrn),
  funcorrPFJetName        ("AKt5PFJets"),
  fGenJetName             ("AKT5GenJets"),
  fLeptonTagElectronsName ("HggLeptonTagElectrons"),
  fLeptonTagMuonsName     ("HggLeptonTagMuons"),

  fIsData                 (false),
  fPhotonsFromBranch      (kTRUE),  
  fPVFromBranch           (kTRUE),
  fGoodElectronsFromBranch(kTRUE),
  fPFJetsFromBranch       (kTRUE),
  // ----------------------------------------
  // collections....
  fPhotons                (0),
  fPFPhotons              (0),
  fElectrons              (0),
  fConversions            (0),
  fTracks                 (0),
  fPileUpDen              (0),
  fPV                     (0),
  fBeamspot               (0),
  fPFCands                (0),
  fMCParticles            (0),
  fPileUp                 (0),
  fSuperClusters          (0),
  fPFJets                 (0),
  fGenJets                (0),
  funcorrPFJets           (0),

  fLeptonTagElectrons     (0),
  fLeptonTagMuons         (0),

  fLoopOnGoodElectrons    (kFALSE),
  fApplyElectronVeto      (kTRUE),  
  fWriteDiphotonTree      (kTRUE),
  fWriteSingleTree        (kTRUE),
  fEnablePFPhotons        (kTRUE),
  fExcludeSinglePrompt    (kFALSE),
  fExcludeDoublePrompt    (kFALSE),
  fEnableJets             (kFALSE),
  fApplyJetId             (kFALSE),
  fApplyLeptonTag         (kFALSE),
  fApplyBTag              (kFALSE),
  fApplyPFMetCorrections  (kFALSE),
  fFillClusterArrays      (kFALSE),
  fFillVertexTree         (kFALSE),
  fPhFixDataFile          (gSystem->Getenv("CMSSW_BASE") +
		           TString("/src/MitPhysics/data/PhotonFixSTART42V13.dat")),
  fTupleName              ("hPhotonTree")
{
  // Constructor
}

PhotonTreeWriter::~PhotonTreeWriter()
{
  // Destructor
}

//--------------------------------------------------------------------------------------------------
void PhotonTreeWriter::Process()
{
  // ------------------------------------------------------------  
  // Process entries of the tree. 
  LoadEventObject(fPhotonBranchName,   fPhotons);
  LoadEventObject(fGoodElectronName,   fGoodElectrons);

  // lepton tag collections
  if( fApplyLeptonTag ) {
    LoadEventObject(fLeptonTagElectronsName, fLeptonTagElectrons);
    LoadEventObject(fLeptonTagMuonsName,     fLeptonTagMuons);
  }

  const BaseCollection *egcol = 0;
  if (fLoopOnGoodElectrons)
    egcol = fGoodElectrons;
  else
    egcol = fPhotons;
  if (egcol->GetEntries()<1)
    return;
  
  if (fEnablePFPhotons) LoadEventObject(fPFPhotonName,   fPFPhotons);
  LoadEventObject(fElectronName,       fElectrons);
  LoadEventObject(fConversionName,     fConversions);
  LoadEventObject(fTrackBranchName,    fTracks);
  LoadEventObject(fPileUpDenName,      fPileUpDen);
  LoadEventObject(fPVName,             fPV);    
  LoadEventObject(fBeamspotName,       fBeamspot);
  LoadEventObject(fPFCandName,         fPFCands);
  LoadEventObject(fSuperClusterName,   fSuperClusters);
  LoadEventObject(fPFMetName,          fPFMet);  
  if (fEnableJets){
    LoadEventObject(fPFJetName,        fPFJets);  
    //LoadEventObject(funcorrPFJetName,  funcorrPFJets);
    LoadBranch(funcorrPFJetName);
    //   if(!fIsData) LoadEventObject(fGenJetName,        fGenJets);
  }
  // ------------------------------------------------------------  
  // load event based information
  Int_t _numPU      = -1.;        // some sensible default values....
  Int_t _numPUminus = -1.;        // some sensible default values....
  Int_t _numPUplus  = -1.;        // some sensible default values....

  Float_t rho  = -99.;
  if( fPileUpDen->GetEntries() > 0 )
    rho  = (Double_t) fPileUpDen->At(0)->RhoRandomLowEta();
  
  const BaseVertex *bsp = dynamic_cast<const BaseVertex*>(fBeamspot->At(0));
    
  if( !fIsData ) {
    LoadBranch(fMCParticleName);
    LoadBranch(fPileUpName);
    if (fEnableJets) LoadEventObject(fGenJetName,        fGenJets);
  }  else fGenJets = NULL;
  
  if( !fIsData ) {
    for (UInt_t i=0; i<fPileUp->GetEntries(); ++i) {
      const PileupInfo *puinfo = fPileUp->At(i);
      if (puinfo->GetBunchCrossing()==0) _numPU = puinfo->GetPU_NumInteractions();
      else if (puinfo->GetBunchCrossing() == -1) _numPUminus = puinfo->GetPU_NumInteractions();
      else if (puinfo->GetBunchCrossing() ==  1) _numPUplus  = puinfo->GetPU_NumInteractions();
    }
  }

  // in case of a MC event, try to find Higgs and Higgs decay Z poisition
  Float_t _pth     = -100.;
  Float_t _decayZ  = -100.;
  Float_t _genmass = -100.;
  if (!fIsData)
    FindHiggsPtAndZ(_pth, _decayZ, _genmass);


  Double_t _spfMet = fPFMet->At(0)->SumEt();
  fDiphotonEvent->leptonTag = -1; // disabled

  fDiphotonEvent->rho = fPileUpDen->At(0)->RhoKt6PFJets();
  fDiphotonEvent->rho25 = fPileUpDen->At(0)->RhoRandomLowEta();
  fDiphotonEvent->rhoold = fPileUpDen->At(0)->Rho();
  fDiphotonEvent->genHiggspt = _pth;
  fDiphotonEvent->genHiggsZ = _decayZ;
  fDiphotonEvent->genmass = _genmass;  
  fDiphotonEvent->gencostheta = -99.;
  fDiphotonEvent->nVtx = fPV->GetEntries();
  fDiphotonEvent->bsX = fBeamspot->At(0)->X();
  fDiphotonEvent->bsY = fBeamspot->At(0)->Y();
  fDiphotonEvent->bsZ = fBeamspot->At(0)->Z();
  fDiphotonEvent->bsSigmaZ = fBeamspot->At(0)->SigmaZ();
  fDiphotonEvent->vtxX = (fDiphotonEvent->nVtx>0) ? fPV->At(0)->X() : -99.;
  fDiphotonEvent->vtxY = (fDiphotonEvent->nVtx>0) ? fPV->At(0)->Y() : -99.;  
  fDiphotonEvent->vtxZ = (fDiphotonEvent->nVtx>0) ? fPV->At(0)->Z() : -99.;
  fDiphotonEvent->numPU = _numPU;
  fDiphotonEvent->numPUminus = _numPUminus;
  fDiphotonEvent->numPUplus = _numPUplus;
  fDiphotonEvent->mass = -99.;
  fDiphotonEvent->ptgg = -99.;
  fDiphotonEvent->costheta = -99.;
  fDiphotonEvent->mt = -99.;
  fDiphotonEvent->cosphimet = -99.;
  fDiphotonEvent->mtele = -99.;
  fDiphotonEvent->cosphimetele = -99.;
  fDiphotonEvent->evt = GetEventHeader()->EvtNum();
  fDiphotonEvent->run = GetEventHeader()->RunNum();
  fDiphotonEvent->lumi = GetEventHeader()->LumiSec();
  fDiphotonEvent->evtcat = -99;
  fDiphotonEvent->nobj = fPhotons->GetEntries();
  fDiphotonEvent->pfmet = fPFMet->At(0)->Pt();
  fDiphotonEvent->pfmetphi = fPFMet->At(0)->Phi();
  fDiphotonEvent->pfmetx = fPFMet->At(0)->Px();
  fDiphotonEvent->pfmety = fPFMet->At(0)->Py();
  fDiphotonEvent->spfMet = _spfMet;
  fDiphotonEvent->masscor = -99.;
  fDiphotonEvent->masscorerr = -99.;
  fDiphotonEvent->masscorele = -99.;
  fDiphotonEvent->masscoreleerr = -99.;
  fDiphotonEvent->ismc = GetEventHeader()->IsMC();
  
  //jets
  const Jet *jet1 = 0;
  const Jet *jet2 = 0;
  const Jet *jetcentral = 0;

  fDiphotonEvent->jet1pt   = -99.;
  fDiphotonEvent->jet1eta  = -99.;
  fDiphotonEvent->jet1phi  = -99.;
  fDiphotonEvent->jet1mass = -99.;
  fDiphotonEvent->jet2pt   = -99.;
  fDiphotonEvent->jet2eta  = -99.;
  fDiphotonEvent->jet2phi  = -99.;
  fDiphotonEvent->jet2mass = -99.;
  fDiphotonEvent->jetcentralpt   = -99.;
  fDiphotonEvent->jetcentraleta  = -99.;
  fDiphotonEvent->jetcentralphi  = -99.;
  fDiphotonEvent->jetcentralmass = -99.;
  fDiphotonEvent->dijetpt = -99.;
  fDiphotonEvent->dijeteta = -99.;
  fDiphotonEvent->dijetphi = -99.;
  fDiphotonEvent->dijetmass = -99.;   
  fDiphotonEvent->jetetaplus = -99.;
  fDiphotonEvent->jetetaminus = -99.;
  fDiphotonEvent->zeppenfeld = -99.;
  fDiphotonEvent->dphidijetgg = -99.;
  
  //uncorrected jets
//   const Jet *uncorrjet1 = 0;
//   const Jet *uncorrjet2 = 0;
//   const Jet *uncorrjetcentral = 0;

//   fDiphotonEvent->uncorrjet1pt   = -99.;
//   fDiphotonEvent->uncorrjet1eta  = -99.;
//   fDiphotonEvent->uncorrjet1phi  = -99.;
//   fDiphotonEvent->uncorrjet1mass = -99.;
//   fDiphotonEvent->uncorrjet2pt   = -99.;
//   fDiphotonEvent->uncorrjet2eta  = -99.;
//   fDiphotonEvent->uncorrjet2phi  = -99.;
//   fDiphotonEvent->uncorrjet2mass = -99.;
//   fDiphotonEvent->uncorrjetcentralpt   = -99.;
//   fDiphotonEvent->uncorrjetcentraleta  = -99.;
//   fDiphotonEvent->uncorrjetcentralphi  = -99.;
//   fDiphotonEvent->uncorrjetcentralmass = -99.;
//   fDiphotonEvent->diuncorrjetpt = -99.;
//   fDiphotonEvent->diuncorrjeteta = -99.;
//   fDiphotonEvent->diuncorrjetphi = -99.;
//   fDiphotonEvent->diuncorrjetmass = -99.; 



  Int_t nhitsbeforevtxmax = 1;
  if (!fApplyElectronVeto)
    nhitsbeforevtxmax = 999;  
  
  if (egcol->GetEntries()>=2) {
    const Particle *p1 = 0;
    const Particle *p2 = 0;
    const Photon *phHard = 0;
    const Photon *phSoft = 0;
    const Electron *ele1 = 0;
    const Electron *ele2 = 0;
    const SuperCluster *sc1 = 0;
    const SuperCluster *sc2 = 0;
    
    if (fLoopOnGoodElectrons) {
      ele1 = fGoodElectrons->At(0);
      ele2 = fGoodElectrons->At(1);
      p1 = ele1;
      p2 = ele2;
      sc1 = ele1->SCluster();
      sc2 = ele2->SCluster();
      phHard = PhotonTools::MatchedPhoton(ele1,fPhotons);
      phSoft = PhotonTools::MatchedPhoton(ele2,fPhotons);
    }
    else {
      phHard = fPhotons->At(0);
      phSoft = fPhotons->At(1);
      p1 = phHard;
      p2 = phSoft;
      sc1 = phHard->SCluster();
      sc2 = phSoft->SCluster();      
      ele1 = PhotonTools::MatchedElectron(phHard,fGoodElectrons);
      ele2 = PhotonTools::MatchedElectron(phSoft,fGoodElectrons);
    }
    
    const DecayParticle *conv1 = PhotonTools::MatchedConversion(sc1,fConversions,bsp,
								nhitsbeforevtxmax);
    const DecayParticle *conv2 = PhotonTools::MatchedConversion(sc2,fConversions,bsp,
								nhitsbeforevtxmax);
    
    const SuperCluster *pfsc1 = PhotonTools::MatchedPFSC(sc1,fPFPhotons, fElectrons);
    const SuperCluster *pfsc2 = PhotonTools::MatchedPFSC(sc2,fPFPhotons, fElectrons);
    
    const MCParticle *phgen1 = 0;
    const MCParticle *phgen2 = 0;
    if (!fIsData) {
      phgen1 = PhotonTools::MatchMC(p1,fMCParticles,!fApplyElectronVeto);
      phgen2 = PhotonTools::MatchMC(p2,fMCParticles,!fApplyElectronVeto);
    }
    
    if (fExcludeSinglePrompt && (phgen1 || phgen2))
      return;
    if (fExcludeDoublePrompt && (phgen1 && phgen2))
      return;
    
    if (!fLoopOnGoodElectrons && phHard->HasPV()) {
      fDiphotonEvent->vtxX = phHard->PV()->X();
      fDiphotonEvent->vtxY = phHard->PV()->Y();
      fDiphotonEvent->vtxZ = phHard->PV()->Z();
      fDiphotonEvent->vtxprob = phHard->VtxProb();
    }
    
    // fill Btag information... set to true if One jet fullfills 
    fDiphotonEvent -> btagJet1       = -1.;
    fDiphotonEvent -> btagJet1Pt     = -99.;
    fDiphotonEvent -> btagJet1Eta    = -99.;

    fDiphotonEvent -> btagJet2       = -1.;
    fDiphotonEvent -> btagJet2Pt     = -99.;
    fDiphotonEvent -> btagJet2Eta    = -99.;

    if( fApplyBTag && fEnableJets ) {
      float highTag     = 0.;
      float highJetPt   = 0.;
      float highJetEta  = 0.;

      float trailTag     = 0.;
      float trailJetPt   = 0.;
      float trailJetEta  = 0.;

      for (UInt_t ijet=0; ijet<fPFJets->GetEntries();++ijet) {
	const Jet *jet = fPFJets->At(ijet);
	if( jet->Pt() < 20. || jet->AbsEta() > 2.4 ) continue;
	if ( jet->CombinedSecondaryVertexBJetTagsDisc() > highTag ) {
	  
	  trailTag = highTag;
	  trailJetPt = highJetPt;
	  trailJetEta = highJetEta;

	  highTag    = jet->CombinedSecondaryVertexBJetTagsDisc();
	  highJetPt  = jet->Pt();
	  highJetEta = jet->Eta();

	} else if ( jet->CombinedSecondaryVertexBJetTagsDisc() > trailTag ) {

	  trailTag    = jet->CombinedSecondaryVertexBJetTagsDisc();
	  trailJetPt  = jet->Pt();
	  trailJetEta = jet->Eta();
	}
      }
      fDiphotonEvent -> btagJet1    = highTag;
      fDiphotonEvent -> btagJet1Pt  = highJetPt;
      fDiphotonEvent -> btagJet1Eta = highJetEta;      

      fDiphotonEvent -> btagJet1    = trailTag;
      fDiphotonEvent -> btagJet1Pt  = trailJetPt;
      fDiphotonEvent -> btagJet1Eta = trailJetEta;      

    }
    
    //fill jet variables
    const Vertex *selvtx = fPV->At(0);
    if (!fLoopOnGoodElectrons && phHard->HasPV()) selvtx = phHard->PV();
    if (fEnableJets) {
      for (UInt_t ijet=0; ijet<fPFJets->GetEntries();++ijet) {
        const Jet *jet = fPFJets->At(ijet);
        if (jet->AbsEta()<4.7 && MathUtils::DeltaR(jet,p1)>0.5 && MathUtils::DeltaR(jet,p2)>0.5) {
          const PFJet *pfjet = dynamic_cast<const PFJet*>(jet);
          if (!pfjet) continue;
          if (fApplyJetId && !fJetId.passCut(pfjet,selvtx,fPV)) continue;
          if (!jet1) jet1 = jet;
          else if (!jet2) jet2 = jet;
          else if (!jetcentral && 0) jetcentral = jet;
        }
        if (jet1&&jet2&&jetcentral) break;
      }      
    }
    
    if (jet1) {
      fDiphotonEvent->jet1pt   = jet1->Pt();
      fDiphotonEvent->jet1eta  = jet1->Eta();
      fDiphotonEvent->jet1phi  = jet1->Phi();
      fDiphotonEvent->jet1mass = jet1->Mass();
    }
    
    if (jet2) {
      fDiphotonEvent->jet2pt   = jet2->Pt();
      fDiphotonEvent->jet2eta  = jet2->Eta();
      fDiphotonEvent->jet2phi  = jet2->Phi();
      fDiphotonEvent->jet2mass = jet2->Mass();
    }
  
    if (jetcentral) {
      fDiphotonEvent->jetcentralpt   = jetcentral->Pt();
      fDiphotonEvent->jetcentraleta  = jetcentral->Eta();
      fDiphotonEvent->jetcentralphi  = jetcentral->Phi();
      fDiphotonEvent->jetcentralmass = jetcentral->Mass();
    }
    
    if (jet1&&jet2){
      FourVectorM momjj = (jet1->Mom() + jet2->Mom());
      
      fDiphotonEvent->dijetpt =  momjj.Pt();
      fDiphotonEvent->dijeteta = momjj.Eta();
      fDiphotonEvent->dijetphi = momjj.Phi();
      fDiphotonEvent->dijetmass = momjj.M();    
      
      if (jet1->Eta()>jet2->Eta()) {
        fDiphotonEvent->jetetaplus = jet1->Eta();
        fDiphotonEvent->jetetaminus = jet2->Eta();
      }
      else {
        fDiphotonEvent->jetetaplus = jet2->Eta();
        fDiphotonEvent->jetetaminus = jet1->Eta();      
      }
    }
    
    //added gen. info of whether a lep. or nutrino is from W or Z --Heng 02/14/2012 12:30 EST
    Double_t _fromZ = -99;
    Double_t _fromW = -99;
    Float_t _zpt = -99;
    Float_t _zEta = -99;
    Float_t _allZpt = -99;
    Float_t _allZEta = -99;

    if( !fIsData ){
      
      // loop over all GEN particles and look for nutrinos whoes mother is Z
      for(UInt_t j=0; j<fMCParticles->GetEntries(); ++j) {
	const MCParticle* p = fMCParticles->At(j);
	if( p->AbsPdgId()==23 ||p->AbsPdgId()==32 || p->AbsPdgId()==33 ) {
	  _allZpt=p->Pt();
	  _allZEta=p->Eta();
	  if (p->HasDaughter(12,kFALSE) || p->HasDaughter(14,kFALSE) || p->HasDaughter(16,kFALSE) ||p->HasDaughter(18,kFALSE) ) {
	  _fromZ=1;
	  _zpt=p->Pt();
	  _zEta=p->Eta();
	  }
	}
	else _fromW=1;
      }
    }
    
    fDiphotonEvent->fromZ = _fromZ;
    fDiphotonEvent->fromW = _fromW;
    fDiphotonEvent->zpt = _zpt;
    fDiphotonEvent->zEta = _zEta;
    fDiphotonEvent->allZpt = _allZpt;
    fDiphotonEvent->allZEta = _allZEta;

    Double_t _dphiMetgg = -99;
    Double_t _cosdphiMetgg = -99;
    Double_t _dphiPhPh = -99;
    Double_t _mass = -99.;
    Double_t _masserr = -99.;
    Double_t _masserrsmeared = -99.;
    Double_t _masserrwrongvtx = -99.;
    Double_t _masserrsmearedwrongvtx = -99.;
    Double_t _ptgg = -99.;
    Double_t _etagg = -99.;
    Double_t _phigg = -99.;    
    Double_t _costheta = -99.;
    PhotonTools::DiphotonR9EtaPtCats _evtcat = PhotonTools::kOctCat0;
    if (phHard && phSoft) {
      _dphiMetgg = MathUtils::DeltaPhi((phHard->Mom()+phSoft->Mom()).Phi(),fPFMet->At(0)->Phi());
      _cosdphiMetgg = TMath::Cos(_dphiMetgg);
      _dphiPhPh = MathUtils::DeltaPhi((phHard->Mom()).Phi(),(phSoft->Mom()).Phi());
      _mass = (phHard->Mom()+phSoft->Mom()).M();
      _masserr = 0.5*_mass*TMath::Sqrt(phHard->EnergyErr()*phHard->EnergyErr()/phHard->E()/phHard->E() + phSoft->EnergyErr()*phSoft->EnergyErr()/phSoft->E()/phSoft->E());
      _masserrsmeared = 0.5*_mass*TMath::Sqrt(phHard->EnergyErrSmeared()*phHard->EnergyErrSmeared()/phHard->E()/phHard->E() + phSoft->EnergyErrSmeared()*phSoft->EnergyErrSmeared()/phSoft->E()/phSoft->E());
      _ptgg = (phHard->Mom()+phSoft->Mom()).Pt();
      _etagg = (phHard->Mom()+phSoft->Mom()).Eta();
      _phigg = (phHard->Mom()+phSoft->Mom()).Phi();
      _costheta = ThreeVector(phHard->Mom()).Unit().Dot(ThreeVector(phSoft->Mom()).Unit());
      _evtcat = PhotonTools::DiphotonR9EtaPtCat(phHard,phSoft);
      
      const Double_t dz = sqrt(2.0)*5.8;
      Double_t deltamvtx = _mass*VertexTools::DeltaMassVtx(phHard->CaloPos().X(),
							   phHard->CaloPos().Y(),
							   phHard->CaloPos().Z(),
							   phSoft->CaloPos().X(),
							   phSoft->CaloPos().Y(),
							   phSoft->CaloPos().Z(),
							   fDiphotonEvent->vtxX,
							   fDiphotonEvent->vtxY,
							   fDiphotonEvent->vtxZ,
							   dz);
      fDiphotonEvent->deltamvtx = deltamvtx;
            
      _masserrwrongvtx        = TMath::Sqrt(_masserr*_masserr + deltamvtx*deltamvtx);
      _masserrsmearedwrongvtx = TMath::Sqrt(_masserrsmeared*_masserrsmeared + deltamvtx*deltamvtx);
            
      if (jet1 && jet2) {
        fDiphotonEvent->zeppenfeld = TMath::Abs(_etagg - 0.5*(jet1->Eta()+jet2->Eta()));
        fDiphotonEvent->dphidijetgg = MathUtils::DeltaPhi( (jet1->Mom()+jet2->Mom()).Phi(), _phigg );
      }
      
      PFJetOArr pfjets;
      for (UInt_t ijet=0; ijet<fPFJets->GetEntries(); ++ijet) {
        const PFJet *pfjet = dynamic_cast<const PFJet*>(fPFJets->At(ijet));
        if (pfjet && MathUtils::DeltaR(*pfjet,*phHard)>0.3 && MathUtils::DeltaR(*pfjet,*phSoft)>0.3) pfjets.Add(pfjet);
      }
      
      PFCandidateOArr pfcands;
      for (UInt_t icand=0; icand<fPFCands->GetEntries(); ++icand) {
        const PFCandidate *pfcand = fPFCands->At(icand);
        if (MathUtils::DeltaR(*pfcand,*phHard)>0.1 && MathUtils::DeltaR(*pfcand,*phSoft)>0.1) pfcands.Add(pfcand);
      }      
      
      const Vertex *firstvtx = fPV->At(0);
      const Vertex *selvtx = fPV->At(0);
      
      if (!fLoopOnGoodElectrons && phHard->HasPV()) {
        selvtx = phHard->PV();
      }
     
      if (0) //disable for now for performance reasons
      {
        Met mmet = fMVAMet.GetMet(  false,
                                  0.,0.,0.,
                                  fPFMet->At(0),
                                  &pfcands,selvtx,fPV, fPileUpDen->At(0)->Rho(),
                                  &pfjets,
                                  int(fPV->GetEntries()),
                                  kFALSE);      
                                  
        TMatrixD *metcov = fMVAMet.GetMetCovariance();
        
        ThreeVector fullmet(mmet.Px() - phHard->Px() - phSoft->Px(),
                            mmet.Py() - phHard->Py() - phSoft->Py(),
                            0.);
        
        ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> > mcov;
        mcov(0,0) = (*metcov)(0,0);
        mcov(0,1) = (*metcov)(0,1);
        mcov(1,0) = (*metcov)(1,0);
        mcov(1,1) = (*metcov)(1,1);
        ROOT::Math::SVector<double,2> vmet;
        vmet(0) = fullmet.X();
        vmet(1) = fullmet.Y();
        mcov.Invert();
        Double_t metsig = sqrt(ROOT::Math::Similarity(mcov,vmet));
                                  
        fDiphotonEvent->mvametsel = fullmet.Rho();
        fDiphotonEvent->mvametselphi = fullmet.Phi();
        fDiphotonEvent->mvametselx = fullmet.X();
        fDiphotonEvent->mvametsely = fullmet.Y();
        fDiphotonEvent->mvametselsig = metsig;
      }
      
      if (0) //disable for now for performance reasons
      {
        Met mmet = fMVAMet.GetMet(  false,
                                  0.,0.,0.,
                                  fPFMet->At(0),
                                  &pfcands,firstvtx,fPV, fPileUpDen->At(0)->Rho(),
                                  &pfjets,
                                  int(fPV->GetEntries()),
                                  kFALSE);      
                                  
        TMatrixD *metcov = fMVAMet.GetMetCovariance();
        
        ThreeVector fullmet(mmet.Px() - phHard->Px() - phSoft->Px(),
                            mmet.Py() - phHard->Py() - phSoft->Py(),
                            0.);
        
        ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> > mcov;
        mcov(0,0) = (*metcov)(0,0);
        mcov(0,1) = (*metcov)(0,1);
        mcov(1,0) = (*metcov)(1,0);
        mcov(1,1) = (*metcov)(1,1);
        ROOT::Math::SVector<double,2> vmet;
        vmet(0) = fullmet.X();
        vmet(1) = fullmet.Y();
        mcov.Invert();
        Double_t metsig = sqrt(ROOT::Math::Similarity(mcov,vmet));
                                  
        fDiphotonEvent->mvametfirst = fullmet.Rho();
        fDiphotonEvent->mvametfirstphi = fullmet.Phi();
        fDiphotonEvent->mvametfirstx = fullmet.X();
        fDiphotonEvent->mvametfirsty = fullmet.Y();
        fDiphotonEvent->mvametfirstsig = metsig;
      }      
      
    }

    fDiphotonEvent->corrpfmet = -99.;
    fDiphotonEvent->corrpfmetphi = -99.;
    fDiphotonEvent->corrpfmetx = -99.;
    fDiphotonEvent->corrpfmety = -99.;

    Met *corrMet =NULL;
    
    if (fApplyPFMetCorrections){
    corrMet = new Met(fPFMet->At(0)->Px(),fPFMet->At(0)->Py());
    
    if (!fIsData){
      PFMetCorrectionTools::correctMet(corrMet,phHard,phSoft,1,1,funcorrPFJets,fGenJets,fPFJets);
      PFMetCorrectionTools::shiftMet(corrMet,fIsData,_spfMet);
    }
    else {
      PFMetCorrectionTools::shiftMet(corrMet,fIsData,_spfMet);
       PFMetCorrectionTools::correctMet(corrMet,phHard,phSoft,0,1,funcorrPFJets,fGenJets,fPFJets);
     }    
 
    fDiphotonEvent->corrpfmet = corrMet->Pt();
    fDiphotonEvent->corrpfmetphi = corrMet->Phi();
    fDiphotonEvent->corrpfmetx = corrMet->Px();
    fDiphotonEvent->corrpfmety = corrMet->Py();

    delete corrMet;
    }

    Float_t _massele = -99.;
    Float_t _ptee = -99.;
    Float_t _costhetaele = -99.;
    if (ele1 && ele2) {
      _massele = (ele1->Mom()+ele2->Mom()).M();
      _ptee = (ele1->Mom()+ele2->Mom()).Pt();
      _costhetaele = ThreeVector(ele1->Mom()).Unit().Dot(ThreeVector(ele2->Mom()).Unit());      
    }
    
    Float_t _gencostheta = -99.;
    if (phgen1 && phgen2) {
      _gencostheta = ThreeVector(phgen1->Mom()).Unit().Dot(ThreeVector(phgen2->Mom()).Unit());
    }
  
    fDiphotonEvent->gencostheta = _gencostheta;
    fDiphotonEvent->dphiMetgg = _dphiMetgg;
    fDiphotonEvent->cosdphiMetgg = _cosdphiMetgg;
    fDiphotonEvent->dphiPhPh = _dphiPhPh;
    fDiphotonEvent->mass = _mass;
    fDiphotonEvent->masserr = _masserr;
    fDiphotonEvent->masserrsmeared = _masserrsmeared;
    fDiphotonEvent->masserrwrongvtx = _masserrwrongvtx;
    fDiphotonEvent->masserrsmearedwrongvtx = _masserrsmearedwrongvtx;    
    fDiphotonEvent->ptgg = _ptgg;
    fDiphotonEvent->etagg = _etagg;
    fDiphotonEvent->phigg = _phigg;
    fDiphotonEvent->costheta =  _costheta;;
    fDiphotonEvent->massele = _massele;
    fDiphotonEvent->ptee = _ptee;
    fDiphotonEvent->costhetaele =  _costhetaele;    
    fDiphotonEvent->evtcat = _evtcat;

    fDiphotonEvent->photons[0].SetVars(phHard,conv1,ele1,pfsc1,phgen1,fPhfixph,fPhfixele,fTracks,fPV,rho,fFillClusterArrays,fElectrons,fApplyElectronVeto);
    fDiphotonEvent->photons[1].SetVars(phSoft,conv2,ele2,pfsc2,phgen2,fPhfixph,fPhfixele,fTracks,fPV,rho,fFillClusterArrays,fElectrons,fApplyElectronVeto);
    
    Float_t ph1ecor    = fDiphotonEvent->photons[0].Ecor();
    Float_t ph1ecorerr = fDiphotonEvent->photons[0].Ecorerr();
    Float_t ph2ecor    = fDiphotonEvent->photons[1].Ecor();
    Float_t ph2ecorerr = fDiphotonEvent->photons[1].Ecorerr();

    Float_t ph1ecorele    = fDiphotonEvent->photons[0].Ecorele();
    Float_t ph1ecoreleerr = fDiphotonEvent->photons[0].Ecoreleerr();
    Float_t ph2ecorele    = fDiphotonEvent->photons[1].Ecorele();
    Float_t ph2ecoreleerr = fDiphotonEvent->photons[1].Ecoreleerr();
    
    
    fDiphotonEvent->masscor    = TMath::Sqrt(2.0*ph1ecor*ph2ecor*(1.0-fDiphotonEvent->costheta));
    fDiphotonEvent->masscorerr = 0.5*fDiphotonEvent->masscor*
      TMath::Sqrt(ph1ecorerr*ph1ecorerr/ph1ecor/ph1ecor + ph2ecorerr*ph2ecorerr/ph2ecor/ph2ecor);
    
    fDiphotonEvent->masscorele = TMath::Sqrt(2.0*ph1ecorele*ph2ecorele*
					     (1.0-fDiphotonEvent->costheta));
    fDiphotonEvent->masscoreleerr = 0.5*fDiphotonEvent->masscorele*
      TMath::Sqrt(ph1ecoreleerr*ph1ecoreleerr/ph1ecorele/ph1ecorele +
		  ph2ecoreleerr*ph2ecoreleerr/ph2ecorele/ph2ecorele);    
    
    //printf("r9 = %5f, photon sigieie = %5f, seed sigieie = %5f\n",phHard->R9(),
    //       phHard->CoviEtaiEta(),sqrt(phHard->SCluster()->Seed()->CoviEtaiEta()));

    // MuonStuff
    fDiphotonEvent-> muonPt  = -99.;
    fDiphotonEvent-> muonEta  = -99.;
    fDiphotonEvent-> muDR1  = -99.;
    fDiphotonEvent-> muDR2  = -99.;
    fDiphotonEvent-> muIso1  = -99.;
    fDiphotonEvent-> muIso2  = -99.;
    fDiphotonEvent-> muIso3  = -99.;
    fDiphotonEvent-> muIso4  = -99.;
    fDiphotonEvent-> muD0  = -99.;
    fDiphotonEvent-> muDZ  = -99.;
    fDiphotonEvent-> muChi2  = -99.;   
    fDiphotonEvent-> muNhits = -99;
    fDiphotonEvent-> muNpixhits = -99;
    fDiphotonEvent-> muNegs = -99;
    fDiphotonEvent-> muNMatch = -99;

    // Electron Stuff
    fDiphotonEvent-> elePt = -99.;
    fDiphotonEvent-> eleEta = -99.;
    fDiphotonEvent-> eleSCEta = -99.;     
    fDiphotonEvent-> eleIso1 = -99.;
    fDiphotonEvent-> eleIso2 = -99.;
    fDiphotonEvent-> eleIso3 = -99.;
    fDiphotonEvent-> eleIso4 = -99.;
    fDiphotonEvent-> eleDist = -99.;
    fDiphotonEvent-> eleDcot = -99.;
    fDiphotonEvent-> eleCoviee = -99.;
    fDiphotonEvent-> eleDphiin = -99.;
    fDiphotonEvent-> eleDetain = -99.;
    fDiphotonEvent-> eleDR1 = -99.;
    fDiphotonEvent-> eleDR2 = -99.;
    fDiphotonEvent-> eleMass1 = -99.;
    fDiphotonEvent-> eleMass2 = -99.;
    fDiphotonEvent-> eleNinnerHits = -99;     

    if( fApplyLeptonTag ) {
      // perform lepton tagging
      // the diphoton event record will have one more entry; i.e. leptonTag
      // leptonTag = -1   -> lepton-taggng was swicthed off
      //           =  0   -> event tagged as 'non-lepton-event'
      //           = +1   -> event tagged as muon-event
      //           = +2   -> event tagged as electron-event
      fDiphotonEvent->leptonTag = 0;

      if ( fLeptonTagMuons->GetEntries() > 0 ) {

	// need to have dR > 1 for with respect to both photons
	for(UInt_t iMuon = 0; iMuon <fLeptonTagMuons->GetEntries(); ++iMuon) {
	  if(MathUtils::DeltaR(fLeptonTagMuons->At(iMuon),phHard) < 1.) continue;
	  if(MathUtils::DeltaR(fLeptonTagMuons->At(iMuon),phSoft) < 1.) continue;

	  fDiphotonEvent->leptonTag = 1;

	  fDiphotonEvent-> muonPt  = fLeptonTagMuons->At(iMuon)->Pt();
	  fDiphotonEvent-> muonEta = fLeptonTagMuons->At(iMuon)->Eta();
	  fDiphotonEvent-> muDR1   = MathUtils::DeltaR(fLeptonTagMuons->At(iMuon),phHard);
	  fDiphotonEvent-> muDR2   = MathUtils::DeltaR(fLeptonTagMuons->At(iMuon),phSoft);
	  fDiphotonEvent-> muIso1   = (fLeptonTagMuons->At(iMuon)->IsoR03SumPt() + fLeptonTagMuons->At(iMuon)->IsoR03EmEt() + fLeptonTagMuons->At(iMuon)->IsoR03HadEt() - fPileUpDen->At(0)->RhoRandomLowEta() * TMath::Pi() * 0.3 * 0.3)/ fLeptonTagMuons->At(iMuon)->Pt();
	  fDiphotonEvent-> muIso2   = (fLeptonTagMuons->At(iMuon)->IsoR03SumPt() + fLeptonTagMuons->At(iMuon)->IsoR03EmEt() + fLeptonTagMuons->At(iMuon)->IsoR03HadEt() - fPileUpDen->At(0)->RhoRandom() * TMath::Pi() * 0.3 * 0.3)/ fLeptonTagMuons->At(iMuon)->Pt();
	  fDiphotonEvent-> muIso3   = (fLeptonTagMuons->At(iMuon)->IsoR03SumPt() + fLeptonTagMuons->At(iMuon)->IsoR03EmEt() + fLeptonTagMuons->At(iMuon)->IsoR03HadEt() - fPileUpDen->At(0)->RhoLowEta() * TMath::Pi() * 0.3 * 0.3)/ fLeptonTagMuons->At(iMuon)->Pt();
	  fDiphotonEvent-> muIso4   = (fLeptonTagMuons->At(iMuon)->IsoR03SumPt() + fLeptonTagMuons->At(iMuon)->IsoR03EmEt() + fLeptonTagMuons->At(iMuon)->IsoR03HadEt() - fPileUpDen->At(0)->Rho() * TMath::Pi() * 0.3 * 0.3)/ fLeptonTagMuons->At(iMuon)->Pt();
	  fDiphotonEvent-> muD0  = TMath::Abs(fLeptonTagMuons->At(iMuon)->BestTrk()->D0Corrected(*fPV->At(0)));
	  fDiphotonEvent-> muDZ  = TMath::Abs(fLeptonTagMuons->At(iMuon)->BestTrk()->DzCorrected(*fPV->At(0)));
	  fDiphotonEvent-> muChi2  = fLeptonTagMuons->At(iMuon)->GlobalTrk()->Chi2()/fLeptonTagMuons->At(iMuon)->GlobalTrk()->Ndof();
	  	  
	  fDiphotonEvent-> muNhits = fLeptonTagMuons->At(iMuon)->BestTrk()->NHits();
	  fDiphotonEvent-> muNpixhits = fLeptonTagMuons->At(iMuon)->BestTrk()->NPixelHits();
	  fDiphotonEvent-> muNegs = fLeptonTagMuons->At(iMuon)->NSegments();
	  fDiphotonEvent-> muNMatch = fLeptonTagMuons->At(iMuon)->NMatches();

	  break;
	}
      }
      if ( fDiphotonEvent->leptonTag < 1 && fLeptonTagElectrons->GetEntries() > 0 ) {
	for(UInt_t iEle = 0; iEle < fLeptonTagElectrons->GetEntries(); ++iEle) {
	  if(MathUtils::DeltaR(fLeptonTagElectrons->At(iEle),phHard) < 1.) continue;
	  if(MathUtils::DeltaR(fLeptonTagElectrons->At(iEle),phSoft) < 1.) continue;

	  // here we also check the mass ....
	  if ( TMath::Abs( (phHard->Mom()+fLeptonTagElectrons->At(iEle)->Mom()).M()-91.19 ) < 5. ) continue;
	  if ( TMath::Abs( (phSoft->Mom()+fLeptonTagElectrons->At(iEle)->Mom()).M()-91.19 ) < 5. ) continue;
	  
	  fDiphotonEvent->leptonTag = 2;

	  fDiphotonEvent-> elePt = fLeptonTagElectrons->At(iEle)->Pt();
	  fDiphotonEvent-> eleEta = fLeptonTagElectrons->At(iEle)->Eta();
	  fDiphotonEvent-> eleSCEta = fLeptonTagElectrons->At(iEle)->SCluster()->Eta();
	  fDiphotonEvent-> eleIso1 = (fLeptonTagElectrons->At(iEle)->TrackIsolationDr03() + fLeptonTagElectrons->At(iEle)->EcalRecHitIsoDr03() + fLeptonTagElectrons->At(iEle)->HcalTowerSumEtDr03() - fPileUpDen->At(0)->RhoRandomLowEta() * TMath::Pi() * 0.3 * 0.3)/fDiphotonEvent-> elePt;
	  fDiphotonEvent-> eleIso2 = (fLeptonTagElectrons->At(iEle)->TrackIsolationDr03() + fLeptonTagElectrons->At(iEle)->EcalRecHitIsoDr03() + fLeptonTagElectrons->At(iEle)->HcalTowerSumEtDr03() - fPileUpDen->At(0)->RhoRandom() * TMath::Pi() * 0.3 * 0.3)/fDiphotonEvent-> elePt;
	  fDiphotonEvent-> eleIso3 = (fLeptonTagElectrons->At(iEle)->TrackIsolationDr03() + fLeptonTagElectrons->At(iEle)->EcalRecHitIsoDr03() + fLeptonTagElectrons->At(iEle)->HcalTowerSumEtDr03() - fPileUpDen->At(0)->RhoLowEta() * TMath::Pi() * 0.3 * 0.3)/fDiphotonEvent-> elePt;
	  fDiphotonEvent-> eleIso4 = (fLeptonTagElectrons->At(iEle)->TrackIsolationDr03() + fLeptonTagElectrons->At(iEle)->EcalRecHitIsoDr03() + fLeptonTagElectrons->At(iEle)->HcalTowerSumEtDr03() - fPileUpDen->At(0)->Rho() * TMath::Pi() * 0.3 * 0.3)/fDiphotonEvent-> elePt;
	  fDiphotonEvent-> eleDist = fLeptonTagElectrons->At(iEle)->ConvPartnerDist();
	  fDiphotonEvent-> eleDcot = fLeptonTagElectrons->At(iEle)->ConvPartnerDCotTheta();
	  fDiphotonEvent-> eleCoviee = fLeptonTagElectrons->At(iEle)->CoviEtaiEta();
	  fDiphotonEvent-> eleDphiin = TMath::Abs(fLeptonTagElectrons->At(iEle)->DeltaPhiSuperClusterTrackAtVtx());
	  fDiphotonEvent-> eleDetain = TMath::Abs(fLeptonTagElectrons->At(iEle)->DeltaEtaSuperClusterTrackAtVtx());
	  fDiphotonEvent-> eleDR1 = MathUtils::DeltaR(fLeptonTagElectrons->At(iEle),phHard);
	  fDiphotonEvent-> eleDR2 = MathUtils::DeltaR(fLeptonTagElectrons->At(iEle),phSoft);
	  fDiphotonEvent-> eleMass1 = (phHard->Mom()+fLeptonTagElectrons->At(iEle)->Mom()).M();
	  fDiphotonEvent-> eleMass2 = (phSoft->Mom()+fLeptonTagElectrons->At(iEle)->Mom()).M();
	  fDiphotonEvent-> eleNinnerHits =      fLeptonTagElectrons->At(iEle)->Trk()->NExpectedHitsInner();

	  break;
	}
      }
    }

    if (fWriteDiphotonTree)
      hCiCTuple->Fill();  
    
    if (fFillVertexTree && phHard && phSoft) {
      for (UInt_t ivtx = 0; ivtx<fPV->GetEntries(); ++ivtx) {
	const Vertex *v = fPV->At(ivtx);
	fDiphotonVtx->SetVars(v,phHard,phSoft,fPFCands,ivtx,fPV->GetEntries(),_decayZ);
	hVtxTree->Fill();
      }
    }
    
  }

  if (!fWriteSingleTree)
    return;

  for (UInt_t iph = 0; iph<egcol->GetEntries(); ++iph) {
    const Particle *p = 0;
    const Photon *ph = 0;
    const Electron *ele = 0;
    const SuperCluster *sc = 0;
    if (fLoopOnGoodElectrons) {
      ele = fGoodElectrons->At(iph);
      p = ele;
      sc = ele->SCluster();
      ph = PhotonTools::MatchedPhoton(ele,fPhotons);
    }
    else {
      ph = fPhotons->At(iph);
      p = ph;
      sc = ph->SCluster();
      ele = PhotonTools::MatchedElectron(ph,fGoodElectrons);    
    }
    
    const DecayParticle *conv = PhotonTools::MatchedConversion(sc,fConversions,bsp,
							       nhitsbeforevtxmax);
    const SuperCluster  *pfsc = PhotonTools::MatchedPFSC(sc,fPFPhotons, fElectrons);
        
    if (!fLoopOnGoodElectrons && ph->HasPV()) {
      fDiphotonEvent->vtxZ = ph->PV()->Z();
    }

    const MCParticle *phgen = 0;
    if( !fIsData ) {
      phgen = PhotonTools::MatchMC(p,fMCParticles,!fApplyElectronVeto);
    }

    if (fExcludeSinglePrompt && phgen) return;

    fDiphotonEvent->mt = -99.;
    fDiphotonEvent->cosphimet = -99.;
    fDiphotonEvent->mtele = -99.;
    fDiphotonEvent->cosphimetele = -99.;

    if (ph) {
      fDiphotonEvent->cosphimet = TMath::Cos(ph->Phi()-fPFMet->At(0)->Phi());
      fDiphotonEvent->mt = TMath::Sqrt(2.0*fPFMet->At(0)->Pt()*ph->Pt()*
				       (1.0-fDiphotonEvent->cosphimet));
    }
    
    if (ele) {
      fDiphotonEvent->cosphimetele = TMath::Cos(ele->Phi()-fPFMet->At(0)->Phi());
      fDiphotonEvent->mtele        = TMath::Sqrt(2.0*fPFMet->At(0)->Pt()*ele->Pt()*
						 (1.0-fDiphotonEvent->cosphimetele));      
    }
    
    fSinglePhoton->SetVars(ph,conv,ele,pfsc,phgen,fPhfixph,fPhfixele,fTracks,fPV,rho,fFillClusterArrays,
			   fElectrons,fApplyElectronVeto);
    hCiCTupleSingle->Fill();
  }



    
  return;
}

//--------------------------------------------------------------------------------------------------
void PhotonTreeWriter::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the photon collection branch.

  if( fApplyLeptonTag ) {
    ReqEventObject(fLeptonTagElectronsName,    fLeptonTagElectrons,    false);  
    ReqEventObject(fLeptonTagMuonsName,        fLeptonTagMuons,        false);  
  }

  ReqEventObject(fPhotonBranchName,fPhotons,      fPhotonsFromBranch);
  if (fEnablePFPhotons) ReqEventObject(fPFPhotonName,fPFPhotons,      true);
  ReqEventObject(fTrackBranchName, fTracks,       true);
  ReqEventObject(fElectronName,    fElectrons,    true);  
  ReqEventObject(fGoodElectronName,fGoodElectrons,fGoodElectronsFromBranch);  
  ReqEventObject(fPileUpDenName,   fPileUpDen,    true);
  ReqEventObject(fPVName,          fPV,           fPVFromBranch);
  ReqEventObject(fConversionName,  fConversions,  true);
  ReqEventObject(fBeamspotName,    fBeamspot,     true);
  ReqEventObject(fPFCandName,      fPFCands,      true);
  ReqEventObject(fSuperClusterName,fSuperClusters,true);
  ReqEventObject(fPFMetName,       fPFMet,        true);
  if (fEnableJets){
    ReqEventObject(fPFJetName,       fPFJets,       fPFJetsFromBranch);
    ReqBranch(funcorrPFJetName, funcorrPFJets);
    //   if (!fIsData) ReqEventObject(fGenJetName, fGenJets, true);
  }
  if (!fIsData) {
    ReqBranch(fPileUpName,         fPileUp);
    ReqBranch(fMCParticleName,     fMCParticles);
    if (fEnableJets) ReqEventObject(fGenJetName, fGenJets, true);
  }
  if (fIsData) {
    fPhFixDataFile = gSystem->Getenv("CMSSW_BASE") +
      TString("/src/MitPhysics/data/PhotonFixGRPV22.dat");
  }
  else {
    fPhFixDataFile = gSystem->Getenv("CMSSW_BASE") +
      TString("/src/MitPhysics/data/PhotonFixSTART42V13.dat");
  }

  fPhfixph.initialise("4_2",std::string(fPhFixDataFile));
  fPhfixele.initialise("4_2e",std::string(fPhFixDataFile));
  
//   fMVAMet.Initialize(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_lowpt.weights.xml"))),
//                       TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_highpt.weights.xml"))),
//                       TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/Utils/python/JetIdParams_cfi.py")),
//                       TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmet_42.root"))),
//                       TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetphi_42.root"))),
//                       TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu1_42.root"))),
//                       TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu2_42.root")))
//                       );  
  
  fMVAMet.Initialize(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_lowpt.weights.xml"))),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_highpt.weights.xml"))),
                      TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/Utils/python/JetIdParams_cfi.py")),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmet_52.root"))),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetphi_52.root"))),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu1cov_52.root"))),
                      TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu2cov_52.root")))
                      );                      
                 
  fJetId.Initialize(JetIDMVA::kMedium,
                          TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_lowpt.weights.xml"))),
                          TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_highpt.weights.xml"))),
                          JetIDMVA::kCut,
                          TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/Utils/python/JetIdParams_cfi.py")))
                          );
                      
  fDiphotonEvent = new PhotonTreeWriterDiphotonEvent;
  fSinglePhoton  = new PhotonTreeWriterPhoton<16>;
  
  TFile *ftmp = TFile::Open(TString::Format("%s_tmp.root",GetName()),"RECREATE");
  
  if (fWriteDiphotonTree) {
    hCiCTuple = new TTree(fTupleName.Data(),fTupleName.Data());
    hCiCTuple->SetAutoSave(300e9);
  }
  TString singlename = fTupleName + TString("Single");
  if (fWriteSingleTree) {
    hCiCTupleSingle = new TTree(singlename,singlename);
    hCiCTupleSingle->SetAutoSave(300e9);
  }
    
  //make flattish tree from classes so we don't have to rely on dictionaries for reading later
  TClass *eclass = TClass::GetClass("mithep::PhotonTreeWriterDiphotonEvent");
  TClass *pclass = TClass::GetClass("mithep::PhotonTreeWriterPhoton<16>");
  TList  *elist  = eclass->GetListOfDataMembers();
  TList  *plist  = pclass->GetListOfDataMembers();
    
  for (int i=0; i<elist->GetEntries(); ++i) {
    const TDataMember *tdm = static_cast<const TDataMember*>(elist->At(i));
    if (!(tdm->IsBasic() && tdm->IsPersistent())) continue;
    TString typestring;
    if (TString(tdm->GetTypeName()).BeginsWith("Char_t")) typestring = "B";
    else if (TString(tdm->GetTypeName()).BeginsWith("UChar_t")) typestring = "b";
    else if (TString(tdm->GetTypeName()).BeginsWith("Short_t")) typestring = "S";
    else if (TString(tdm->GetTypeName()).BeginsWith("UShort_t")) typestring = "s";
    else if (TString(tdm->GetTypeName()).BeginsWith("Int_t")) typestring = "I";
    else if (TString(tdm->GetTypeName()).BeginsWith("UInt_t")) typestring = "i";
    else if (TString(tdm->GetTypeName()).BeginsWith("Float_t")) typestring = "F";
    else if (TString(tdm->GetTypeName()).BeginsWith("Double_t")) typestring = "D";
    else if (TString(tdm->GetTypeName()).BeginsWith("Long64_t")) typestring = "L";
    else if (TString(tdm->GetTypeName()).BeginsWith("ULong64_t")) typestring = "l";
    else if (TString(tdm->GetTypeName()).BeginsWith("Bool_t")) typestring = "O";
    else continue;
    //printf("%s %s: %i\n",tdm->GetTypeName(),tdm->GetName(),int(tdm->GetOffset()));
    Char_t *addr = (Char_t*)fDiphotonEvent;
    assert(sizeof(Char_t)==1);
    if (fWriteDiphotonTree) hCiCTuple->Branch(tdm->GetName(),addr + tdm->GetOffset(),TString::Format("%s/%s",tdm->GetName(),typestring.Data()));
    if (fWriteSingleTree) hCiCTupleSingle->Branch(tdm->GetName(),addr + tdm->GetOffset(),TString::Format("%s/%s",tdm->GetName(),typestring.Data()));
  }

  for (int iph=0; iph<2; ++iph) {
    for (int i=0; i<plist->GetEntries(); ++i) {
      const TDataMember *tdm = static_cast<const TDataMember*>(plist->At(i));
      if (!(tdm->IsBasic() && tdm->IsPersistent())) continue;
      TString typestring;
      if (TString(tdm->GetTypeName()).BeginsWith("Char_t")) typestring = "B";
      else if (TString(tdm->GetTypeName()).BeginsWith("UChar_t")) typestring = "b";
      else if (TString(tdm->GetTypeName()).BeginsWith("Short_t")) typestring = "S";
      else if (TString(tdm->GetTypeName()).BeginsWith("UShort_t")) typestring = "s";
      else if (TString(tdm->GetTypeName()).BeginsWith("Int_t")) typestring = "I";
      else if (TString(tdm->GetTypeName()).BeginsWith("UInt_t")) typestring = "i";
      else if (TString(tdm->GetTypeName()).BeginsWith("Float_t")) typestring = "F";
      else if (TString(tdm->GetTypeName()).BeginsWith("Double_t")) typestring = "D";
      else if (TString(tdm->GetTypeName()).BeginsWith("Long64_t")) typestring = "L";
      else if (TString(tdm->GetTypeName()).BeginsWith("ULong64_t")) typestring = "l";
      else if (TString(tdm->GetTypeName()).BeginsWith("Bool_t")) typestring = "O";
      else continue;
      //printf("%s\n",tdm->GetTypeName());
      TString varname = TString::Format("ph%d.%s",iph+1,tdm->GetName());
      if (tdm->GetArrayDim()==1) {
        varname = TString::Format("%s[%i]",varname.Data(),tdm->GetMaxIndex(0));
      }
      
      //printf("typename = %s, arraydim = %i, arraysize = %i,varname = %s\n", tdm->GetTypeName(), tdm->GetArrayDim(), tdm->GetMaxIndex(0), varname.Data());
      
      Char_t *addr = (Char_t*)&fDiphotonEvent->photons[iph];
      assert(sizeof(Char_t)==1);
      if (fWriteDiphotonTree) hCiCTuple->Branch(varname,addr+tdm->GetOffset(),TString::Format("%s/%s",varname.Data(),typestring.Data()));
      
      if (iph==0) {
        TString singlename = TString::Format("ph.%s",tdm->GetName());
        if (tdm->GetArrayDim()==1) {
          singlename = TString::Format("%s[%i]",singlename.Data(),tdm->GetMaxIndex(0));
        }       
        Char_t *addrsingle = (Char_t*)fSinglePhoton;
        if (fWriteSingleTree) hCiCTupleSingle->Branch(singlename,addrsingle+tdm->GetOffset(),TString::Format("%s/%s",singlename.Data(),typestring.Data()));
      }
    }
  }
  
  if (fWriteDiphotonTree)
    AddOutput(hCiCTuple);
  if (fWriteSingleTree)
    AddOutput(hCiCTupleSingle);
  
  if (fFillVertexTree) {
    fDiphotonVtx = new PhotonTreeWriterVtx;
    hVtxTree = new TTree("hVtxTree","hVtxTree");
    hVtxTree->SetAutoSave(300e9);
    AddOutput(hVtxTree);
    
    TClass *vclass = TClass::GetClass("mithep::PhotonTreeWriterVtx");
    TList  *vlist  = vclass->GetListOfDataMembers();

    for (int i=0; i<vlist->GetEntries(); ++i) {
      const TDataMember *tdm = static_cast<const TDataMember*>(vlist->At(i));
      if (!(tdm->IsBasic() && tdm->IsPersistent())) continue;
      TString typestring;
      if (TString(tdm->GetTypeName()).BeginsWith("Char_t")) typestring = "B";
      else if (TString(tdm->GetTypeName()).BeginsWith("UChar_t")) typestring = "b";
      else if (TString(tdm->GetTypeName()).BeginsWith("Short_t")) typestring = "S";
      else if (TString(tdm->GetTypeName()).BeginsWith("UShort_t")) typestring = "s";
      else if (TString(tdm->GetTypeName()).BeginsWith("Int_t")) typestring = "I";
      else if (TString(tdm->GetTypeName()).BeginsWith("UInt_t")) typestring = "i";
      else if (TString(tdm->GetTypeName()).BeginsWith("Float_t")) typestring = "F";
      else if (TString(tdm->GetTypeName()).BeginsWith("Double_t")) typestring = "D";
      else if (TString(tdm->GetTypeName()).BeginsWith("Long64_t")) typestring = "L";
      else if (TString(tdm->GetTypeName()).BeginsWith("ULong64_t")) typestring = "l";
      else if (TString(tdm->GetTypeName()).BeginsWith("Bool_t")) typestring = "O";
      else continue;
      //printf("%s %s: %i\n",tdm->GetTypeName(),tdm->GetName(),int(tdm->GetOffset()));
      Char_t *addr = (Char_t*)fDiphotonVtx;
      assert(sizeof(Char_t)==1);
      hVtxTree->Branch(tdm->GetName(),addr + tdm->GetOffset(),TString::Format("%s/%s",tdm->GetName(),typestring.Data()));
    }
  }
  
}

// ----------------------------------------------------------------------------------------
// some helpfer functions....
void PhotonTreeWriter::FindHiggsPtAndZ(Float_t& pt, Float_t& decayZ, Float_t& mass)
{
  pt     = -999.;
  decayZ = -999.;
  mass   = -999.;

  // loop over all GEN particles and look for status 1 photons
  for(UInt_t i=0; i<fMCParticles->GetEntries(); ++i) {
    const MCParticle* p = fMCParticles->At(i);
    if (p->Is(MCParticle::kH) || (!fApplyElectronVeto &&
				  (p->AbsPdgId()==23 || p->AbsPdgId()==24))) {
      pt     = p->Pt();
      decayZ = p->DecayVertex().Z();
      mass   = p->Mass();
      break;
    }
  }
  
  return;
 }


Float_t PhotonTreeWriter::GetEventCat(PhotonTools::CiCBaseLineCats cat1,
                                      PhotonTools::CiCBaseLineCats cat2) {
  
  bool ph1IsEB  = (cat1 ==  PhotonTools::kCiCCat1 || cat1 == PhotonTools::kCiCCat2);
  bool ph2IsEB  = (cat2 ==  PhotonTools::kCiCCat1 || cat2 == PhotonTools::kCiCCat2);

  bool ph1IsHR9 = (cat1 ==  PhotonTools::kCiCCat1 || cat1 == PhotonTools::kCiCCat3);
  bool ph2IsHR9 = (cat2 ==  PhotonTools::kCiCCat1 || cat2 == PhotonTools::kCiCCat3);
  
  if( ph1IsEB && ph2IsEB )
    return ( ph1IsHR9 && ph2IsHR9 ? 0. : 1.);
  
  return ( ph1IsHR9 && ph2IsHR9 ? 2. : 3.);
}

template <int NClus>
void PhotonTreeWriterPhoton<NClus>::SetVars(const Photon *p, const DecayParticle *c, const Electron *ele,
				     const SuperCluster *pfsc, const MCParticle *m,
				     PhotonFix &phfixph, PhotonFix &phfixele,
				     const TrackCol* trackCol,const VertexCol* vtxCol,Double_t rho,
				     Bool_t fillclusterarrays, 
				     const ElectronCol* els, Bool_t applyElectronVeto) {
  
  const SuperCluster *s = 0;
  if (p)
    s = p->SCluster();
  else
    s = ele->SCluster();
  const BasicCluster *b = s->Seed();
  const BasicCluster *b2 = 0;
  Double_t ebcmax = -99.;
  for (UInt_t i=0; i<s->ClusterSize(); ++i) {
    const BasicCluster *bc = s->Cluster(i);
    if (bc->Energy() > ebcmax && bc !=b) {
      b2 = bc;
      ebcmax = bc->Energy();
    }
  }
  
  const BasicCluster *bclast = 0;
  Double_t ebcmin = 1e6;
  for (UInt_t i=0; i<s->ClusterSize(); ++i) {
    const BasicCluster *bc = s->Cluster(i);
    if (bc->Energy() < ebcmin && bc !=b) {
      bclast = bc;
      ebcmin = bc->Energy();
    }
  }

  const BasicCluster *bclast2 = 0;
  ebcmin = 1e6;
  for (UInt_t i=0; i<s->ClusterSize(); ++i) {
    const BasicCluster *bc = s->Cluster(i);
    if (bc->Energy() < ebcmin && bc !=b && bc!=bclast) {
      bclast2 = bc;
      ebcmin = bc->Energy();
    }
  }
  
  if (p) {
    hasphoton = kTRUE;
    e = p->E();
    pt = p->Pt();
    eta = p->Eta();
    phi = p->Phi();
    r9 = p->R9();
    e3x3 = p->E33();
    e5x5 = p->E55();
    hovere = p->HadOverEm();
    hoveretower = p->HadOverEmTow();
    sigietaieta = p->CoviEtaiEta();      
    phcat = PhotonTools::CiCBaseLineCat(p);
    eerr = p->EnergyErr();
    eerrsmeared = p->EnergyErrSmeared();
    esmearing = p->EnergySmearing();
    idmva = p->IdMva();
    hcalisodr03 = p->HcalTowerSumEtDr03();
    ecalisodr03 = p->EcalRecHitIsoDr03();        
    trkisohollowdr03 = p->HollowConeTrkIsoDr03();
    hcalisodr04 = p->HcalTowerSumEtDr04();
    ecalisodr04 = p->EcalRecHitIsoDr04();        
    trkisohollowdr04 = p->HollowConeTrkIsoDr04();
    
    
    const Vertex *vtx = vtxCol->At(0);
    if (p->HasPV()) vtx = p->PV();
    
    UInt_t wVtxInd = 0;
    
    trackiso1 = IsolationTools::CiCTrackIsolation(p,vtx, 0.3, 0.02, 0.0, 0.0, 0.1, 1.0,
  						      trackCol, NULL, NULL,
  						      (!applyElectronVeto ? els : NULL) );
  	//Question Ming:whyfPV->At(0) instead of selected vertex using ranking method?
      
    // track iso worst vtx
    trackiso2 = IsolationTools::CiCTrackIsolation(p,vtx, 0.4, 0.02, 0.0, 0.0, 0.1, 1.0,
  						      trackCol, &wVtxInd,vtxCol,
  						      (!applyElectronVeto ? els : NULL) );
    combiso1 = ecalisodr03+hcalisodr04+trackiso1 - 0.17*rho;
    combiso2 = ecalisodr04+hcalisodr04+trackiso2 - 0.52*rho;        
  }
  else {
    hasphoton = kFALSE;
    e = -99.;
    pt = -99.;
    eta = -99.;
    phi = -99.;
    r9 = b->E3x3()/s->RawEnergy();
    e3x3 = b->E3x3();
    e5x5 = b->E5x5();
    hovere = ele->HadronicOverEm();
    sigietaieta = ele->CoviEtaiEta();      
    phcat = -99;    
    eerr = -99.;
    eerrsmeared = -99.;
    esmearing = 0.;
    idmva = -99.;
    hcalisodr03 = -99;
    ecalisodr03 = -99;
    trkisohollowdr03 = -99;
    hcalisodr04 = -99;
    ecalisodr04 = -99;
    trkisohollowdr04 = -99;        
    trackiso1 = -99.;
    trackiso2 = -99.;
    combiso1 = -99.;
    combiso2 = -99.;
  }
  
  sce = s->Energy();
  scrawe = s->RawEnergy();
  scpse = s->PreshowerEnergy();
  scpssigmaxx = s->PsEffWidthSigmaXX();
  scpssigmayy = s->PsEffWidthSigmaYY();
  sceta = s->Eta();
  scphi = s->Phi();
  scnclusters = s->ClusterSize();
  scnhits = s->NHits();
  scetawidth = -99.;
  scphiwidth = -99.;
  if (p) {
    scetawidth = p->EtaWidth();
    scphiwidth = p->PhiWidth();
  }
  else {
    scetawidth = s->EtaWidth();
    scphiwidth = s->PhiWidth();        
  }
  isbarrel = (s->AbsEta()<1.5);
  isr9reco = (isbarrel && r9>0.94) || (!isbarrel && r9>0.95);
  isr9cat = (r9>0.94);
  
  eseed = b->Energy();      
  etaseed = b->Eta();
  phiseed = b->Phi();
  ietaseed = b->IEta();
  iphiseed = b->IPhi();
  ixseed = b->IX();
  iyseed = b->IY();
  etacryseed = b->EtaCry();
  phicryseed = b->PhiCry();
  xcryseed = b->XCry();
  ycryseed = b->YCry();
  thetaaxisseed = b->ThetaAxis();
  phiaxisseed = b->PhiAxis();
  sigietaietaseed = TMath::Sqrt(b->CoviEtaiEta());
  sigiphiphiseed = TMath::Sqrt(b->CoviPhiiPhi());
  if (isnan(sigiphiphiseed)) sigiphiphiseed = -99.;
  covietaiphiseed = b->CoviEtaiPhi();
  if (isnan(covietaiphiseed)) covietaiphiseed = -99.;
  e3x3seed = b->E3x3();
  e5x5seed = b->E5x5();
  emaxseed = b->EMax();
  e2ndseed = b->E2nd();
  etopseed = b->ETop();
  ebottomseed = b->EBottom();
  eleftseed = b->ELeft();
  erightseed = b->ERight();
  e1x3seed = b->E1x3();
  e3x1seed = b->E3x1();
  e1x5seed = b->E1x5();
  e2x2seed = b->E2x2();
  e4x4seed = b->E4x4();
  e2x5maxseed = b->E2x5Max();
  e2x5topseed = b->E2x5Top();
  e2x5bottomseed = b->E2x5Bottom();
  e2x5leftseed = b->E2x5Left();
  e2x5rightseed = b->E2x5Right();
  xseedseed = b->Pos().X();
  yseedseed = b->Pos().Y();
  zseedseed = b->Pos().Z();
  nhitsseed = b->NHits();

  if (b2) {
    ebc2 = b2->Energy();      
    etabc2 = b2->Eta();
    phibc2 = b2->Phi();
    ietabc2 = b2->IEta();
    iphibc2 = b2->IPhi();
    ixbc2 = b2->IX();
    iybc2 = b2->IY();
    etacrybc2 = b2->EtaCry();
    phicrybc2 = b2->PhiCry();
    xcrybc2 = b2->XCry();
    ycrybc2 = b2->YCry();
    thetaaxisbc2 = b2->ThetaAxis();
    phiaxisbc2 = b2->PhiAxis();
    sigietaietabc2 = TMath::Sqrt(b2->CoviEtaiEta());
    sigiphiphibc2 = TMath::Sqrt(b2->CoviPhiiPhi());
    if (isnan(sigiphiphibc2)) sigiphiphibc2 = -99.;
    covietaiphibc2 = b2->CoviEtaiPhi();
    if (isnan(covietaiphibc2)) covietaiphibc2 = -99.;
    e3x3bc2 = b2->E3x3();
    e5x5bc2 = b2->E5x5();
    emaxbc2 = b2->EMax();
    e2ndbc2 = b2->E2nd();
    etopbc2 = b2->ETop();
    ebottombc2 = b2->EBottom();
    eleftbc2 = b2->ELeft();
    erightbc2 = b2->ERight();
    e1x3bc2 = b2->E1x3();
    e3x1bc2 = b2->E3x1();
    e1x5bc2 = b2->E1x5();
    e2x2bc2 = b2->E2x2();
    e4x4bc2 = b2->E4x4();
    e2x5maxbc2 = b2->E2x5Max();
    e2x5topbc2 = b2->E2x5Top();
    e2x5bottombc2 = b2->E2x5Bottom();
    e2x5leftbc2 = b2->E2x5Left();
    e2x5rightbc2 = b2->E2x5Right();
    xbc2bc2 = b2->Pos().X();
    ybc2bc2 = b2->Pos().Y();
    zbc2bc2 = b2->Pos().Z();      
    nhitsbc2= b2->NHits();
  }
  else {
    ebc2 = 0;
    etabc2 = 0;
    phibc2 = 0;
    ietabc2 = 0;
    iphibc2 = 0;
    ixbc2 = 0;
    iybc2 = 0;
    etacrybc2 = 0;
    phicrybc2 = 0;
    xcrybc2 = 0;
    ycrybc2 = 0;
    thetaaxisbc2 = 0;
    phiaxisbc2 = 0;
    sigietaietabc2 = 0;
    sigiphiphibc2 = 0;
    covietaiphibc2 = 0;
    e3x3bc2 = 0;
    e5x5bc2 = 0;
    emaxbc2 = 0;
    e2ndbc2 = 0;
    etopbc2 = 0;
    ebottombc2 = 0;
    eleftbc2 = 0;
    erightbc2 = 0;
    e1x3bc2 = 0;
    e3x1bc2 = 0;
    e1x5bc2 = 0;
    e2x2bc2 = 0;
    e4x4bc2 = 0;
    e2x5maxbc2 = 0;
    e2x5topbc2 = 0;
    e2x5bottombc2 = 0;
    e2x5leftbc2 = 0;
    e2x5rightbc2 = 0;
    xbc2bc2 = 0;
    ybc2bc2 = 0;
    zbc2bc2 = 0;
    nhitsbc2 = 0;
  }

  if (bclast) {
    ebclast = bclast->Energy();      
    etabclast = bclast->Eta();
    phibclast = bclast->Phi();
    ietabclast = bclast->IEta();
    iphibclast = bclast->IPhi();
    ixbclast = bclast->IX();
    iybclast = bclast->IY();
    etacrybclast = bclast->EtaCry();
    phicrybclast = bclast->PhiCry();
    xcrybclast = bclast->XCry();
    ycrybclast = bclast->YCry();
    thetaaxisbclast = bclast->ThetaAxis();
    phiaxisbclast = bclast->PhiAxis();
    sigietaietabclast = TMath::Sqrt(bclast->CoviEtaiEta());
    sigiphiphibclast = TMath::Sqrt(bclast->CoviPhiiPhi());
    if (isnan(sigiphiphibclast)) sigiphiphibclast = -99.;
    covietaiphibclast = bclast->CoviEtaiPhi();
    if (isnan(covietaiphibclast)) covietaiphibclast = -99.;
    e3x3bclast = bclast->E3x3();
    nhitsbclast = bclast->NHits();
  }
  else {
    ebclast = 0;
    etabclast = 0;
    phibclast = 0;
    ietabclast = 0;
    iphibclast = 0;
    ixbclast = 0;
    iybclast = 0;
    etacrybclast = 0;
    phicrybclast = 0;
    xcrybclast = 0;
    ycrybclast = 0;
    thetaaxisbclast = 0;
    phiaxisbclast = 0;
    sigietaietabclast = 0;
    sigiphiphibclast = 0;
    covietaiphibclast = 0;
    e3x3bclast = 0;
    nhitsbclast = 0;
  }

  if (bclast2) {
    ebclast2 = bclast2->Energy();      
    etabclast2 = bclast2->Eta();
    phibclast2 = bclast2->Phi();
    ietabclast2 = bclast2->IEta();
    iphibclast2 = bclast2->IPhi();
    ixbclast2 = bclast2->IX();
    iybclast2 = bclast2->IY();
    etacrybclast2 = bclast2->EtaCry();
    phicrybclast2 = bclast2->PhiCry();
    xcrybclast2 = bclast2->XCry();
    ycrybclast2 = bclast2->YCry();
    thetaaxisbclast2 = bclast2->ThetaAxis();
    phiaxisbclast2 = bclast2->PhiAxis();
    sigietaietabclast2 = TMath::Sqrt(bclast2->CoviEtaiEta());
    sigiphiphibclast2 = TMath::Sqrt(bclast2->CoviPhiiPhi());
    if (isnan(sigiphiphibclast2)) sigiphiphibclast2 = -99.;
    covietaiphibclast2 = bclast2->CoviEtaiPhi();
    if (isnan(covietaiphibclast2)) covietaiphibclast2 = -99.;
    e3x3bclast2 = bclast2->E3x3();
    nhitsbclast2 = bclast2->NHits();
  }
  else {
    ebclast2 = 0;
    etabclast2 = 0;
    phibclast2 = 0;
    ietabclast2 = 0;
    iphibclast2 = 0;
    ixbclast2 = 0;
    iybclast2 = 0;
    etacrybclast2 = 0;
    phicrybclast2 = 0;
    xcrybclast2 = 0;
    ycrybclast2 = 0;
    thetaaxisbclast2 = 0;
    phiaxisbclast2 = 0;
    sigietaietabclast2 = 0;
    sigiphiphibclast2 = 0;
    covietaiphibclast2 = 0;
    e3x3bclast2 = 0;
    nhitsbclast2 = 0;
  }

  for (UInt_t iclus=0; iclus<NClus; ++iclus) {
    if (fillclusterarrays && iclus < s->ClusterSize() ) {
      const BasicCluster *ib =s->Cluster(iclus);
      
      ebcs[iclus] = ib->Energy();      
      etabcs[iclus] = ib->Eta();
      phibcs[iclus] = ib->Phi();
      ietabcs[iclus] = ib->IEta();
      iphibcs[iclus] = ib->IPhi();
      ixbcs[iclus] = ib->IX();
      iybcs[iclus] = ib->IY();
      etacrybcs[iclus] = ib->EtaCry();
      phicrybcs[iclus] = ib->PhiCry();
      xcrybcs[iclus] = ib->XCry();
      ycrybcs[iclus] = ib->YCry();
      sigietaietabcs[iclus] = TMath::Sqrt(ib->CoviEtaiEta());
      sigiphiphibcs[iclus] = TMath::Sqrt(ib->CoviPhiiPhi());
      covietaiphibcs[iclus] = ib->CoviEtaiPhi();
      sigetaetabcs[iclus] = TMath::Sqrt(ib->CovEtaEta());
      sigphiphibcs[iclus] = TMath::Sqrt(ib->CovPhiPhi());
      covetaphibcs[iclus] = ib->CovEtaPhi();
      e3x3bcs[iclus] = ib->E3x3();
      e5x5bcs[iclus] = ib->E5x5();
      emaxbcs[iclus] = ib->EMax();
      e2ndbcs[iclus] = ib->E2nd();
      etopbcs[iclus] = ib->ETop();
      ebottombcs[iclus] = ib->EBottom();
      eleftbcs[iclus] = ib->ELeft();
      erightbcs[iclus] = ib->ERight();
      e1x3bcs[iclus] = ib->E1x3();
      e3x1bcs[iclus] = ib->E3x1();
      e1x5bcs[iclus] = ib->E1x5();
      e2x2bcs[iclus] = ib->E2x2();
      e4x4bcs[iclus] = ib->E4x4();
      e2x5maxbcs[iclus] = ib->E2x5Max();
      e2x5topbcs[iclus] = ib->E2x5Top();
      e2x5bottombcs[iclus] = ib->E2x5Bottom();
      e2x5leftbcs[iclus] = ib->E2x5Left();
      e2x5rightbcs[iclus] = ib->E2x5Right();  
      nhitsbcs[iclus]= ib->NHits();
    }
    else {
      ebcs[iclus] = -999;
      etabcs[iclus] = -999;
      phibcs[iclus] = -999;
      ietabcs[iclus] = -999;
      iphibcs[iclus] = -999;
      ixbcs[iclus] = -999;
      iybcs[iclus] = -999;
      etacrybcs[iclus] = -999;
      phicrybcs[iclus] = -999;
      xcrybcs[iclus] = -999;
      ycrybcs[iclus] = -999;
      sigietaietabcs[iclus] = -999;
      sigiphiphibcs[iclus] = -999;
      covietaiphibcs[iclus] = -999;
      sigetaetabcs[iclus] = -999;
      sigphiphibcs[iclus] = -999;
      covetaphibcs[iclus] = -999;          
      e3x3bcs[iclus] = -999;
      e5x5bcs[iclus] = -999;
      emaxbcs[iclus] = -999;
      e2ndbcs[iclus] = -999;
      etopbcs[iclus] = -999;
      ebottombcs[iclus] = -999;
      eleftbcs[iclus] = -999;
      erightbcs[iclus] = -999;
      e1x3bcs[iclus] = -999;
      e3x1bcs[iclus] = -999;
      e1x5bcs[iclus] = -999;
      e2x2bcs[iclus] = -999;
      e4x4bcs[iclus] = -999;
      e2x5maxbcs[iclus] = -999;
      e2x5topbcs[iclus] = -999;
      e2x5bottombcs[iclus] = -999;
      e2x5leftbcs[iclus] = -999;
      e2x5rightbcs[iclus] = -999;
      nhitsbcs[iclus] = -999;
    }
  }

  for (UInt_t iclus=0; iclus<NClus; ++iclus) {
    if (fillclusterarrays && pfsc && iclus < pfsc->ClusterSize() ) {
      const BasicCluster *ib =pfsc->Cluster(iclus);
      
      epfbcs[iclus] = ib->Energy();      
      etapfbcs[iclus] = ib->Eta();
      phipfbcs[iclus] = ib->Phi();
      ietapfbcs[iclus] = ib->IEta();
      iphipfbcs[iclus] = ib->IPhi();
      ixpfbcs[iclus] = ib->IX();
      iypfbcs[iclus] = ib->IY();
      etacrypfbcs[iclus] = ib->EtaCry();
      phicrypfbcs[iclus] = ib->PhiCry();
      xcrypfbcs[iclus] = ib->XCry();
      ycrypfbcs[iclus] = ib->YCry();
      sigietaietapfbcs[iclus] = TMath::Sqrt(ib->CoviEtaiEta());
      sigiphiphipfbcs[iclus] = TMath::Sqrt(ib->CoviPhiiPhi());
      covietaiphipfbcs[iclus] = ib->CoviEtaiPhi();
      sigetaetapfbcs[iclus] = TMath::Sqrt(ib->CovEtaEta());
      sigphiphipfbcs[iclus] = TMath::Sqrt(ib->CovPhiPhi());
      covetaphipfbcs[iclus] = ib->CovEtaPhi();
      e3x3pfbcs[iclus] = ib->E3x3();
      e5x5pfbcs[iclus] = ib->E5x5();
      emaxpfbcs[iclus] = ib->EMax();
      e2ndpfbcs[iclus] = ib->E2nd();
      etoppfbcs[iclus] = ib->ETop();
      ebottompfbcs[iclus] = ib->EBottom();
      eleftpfbcs[iclus] = ib->ELeft();
      erightpfbcs[iclus] = ib->ERight();
      e1x3pfbcs[iclus] = ib->E1x3();
      e3x1pfbcs[iclus] = ib->E3x1();
      e1x5pfbcs[iclus] = ib->E1x5();
      e2x2pfbcs[iclus] = ib->E2x2();
      e4x4pfbcs[iclus] = ib->E4x4();
      e2x5maxpfbcs[iclus] = ib->E2x5Max();
      e2x5toppfbcs[iclus] = ib->E2x5Top();
      e2x5bottompfbcs[iclus] = ib->E2x5Bottom();
      e2x5leftpfbcs[iclus] = ib->E2x5Left();
      e2x5rightpfbcs[iclus] = ib->E2x5Right();  
      nhitspfbcs[iclus]= ib->NHits();
    }
    else {
      epfbcs[iclus] = -999;
      etapfbcs[iclus] = -999;
      phipfbcs[iclus] = -999;
      ietapfbcs[iclus] = -999;
      iphipfbcs[iclus] = -999;
      ixpfbcs[iclus] = -999;
      iypfbcs[iclus] = -999;
      etacrypfbcs[iclus] = -999;
      phicrypfbcs[iclus] = -999;
      xcrypfbcs[iclus] = -999;
      ycrypfbcs[iclus] = -999;
      sigietaietapfbcs[iclus] = -999;
      sigiphiphipfbcs[iclus] = -999;
      covietaiphipfbcs[iclus] = -999;
      sigetaetapfbcs[iclus] = -999;
      sigphiphipfbcs[iclus] = -999;
      covetaphipfbcs[iclus] = -999;          
      e3x3pfbcs[iclus] = -999;
      e5x5pfbcs[iclus] = -999;
      emaxpfbcs[iclus] = -999;
      e2ndpfbcs[iclus] = -999;
      etoppfbcs[iclus] = -999;
      ebottompfbcs[iclus] = -999;
      eleftpfbcs[iclus] = -999;
      erightpfbcs[iclus] = -999;
      e1x3pfbcs[iclus] = -999;
      e3x1pfbcs[iclus] = -999;
      e1x5pfbcs[iclus] = -999;
      e2x2pfbcs[iclus] = -999;
      e4x4pfbcs[iclus] = -999;
      e2x5maxpfbcs[iclus] = -999;
      e2x5toppfbcs[iclus] = -999;
      e2x5bottompfbcs[iclus] = -999;
      e2x5leftpfbcs[iclus] = -999;
      e2x5rightpfbcs[iclus] = -999;
      nhitspfbcs[iclus] = -999;
    }
  }

  for (UInt_t iclus=0; iclus<100; ++iclus) {
    if (fillclusterarrays && pfsc && iclus < pfsc->NPsClusts() ) {
      const PsCluster *ib = pfsc->PsClust(iclus);
      
      epsc[iclus] = ib->Energy();
      etapsc[iclus] = ib->Eta();
      phipsc[iclus] = ib->Phi();
      planepsc[iclus] = ib->PsPlane();
    }
    else {
      epsc[iclus] = -999;
      etapsc[iclus] = -999;
      phipsc[iclus] = -999;
      planepsc[iclus] = 0;      
    }
  }

  //initialize photon energy corrections if needed
  /*if (!PhotonFix::initialised()) {
    PhotonFix::initialise("4_2",std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFix.dat")).Data()));  
  }*/ 

  phfixph.setup(e,sceta,scphi,r9);
  phfixele.setup(e,sceta,scphi,r9);      
  
  const Float_t dval = -99.;
  ecor = phfixph.fixedEnergy();
  ecorerr = phfixph.sigmaEnergy();
  ecorele = phfixele.fixedEnergy();
  ecoreleerr = phfixele.sigmaEnergy();
  if (phfixph.isbarrel()) {
   etac = phfixph.etaC();
   etas = phfixph.etaS();
   etam = phfixph.etaM();
   phic = phfixph.phiC();
   phis = phfixph.phiS();
   phim = phfixph.phiM();
   xz = dval;
   xc = dval;
   xs = dval;
   xm = dval;
   yz = dval;
   yc = dval;
   ys = dval;
   ym = dval;
  }
  else {
   etac = dval;
   etas = dval;
   etam = dval;
   phic = dval;
   phis = dval;
   phim = dval;
   xz = phfixph.xZ();
   xc = phfixph.xC();
   xs = phfixph.xS();
   xm = phfixph.xM();
   yz = phfixph.yZ();
   yc = phfixph.yC();
   ys = phfixph.yS();
   ym = phfixph.yM();
  }   
  
  if (c) {
    hasconversion = kTRUE;
    convp = c->P();
    convpt = c->Pt();
    conveta = c->Eta();
    convphi = c->Phi();
    ThreeVector dirconvsc = ThreeVector(s->Point()) - c->Position();
    convdeta = c->Eta() - dirconvsc.Eta();
    convdphi = MathUtils::DeltaPhi(c->Phi(),dirconvsc.Phi());
    convvtxrho = c->Position().Rho();
    convvtxz =   c->Position().Z();
    convvtxphi = c->Position().Phi();        
   
    const StableData *leadsd = dynamic_cast<const StableData*>(c->DaughterDat(0));
    const StableData *trailsd = dynamic_cast<const StableData*>(c->DaughterDat(1));
    if (leadsd->Pt()<trailsd->Pt()) {
      const StableData *sdtmp = leadsd;
      leadsd = trailsd;
      trailsd = sdtmp;
    }
    
    const Track *leadtrack = dynamic_cast<const StableParticle*>(leadsd->Original())->Trk();
    const Track *trailtrack = dynamic_cast<const StableParticle*>(trailsd->Original())->Trk();
    
    convleadpt = leadsd->Pt();
    convtrailpt = trailsd->Pt();
    convleadtrackpt = leadtrack->Pt();
    convleadtrackalgo = leadtrack->Algo();
    if (convleadtrackalgo==29) convleadtrackalgos=2; //gsf track
    else if (convleadtrackalgo==15 ||convleadtrackalgo==16) convleadtrackalgos=1; //ecal-seeded track
    else convleadtrackalgos = 0; //std iterative track
    convleadtrackcharge = leadtrack->Charge();
    convtrailtrackpt = trailtrack->Pt();
    convtrailtrackalgo = trailtrack->Algo();
    if (convtrailtrackalgo==29) convtrailtrackalgos=2; //gsf track
    else if (convtrailtrackalgo==15 ||convtrailtrackalgo==16) convtrailtrackalgos=1; //ecal-seeded track
    else convtrailtrackalgos = 0; //std iterative track
    trailtrackcharge = trailtrack->Charge();      
  }
  else {
    hasconversion = kFALSE;
    convp = -99.;
    convpt = -99.;
    conveta = -99.;
    convphi = -99.;
    convdeta = -99.;
    convdphi = -99.;
    convvtxrho = -99.;
    convvtxz = -999.;
    convvtxphi = -99.;        
    convleadpt = -99.;
    convtrailpt = -99.;
    convleadtrackpt = -99.;
    convleadtrackalgo = -99;
    convleadtrackalgos = -99;
    convleadtrackcharge = 0;
    convtrailtrackpt = -99.;
    convtrailtrackalgo = -99;
    convtrailtrackalgos = -99;      
    trailtrackcharge = 0;      
  }
  
  //electron quantities
  if (ele) {
    haselectron = kTRUE;
    eleisecaldriven = ele->IsEcalDriven();
    eleistrackerdriven = ele->IsTrackerDriven();
    elee = ele->E();
    elept = ele->Pt();
    eleeta = ele->Eta();
    elephi = ele->Phi();
    elecharge = ele->Charge();
    elefbrem = ele->FBrem();
    eledeta = ele->DeltaEtaSuperClusterTrackAtVtx();
    eledphi = ele->DeltaPhiSuperClusterTrackAtVtx();
    elep = s->Energy()/ele->ESuperClusterOverP();
    elepin = ele->PIn();
    elepout = ele->POut();        
  }
  else {
    haselectron = kFALSE;
    eleisecaldriven = kFALSE;
    eleistrackerdriven = kFALSE;
    elee = -99.;
    elept = -99.;
    eleeta = -99.;
    elephi = -99.;
    elecharge = -99;
    elefbrem = -99.;
    eledeta = -99.;
    eledphi = -99.;
    elep = -99.;
    elepin = -99.;
    elepout = -99.;
  }
  
  //pf supercluster quantities
  if (pfsc) {
    haspfsc = kTRUE;
    pfsce = pfsc->Energy();
    pfscrawe = pfsc->RawEnergy();
    pfsceta = pfsc->Eta();
    pfscphi = pfsc->Phi();        
    pfscnclusters = pfsc->NClusters();
    pfscnhits = pfsc->NHits();
    pfscetawidth = pfsc->EtaWidth();
    pfscphiwidth = pfsc->PhiWidth();
    pfscnpsclusters = pfsc->NPsClusts();   
  }
  else {
    haspfsc = kFALSE;
    pfsce = -99.;
    pfscrawe = -99.;
    pfsceta = -99.;
    pfscphi = -99.;
    pfscnclusters = 0;
    pfscnhits = 0;
    pfscetawidth = -99.;
    pfscphiwidth = -99.;
    pfscnpsclusters = 0;
  }
  
  genz = -99.;
  if (m) {
    ispromptgen = kTRUE;
    gene = m->E();
    genpt = m->Pt();
    geneta = m->Eta();
    genphi = m->Phi();
    const MCParticle *mm = m->DistinctMother();
    if (mm) genz = mm->DecayVertex().Z();
    pdgid = m->PdgId();
    if (mm) motherpdgid = mm->PdgId();
    else motherpdgid = -99;
  }
  else {
    ispromptgen = kFALSE;
    gene = -99.;
    genpt = -99.;
    geneta = -99.;
    genphi = -99.;
    pdgid = -99;
    motherpdgid = -99;
  }
}

void PhotonTreeWriterVtx::SetVars(const Vertex *v, const Photon *p1, const Photon *p2, const PFCandidateCol *pfcands, Int_t idx, Int_t numvtx, Float_t genvtxz)
{
 
  //printf("start\n");
  
  n = idx;
  nvtx = numvtx;
  zgen = genvtxz;
  
  x = v->X();
  y = v->Y();
  z = v->Z();
  
  Double_t dsumpt = 0.;
  Double_t dsumptsq = 0.;
  
  nchalltoward = 0;
  nchalltransverse = 0;
  nchallaway = 0;
  nchcuttoward = 0;
  nchcuttransverse = 0;
  nchcutaway = 0;  
  
  ThreeVector vtxmom;
  
  //printf("mom\n");
  FourVectorM diphoton = p1->MomVtx(v->Position()) + p2->MomVtx(v->Position());
  //printf("done mom\n");
  ptgg = diphoton.Pt();
  phigg = diphoton.Phi();
  etagg = diphoton.Eta();
  mgg = diphoton.M();
  pxgg = diphoton.Px();
  pygg = diphoton.Py();
  pzgg = diphoton.Pz();
  
  //printf("loop\n");
  
  for (UInt_t i = 0; i<pfcands->GetEntries(); ++i) {
    const PFCandidate *pfc = pfcands->At(i);
    if (pfc->PFType()!=PFCandidate::eHadron || !pfc->HasTrackerTrk()) continue;
    if (TMath::Abs( pfc->TrackerTrk()->DzCorrected(*v) ) > 0.2) continue;
    if (TMath::Abs( pfc->TrackerTrk()->D0Corrected(*v) ) > 0.1) continue;
    
    vtxmom += ThreeVector(pfc->Px(),pfc->Py(),pfc->Pz());
    
    dsumpt += pfc->Pt();
    dsumptsq += pfc->Pt()*pfc->Pt();
    
    Double_t dphi = TMath::Abs(MathUtils::DeltaPhi(*pfc,diphoton));
    if (dphi<(TMath::Pi()/3.0)) {
      ++nchalltoward;
      if (pfc->Pt()>0.5) ++nchcuttoward;
    }
    else if (dphi>(2.0*TMath::Pi()/3.0)) {
      ++nchallaway;
      if (pfc->Pt()>0.5) ++nchcutaway;      
    }
    else {
      ++nchalltransverse;
      if (pfc->Pt()>0.5) ++nchcuttransverse;            
    }
     
  }
  
  //printf("doneloop\n");

    
  sumpt = dsumpt;
  sumptsq = dsumptsq;
  
  pt = vtxmom.Rho();
  phi = vtxmom.Phi();
  eta = vtxmom.Eta();
  px = vtxmom.X();
  py = vtxmom.Y();
  pz = vtxmom.Z();
  
  //printf("done\n");
  
  return;
  
}
