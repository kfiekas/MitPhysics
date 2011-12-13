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
#include "TDataMember.h"
#include <TNtuple.h>
#include <TRandom3.h>
#include <TSystem.h>

using namespace mithep;

ClassImp(mithep::PhotonTreeWriter)
ClassImp(mithep::PhotonTreeWriterPhoton)
ClassImp(mithep::PhotonTreeWriterDiphotonEvent)


//--------------------------------------------------------------------------------------------------
PhotonTreeWriter::PhotonTreeWriter(const char *name, const char *title) : 
  // Base Module...
  BaseMod            (name,title),

  // define all the Branches to load
  fPhotonBranchName  (Names::gkPhotonBrn),
  fElectronName      (Names::gkElectronBrn),
  fGoodElectronName  (Names::gkElectronBrn),  
  fConversionName    (Names::gkMvfConversionBrn),  
  fTrackBranchName   (Names::gkTrackBrn),
  fPileUpDenName     (Names::gkPileupEnergyDensityBrn),
  fPVName            (Names::gkPVBeamSpotBrn),
  fBeamspotName      (Names::gkBeamSpotBrn),
  fPFCandName        (Names::gkPFCandidatesBrn),
  fMCParticleName    (Names::gkMCPartBrn),
  fPileUpName        (Names::gkPileupInfoBrn),  
  fSuperClusterName  ("PFSuperClusters"),
  fPFMetName         ("PFMet"),
  fPFJetName         (Names::gkPFJetBrn),

  
  fIsData            (false),
  fPhotonsFromBranch (true),  
  fPVFromBranch      (true),
  fGoodElectronsFromBranch (kTRUE),
  fPFJetsFromBranch  (kTRUE),

  // ----------------------------------------
  // collections....
  fPhotons           (0),
  fElectrons         (0),
  fConversions       (0),
  fTracks            (0),
  fPileUpDen         (0),
  fPV                (0),
  fBeamspot          (0),
  fPFCands           (0),
  fMCParticles       (0),
  fPileUp            (0),
  fSuperClusters     (0),
  fPFJets            (0),

  fLoopOnGoodElectrons(kFALSE),
  fApplyElectronVeto(kTRUE),  

  fWriteDiphotonTree(kTRUE),
  fWriteSingleTree(kTRUE),
  fExcludeSinglePrompt(kFALSE),
  fExcludeDoublePrompt(kFALSE),
  fEnableJets(kFALSE),
  fPhFixDataFile(gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFixSTART42V13.dat")),
  fTupleName         ("hPhotonTree")


{
  // Constructor.
}

PhotonTreeWriter::~PhotonTreeWriter(){
;
}

//--------------------------------------------------------------------------------------------------
void PhotonTreeWriter::Process()
{
  // ------------------------------------------------------------  
  // Process entries of the tree. 
  
  LoadEventObject(fPhotonBranchName,   fPhotons);
  LoadEventObject(fGoodElectronName,       fGoodElectrons);
  
  const BaseCollection *egcol = 0;
  if (fLoopOnGoodElectrons) {
    egcol = fGoodElectrons;
  }
  else {
    egcol = fPhotons;
  }
  
  if (egcol->GetEntries()<1) return;
  
  LoadEventObject(fElectronName,       fElectrons);
  LoadEventObject(fConversionName,     fConversions);
  LoadEventObject(fTrackBranchName,    fTracks);
  LoadEventObject(fPileUpDenName,      fPileUpDen);
  LoadEventObject(fPVName,             fPV);    
  LoadEventObject(fBeamspotName,       fBeamspot);
  LoadEventObject(fPFCandName,         fPFCands);
  LoadEventObject(fSuperClusterName,   fSuperClusters);
  LoadEventObject(fPFMetName,   fPFMet);  
  if (fEnableJets) LoadEventObject(fPFJetName,   fPFJets);  

  // ------------------------------------------------------------  
  // load event based information
  Int_t _numPU = -1.;        // some sensible default values....
  Int_t _numPUminus = -1.;        // some sensible default values....
  Int_t _numPUplus = -1.;        // some sensible default values....

  Float_t _tRho  = -99.;
  if( fPileUpDen->GetEntries() > 0 )
    _tRho  = (Double_t) fPileUpDen->At(0)->RhoRandomLowEta();
  
  const BaseVertex *bsp = dynamic_cast<const BaseVertex*>(fBeamspot->At(0));
  
  if( !fIsData ) {
    LoadBranch(fMCParticleName);
    LoadBranch(fPileUpName);
  }  
  
  if( !fIsData ) {
    for (UInt_t i=0; i<fPileUp->GetEntries(); ++i) {
      const PileupInfo *puinfo = fPileUp->At(i);
      if (puinfo->GetBunchCrossing()==0) _numPU = puinfo->GetPU_NumInteractions();
      else if (puinfo->GetBunchCrossing() == -1) _numPUminus = puinfo->GetPU_NumInteractions();
      else if (puinfo->GetBunchCrossing() ==  1) _numPUplus  = puinfo->GetPU_NumInteractions();
    }
  }

  // in case of a MC event, try to find Higgs and Higgs decay Z poisition
  Float_t _pth    = -100.;
  Float_t _decayZ = -100.;
  Float_t _genmass = -100.;
  if( !fIsData ) FindHiggsPtAndZ(_pth, _decayZ, _genmass);

  fDiphotonEvent->rho = _tRho;
  fDiphotonEvent->genHiggspt = _pth;
  fDiphotonEvent->genHiggsZ = _decayZ;
  fDiphotonEvent->genmass = _genmass;  
  fDiphotonEvent->gencostheta = -99.;
  fDiphotonEvent->nVtx = fPV->GetEntries();
  fDiphotonEvent->bsX = fBeamspot->At(0)->X();
  fDiphotonEvent->bsY = fBeamspot->At(0)->Y();
  fDiphotonEvent->bsZ = fBeamspot->At(0)->Z();
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
  
  Int_t nhitsbeforevtxmax = 1;
  if (!fApplyElectronVeto) nhitsbeforevtxmax = 999;  
  
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
    
    const DecayParticle *conv1 = PhotonTools::MatchedConversion(sc1,fConversions,bsp,nhitsbeforevtxmax);
    const DecayParticle *conv2 = PhotonTools::MatchedConversion(sc2,fConversions,bsp,nhitsbeforevtxmax);
    
    const SuperCluster *pfsc1 = PhotonTools::MatchedSC(sc1,fSuperClusters);
    const SuperCluster *pfsc2 = PhotonTools::MatchedSC(sc2,fSuperClusters);
    
    const MCParticle *phgen1 = 0;
    const MCParticle *phgen2 = 0;
    if( !fIsData ) {
      phgen1 = PhotonTools::MatchMC(p1,fMCParticles,!fApplyElectronVeto);
      phgen2 = PhotonTools::MatchMC(p2,fMCParticles,!fApplyElectronVeto);
    }
    
/*    if (phgen1 && phgen2) {
      printf("p1     pt = %5f, eta = %5f, phi = %5f\n",p1->Pt(),p1->Eta(),p1->Phi());
      printf("p2     pt = %5f, eta = %5f, phi = %5f\n",p2->Pt(),p2->Eta(),p2->Phi());
      printf("phgen1 pt = %5f, eta = %5f, phi = %5f, pdg = %i, motherpdg = %i, distinctmotherpdg = %i, distinctmotherstatus = %i\n",phgen1->Pt(),phgen1->Eta(),phgen1->Phi(),phgen1->PdgId(),phgen1->Mother()->PdgId(),phgen1->DistinctMother()->PdgId(),phgen1->DistinctMother()->Status());
      printf("phgen2 pt = %5f, eta = %5f, phi = %5f, pdg = %i, motherpdg = %i, distinctmotherpdg = %i, distinctmotherstatus = %i\n",phgen2->Pt(),phgen2->Eta(),phgen2->Phi(),phgen2->PdgId(),phgen2->Mother()->PdgId(),phgen2->DistinctMother()->PdgId(),phgen2->DistinctMother()->Status());
    }  */  
    
    
    if (fExcludeSinglePrompt && (phgen1 || phgen2) ) return;
    if (fExcludeDoublePrompt && (phgen1 && phgen2) ) return;
    
    if (!fLoopOnGoodElectrons && phHard->HasPV()) {
      fDiphotonEvent->vtxX = phHard->PV()->X();
      fDiphotonEvent->vtxY = phHard->PV()->Y();
      fDiphotonEvent->vtxZ = phHard->PV()->Z();
      fDiphotonEvent->vtxprob = phHard->VtxProb();
    }
    
    //fill jet variables
    if (fEnableJets) {
      for (UInt_t ijet=0; ijet<fPFJets->GetEntries();++ijet) {
        const Jet *jet = fPFJets->At(ijet);
        if (jet->AbsEta()<4.7 && MathUtils::DeltaR(jet,p1)>0.3 && MathUtils::DeltaR(jet,p2)>0.3) {
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
      _mass = (phHard->Mom()+phSoft->Mom()).M();
      _masserr = 0.5*_mass*TMath::Sqrt(phHard->EnergyErr()*phHard->EnergyErr()/phHard->E()/phHard->E() + phSoft->EnergyErr()*phSoft->EnergyErr()/phSoft->E()/phSoft->E());
      _masserrsmeared = 0.5*_mass*TMath::Sqrt(phHard->EnergyErrSmeared()*phHard->EnergyErrSmeared()/phHard->E()/phHard->E() + phSoft->EnergyErrSmeared()*phSoft->EnergyErrSmeared()/phSoft->E()/phSoft->E());
      _ptgg = (phHard->Mom()+phSoft->Mom()).Pt();
      _etagg = (phHard->Mom()+phSoft->Mom()).Eta();
      _phigg = (phHard->Mom()+phSoft->Mom()).Phi();
      _costheta = ThreeVector(phHard->Mom()).Unit().Dot(ThreeVector(phSoft->Mom()).Unit());
      _evtcat = PhotonTools::DiphotonR9EtaPtCat(phHard,phSoft);
      
      const Double_t dz = sqrt(2.0)*5.8;
      Double_t deltamvtx = _mass*VertexTools::DeltaMassVtx(phHard->CaloPos().X(), phHard->CaloPos().Y(), phHard->CaloPos().Z(),
            phSoft->CaloPos().X(), phSoft->CaloPos().Y(), phSoft->CaloPos().Z(),
            dz);
            
      fDiphotonEvent->deltamvtx = deltamvtx;
            
      _masserrwrongvtx = TMath::Sqrt(_masserr*_masserr + deltamvtx*deltamvtx);
      _masserrsmearedwrongvtx = TMath::Sqrt(_masserrsmeared*_masserrsmeared + deltamvtx*deltamvtx);
            
      if (jet1 && jet2) {
        fDiphotonEvent->zeppenfeld = TMath::Abs(_etagg - 0.5*(jet1->Eta()+jet2->Eta()));
        fDiphotonEvent->dphidijetgg = MathUtils::DeltaPhi( (jet1->Mom()+jet2->Mom()).Phi(), _phigg );
      }
      
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

    fDiphotonEvent->photons[0].SetVars(phHard,conv1,ele1,pfsc1,phgen1,fPhfixph,fPhfixele);
    fDiphotonEvent->photons[1].SetVars(phSoft,conv2,ele2,pfsc2,phgen2,fPhfixph,fPhfixele);
    
    Float_t ph1ecor    = fDiphotonEvent->photons[0].Ecor();
    Float_t ph1ecorerr = fDiphotonEvent->photons[0].Ecorerr();
    Float_t ph2ecor    = fDiphotonEvent->photons[1].Ecor();
    Float_t ph2ecorerr = fDiphotonEvent->photons[1].Ecorerr();

    Float_t ph1ecorele    = fDiphotonEvent->photons[0].Ecorele();
    Float_t ph1ecoreleerr = fDiphotonEvent->photons[0].Ecoreleerr();
    Float_t ph2ecorele    = fDiphotonEvent->photons[1].Ecorele();
    Float_t ph2ecoreleerr = fDiphotonEvent->photons[1].Ecoreleerr();
    
    
    fDiphotonEvent->masscor = TMath::Sqrt(2.0*ph1ecor*ph2ecor*(1.0-fDiphotonEvent->costheta));
    fDiphotonEvent->masscorerr = 0.5*fDiphotonEvent->masscor*TMath::Sqrt(ph1ecorerr*ph1ecorerr/ph1ecor/ph1ecor + ph2ecorerr*ph2ecorerr/ph2ecor/ph2ecor);
    
    fDiphotonEvent->masscorele = TMath::Sqrt(2.0*ph1ecorele*ph2ecorele*(1.0-fDiphotonEvent->costheta));
    fDiphotonEvent->masscoreleerr = 0.5*fDiphotonEvent->masscorele*TMath::Sqrt(ph1ecoreleerr*ph1ecoreleerr/ph1ecorele/ph1ecorele + ph2ecoreleerr*ph2ecoreleerr/ph2ecorele/ph2ecorele);    
    
    //printf("r9 = %5f, photon sigieie = %5f, seed sigieie = %5f\n",phHard->R9(),phHard->CoviEtaiEta(),sqrt(phHard->SCluster()->Seed()->CoviEtaiEta()));
    
    if (fWriteDiphotonTree) hCiCTuple->Fill();  
        
  }

  if (!fWriteSingleTree) return;


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
    
    const DecayParticle *conv = PhotonTools::MatchedConversion(sc,fConversions,bsp,nhitsbeforevtxmax);
    const SuperCluster *pfsc = PhotonTools::MatchedSC(sc,fSuperClusters);

    
    
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
      fDiphotonEvent->mt = TMath::Sqrt(2.0*fPFMet->At(0)->Pt()*ph->Pt()*(1.0-fDiphotonEvent->cosphimet));
    }
    
    if (ele) {
      fDiphotonEvent->cosphimetele = TMath::Cos(ele->Phi()-fPFMet->At(0)->Phi());
      fDiphotonEvent->mtele = TMath::Sqrt(2.0*fPFMet->At(0)->Pt()*ele->Pt()*(1.0-fDiphotonEvent->cosphimetele));      
    }
    
    fSinglePhoton->SetVars(ph,conv,ele,pfsc,phgen,fPhfixph,fPhfixele);
    hCiCTupleSingle->Fill();
    
  }

    
  return;

}

//--------------------------------------------------------------------------------------------------
void PhotonTreeWriter::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the photon collection branch.

  ReqEventObject(fPhotonBranchName,   fPhotons,    fPhotonsFromBranch);
  ReqEventObject(fTrackBranchName,    fTracks,     true);
  ReqEventObject(fElectronName,       fElectrons,  true);  
  ReqEventObject(fGoodElectronName,       fGoodElectrons,   fGoodElectronsFromBranch);  
  ReqEventObject(fPileUpDenName,      fPileUpDen,  true);
  ReqEventObject(fPVName,             fPV,         fPVFromBranch);
  ReqEventObject(fConversionName,     fConversions,true);
  ReqEventObject(fBeamspotName,       fBeamspot,   true);
  ReqEventObject(fPFCandName,         fPFCands,    true);
  ReqEventObject(fSuperClusterName,   fSuperClusters, true);
  ReqEventObject(fPFMetName,   fPFMet, true);
  if (fEnableJets) ReqEventObject(fPFJetName,   fPFJets, fPFJetsFromBranch);
  
  if (!fIsData) {
    ReqBranch(fPileUpName,            fPileUp);
    ReqBranch(fMCParticleName,        fMCParticles);
  }
  

  if (fIsData) {
    fPhFixDataFile = gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFixGRPV22.dat");
  }
  else {
    fPhFixDataFile = gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFixSTART42V13.dat");
  }
  
  //initialize photon energy corrections
  //PhotonFix::initialise("4_2",std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFix.dat")).Data()));  

  fPhfixph.initialise("4_2",std::string(fPhFixDataFile));
  fPhfixele.initialise("4_2e",std::string(fPhFixDataFile));
  
  fDiphotonEvent = new PhotonTreeWriterDiphotonEvent;
  fSinglePhoton = new PhotonTreeWriterPhoton;

  
  if (fWriteDiphotonTree) hCiCTuple = new TTree(fTupleName.Data(),fTupleName.Data());
  TString singlename = fTupleName + TString("Single");
  if (fWriteSingleTree) hCiCTupleSingle = new TTree(singlename,singlename);
  
  //make flattish tree from classes so we don't have to rely on dictionaries for reading later
  TClass *eclass = TClass::GetClass("mithep::PhotonTreeWriterDiphotonEvent");
  TClass *pclass = TClass::GetClass("mithep::PhotonTreeWriterPhoton");
  TList *elist = eclass->GetListOfDataMembers();
  TList *plist = pclass->GetListOfDataMembers();
    
  for (int i=0; i<elist->GetEntries(); ++i) {
    const TDataMember *tdm = static_cast<const TDataMember*>(elist->At(i));
    if (!(tdm->IsBasic() && tdm->IsPersistent())) continue;
    TString typestring;
    if (TString(tdm->GetTypeName())=="Char_t") typestring = "B";
    else if (TString(tdm->GetTypeName())=="UChar_t") typestring = "b";
    else if (TString(tdm->GetTypeName())=="Short_t") typestring = "S";
    else if (TString(tdm->GetTypeName())=="UShort_t") typestring = "s";
    else if (TString(tdm->GetTypeName())=="Int_t") typestring = "I";
    else if (TString(tdm->GetTypeName())=="UInt_t") typestring = "i";
    else if (TString(tdm->GetTypeName())=="Float_t") typestring = "F";
    else if (TString(tdm->GetTypeName())=="Double_t") typestring = "D";
    else if (TString(tdm->GetTypeName())=="Long64_t") typestring = "L";
    else if (TString(tdm->GetTypeName())=="ULong64_t") typestring = "l";
    else if (TString(tdm->GetTypeName())=="Bool_t") typestring = "O";
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
      if (TString(tdm->GetTypeName())=="Char_t") typestring = "B";
      else if (TString(tdm->GetTypeName())=="UChar_t") typestring = "b";
      else if (TString(tdm->GetTypeName())=="Short_t") typestring = "S";
      else if (TString(tdm->GetTypeName())=="UShort_t") typestring = "s";
      else if (TString(tdm->GetTypeName())=="Int_t") typestring = "I";
      else if (TString(tdm->GetTypeName())=="UInt_t") typestring = "i";
      else if (TString(tdm->GetTypeName())=="Float_t") typestring = "F";
      else if (TString(tdm->GetTypeName())=="Double_t") typestring = "D";
      else if (TString(tdm->GetTypeName())=="Long64_t") typestring = "L";
      else if (TString(tdm->GetTypeName())=="ULong64_t") typestring = "l";
      else if (TString(tdm->GetTypeName())=="Bool_t") typestring = "O";
      else continue;
      //printf("%s\n",tdm->GetTypeName());
      TString varname = TString::Format("ph%d.%s",iph+1,tdm->GetName());
      
      Char_t *addr = (Char_t*)&fDiphotonEvent->photons[iph];
      assert(sizeof(Char_t)==1);
      if (fWriteDiphotonTree) hCiCTuple->Branch(varname,addr+tdm->GetOffset(),TString::Format("%s/%s",varname.Data(),typestring.Data()));
      
      if (iph==0) {
        TString singlename = TString::Format("ph.%s",tdm->GetName());
        Char_t *addrsingle = (Char_t*)fSinglePhoton;
        if (fWriteSingleTree) hCiCTupleSingle->Branch(singlename,addrsingle+tdm->GetOffset(),TString::Format("%s/%s",singlename.Data(),typestring.Data()));
      }
    }
  }


  if (fWriteDiphotonTree) AddOutput(hCiCTuple);
  if (fWriteSingleTree) AddOutput(hCiCTupleSingle);


}

// ----------------------------------------------------------------------------------------
// some helpfer functions....
void PhotonTreeWriter::FindHiggsPtAndZ(Float_t& pt, Float_t& decayZ, Float_t& mass) {

  pt = -999.;
  decayZ = -999.;
  mass = -999.;

  // loop over all GEN particles and look for status 1 photons
  for(UInt_t i=0; i<fMCParticles->GetEntries(); ++i) {
    const MCParticle* p = fMCParticles->At(i);
    if( p->Is(MCParticle::kH) || (!fApplyElectronVeto && (p->AbsPdgId()==23||p->AbsPdgId()==24) ) ) {
      pt=p->Pt();
      decayZ = p->DecayVertex().Z();
      mass = p->Mass();
      break;
    }
  }
  
  return;
 }


Float_t PhotonTreeWriter::GetEventCat(PhotonTools::CiCBaseLineCats cat1, PhotonTools::CiCBaseLineCats cat2) {
  
  bool ph1IsEB = (cat1 ==  PhotonTools::kCiCCat1 || cat1 == PhotonTools::kCiCCat2);
  bool ph2IsEB = (cat2 ==  PhotonTools::kCiCCat1 || cat2 == PhotonTools::kCiCCat2);

  bool ph1IsHR9 = (cat1 ==  PhotonTools::kCiCCat1 || cat1 == PhotonTools::kCiCCat3);
  bool ph2IsHR9 = (cat2 ==  PhotonTools::kCiCCat1 || cat2 == PhotonTools::kCiCCat3);
  
  if( ph1IsEB && ph2IsEB )
    return ( ph1IsHR9 && ph2IsHR9 ? 0. : 1.);
  
  return ( ph1IsHR9 && ph2IsHR9 ? 2. : 3.);
}

void PhotonTreeWriterPhoton::SetVars(const Photon *p, const DecayParticle *c, const Electron *ele, const SuperCluster *pfsc, const MCParticle *m, PhotonFix &phfixph, PhotonFix &phfixele) {
  
      const SuperCluster *s = 0;
      if (p) {
        s = p->SCluster();
      }
      else {
        s = ele->SCluster();
      }
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
        sigietaieta = p->CoviEtaiEta();      
        phcat = PhotonTools::CiCBaseLineCat(p);
        eerr = p->EnergyErr();
        eerrsmeared = p->EnergyErrSmeared();
        esmearing = p->EnergySmearing();
        idmva = p->IdMva();
        hcalisodr03 = p->HcalTowerSumEtDr03();
        ecalisodr03 = p->EcalRecHitIsoDr03();
        trkisohollowdr03 = p->HollowConeTrkIsoDr03();
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
      }
       
      
      sce = s->Energy();
      scrawe = s->RawEnergy();
      scpse = s->PreshowerEnergy();
      sceta = s->Eta();
      scphi = s->Phi();
      scnclusters = s->ClusterSize();
      scnhits = s->NHits();
      scetawidth = s->EtaWidth();
      scphiwidth = s->PhiWidth();
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
      }
      else {
        haspfsc = kFALSE;
        pfsce = -99.;
        pfscrawe = -99.;
        pfsceta = -99.;
        pfscphi = -99.;
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

