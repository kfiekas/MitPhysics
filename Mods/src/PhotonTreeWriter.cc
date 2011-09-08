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

  
  fIsData            (false),
  fPhotonsFromBranch (true),  
  fPVFromBranch      (true),
  fGoodElectronsFromBranch (kTRUE),

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

  fLoopOnGoodElectrons(kFALSE),
  fInvertElectronVeto(kFALSE),  

  fWriteDiphotonTree(kTRUE),
  fWriteSingleTree(kTRUE),
  fExcludeSinglePrompt(kFALSE),
  fExcludeDoublePrompt(kFALSE),
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
  
  Int_t nhitsbeforevtxmax = 1;
  if (fInvertElectronVeto) nhitsbeforevtxmax = 999;  
  
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
      phgen1 = PhotonTools::MatchMC(p1,fMCParticles,fInvertElectronVeto);
      phgen2 = PhotonTools::MatchMC(p2,fMCParticles,fInvertElectronVeto);
    }
    
    if (fExcludeSinglePrompt && (phgen1 || phgen2) ) return;
    if (fExcludeDoublePrompt && (phgen1 && phgen2) ) return;
    
    if (!fLoopOnGoodElectrons && phHard->HasPV()) {
      fDiphotonEvent->vtxX = phHard->PV()->X();
      fDiphotonEvent->vtxY = phHard->PV()->Y();
      fDiphotonEvent->vtxZ = phHard->PV()->Z();
    }
    
    Float_t _mass = -99.;
    Float_t _masserr = -99.;
    Float_t _masserrsmeared = -99.;
    Float_t _ptgg = -99.;
    Float_t _costheta = -99.;
    PhotonTools::DiphotonR9EtaPtCats _evtcat = PhotonTools::kOctCat0;
    if (phHard && phSoft) {
      _mass = (phHard->Mom()+phSoft->Mom()).M();
      _masserr = 0.5*_mass*TMath::Sqrt(phHard->EnergyErr()*phHard->EnergyErr()/phHard->E()/phHard->E() + phSoft->EnergyErr()*phSoft->EnergyErr()/phSoft->E()/phSoft->E());
      _masserrsmeared = 0.5*_mass*TMath::Sqrt(phHard->EnergyErrSmeared()*phHard->EnergyErrSmeared()/phHard->E()/phHard->E() + phSoft->EnergyErrSmeared()*phSoft->EnergyErrSmeared()/phSoft->E()/phSoft->E());
      _ptgg = (phHard->Mom()+phSoft->Mom()).Pt();
      _costheta = ThreeVector(phHard->Mom()).Unit().Dot(ThreeVector(phSoft->Mom()).Unit());
      _evtcat = PhotonTools::DiphotonR9EtaPtCat(phHard,phSoft);
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
    fDiphotonEvent->ptgg = _ptgg;
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
      phgen = PhotonTools::MatchMC(p,fMCParticles,fInvertElectronVeto);
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
  
  if (!fIsData) {
    ReqBranch(fPileUpName,            fPileUp);
    ReqBranch(fMCParticleName,        fMCParticles);
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
    if( p->Is(MCParticle::kH) || (fInvertElectronVeto && (p->AbsPdgId()==23||p->AbsPdgId()==24) ) ) {
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
  
//       if (p && ele) {
//         printf("p    : r9 = %5f, e33 = %5f, e55 = %5f, hovere = %5f, sigieie = %5f\n",p->R9(),p->E33(),p->E55(),p->HadOverEm(),p->CoviEtaiEta());
//         printf("ele/b: r9 = %5f, e33 = %5f, e55 = %5f, hovere = %5f, sigieie = %5f\n",b->E3x3()/s->RawEnergy(),b->E3x3(),b->E5x5(),ele->HadronicOverEm(),ele->CoviEtaiEta());
//       }
      
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
      
      
      
      
      sigiphiphi = TMath::Sqrt(b->CoviPhiiPhi());
      if (isnan(sigiphiphi)) sigiphiphi = -99.;
      covietaiphi = b->CoviEtaiPhi();
      if (isnan(covietaiphi)) covietaiphi = -99.;
      emax = b->EMax();
      e2nd = b->E2nd();
      etop = b->ETop();
      ebottom = b->EBottom();
      eleft = b->ELeft();
      eright = b->ERight();
      e1x3 = b->E1x3();
      e3x1 = b->E3x1();
      e1x5 = b->E1x5();
      e2x2 = b->E2x2();
      e4x4 = b->E4x4();
      e2x5max = b->E2x5Max();
      e2x5top = b->E2x5Top();
      e2x5bottom = b->E2x5Bottom();
      e2x5left = b->E2x5Left();
      e2x5right = b->E2x5Right();
      eseed = b->Energy();
      
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
      }
      else {
        ispromptgen = kFALSE;
        gene = -99.;
        genpt = -99.;
        geneta = -99.;
        genphi = -99.;
      }
            
}

