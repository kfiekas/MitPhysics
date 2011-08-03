#include "MitPhysics/Mods/interface/PhotonTreeWriter.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/StableParticle.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/PhotonFix.h"
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
  // MC specific stuff...
  fMCParticleName    (Names::gkMCPartBrn),
  fPileUpName        (Names::gkPileupInfoBrn),

  
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

  fInvertElectronVeto(kFALSE),  

  fWriteDiphotonTree(kTRUE),
  fWriteSingleTree(kTRUE),
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
  
  if (fPhotons->GetEntries()<1) return;
  
  LoadEventObject(fElectronName,       fElectrons);
  LoadEventObject(fGoodElectronName,       fGoodElectrons);
  LoadEventObject(fConversionName,     fConversions);
  LoadEventObject(fTrackBranchName,    fTracks);
  LoadEventObject(fPileUpDenName,      fPileUpDen);
  LoadEventObject(fPVName,             fPV);    
  LoadEventObject(fBeamspotName,       fBeamspot);
  LoadEventObject(fPFCandName,         fPFCands);

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
  fDiphotonEvent->vtxZ = (fDiphotonEvent->nVtx>0) ? fPV->At(0)->Z() : -99.;
  fDiphotonEvent->numPU = _numPU;
  fDiphotonEvent->numPUminus = _numPUminus;
  fDiphotonEvent->numPUplus = _numPUplus;
  fDiphotonEvent->mass = -99.;
  fDiphotonEvent->ptgg = -99.;
  fDiphotonEvent->costheta = -99.;
  fDiphotonEvent->evt = GetEventHeader()->EvtNum();
  fDiphotonEvent->run = GetEventHeader()->RunNum();
  fDiphotonEvent->lumi = GetEventHeader()->LumiSec();
  fDiphotonEvent->evtcat = -99;
  fDiphotonEvent->masscor = -99.;
  fDiphotonEvent->masscorerr = -99.;

  Int_t nhitsbeforevtxmax = 1;
  if (fInvertElectronVeto) nhitsbeforevtxmax = 999;  
  
  if (fPhotons->GetEntries()>=2) {
    Bool_t doFill = kTRUE;
    
    const Photon *phHard = fPhotons->At(0);
    const Photon *phSoft = fPhotons->At(1);

    if (phHard->HasPV()) {
      fDiphotonEvent->vtxZ = phHard->PV()->Z();
    }
    
    Float_t _mass = ( doFill ? (phHard->Mom()+phSoft->Mom()).M()  : -100.);
    Float_t _ptgg = ( doFill ? (phHard->Mom()+phSoft->Mom()).Pt() : -100.);
    

    const DecayParticle *conv1 = PhotonTools::MatchedConversion(phHard,fConversions,bsp,nhitsbeforevtxmax);
    const DecayParticle *conv2 = PhotonTools::MatchedConversion(phSoft,fConversions,bsp,nhitsbeforevtxmax);

    const MCParticle *phgen1 = 0;
    const MCParticle *phgen2 = 0;
    if( !fIsData ) {
      phgen1 = PhotonTools::MatchMC(phHard,fMCParticles,!fInvertElectronVeto);
      phgen2 = PhotonTools::MatchMC(phSoft,fMCParticles,!fInvertElectronVeto);
    }
    
    Float_t _gencostheta = -99.;
    if (phgen1 && phgen2) {
      _gencostheta = ThreeVector(phgen1->Mom()).Unit().Dot(ThreeVector(phgen2->Mom()).Unit());
    }
  
    fDiphotonEvent->gencostheta = _gencostheta;
    fDiphotonEvent->mass = _mass;
    fDiphotonEvent->ptgg = _ptgg;
    fDiphotonEvent->costheta =  ThreeVector(phHard->Mom()).Unit().Dot(ThreeVector(phSoft->Mom()).Unit());
    fDiphotonEvent->evtcat = PhotonTools::DiphotonR9EtaPtCat(phHard,phSoft);

    fDiphotonEvent->photons[0].SetVars(phHard,conv1,phgen1);
    fDiphotonEvent->photons[1].SetVars(phSoft,conv2,phgen2);
    
    Float_t ph1ecor    = fDiphotonEvent->photons[0].Ecor();
    Float_t ph1ecorerr = fDiphotonEvent->photons[0].Ecorerr();
    Float_t ph2ecor    = fDiphotonEvent->photons[1].Ecor();
    Float_t ph2ecorerr = fDiphotonEvent->photons[1].Ecorerr();
    
    fDiphotonEvent->masscor = TMath::Sqrt(2.0*ph1ecor*ph2ecor*(1.0-fDiphotonEvent->costheta));
    fDiphotonEvent->masscorerr = 0.5*fDiphotonEvent->masscor*TMath::Sqrt(ph1ecorerr*ph1ecorerr/ph1ecor/ph1ecor + ph2ecorerr*ph2ecorerr/ph2ecor/ph2ecor);
    
    //printf("r9 = %5f, photon sigieie = %5f, seed sigieie = %5f\n",phHard->R9(),phHard->CoviEtaiEta(),sqrt(phHard->SCluster()->Seed()->CoviEtaiEta()));
    
    if (fWriteDiphotonTree) hCiCTuple->Fill();  
        
  }

  if (!fWriteSingleTree) return;


  for (UInt_t iph = 0; iph<fPhotons->GetEntries(); ++iph) {
    const Photon *ph = fPhotons->At(iph);
    
    if (ph->HasPV()) {
      fDiphotonEvent->vtxZ = ph->PV()->Z();
    }

    const MCParticle *phgen = 0;
    if( !fIsData ) {
      phgen = PhotonTools::MatchMC(ph,fMCParticles,!fInvertElectronVeto);
    }

    const DecayParticle *conv = PhotonTools::MatchedConversion(ph,fConversions,bsp,nhitsbeforevtxmax);
    
    fSinglePhoton->SetVars(ph,conv,phgen);
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
  
  if (!fIsData) {
    ReqBranch(fPileUpName,            fPileUp);
    ReqBranch(fMCParticleName,        fMCParticles);
  }
  

  
  //initialize photon energy corrections
  //PhotonFix::initialise("4_2",std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFix.dat")).Data()));  

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
    if( p->Is(MCParticle::kH) || (!fInvertElectronVeto && p->AbsPdgId()==23) ) {
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

void PhotonTreeWriterPhoton::SetVars(const Photon *p, const DecayParticle *c, const MCParticle *m) {
      e = p->E();
      pt = p->Pt();
      eta = p->Eta();
      phi = p->Phi();
      r9 = p->R9();
      e3x3 = p->E33();
      e5x5 = p->E55();
      const SuperCluster *s = p->SCluster();
      sce = s->Energy();
      scrawe = s->RawEnergy();
      scpse = s->PreshowerEnergy();
      sceta = s->Eta();
      scphi = s->Phi();
      scnclusters = s->ClusterSize();
      scnhits = s->NHits();
      hovere = p->HadOverEm();
      sigietaieta = p->CoviEtaiEta();
      isbarrel = (s->AbsEta()<1.5);
      isr9reco = (isbarrel && r9>0.94) || (!isbarrel && r9>0.95);
      isr9cat = (r9>0.94);
      
      phcat = PhotonTools::CiCBaseLineCat(p);
      
/*      if (isbarrel) {
        if (isr9cat) phcat = 1;
        else phcat = 2;
      }
      else {
        if (isr9cat) phcat = 3;
        else phcat = 4;
      }*/
      
      const BasicCluster *b = s->Seed();
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

      //initialize photon energy corrections if needed
      if (!PhotonFix::initialised()) {
        PhotonFix::initialise("4_2",std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFix.dat")).Data()));  
      }
      
      PhotonFix fix(e,sceta,scphi,r9);
      const Float_t dval = -99.;
      ecor = fix.fixedEnergy();
      ecorerr = fix.sigmaEnergy();
      if (isbarrel) {
       etac = fix.etaC();
       etas = fix.etaS();
       etam = fix.etaM();
       phic = fix.phiC();
       phis = fix.phiS();
       phim = fix.phiM();
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
       xz = fix.xZ();
       xc = fix.xC();
       xs = fix.xS();
       xm = fix.xM();
       yz = fix.yZ();
       yc = fix.yC();
       ys = fix.yS();
       ym = fix.yM();
      }   
      
      if (c) {
        hasconversion = kTRUE;
        convp = c->P();
        convpt = c->Pt();
        conveta = c->Eta();
        convphi = c->Phi();
        ThreeVector dirconvsc = ThreeVector(p->SCluster()->Point()) - c->Position();
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
