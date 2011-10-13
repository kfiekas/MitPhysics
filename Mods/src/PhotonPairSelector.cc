#include "MitPhysics/Mods/interface/PhotonPairSelector.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/StableParticle.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/PhotonFix.h"
#include "MitPhysics/Utils/interface/MVATools.h"
#include "TDataMember.h"
#include <TNtuple.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TH1D.h>

using namespace mithep;

ClassImp(mithep::PhotonPairSelector)

//--------------------------------------------------------------------------------------------------
PhotonPairSelector::PhotonPairSelector(const char *name, const char *title) : 
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

  fGoodPhotonsName   (ModNames::gkGoodPhotonsName),

  // ----------------------------------------
  // Selection Types
  fPhotonSelType     ("NoSelection"),
  fVertexSelType     ("StdSelection"),
  fPhSelType         (kNoPhSelection),
  fVtxSelType        (kStdVtxSelection),
  
  // ----------------------------------------
  fPhotonPtMin       (20.0),
  fPhotonEtaMax      (2.5),

  fLeadingPtMin      (40.0),
  fTrailingPtMin     (30.0),

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

  // ---------------------------------------
  fDataEnCorr_EB_hR9 (0.),
  fDataEnCorr_EB_lR9 (0.),
  fDataEnCorr_EE_hR9 (0.),
  fDataEnCorr_EE_lR9 (0.),

  fRunStart          (0),
  fRunEnd            (0),

  fMCSmear_EB_hR9    (0.),
  fMCSmear_EB_lR9    (0.),
  fMCSmear_EE_hR9    (0.),
  fMCSmear_EE_lR9    (0.),

  // ---------------------------------------
  rng                (new TRandom3()),  
  fDoRegression      (kFALSE),
  fPhFixString       ("4_2"),
  fRegWeights         (gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/gbrph.root")),
  fEtaCorrections    (0),
  // ---------------------------------------
  fDoDataEneCorr     (true),
  fDoMCSmear         (true),
  fDoVtxSelection    (true),
  fApplyEleVeto      (true),
  fInvertElectronVeto(kFALSE),
  //MVA
  fVariableType      (2), 
  fEndcapWeights      (gSystem->Getenv("CMSSW_BASE")+TString("/src/MitPhysics/data/TMVAClassificationPhotonID_NewMotherId_Endcap_PtMin30_IsoCut250_VariableType2_BDTnCuts2000_ApplyElecVeto1_PuWeight_BDT.weights.xml")),
  fBarrelWeights      (gSystem->Getenv("CMSSW_BASE")+TString("/src/MitPhysics/data/TMVAClassificationPhotonID_NewMotherId_Barrel_PtMin30_IsoCut250_VariableType2_BDTnCuts2000_ApplyElecVeto1_PuWeight_BDT.weights.xml")),
  fbdtCutBarrel      (0.0031324),
  fbdtCutEndcap      (0.0086)
{
  // Constructor.
}

PhotonPairSelector::~PhotonPairSelector(){
  if(rng) delete rng;
}

//--------------------------------------------------------------------------------------------------
void PhotonPairSelector::Process()
{
  // ------------------------------------------------------------  
  // Process entries of the tree. 
  LoadEventObject(fPhotonBranchName,   fPhotons);
  
  // -----------------------------------------------------------
  // OUtput Photon Collection. It will ALWAYS conatrin either 0 or 2 Photons
  PhotonOArr *GoodPhotons = new PhotonOArr;
  GoodPhotons->SetName(fGoodPhotonsName);
  GoodPhotons->SetOwner(kTRUE);
  // add to event for other modules to use
  AddObjThisEvt(GoodPhotons);  
  
  if (fPhotons->GetEntries()<2) return;
  
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
  Float_t _tRho  = -99.;
  if( fPileUpDen->GetEntries() > 0 )
    _tRho  = (Double_t) fPileUpDen->At(0)->RhoRandomLowEta();
  
  const BaseVertex *bsp = dynamic_cast<const BaseVertex*>(fBeamspot->At(0));

  // ------------------------------------------------------------  
  // Get Event header for Run info etc.
  const EventHeader* evtHead = this->GetEventHeader();
  unsigned int evtNum = evtHead->EvtNum();
  //Float_t _evtNum1   = (Float_t) ( (int) (evtNum/10000.) );
  //Float_t _evtNum2   = (Float_t) ( (int) (evtNum % 10000)  );
  UInt_t   runNumber = evtHead->RunNum();
  Float_t _runNum    = (Float_t) runNumber;
  Float_t _lumiSec   = (Float_t) evtHead->LumiSec();

  // ------------------------------------------------------------
  // here we'll store the preselected Photons (and which CiCCategory they are...)
  PhotonOArr* preselPh  = new PhotonOArr;
  std::vector<PhotonTools::CiCBaseLineCats> preselCat;
  preselCat.resize(0);
  
  // 1. we do the pre-selection; but keep the non-passing photons in a secont container...
  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {    
    const Photon *ph = fPhotons->At(i);

    if(ph->SCluster()->AbsEta()>= fPhotonEtaMax || (ph->SCluster()->AbsEta()>=1.4442 && ph->SCluster()->AbsEta()<=1.566)) continue;
    if(ph->Et()                <  fPhotonPtMin)     continue;
    if(ph->HadOverEm()         >  0.15)     continue;
    if(ph->IsEB()) {
      if(ph->CoviEtaiEta() > 0.013) continue;      
    } else {
      if(ph->CoviEtaiEta() > 0.03) continue;
    }    
    preselPh->Add(ph);
  }

  if (preselPh->GetEntries()<2) return;

  // Sorry... need the second loop here in order to sort & assign the right Categories..
  preselPh->Sort();
  for(unsigned int iPh = 0; iPh <preselPh->GetEntries(); ++iPh) {
    const Photon* ph = preselPh->At(iPh);
    preselCat.push_back(PhotonTools::CiCBaseLineCat(ph));
  }
  
  // ------------------------------------------------------------
  // compute how many pairs there are ...
  unsigned int numPairs = 0;
  if( preselPh->GetEntries() > 0) numPairs = (preselPh->GetEntries()-1)*preselPh->GetEntries()/2;  
  // ... and create all possible pairs of pre-selected photons
  std::vector<unsigned int> idx1st;
  std::vector<unsigned int> idx2nd;
  std::vector<PhotonTools::CiCBaseLineCats> cat1st;
  std::vector<PhotonTools::CiCBaseLineCats> cat2nd;
  // ... this will be used to store whether a givne pair passes the cuts
  std::vector<bool> pairPasses;
  
  if(numPairs > 0) {
    for(unsigned int i1st = 0; i1st <preselPh->GetEntries() - 1; ++i1st) {
      for(unsigned int i2nd = i1st + 1; i2nd <preselPh->GetEntries(); ++i2nd) {
	idx1st.push_back(i1st);
	idx2nd.push_back(i2nd);
	pairPasses.push_back(true);
      }
    }
  }

  // ------------------------------------------------------------  
  // array to store the index of 'chosen Vtx' for each pair
  const Vertex** theVtx        = new const Vertex*[numPairs];    // holds the 'chosen' Vtx for each Pair
  Photon**       fixPh1st      = new       Photon*[numPairs];    // holds the 1st Photon for each Pair       
  Photon**       fixPh2nd      = new       Photon*[numPairs];    // holds the 2nd photon for each Pair
  

  // store pair-indices for pairs passing the selection
  std::vector<unsigned int> passPairs;
  passPairs.resize(0);
  
  // ------------------------------------------------------------  
  // Loop over all Pairs and to the 'incredible machine' running....
  for(unsigned int iPair = 0; iPair < numPairs; ++iPair) {    

    // first we need a hard copy of the incoming photons
    fixPh1st[iPair] = new Photon(*preselPh->At(idx1st[iPair]));
    fixPh2nd[iPair] = new Photon(*preselPh->At(idx2nd[iPair]));
    // we also store the category, so we don't have to ask all the time...
    cat1st.push_back(preselCat[idx1st[iPair]]);
    cat2nd.push_back(preselCat[idx2nd[iPair]]);

    if (fDoRegression) {
      if (!egcor.IsInitialized()) {
        //egcor.Initialize(!fIsData,"4_2",gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFixGRPV22.dat"),"/scratch/bendavid/root/weights-Base/TMVARegressionebph_BDTG.weights.xml","/scratch/bendavid/root/weights-Base/TMVARegressionVarianceebph_BDTG.weights.xml","/scratch/bendavid/root/weights-Base/TMVARegressioneeph_BDTG.weights.xml","/scratch/bendavid/root/weights-Base/TMVARegressionVarianceeeph_BDTG.weights.xml");
        egcor.Initialize(!fIsData,fPhFixString,fPhFixFile,fRegWeights);
      }
    
      egcor.CorrectEnergyWithError(fixPh1st[iPair]);
      egcor.CorrectEnergyWithError(fixPh2nd[iPair]);
      
      ThreeVectorC scpos1 = fixPh1st[iPair]->SCluster()->Point();
      ThreeVectorC scpos2 = fixPh2nd[iPair]->SCluster()->Point();
      
      fixPh1st[iPair]->SetCaloPosXYZ(scpos1.X(),scpos1.Y(),scpos1.Z());
      fixPh2nd[iPair]->SetCaloPosXYZ(scpos2.X(),scpos2.Y(),scpos2.Z());
      
      
    }
    
    // now we dicide if we either scale (Data) or Smear (MC) the Photons
    if (fIsData) {
      if(fDoDataEneCorr) {
	// statring with scale = 1.
	double scaleFac1 = 1.;
	double scaleFac2 = 1.;
        
        //eta-dependent corrections
        if (fEtaCorrections) {
          double etacor1 = fEtaCorrections->GetBinContent(fEtaCorrections->GetXaxis()->FindFixBin(fixPh1st[iPair]->SCluster()->Eta()));
          double etacor2 = fEtaCorrections->GetBinContent(fEtaCorrections->GetXaxis()->FindFixBin(fixPh2nd[iPair]->SCluster()->Eta()));
          
          if (fixPh1st[iPair]->SCluster()->AbsEta()>1.5) scaleFac1 *= (etacor1*etacor1);
          if (fixPh2nd[iPair]->SCluster()->AbsEta()>1.5) scaleFac2 *= (etacor2*etacor2);
        }
        
	// checking the run Rangees ...
	Int_t runRange = FindRunRangeIdx(runNumber);
	if(runRange > -1) { 
	  scaleFac1 /= (1.0+GetDataEnCorr(runRange, cat1st[iPair]));
	  scaleFac2 /= (1.0+GetDataEnCorr(runRange, cat2nd[iPair]));
	}      
	PhotonTools::ScalePhoton(fixPh1st[iPair], scaleFac1);
	PhotonTools::ScalePhoton(fixPh2nd[iPair], scaleFac2);
      }
    } 
    
    if(fDoMCSmear) {      
      
      double width1 = GetMCSmearFac(cat1st[iPair]);
      double width2 = GetMCSmearFac(cat2nd[iPair]);

      if (!fIsData) {
        // get the seed to do deterministic smearing...
        UInt_t seedBase = (UInt_t) evtNum + (UInt_t) _runNum + (UInt_t) _lumiSec;
        UInt_t seed1    = seedBase + (UInt_t) fixPh1st[iPair]->E() + (UInt_t) (TMath::Abs(10.*fixPh1st[iPair]->SCluster()->Eta()));
        UInt_t seed2    = seedBase + (UInt_t) fixPh2nd[iPair]->E() + (UInt_t) (TMath::Abs(10.*fixPh2nd[iPair]->SCluster()->Eta()));
        // get the smearing for MC photons..

        PhotonTools::SmearPhoton(fixPh1st[iPair], rng, width1, seed1);
        PhotonTools::SmearPhoton(fixPh2nd[iPair], rng, width2, seed2);
      }

    
      PhotonTools::SmearPhotonError(fixPh1st[iPair], width1);
      PhotonTools::SmearPhotonError(fixPh2nd[iPair], width2);
    
    }
    

    // store the vertex for this pair
    switch( fVtxSelType ){
    case kStdVtxSelection:
      theVtx[iPair] = fPV->At(0);
      break;
    case kCiCVtxSelection:
      theVtx[iPair] = VertexTools::findVtxBasicRanking(fixPh1st[iPair],fixPh2nd[iPair], bsp, fPV, fConversions);
      break;
    case kMITVtxSelection:
      // need PFCandidate Collection
      theVtx[iPair] = VertexTools::BestVtx(fPFCands, fPV, bsp, mithep::FourVector((fixPh1st[iPair]->Mom()+fixPh2nd[iPair]->Mom()))); 
      break;
    default:
      theVtx[iPair] = fPV->At(0);

    }
    
    //set PV ref in photons
    fixPh1st[iPair]->SetPV(theVtx[iPair]);
    fixPh2nd[iPair]->SetPV(theVtx[iPair]);

    // fix the kinematics for both events
    FourVectorM newMom1st = fixPh1st[iPair]->MomVtx(theVtx[iPair]->Position());
    FourVectorM newMom2nd = fixPh2nd[iPair]->MomVtx(theVtx[iPair]->Position());
    fixPh1st[iPair]->SetMom(newMom1st.X(), newMom1st.Y(), newMom1st.Z(), newMom1st.E());
    fixPh2nd[iPair]->SetMom(newMom2nd.X(), newMom2nd.Y(), newMom2nd.Z(), newMom2nd.E());

    /* Float_t bdt1=-99;
    Float_t bdt2=-99;

    if(fixPh1st[iPair]->HasPV() && fixPh2nd[iPair]->HasPV()){
      bdt1 = fTool.GetMVAbdtValue(fixPh1st[iPair],fixPh1st[iPair]->PV(),fTracks,fPV,_tRho,fElectrons);
      bdt2 = fTool.GetMVAbdtValue(fixPh2nd[iPair],fixPh2nd[iPair]->PV(),fTracks,fPV,_tRho,fElectrons);
    }
    else{
      bdt1 = fTool.GetMVAbdtValue(fixPh1st[iPair],fPV->At(0),fTracks,fPV,_tRho,fElectrons);
      bdt2 = fTool.GetMVAbdtValue(fixPh2nd[iPair],fPV->At(0),fTracks,fPV,_tRho,fElectrons);
    }
 
    fixPh1st[iPair]->SetBDT(bdt1);
    fixPh2nd[iPair]->SetBDT(bdt2); */

    // check if both photons pass the CiC selection
    // FIX-ME: Add other possibilities....
    bool pass1 = false;
    bool pass2 = false;
 
    switch( fPhSelType ){
    case kNoPhSelection:
      pass1 = ( fixPh1st[iPair]->Pt() > fLeadingPtMin  );
      pass2 = ( fixPh2nd[iPair]->Pt() > fTrailingPtMin );
      break;
    case kCiCPhSelection:


      pass1 = PhotonTools::PassCiCSelection(fixPh1st[iPair], theVtx[iPair], fTracks, fElectrons, fPV, _tRho, fLeadingPtMin, fApplyEleVeto);
      if(pass1) pass2 = PhotonTools::PassCiCSelection(fixPh2nd[iPair], theVtx[iPair], fTracks, fElectrons, fPV, _tRho, fTrailingPtMin, fApplyEleVeto);

      break;
    case kMVAPhSelection://MVA
      pass1 = fTool.PassMVASelection(fixPh1st[iPair],theVtx[iPair],fTracks,fPV,_tRho,fElectrons,fLeadingPtMin,fbdtCutBarrel,fbdtCutEndcap);
      if(pass1) pass2 = fTool.PassMVASelection(fixPh2nd[iPair],theVtx[iPair],fTracks,fPV,_tRho,fElectrons,fTrailingPtMin,fbdtCutBarrel,fbdtCutEndcap);
      
      break;
    case kMITPhSelection:
      // FIX-ME: This is a place-holder.. MIT guys: Please work hard... ;)
      pass1 = ( fixPh1st[iPair]->Pt() > fLeadingPtMin  );
      pass2 = ( fixPh2nd[iPair]->Pt() > fTrailingPtMin );
      break;
    default:
      pass1 = true;
      pass2 = true;
    }
    
    //match to good electrons if requested
    if (fInvertElectronVeto) {
      pass1 &= !PhotonTools::PassElectronVeto(fixPh1st[iPair],fGoodElectrons);
      pass2 &= !PhotonTools::PassElectronVeto(fixPh2nd[iPair],fGoodElectrons);
    }
    // finally, if both Photons pass the selections, add the pair to the 'passing Pairs)
    if( pass1 && pass2 ) passPairs.push_back(iPair);
  }
  
  
  // ---------------------------------------------------------------
  // ... we're almost done, stau focused...
  // loop over all passing pairs and find the one with the highest sum Et
  const Vertex* _theVtx  = NULL;
  Photon*        phHard  = NULL;
  Photon*        phSoft  = NULL;

  PhotonTools::CiCBaseLineCats catPh1 = PhotonTools::kCiCNoCat;
  PhotonTools::CiCBaseLineCats catPh2 = PhotonTools::kCiCNoCat;
  
  double maxSumEt = 0.;
  for(unsigned int iPair=0; iPair<passPairs.size(); ++iPair){
    double sumEt = fixPh1st[passPairs[iPair]]->Et();
    sumEt += fixPh2nd[passPairs[iPair]]->Et();
    if( sumEt > maxSumEt ) {
      maxSumEt = sumEt;
      phHard = fixPh1st[passPairs[iPair]];
      phSoft = fixPh2nd[passPairs[iPair]];
      catPh1 = cat1st[passPairs[iPair]];
      catPh2 = cat2nd[passPairs[iPair]];
      _theVtx = theVtx[iPair];
    }
  }
  
  // ---------------------------------------------------------------
  // we have the Photons (*PARTY*)... compute some useful qunatities

  

  if(phHard && phSoft) {
    GoodPhotons->AddOwned(phHard);
    GoodPhotons->AddOwned(phSoft);
  }

  
  // sort according to pt
  GoodPhotons->Sort();
  
  // delete auxiliary photon collection...
  delete preselPh;
  delete[] theVtx;
    
  return;

}

//--------------------------------------------------------------------------------------------------
void PhotonPairSelector::SlaveBegin()
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
  
  if      (fPhotonSelType.CompareTo("CiCSelection") == 0) 
    fPhSelType =       kCiCPhSelection;
  else if (fPhotonSelType.CompareTo("MVASelection") == 0) //MVA
    fPhSelType =       kMVAPhSelection;
  else if (fPhotonSelType.CompareTo("MITSelection") == 0) 
    fPhSelType =       kMITPhSelection;
  else 
    fPhSelType =       kNoPhSelection;

  if      (fVertexSelType.CompareTo("CiCSelection") == 0) 
    fVtxSelType =       kCiCVtxSelection;
  else if (fVertexSelType.CompareTo("MITSelection") == 0) 
    fVtxSelType =       kMITVtxSelection;
  else 
    fVtxSelType =       kStdVtxSelection;  

  if (fIsData) {
    fPhFixFile = gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFixGRPV22.dat");
  }
  else {
    fPhFixFile = gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFixSTART42V13.dat");
  }

  printf("initialize pairselc\n");

  fTool.InitializeMVA(fVariableType,fEndcapWeights,fBarrelWeights);

}

// ----------------------------------------------------------------------------------------
// some helpfer functions....
void PhotonPairSelector::FindHiggsPtAndZ(Float_t& pt, Float_t& decayZ, Float_t& mass) {

  pt = -999.;
  decayZ = -999.;
  mass = -999.;

  // loop over all GEN particles and look for status 1 photons
  for(UInt_t i=0; i<fMCParticles->GetEntries(); ++i) {
    const MCParticle* p = fMCParticles->At(i);
    if( p->Is(MCParticle::kH) || (!fApplyEleVeto && p->AbsPdgId()==23) ) {
      pt=p->Pt();
      decayZ = p->DecayVertex().Z();
      mass = p->Mass();
      break;
    }
  }
  
  return;
 }

// this routine looks for the idx of the run-range
Int_t PhotonPairSelector::FindRunRangeIdx(UInt_t run) {
  Int_t runRange=-1;
  for(UInt_t iRun = 0; iRun<fRunStart.size(); ++iRun) {
    if( run >= fRunStart[iRun] && run <= fRunEnd[iRun]) {
      runRange = (Int_t) iRun;
      return runRange;
    }
  }
  return runRange;
}


Double_t PhotonPairSelector::GetDataEnCorr(Int_t runRange, PhotonTools::CiCBaseLineCats cat) {
  switch( cat ) {
  case PhotonTools::kCiCCat1:
    return fDataEnCorr_EB_hR9[runRange];
  case PhotonTools::kCiCCat2:
    return fDataEnCorr_EB_lR9[runRange];
  case PhotonTools::kCiCCat3:
    return fDataEnCorr_EE_hR9[runRange];
  case PhotonTools::kCiCCat4:
    return fDataEnCorr_EE_lR9[runRange];
  default:
    return 1.;
  }
}


Double_t PhotonPairSelector::GetMCSmearFac(PhotonTools::CiCBaseLineCats cat) {
  switch( cat ) {
  case PhotonTools::kCiCCat1:
    return fMCSmear_EB_hR9;
  case PhotonTools::kCiCCat2:
    return fMCSmear_EB_lR9;
  case PhotonTools::kCiCCat3:
    return fMCSmear_EE_hR9;
  case PhotonTools::kCiCCat4:
    return fMCSmear_EE_lR9;
  default:
    return 1.;
  }
}

Float_t PhotonPairSelector::GetEventCat(PhotonTools::CiCBaseLineCats cat1, PhotonTools::CiCBaseLineCats cat2) {
  
  bool ph1IsEB = (cat1 ==  PhotonTools::kCiCCat1 || cat1 == PhotonTools::kCiCCat2);
  bool ph2IsEB = (cat2 ==  PhotonTools::kCiCCat1 || cat2 == PhotonTools::kCiCCat2);

  bool ph1IsHR9 = (cat1 ==  PhotonTools::kCiCCat1 || cat1 == PhotonTools::kCiCCat3);
  bool ph2IsHR9 = (cat2 ==  PhotonTools::kCiCCat1 || cat2 == PhotonTools::kCiCCat3);
  
  if( ph1IsEB && ph2IsEB )
    return ( ph1IsHR9 && ph2IsHR9 ? 0. : 1.);
  
  return ( ph1IsHR9 && ph2IsHR9 ? 2. : 3.);
}

