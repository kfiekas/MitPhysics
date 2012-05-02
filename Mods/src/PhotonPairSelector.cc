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
  BaseMod                        (name,title),
  // define all the Branches to load
  fPhotonBranchName              (Names::gkPhotonBrn),
  fElectronName                  (Names::gkElectronBrn),
  fGoodElectronName              (Names::gkElectronBrn),
  fConversionName                (Names::gkMvfConversionBrn),
  fTrackBranchName               (Names::gkTrackBrn),
  fPileUpDenName                 (Names::gkPileupEnergyDensityBrn),
  fPVName                        (Names::gkPVBeamSpotBrn),
  fBeamspotName                  (Names::gkBeamSpotBrn),
  fPFCandName                    (Names::gkPFCandidatesBrn),
  // MC specific stuff...
  fMCParticleName                (Names::gkMCPartBrn),
  fPileUpName                    (Names::gkPileupInfoBrn),
  fGoodPhotonsName               (ModNames::gkGoodPhotonsName),
  // ----------------------------------------
  // Selection Types
  fPhotonSelType                 ("NoSelection"),
  fVertexSelType                 ("StdSelection"),
  fPhSelType                     (kNoPhSelection),
  fVtxSelType                    (kStdVtxSelection),
  // ----------------------------------------
  fPhotonPtMin                   (20.0),
  fPhotonEtaMax                  (2.5),
  fLeadingPtMin                  (100.0/3.0),
  fTrailingPtMin                 (100.0/4.0),
  fIsData                        (false),
  fPhotonsFromBranch             (true),
  fPVFromBranch                  (true),
  fGoodElectronsFromBranch       (kTRUE),
  // ----------------------------------------
  // collections....
  fPhotons                       (0),
  fElectrons                     (0),
  fConversions                   (0),
  fTracks                        (0),
  fPileUpDen                     (0),
  fPV                            (0),
  fBeamspot                      (0),
  fPFCands                       (0),
  fMCParticles                   (0),
  fPileUp                        (0),
  // ---------------------------------------
  fDataEnCorr_EBlowEta_hR9central(0.),
  fDataEnCorr_EBlowEta_hR9gap    (0.),
  fDataEnCorr_EBlowEta_lR9       (0.),
  fDataEnCorr_EBhighEta_hR9      (0.),
  fDataEnCorr_EBhighEta_lR9      (0.),
  fDataEnCorr_EElowEta_hR9       (0.),
  fDataEnCorr_EElowEta_lR9       (0.),
  fDataEnCorr_EEhighEta_hR9      (0.),
  fDataEnCorr_EEhighEta_lR9      (0.),
  fRunStart                      (0),
  fRunEnd                        (0),
  fMCSmear_EBlowEta_hR9central   (0.),
  fMCSmear_EBlowEta_hR9gap       (0.),
  fMCSmear_EBlowEta_lR9          (0.),
  fMCSmear_EBhighEta_hR9         (0.),
  fMCSmear_EBhighEta_lR9         (0.),
  fMCSmear_EElowEta_hR9          (0.),
  fMCSmear_EElowEta_lR9          (0.),
  fMCSmear_EEhighEta_hR9         (0.),
  fMCSmear_EEhighEta_lR9         (0.),
  // ---------------------------------------
  fRng                           (new TRandom3()),
  fPhFixString                   ("4_2"),
  fEtaCorrections                (0),
  // ---------------------------------------
  fDoDataEneCorr                 (true),
  fDoMCSmear                     (true),
  fDoVtxSelection                (true),
  fApplyEleVeto                  (true),
  fInvertElectronVeto            (kFALSE),
  //MVA
  fVariableType                  (10), //please use 4 which is the correct type
  fEndcapWeights                 (gSystem->Getenv("CMSSW_BASE")+
                                  TString("/src/MitPhysics/data/TMVAClassificationPhotonID_")+
                                  TString("Endcap_PassPreSel_Variable_10_BDTnCuts2000_BDT.")+
                                  TString("weights.xml")),
  fBarrelWeights                 (gSystem->Getenv("CMSSW_BASE")+
                                  TString("/src/MitPhysics/data/TMVAClassificationPhotonID_")+
                                  TString("Barrel_PassPreSel_Variable_10_BDTnCuts2000_BDT.")+
                                  TString("weights.xml")),
  fbdtCutBarrel                  (0.0744), //cuts give same eff (relative to presel) with cic
  fbdtCutEndcap                  (0.0959), //cuts give same eff (relative to presel) with cic
  fDoMCR9Scaling                 (kFALSE),
  fMCR9ScaleEB                   (1.0),
  fMCR9ScaleEE                   (1.0),
  fDoMCSigIEtaIEtaScaling        (kFALSE),
  fDoMCWidthScaling              (kFALSE),
  fDoMCErrScaling                (kFALSE),
  fMCErrScaleEB                  (1.0),
  fMCErrScaleEE                  (1.0),
  fRelativePtCuts                (kFALSE)
{
  // Constructor.
}

PhotonPairSelector::~PhotonPairSelector()
{
  if (fRng)
    delete fRng;
}

//--------------------------------------------------------------------------------------------------
void PhotonPairSelector::Process()
{
  // ------------------------------------------------------------
  // Process entries of the tree.
  LoadEventObject(fPhotonBranchName,   fPhotons);

  // -----------------------------------------------------------
  // Output Photon Collection. It will ALWAYS contain either 0 or 2 photons
  PhotonOArr  *GoodPhotons = new PhotonOArr;
  GoodPhotons->SetName(fGoodPhotonsName);
  GoodPhotons->SetOwner(kTRUE);
  // add to event for other modules to use
  AddObjThisEvt(GoodPhotons);

  if (fPhotons->GetEntries()<2)
    return;

  LoadEventObject(fElectronName,       fElectrons);
  LoadEventObject(fGoodElectronName,   fGoodElectrons);
  LoadEventObject(fConversionName,     fConversions);
  LoadEventObject(fTrackBranchName,    fTracks);
  LoadEventObject(fPileUpDenName,      fPileUpDen);
  LoadEventObject(fPVName,             fPV);
  LoadEventObject(fBeamspotName,       fBeamspot);
  LoadEventObject(fPFCandName,         fPFCands);

  // ------------------------------------------------------------
  // load event based information
  Float_t rho = -99.;
  if (fPileUpDen->GetEntries() > 0)
    rho = (Double_t) fPileUpDen->At(0)->RhoRandomLowEta();
  const BaseVertex *bsp = dynamic_cast<const BaseVertex*>(fBeamspot->At(0));

  // ------------------------------------------------------------
  // Get Event header for Run info etc.
  const EventHeader* evtHead   = this->GetEventHeader();
  unsigned int       evtNum    = evtHead->EvtNum();
  UInt_t             runNumber = evtHead->RunNum();
  Float_t            _runNum   = (Float_t) runNumber;
  Float_t            _lumiSec  = (Float_t) evtHead->LumiSec();

  // ------------------------------------------------------------
  // here we'll store the preselected Photons (and which CiCCategory they are...)
  PhotonOArr* preselPh  = new PhotonOArr;
  std::vector<PhotonTools::CiCBaseLineCats> preselCat;
  preselCat.resize(0);

  // 1. we do the pre-selection; but keep the non-passing photons in a second container...
  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {
    const Photon *ph = fPhotons->At(i);

    if (ph->SCluster()->AbsEta()>= fPhotonEtaMax ||
        (ph->SCluster()->AbsEta()>=1.4442 && ph->SCluster()->AbsEta()<=1.566))
      continue;
    if (ph->Et()                <  fPhotonPtMin)
      continue;
    if (ph->HadOverEm()         >  0.15)
      continue;

    if (ph->IsEB()) {
      if (ph->CoviEtaiEta()     > 0.015)
        continue;
    }
    else {
      if (ph->CoviEtaiEta()     > 0.035)
        continue;
    }
    // photon passes the preselection
    ph->Mark();        // mark for later skimming
    preselPh->Add(ph);
  }

  if (preselPh->GetEntries()<2) {
    this->SkipEvent();
    delete preselPh;
    return;
  }

  // second loop: sort & assign the right categories..
  preselPh->Sort();
  for (unsigned int iPh = 0; iPh <preselPh->GetEntries(); ++iPh) {
    const Photon* ph = preselPh->At(iPh);
    preselCat.push_back(PhotonTools::CiCBaseLineCat(ph));
  }

  // ------------------------------------------------------------
  // compute how many pairs there are ...
  unsigned int numPairs = 0;
  if (preselPh->GetEntries() > 0)
    numPairs = (preselPh->GetEntries()-1)*preselPh->GetEntries()/2;

  // ... and create all possible pairs of pre-selected photons
  std::vector<unsigned int>                 idx1st;
  std::vector<unsigned int>                 idx2nd;
  std::vector<PhotonTools::CiCBaseLineCats> cat1st;
  std::vector<PhotonTools::CiCBaseLineCats> cat2nd;

  // ... this will be used to store whether a given pair passes the cuts
  std::vector<bool> pairPasses;

  if (numPairs > 0) {
    for (unsigned int i1st = 0; i1st <preselPh->GetEntries() - 1; ++i1st) {
      for (unsigned int i2nd = i1st + 1; i2nd <preselPh->GetEntries(); ++i2nd) {
        idx1st.push_back(i1st);
        idx2nd.push_back(i2nd);
        pairPasses.push_back(true);
      }
    }
  }

  // ------------------------------------------------------------
  // array to store the index of 'chosen Vtx' for each pair
//   const Vertex** theVtx        = new const Vertex*[numPairs];    // holds the 'chosen' Vtx for each Pair
//   Photon**       fixPh1st      = new       Photon*[numPairs];    // holds the 1st Photon for each Pair
//   Photon**       fixPh2nd      = new       Photon*[numPairs];    // holds the 2nd photon for each Pair

  std::vector<const Vertex*> theVtx; // holds the 'chosen' Vtx for each Pair
  std::vector<Photon*> fixPh1st;     // holds the 1st Photon for each Pair
  std::vector<Photon*> fixPh2nd;     // holds the 2nd photon for each Pair

  theVtx.reserve(numPairs);
  fixPh1st.reserve(numPairs);
  fixPh2nd.reserve(numPairs);

  // store pair-indices for pairs passing the selection
  std::vector<unsigned int> passPairs;
  passPairs.resize(0);

  // ------------------------------------------------------------
  // Loop over all Pairs and do the 'incredible machine' running....
  for (unsigned int iPair = 0; iPair < numPairs; ++iPair) {
    // first we need a hard copy of the incoming photons
    fixPh1st.push_back(new Photon(*preselPh->At(idx1st[iPair])));
    fixPh2nd.push_back(new Photon(*preselPh->At(idx2nd[iPair])));
    // we also store the category, so we don't have to ask all the time...
    cat1st.push_back(preselCat[idx1st[iPair]]);
    cat2nd.push_back(preselCat[idx2nd[iPair]]);

    //scale regression sigmaE in MC if activated
    if (fDoMCErrScaling && !fIsData) {
      if (fixPh1st[iPair]->SCluster()->AbsEta()<1.5)
        PhotonTools::ScalePhotonError(fixPh1st[iPair],fMCErrScaleEB);
      else
        PhotonTools::ScalePhotonError(fixPh1st[iPair],fMCErrScaleEE);

      if (fixPh2nd[iPair]->SCluster()->AbsEta()<1.5)
        PhotonTools::ScalePhotonError(fixPh2nd[iPair],fMCErrScaleEB);
      else
        PhotonTools::ScalePhotonError(fixPh2nd[iPair],fMCErrScaleEE);
    }

    //scale R9 in Monte Carlo if activated
    if (fDoMCR9Scaling && !fIsData) {
      if (fixPh1st[iPair]->SCluster()->AbsEta()<1.5)
        PhotonTools::ScalePhotonR9(fixPh1st[iPair],fMCR9ScaleEB);
      else
        PhotonTools::ScalePhotonR9(fixPh1st[iPair],fMCR9ScaleEE);

      if (fixPh2nd[iPair]->SCluster()->AbsEta()<1.5)
        PhotonTools::ScalePhotonR9(fixPh2nd[iPair],fMCR9ScaleEB);
      else
        PhotonTools::ScalePhotonR9(fixPh2nd[iPair],fMCR9ScaleEE);
    }

    if (fDoMCSigIEtaIEtaScaling && !fIsData) {
      if (fixPh1st[iPair]->SCluster()->AbsEta()<1.5) fixPh1st[iPair]->SetCoviEtaiEta(0.87*fixPh1st[iPair]->CoviEtaiEta() + 0.0011);
      else fixPh1st[iPair]->SetCoviEtaiEta(0.99*fixPh1st[iPair]->CoviEtaiEta());

      if (fixPh2nd[iPair]->SCluster()->AbsEta()<1.5) fixPh2nd[iPair]->SetCoviEtaiEta(0.87*fixPh2nd[iPair]->CoviEtaiEta() + 0.0011);
      else fixPh2nd[iPair]->SetCoviEtaiEta(0.99*fixPh2nd[iPair]->CoviEtaiEta());
    }


    if (fDoMCWidthScaling && !fIsData) {
      fixPh1st[iPair]->SetEtaWidth(0.99*fixPh1st[iPair]->EtaWidth());
      fixPh1st[iPair]->SetPhiWidth(0.99*fixPh1st[iPair]->PhiWidth());

      fixPh2nd[iPair]->SetEtaWidth(0.99*fixPh2nd[iPair]->EtaWidth());
      fixPh2nd[iPair]->SetPhiWidth(0.99*fixPh2nd[iPair]->PhiWidth());
    }

    PhotonTools::eScaleCats escalecat1 = PhotonTools::EScaleCat(fixPh1st[iPair]);
    PhotonTools::eScaleCats escalecat2 = PhotonTools::EScaleCat(fixPh2nd[iPair]);

    // now we dicide if we either scale (Data) or Smear (MC) the Photons
    if (fIsData) {
      if (fDoDataEneCorr) {
        // starting with scale = 1.
        double scaleFac1 = 1.;
        double scaleFac2 = 1.;

        //eta-dependent corrections

        // checking the run Rangees ...
        Int_t runRange = FindRunRangeIdx(runNumber);
        if(runRange > -1) {
          scaleFac1 *= GetDataEnCorr(runRange, escalecat1);
          scaleFac2 *= GetDataEnCorr(runRange, escalecat2);
        }
        PhotonTools::ScalePhoton(fixPh1st[iPair], scaleFac1);
        PhotonTools::ScalePhoton(fixPh2nd[iPair], scaleFac2);
      }
    }

    if (fDoMCSmear) {
      double width1 = GetMCSmearFac(escalecat1);
      double width2 = GetMCSmearFac(escalecat2);
      if (!fIsData) {
        // get the seed to do deterministic smearing...
        UInt_t seedBase = (UInt_t) evtNum + (UInt_t) _runNum + (UInt_t) _lumiSec;
        UInt_t seed1    = seedBase + (UInt_t) fixPh1st[iPair]->E() +
          (UInt_t) (TMath::Abs(10.*fixPh1st[iPair]->SCluster()->Eta()));
        UInt_t seed2    = seedBase + (UInt_t) fixPh2nd[iPair]->E() +
          (UInt_t) (TMath::Abs(10.*fixPh2nd[iPair]->SCluster()->Eta()));
        // get the smearing for MC photons..
        PhotonTools::SmearPhoton(fixPh1st[iPair], fRng, width1, seed1);
        PhotonTools::SmearPhoton(fixPh2nd[iPair], fRng, width2, seed2);
      }
      PhotonTools::SmearPhotonError(fixPh1st[iPair], width1);
      PhotonTools::SmearPhotonError(fixPh2nd[iPair], width2);
    }

    //probability that selected vertex is the correct one
    Double_t vtxProb = 1.0;

    // store the vertex for this pair
    switch (fVtxSelType) {

    case kStdVtxSelection:
      theVtx[iPair] = fPV->At(0);
      break;

    case kCiCVtxSelection:
      theVtx[iPair] = fVtxTools.findVtxBasicRanking(fixPh1st[iPair],fixPh2nd[iPair], bsp, fPV,
                                                    fConversions,kFALSE,vtxProb);
      break;

    case kCiCMVAVtxSelection:
      theVtx[iPair] = fVtxTools.findVtxBasicRanking(fixPh1st[iPair],fixPh2nd[iPair], bsp, fPV,
                                                    fConversions,kTRUE,vtxProb);
      break;

    case kMITVtxSelection:
      // need PFCandidate Collection
      theVtx[iPair] = VertexTools::BestVtx(fPFCands, fPV, bsp,
                                           mithep::FourVector((fixPh1st[iPair]->Mom()+
                                                               fixPh2nd[iPair]->Mom())));
      break;
    default:
      theVtx[iPair] = fPV->At(0);
    }

    //set PV ref in photons
    fixPh1st[iPair]->SetPV(theVtx[iPair]);
    fixPh2nd[iPair]->SetPV(theVtx[iPair]);
    fixPh1st[iPair]->SetVtxProb(vtxProb);
    fixPh2nd[iPair]->SetVtxProb(vtxProb);

    // fix the kinematics for both events
    FourVectorM newMom1st = fixPh1st[iPair]->MomVtx(theVtx[iPair]->Position());
    FourVectorM newMom2nd = fixPh2nd[iPair]->MomVtx(theVtx[iPair]->Position());
    fixPh1st[iPair]->SetMom(newMom1st.X(), newMom1st.Y(), newMom1st.Z(), newMom1st.E());
    fixPh2nd[iPair]->SetMom(newMom2nd.X(), newMom2nd.Y(), newMom2nd.Z(), newMom2nd.E());

    double pairmass = (fixPh1st[iPair]->Mom() + fixPh2nd[iPair]->Mom()).M();

    double leadptcut = fLeadingPtMin;
    double trailptcut = fTrailingPtMin;

    if (fixPh2nd[iPair]->Pt() > fixPh1st[iPair]->Pt()) {
      leadptcut = fTrailingPtMin;
      trailptcut = fLeadingPtMin;
    }


    if (fRelativePtCuts) {
      leadptcut = leadptcut*pairmass;
      trailptcut = trailptcut*pairmass;
    }


    //compute id bdt values
    Double_t bdt1 = fTool.GetMVAbdtValue(fixPh1st[iPair],theVtx[iPair],fTracks,fPV,rho, fElectrons, fApplyEleVeto);
    Double_t bdt2 = fTool.GetMVAbdtValue(fixPh2nd[iPair],theVtx[iPair],fTracks,fPV,rho, fElectrons, fApplyEleVeto);

    fixPh1st[iPair]->SetIdMva(bdt1);
    fixPh2nd[iPair]->SetIdMva(bdt2);


    //printf("applying id\n");

    // check if both photons pass the CiC selection
    // FIX-ME: Add other possibilities....
    bool pass1 = false;
    bool pass2 = false;

    switch (fPhSelType) {
    case kNoPhSelection:
      pass1 = (fixPh1st[iPair]->Pt() > leadptcut );
      pass2 = (fixPh2nd[iPair]->Pt() > trailptcut);
      break;
    case kCiCPhSelection:
      pass1 = PhotonTools::PassCiCSelection(fixPh1st[iPair], theVtx[iPair], fTracks,
                                            fElectrons, fPV, rho, leadptcut, fApplyEleVeto);
      if (pass1)
        pass2 = PhotonTools::PassCiCSelection(fixPh2nd[iPair], theVtx[iPair], fTracks,
                                              fElectrons, fPV, rho, trailptcut, fApplyEleVeto);
      break;
    case kMVAPhSelection://MVA
      pass1 = fixPh1st[iPair]->Pt()>leadptcut                                              &&
	PhotonTools::PassSinglePhotonPresel(fixPh1st[iPair],fElectrons,fConversions,bsp,
					    fTracks,theVtx[iPair],rho,fApplyEleVeto)       &&
	fTool.PassMVASelection(fixPh1st[iPair],theVtx[iPair],fTracks,fPV,rho,
			       fbdtCutBarrel,fbdtCutEndcap, fElectrons, fApplyEleVeto);
      if (pass1)
	pass2 = fixPh2nd[iPair]->Pt() > trailptcut                                         &&
	  PhotonTools::PassSinglePhotonPresel(fixPh2nd[iPair],fElectrons,fConversions,bsp,
					      fTracks,theVtx[iPair],rho,fApplyEleVeto)     &&
	  fTool.PassMVASelection(fixPh2nd[iPair],theVtx[iPair],fTracks,fPV,rho,
				 fbdtCutBarrel,fbdtCutEndcap, fElectrons, fApplyEleVeto);

      break;
    case kMITPhSelection:
      // loose preselection for mva
      pass1 = fixPh1st[iPair]->Pt() > leadptcut &&
	PhotonTools::PassSinglePhotonPresel(fixPh1st[iPair],fElectrons,fConversions,bsp,
					    fTracks,theVtx[iPair],rho,fApplyEleVeto,
					    fInvertElectronVeto);
      if (pass1)
	pass2 = fixPh2nd[iPair]->Pt() > trailptcut &&
	  PhotonTools::PassSinglePhotonPresel(fixPh2nd[iPair],fElectrons,fConversions,bsp,
					      fTracks,theVtx[iPair],rho,fApplyEleVeto,
					      fInvertElectronVeto);

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
    // finally, if both Photons pass the selections, add the pair to the passing pairs
    if (pass1 && pass2)
      passPairs.push_back(iPair);
  }


  // ---------------------------------------------------------------
  // ... we're almost done, stay focused...
  // loop over all passing pairs and find the one with the highest sum Et
  const Vertex*  vtx     = NULL;
  Photon*        phHard  = NULL;
  Photon*        phSoft  = NULL;

  PhotonTools::CiCBaseLineCats catPh1 = PhotonTools::kCiCNoCat;
  PhotonTools::CiCBaseLineCats catPh2 = PhotonTools::kCiCNoCat;

  double maxSumEt = 0.;
  for (unsigned int iPair=0; iPair<passPairs.size(); ++iPair) {
    double sumEt  = fixPh1st[passPairs[iPair]]->Et()
      +             fixPh2nd[passPairs[iPair]]->Et();
    if (sumEt > maxSumEt) {
      maxSumEt = sumEt;
      phHard   = fixPh1st[passPairs[iPair]];
      phSoft   = fixPh2nd[passPairs[iPair]];
      catPh1   = cat1st  [passPairs[iPair]];
      catPh2   = cat2nd  [passPairs[iPair]];
      vtx      = theVtx[iPair];
    }
  }

  for(unsigned int iPair = 0; iPair < numPairs; ++iPair) {
    if (fixPh1st[iPair]!=phHard) delete fixPh1st[iPair];
    if (fixPh2nd[iPair]!=phSoft) delete fixPh2nd[iPair];
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
  //delete[] theVtx;

  return;

  return;
}

//---------------------------------------------------------------------------------------------------
void PhotonPairSelector::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here, we just request the
  // photon collection branch.

  // load all branches
  ReqEventObject(fPhotonBranchName,   fPhotons,      fPhotonsFromBranch);
  ReqEventObject(fTrackBranchName,    fTracks,       true);
  ReqEventObject(fElectronName,       fElectrons,    true);
  ReqEventObject(fGoodElectronName,   fGoodElectrons,fGoodElectronsFromBranch);
  ReqEventObject(fPileUpDenName,      fPileUpDen,    true);
  ReqEventObject(fPVName,             fPV,           fPVFromBranch);
  ReqEventObject(fConversionName,     fConversions,  true);
  ReqEventObject(fBeamspotName,       fBeamspot,     true);
  ReqEventObject(fPFCandName,         fPFCands,      true);
  if (!fIsData) {
    ReqBranch(fPileUpName,            fPileUp);
    ReqBranch(fMCParticleName,        fMCParticles);
  }
  // determine photon selection type
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
  else if (fVertexSelType.CompareTo("CiCMVASelection") == 0)
    fVtxSelType =       kCiCMVAVtxSelection;
  else if (fVertexSelType.CompareTo("ZeroVtxSelection") == 0)
    fVtxSelType =       kStdVtxSelection;
  else {
    std::cerr<<" Vertex Seclection "<<fVertexSelType<<" not implemented."<<std::endl;
    return;
  }

  if (fIsData)
    fPhFixFile = gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFixGRPV22.dat");
  else
    fPhFixFile = gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFixSTART42V13.dat");

  printf("initialize photon pair selector\n");

  fTool.InitializeMVA(fVariableType,fEndcapWeights,fBarrelWeights);
  fVtxTools.InitP();

}

// ----------------------------------------------------------------------------------------
// some helpfer functions....
void PhotonPairSelector::FindHiggsPtAndZ(Float_t& pt, Float_t& decayZ, Float_t& mass)
{
  pt     = -999.;
  decayZ = -999.;
  mass   = -999.;

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

//---------------------------------------------------------------------------------------------------
Int_t PhotonPairSelector::FindRunRangeIdx(UInt_t run)
{
  // this routine looks for the idx of the run-range
  Int_t runRange=-1;
  for (UInt_t iRun = 0; iRun<fRunStart.size(); ++iRun) {
    if (run >= fRunStart[iRun] && run <= fRunEnd[iRun]) {
      runRange = (Int_t) iRun;
      return runRange;
    }
  }
  return runRange;
}

//---------------------------------------------------------------------------------------------------
Double_t PhotonPairSelector::GetDataEnCorr(Int_t runRange, PhotonTools::eScaleCats cat)
{
  switch (cat) {
  case PhotonTools::kEBhighEtaGold:
    return fDataEnCorr_EBhighEta_hR9[runRange];
  case PhotonTools::kEBhighEtaBad:
    return fDataEnCorr_EBhighEta_lR9[runRange];
  case PhotonTools::kEBlowEtaGoldCenter:
    return fDataEnCorr_EBlowEta_hR9central[runRange];
  case PhotonTools::kEBlowEtaGoldGap:
    return fDataEnCorr_EBlowEta_hR9gap[runRange];
  case PhotonTools::kEBlowEtaBad:
    return fDataEnCorr_EBlowEta_lR9[runRange];
  case PhotonTools::kEEhighEtaGold:
    return fDataEnCorr_EEhighEta_hR9[runRange];
  case PhotonTools::kEEhighEtaBad:
    return fDataEnCorr_EEhighEta_lR9[runRange];
  case PhotonTools::kEElowEtaGold:
    return fDataEnCorr_EElowEta_hR9[runRange];
  case PhotonTools::kEElowEtaBad:
    return fDataEnCorr_EElowEta_lR9[runRange];
  default:
    return 1.;
  }
}

//---------------------------------------------------------------------------------------------------
Double_t PhotonPairSelector::GetMCSmearFac(PhotonTools::eScaleCats cat)
{
  switch (cat) {
  case PhotonTools::kEBhighEtaGold:
    return fMCSmear_EBhighEta_hR9;
  case PhotonTools::kEBhighEtaBad:
    return fMCSmear_EBhighEta_lR9;
  case PhotonTools::kEBlowEtaGoldCenter:
    return fMCSmear_EBlowEta_hR9central;
  case PhotonTools::kEBlowEtaGoldGap:
    return fMCSmear_EBlowEta_hR9gap;
  case PhotonTools::kEBlowEtaBad:
    return fMCSmear_EBlowEta_lR9;
  case PhotonTools::kEEhighEtaGold:
    return fMCSmear_EEhighEta_hR9;
  case PhotonTools::kEEhighEtaBad:
    return fMCSmear_EEhighEta_lR9;
  case PhotonTools::kEElowEtaGold:
    return fMCSmear_EElowEta_hR9;
  case PhotonTools::kEElowEtaBad:
    return fMCSmear_EElowEta_lR9;
  default:
    return 1.;
  }
}

//---------------------------------------------------------------------------------------------------
Float_t PhotonPairSelector::GetEventCat(PhotonTools::CiCBaseLineCats cat1,
                                        PhotonTools::CiCBaseLineCats cat2)
{
  bool ph1IsEB  = (cat1 ==  PhotonTools::kCiCCat1 || cat1 == PhotonTools::kCiCCat2);
  bool ph2IsEB  = (cat2 ==  PhotonTools::kCiCCat1 || cat2 == PhotonTools::kCiCCat2);

  bool ph1IsHR9 = (cat1 ==  PhotonTools::kCiCCat1 || cat1 == PhotonTools::kCiCCat3);
  bool ph2IsHR9 = (cat2 ==  PhotonTools::kCiCCat1 || cat2 == PhotonTools::kCiCCat3);

  if (ph1IsEB && ph2IsEB)
    return (ph1IsHR9 && ph2IsHR9 ? 0. : 1.);

  return (ph1IsHR9 && ph2IsHR9 ? 2. : 3.);
}

