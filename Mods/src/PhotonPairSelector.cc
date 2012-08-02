#include "MitPhysics/Mods/interface/PhotonPairSelector.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/StableParticle.h"
#include "MitAna/DataTree/interface/DecayParticleFwd.h"
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
#include "Math/SMatrix.h"
#include "Math/SVector.h"

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
  fPFConversionName              ("PFPhotonConversions"),  
  fTrackBranchName               (Names::gkTrackBrn),
  fPileUpDenName                 (Names::gkPileupEnergyDensityBrn),
  fPVName                        (Names::gkPVBeamSpotBrn),
  fBeamspotName                  (Names::gkBeamSpotBrn),
  fPFCandName                    (Names::gkPFCandidatesBrn),
  // MC specific stuff...
  fMCParticleName                (Names::gkMCPartBrn),
  fPileUpName                    (Names::gkPileupInfoBrn),
  fJetsName                      (Names::gkPFJetBrn),
  fPFMetName                     ("PFMet"),    
  fGoodPhotonsName               (ModNames::gkGoodPhotonsName),
  fChosenVtxName                 ("HggChosenVtx"),
  // ----------------------------------------
  // Selection Types
  fPhotonSelType                 ("NoSelection"),
  fVertexSelType                 ("StdSelection"),
  fPhSelType                     (kNoPhSelection),
  fVtxSelType                    (kStdVtxSelection),
  //-----------------------------------------
  // Id Types fab: shouldn't we have kNone as default ???
  fIdMVAType                     ("2011IdMVA"),
  fIdType                        (MVATools::k2011IdMVA),
  //-----------------------------------------
  // preselection Type
  fShowerShapeType               ("2011ShowerShape"),
  fSSType                        (PhotonTools::k2011ShowerShape),
  // ----------------------------------------
  fPhotonPtMin                   (20.0),
  fPhotonEtaMax                  (2.5),
  fLeadingPtMin                  (100.0/3.0),
  fTrailingPtMin                 (100.0/4.0),
  fIsData                        (false),
  fPhotonsFromBranch             (true),
  fPVFromBranch                  (true),
  fGoodElectronsFromBranch       (kTRUE),
  fUseSingleLegConversions       (kFALSE),
  // ----------------------------------------
  // collections....
  fPhotons                       (0),
  fElectrons                     (0),
  fConversions                   (0),
  fPFConversions                 (0),
  fTracks                        (0),
  fPileUpDen                     (0),
  fPV                            (0),
  fBeamspot                      (0),
  fPFCands                       (0),
  fMCParticles                   (0),
  fPileUp                        (0),
  fJets                          (0),
  fPFMet                         (0),
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
  // ------------------------------------------------------------------
  // this block should eventually be deleted... 
  fVariableType_2011             (10), 

  fEndcapWeights_2011            (gSystem->Getenv("CMSSW_BASE")+
				  TString("/src/MitPhysics/data/TMVAClassificationPhotonID_")+
				  TString("Endcap_PassPreSel_Variable_10_BDTnCuts2000_BDT.")+
				  TString("weights.xml")),
  fBarrelWeights_2011            (gSystem->Getenv("CMSSW_BASE")+
				  TString("/src/MitPhysics/data/TMVAClassificationPhotonID_")+
				  TString("Barrel_PassPreSel_Variable_10_BDTnCuts2000_BDT.")+
				  TString("weights.xml")),
  fVariableType_2012_globe       (1201),
  //fEndcapWeights_2012_globe      (gSystem->Getenv("CMSSW_BASE")+
  //TString("/src/MitPhysics/data/")+
  //				  TString("TMVA_EEpf_BDT_globe.")+
  //				  TString("weights.xml")),
  //fBarrelWeights_2012_globe      (gSystem->Getenv("CMSSW_BASE")+
  //				  TString("/src/MitPhysics/data/")+
  //				  TString("TMVA_EBpf_BDT_globe.")+
  //				  TString("weights.xml")),
  //fEndcapWeights_2012_globe      (gSystem->Getenv("CMSSW_BASE")+
  //				  TString("/src/MitPhysics/data/")+
  //				  TString("2012ICHEP_PhotonID_Endcap_BDT.")+
  //				  TString("weights.xml")),
  fEndcapWeights_2012_globe      (gSystem->Getenv("CMSSW_BASE")+
  				  TString("/src/MitPhysics/data/")+
  				  TString("2012ICHEP_PhotonID_Endcap_BDT.")+
  				  TString("weights_PSCorr.xml")),
  fBarrelWeights_2012_globe      (gSystem->Getenv("CMSSW_BASE")+
				  TString("/src/MitPhysics/data/")+
				  TString("2012ICHEP_PhotonID_Barrel_BDT.")+
				  TString("weights.xml")),
  // --------------------------------------------------------------------------

  fbdtCutBarrel                  (0.0744), //cuts give same eff (relative to presel) with cic
  fbdtCutEndcap                  (0.0959), //cuts give same eff (relative to presel) with cic

  // --------------------------------------------------------------------------

  fDoShowerShapeScaling          (kFALSE),
  //
  fMCErrScaleEB                  (1.0),
  fMCErrScaleEE                  (1.0),
  fRelativePtCuts                (kFALSE),

  fRhoType                       (RhoUtilities::CMS_RHO_RHOKT6PFJETS)

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


  VertexOArr  *ChosenVtx   = new VertexOArr;
  ChosenVtx->SetName(fChosenVtxName);
  ChosenVtx->SetOwner(kTRUE);

  // add to event for other modules to use
  AddObjThisEvt(GoodPhotons);
  AddObjThisEvt(ChosenVtx);

  //AddObjThisEvt(ChosenVtx);

  if (fPhotons->GetEntries()<2)
    return;

  LoadEventObject(fElectronName,       fElectrons);
  LoadEventObject(fGoodElectronName,   fGoodElectrons);
  LoadEventObject(fConversionName,     fConversions);
  if (fUseSingleLegConversions) LoadEventObject(fPFConversionName,     fPFConversions);
  LoadEventObject(fTrackBranchName,    fTracks);
  LoadEventObject(fPileUpDenName,      fPileUpDen);
  LoadEventObject(fPVName,             fPV);
  LoadEventObject(fBeamspotName,       fBeamspot);
  LoadEventObject(fPFCandName,         fPFCands);
  LoadEventObject(fJetsName,           fJets);
  LoadEventObject(fPFMetName,          fPFMet);

  if (!fIsData) {
    LoadBranch(fMCParticleName);
  }
  
  
  // ------------------------------------------------------------
  // load event based information
  Float_t rho = -99.;
  if (fPileUpDen->GetEntries() > 0)
    rho = (Double_t) fPileUpDen->At(0)->RhoRandomLowEta();
  Float_t rho2012 = -99;
  if (fPileUpDen->At(0)->RhoKt6PFJets()>0.) rho2012 = fPileUpDen->At(0)->RhoKt6PFJets();
  else rho2012 = fPileUpDen->At(0)->Rho();
  const BaseVertex *bsp = dynamic_cast<const BaseVertex*>(fBeamspot->At(0));

  Float_t theRho = rho;
  switch (fRhoType) {
  case RhoUtilities::CMS_RHO_RHOKT6PFJETS:
    theRho = rho2012;
    break;
  case RhoUtilities::MIT_RHO_VORONOI_LOW_ETA:       
    theRho = ( fPileUpDen->GetEntries() ? fPileUpDen->At(0)->RhoLowEta(): rho );
    break;
  case RhoUtilities::MIT_RHO_VORONOI_HIGH_ETA:
    theRho = ( fPileUpDen->GetEntries() ? fPileUpDen->At(0)->Rho() : rho );
    break;    
  case RhoUtilities::MIT_RHO_RANDOM_LOW_ETA:
    theRho = ( fPileUpDen->GetEntries() ? fPileUpDen->At(0)->RhoRandomLowEta() : rho );
    break;    
  case RhoUtilities::MIT_RHO_RANDOM_HIGH_ETA:       
    theRho = ( fPileUpDen->GetEntries() ? fPileUpDen->At(0)->RhoRandom() : rho );
    break;
  default:
    theRho = rho;
  }
    

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
      if (ph->CoviEtaiEta()     >  0.015)
        continue;
    }
    else {
      if (ph->CoviEtaiEta()     >  0.035)
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
  
  //fill conversion collection for vertex selection, adding single leg conversions if needed
  //note that momentum of single leg conversions needs to be recomputed from the track
  //as it is not filled properly
  DecayParticleOArr vtxconversions;
  if (fUseSingleLegConversions) {
    vtxconversions.SetOwner(kTRUE);
    for (UInt_t iconv=0; iconv<fConversions->GetEntries(); ++iconv) {
      DecayParticle *conv = new DecayParticle(*fConversions->At(iconv));
      vtxconversions.AddOwned(conv);
    }
    
    for (UInt_t iconv=0; iconv<fPFConversions->GetEntries(); ++iconv) {
      const DecayParticle *c = fPFConversions->At(iconv);
      if (c->NDaughters()!=1) continue;
      
      DecayParticle *conv = new DecayParticle(*c);
      const Track *trk = static_cast<const StableParticle*>(conv->Daughter(0))->Trk();
      conv->SetMom(trk->Px(), trk->Py(), trk->Pz(), trk->P());
      vtxconversions.AddOwned(conv);
    }    
  }
  else {
    for (UInt_t iconv=0; iconv<fConversions->GetEntries(); ++iconv) {
      const DecayParticle *c = fConversions->At(iconv);
      vtxconversions.Add(c);
    }
  }
  

  float higgspt = -99.;
  float higgsz = -99.;
  float higgsmass = -99.;
  
  if (!fIsData) FindHiggsPtAndZ(higgspt,higgsz,higgsmass);

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

    // do showershape scaling.... (fab: outsourced to PhotonTools)
    if (fDoShowerShapeScaling && !fIsData) {
      PhotonTools::ScalePhotonShowerShapes(fixPh1st[iPair],fSSType);
      PhotonTools::ScalePhotonShowerShapes(fixPh2nd[iPair],fSSType);
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
        else {
	  printf("Error: Run Range not found for data energy scale correction\n");
	  assert(0);
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
                                                    &vtxconversions,kFALSE,vtxProb);
      break;
      
    case kCiCMVAVtxSelection:
      theVtx[iPair] = fVtxTools.findVtxBasicRanking(fixPh1st[iPair],fixPh2nd[iPair], bsp, fPV,
                                                    &vtxconversions,kTRUE,vtxProb);
      break;
      
    case kMITVtxSelection:
      // need PFCandidate Collection
      theVtx[iPair] = VertexTools::BestVtx(fPFCands, fPV, bsp,
                                           mithep::FourVector((fixPh1st[iPair]->Mom()+
                                                               fixPh2nd[iPair]->Mom())));
      break;
      
    case kMetSigVtxSelection: {
      // need PFCandidate Collection, otherwise use 0
      if( !fJets || !fPFCands ) {
	theVtx[iPair] = fPV->At(0);
	break;
      }      
      
      PFJetOArr pfjets;
      for (UInt_t ijet=0; ijet<fJets->GetEntries(); ++ijet) {
	const PFJet *pfjet = dynamic_cast<const PFJet*>(fJets->At(ijet));
	if (pfjet && MathUtils::DeltaR(*pfjet,*fixPh1st[iPair])>0.3 && MathUtils::DeltaR(*pfjet,*fixPh2nd[iPair])>0.3) pfjets.Add(pfjet);
      }
      
      PFCandidateOArr pfcands;
      for (UInt_t icand=0; icand<fPFCands->GetEntries(); ++icand) {
	const PFCandidate *pfcand = fPFCands->At(icand);
	if (MathUtils::DeltaR(*pfcand,*fixPh1st[iPair])>0.05 && MathUtils::DeltaR(*pfcand,*fixPh2nd[iPair])>0.05) pfcands.Add(pfcand);
      }      
      
      double minsig = 1e6;
      for (UInt_t ivtx=0; ivtx<fPV->GetEntries(); ++ivtx) {
	const Vertex *v = fPV->At(ivtx);
	
	//         Met mmet = fMVAMet.GetMet(  false,
	//                                   fixPh1st[iPair]->Pt(),fixPh1st[iPair]->Phi(),fixPh1st[iPair]->Eta(),
	//                                   fixPh2nd[iPair]->Pt(),fixPh2nd[iPair]->Phi(),fixPh2nd[iPair]->Eta(),
	//                                   fPFMet->At(0),
	//                                   &pfcands,fPV->At(ivtx),fPV, fPileUpDen->At(0)->Rho(),
	//                                   &pfjets,
	//                                   int(fPV->GetEntries()),
	//                                   kTRUE);
	
        Met mmet = fMVAMet.GetMet(  false,
				    0.,0.,0.,
				    0.,0.,0.,
				    fPFMet->At(0),
				    &pfcands,fPV->At(ivtx),fPV, fPileUpDen->At(0)->Rho(),
				    &pfjets,
				    int(fPV->GetEntries()),
				    kFALSE);
	
	ThreeVector fullmet(mmet.Px() - fixPh1st[iPair]->Px() - fixPh2nd[iPair]->Px(),
			    mmet.Py() - fixPh1st[iPair]->Py() - fixPh2nd[iPair]->Py(),
			    0.);
	
	// 	ThreeVector fullmet(mmet.Px(),
	// 			    mmet.Py(),
	// 			    0.);	
	
	TMatrixD *metcov = fMVAMet.GetMetCovariance();
	
	//Double_t metsigma =  sqrt(fullmet.X()*fullmet.X()*(*metcov)(0,0) + 2.*fullmet.X()*fullmet.Y()*(*metcov)(1,0)  + fullmet.Y()*fullmet.Y()*(*metcov)(1,1))/fullmet.Rho(); 
	
	//Double_t metsig = fullmet.Rho()/metsigma;
	
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
	
	
	Double_t sumptsq = 0.;
	for (UInt_t ipfc = 0; ipfc<fPFCands->GetEntries(); ++ipfc) {
	  const PFCandidate *pfc = fPFCands->At(ipfc);
	  if (pfc->PFType()!=PFCandidate::eHadron || !pfc->HasTrackerTrk()) continue;
	  if (TMath::Abs( pfc->TrackerTrk()->DzCorrected(*v) ) > 0.2) continue;
	  if (TMath::Abs( pfc->TrackerTrk()->D0Corrected(*v) ) > 0.1) continue;	  
	  sumptsq += pfc->Pt()*pfc->Pt();
	}
	
	if ( sumptsq < 10.0 ) metsig = -99;
	
        if (metsig<minsig && metsig>0.) {
          minsig = metsig;
          theVtx[iPair] = fPV->At(ivtx);
        }
        
        printf("ivtx = %i, sumptsq = %5f, met = %5f, metsig = %5f, dzgen = %5f\n",ivtx,sumptsq, fullmet.Rho(), metsig,fPV->At(ivtx)->Z()-higgsz);
        
      }
      
      double testpfmetx = 0.;
      double testpfmety = 0.;
      for (UInt_t icand=0; icand<pfcands.GetEntries(); ++icand) {
	const PFCandidate *pfcand = pfcands.At(icand);
	testpfmetx -= pfcand->Px();
	testpfmety -= pfcand->Py();
	//pfcands.Add(pfcand);
      }            
      
      testpfmetx -= fixPh1st[iPair]->Px();
      testpfmetx -= fixPh2nd[iPair]->Px();
      
      testpfmety -= fixPh1st[iPair]->Py();
      testpfmety -= fixPh2nd[iPair]->Py();     
      
      double testpfmet = sqrt(testpfmetx*testpfmetx + testpfmety*testpfmety);
      
      printf("bestvtx: metsig = %5f, dzgen = %5f, diphopt = %5f, pfmet = %5f, testpfmet = %5f\n", minsig, theVtx[iPair]->Z()-higgsz, higgspt,fPFMet->At(0)->Pt(),testpfmet);
    
    }
      break;      
    default:
      theVtx[iPair] = fPV->At(0);
    }
    
    //set PV ref in photons
    fixPh1st[iPair]->SetPV(theVtx[iPair]);
    fixPh2nd[iPair]->SetPV(theVtx[iPair]);
    fixPh1st[iPair]->SetVtxProb(vtxProb);
    fixPh2nd[iPair]->SetVtxProb(vtxProb);
        
    ThreeVector vtxpos = theVtx[iPair]->Position();
    
    // fix the kinematics for both events
    FourVectorM newMom1st = fixPh1st[iPair]->MomVtx(vtxpos);
    FourVectorM newMom2nd = fixPh2nd[iPair]->MomVtx(vtxpos);
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
    Double_t bdt1 = -99.;
    Double_t bdt2 = -99.;
    // ---------------------------------------------------------------------------------------------------------------
    // using new interface letting the MVATools handle the type... (fab)
    if( fIdType != MVATools::kNone ) {    // not strictly needed, but cold spped up things slightly if no MVA is needed...
      bdt1 = fTool.GetMVAbdtValue(fixPh1st[iPair],theVtx[iPair],fTracks,fPV,theRho,fPFCands,fElectrons,fApplyEleVeto);
      bdt2 = fTool.GetMVAbdtValue(fixPh2nd[iPair],theVtx[iPair],fTracks,fPV,theRho,fPFCands,fElectrons,fApplyEleVeto);
    }
    fixPh1st[iPair]->SetIdMva(bdt1);
    fixPh2nd[iPair]->SetIdMva(bdt2);
    // ---------------------------------------------------------------------------------------------------------------
    // superseedded by above code... (fab)
    //     switch (fIdType) {
    //     case MVATools::k2011IdMVA:
    //       bdt1 = fTool.GetMVAbdtValue_2011(fixPh1st[iPair],theVtx[iPair],fTracks,fPV,rho,fElectrons,fApplyEleVeto);
    //       bdt2 = fTool.GetMVAbdtValue_2011(fixPh2nd[iPair],theVtx[iPair],fTracks,fPV,rho,fElectrons,fApplyEleVeto);
    //       fixPh1st[iPair]->SetIdMva(bdt1);
    //       fixPh2nd[iPair]->SetIdMva(bdt2);
    //       break;
    
    //     case MVATools::k2012IdMVA_globe:
    //       bdt1 = fTool.GetMVAbdtValue_2012_globe(fixPh1st[iPair],theVtx[iPair],fTracks,fPV,rho2012,fPFCands,fElectrons,fApplyEleVeto);
    //       bdt2 = fTool.GetMVAbdtValue_2012_globe(fixPh2nd[iPair],theVtx[iPair],fTracks,fPV,rho2012,fPFCands,fElectrons,fApplyEleVeto);
    //       fixPh1st[iPair]->SetIdMva(bdt1);
    //       fixPh2nd[iPair]->SetIdMva(bdt2);
    //       break;
    
    //     default:
    //       fixPh1st[iPair]->SetIdMva(-99.);
    //       fixPh2nd[iPair]->SetIdMva(-99.);
    //     }
    // ---------------------------------------------------------------------------------------------------------------    

    //printf("applying id\n");

    // check if both photons pass the CiC selection
    // FIX-ME: Add other possibilities....
    bool pass1 = false;
    bool pass2 = false;
    
    switch (fPhSelType) {

      // --------------------------------------------------------------------
      // trivial (no) selection with only pt-cuts
    case kNoPhSelection:
      pass1 = (fixPh1st[iPair]->Pt() > leadptcut );
      pass2 = (fixPh2nd[iPair]->Pt() > trailptcut);
      break;

      // --------------------------------------------------------------------
      // CiC4 Selection as used for the 2011 7TeV Baseline analysis
    case kCiCPhSelection:
      pass1 = PhotonTools::PassCiCSelection(fixPh1st[iPair], theVtx[iPair], fTracks,
                                            fElectrons, fPV, rho, leadptcut, fApplyEleVeto);
      if (pass1)
        pass2 = PhotonTools::PassCiCSelection(fixPh2nd[iPair], theVtx[iPair], fTracks,
                                              fElectrons, fPV, rho, trailptcut, fApplyEleVeto);
      break;

      // --------------------------------------------------------------------
      // PF-CiC4 Selection, as used in the 2012 8TeV Baseline analysis
    case kCiCPFPhSelection:
      
      // loose preselection for mva
      pass1 = fixPh1st[iPair]->Pt() > leadptcut &&
	PhotonTools::PassSinglePhotonPreselPFISO(fixPh1st[iPair],fElectrons,fConversions,bsp,
					    fTracks,theVtx[iPair],rho2012,fPFCands,fApplyEleVeto,
					    fInvertElectronVeto);
      if (pass1)
	pass2 = fixPh2nd[iPair]->Pt() > trailptcut &&
	  PhotonTools::PassSinglePhotonPreselPFISO(fixPh2nd[iPair],fElectrons,fConversions,bsp,
					      fTracks,theVtx[iPair],rho2012,fPFCands,fApplyEleVeto,
					      fInvertElectronVeto);      
      
      pass1 = pass1 && PhotonTools::PassCiCPFIsoSelection(fixPh1st[iPair], theVtx[iPair], fPFCands,
                                             fPV, rho2012, leadptcut);
      if (pass1)
        pass2 = pass2 && PhotonTools::PassCiCPFIsoSelection(fixPh2nd[iPair], theVtx[iPair], fPFCands,
                                               fPV, rho2012, trailptcut);
      break;      
      
      // --------------------------------------------------------------------
      // MVA selection
    case kMVAPhSelection://MVA

      pass1 = fixPh1st[iPair]->Pt()>leadptcut                                              &&
	PhotonTools::PassSinglePhotonPresel(fixPh1st[iPair],fElectrons,fConversions,bsp,
					    fTracks,theVtx[iPair],rho,fApplyEleVeto)       &&
	( ( fixPh1st[iPair]->SCluster()->AbsEta() < 1.5 && fixPh1st[iPair]->IdMva() > fbdtCutBarrel )
	  || ( fixPh1st[iPair]->IdMva() > fbdtCutEndcap ) );
	  
	// we've already compyted the MVA varible above, not needed to redo...
// 	fTool.PassMVASelection(fixPh1st[iPair],theVtx[iPair],fTracks,fPV,rho,
// 			       fbdtCutBarrel,fbdtCutEndcap, fElectrons, fApplyEleVeto);


      if (pass1)
	pass2 = fixPh2nd[iPair]->Pt() > trailptcut                                         &&
	  PhotonTools::PassSinglePhotonPresel(fixPh2nd[iPair],fElectrons,fConversions,bsp,
					      fTracks,theVtx[iPair],rho,fApplyEleVeto)     &&
	( ( fixPh2nd[iPair]->SCluster()->AbsEta() < 1.5 && fixPh2nd[iPair]->IdMva() > fbdtCutBarrel )
	  || ( fixPh2nd[iPair]->IdMva() > fbdtCutEndcap ) );

// 	  fTool.PassMVASelection(fixPh2nd[iPair],theVtx[iPair],fTracks,fPV,rho,
// 				 fbdtCutBarrel,fbdtCutEndcap, fElectrons, fApplyEleVeto);

      break;

      // --------------------------------------------------------------------
      // MIT proposed pre-selection as used in the 2011 7TeV MVA analyses
    case kMITPhSelection:
      // loose preselection for mva
      pass1 = fixPh1st[iPair]->Pt() > leadptcut &&
	PhotonTools::PassSinglePhotonPresel(fixPh1st[iPair],fElectrons,fConversions,bsp,
					    fTracks,theVtx[iPair],rho2012,fApplyEleVeto,
					    fInvertElectronVeto);
      if (pass1)
	pass2 = fixPh2nd[iPair]->Pt() > trailptcut &&
	  PhotonTools::PassSinglePhotonPresel(fixPh2nd[iPair],fElectrons,fConversions,bsp,
					      fTracks,theVtx[iPair],rho2012,fApplyEleVeto,
					      fInvertElectronVeto);

      break;
     
      // --------------------------------------------------------------------
      // updated (PF absed) pre-selection as used in the 2012 8TeV MVA/Baseline analyses
    case kMITPFPhSelection:
      // loose preselection for mva
      pass1 = fixPh1st[iPair]->Pt() > leadptcut &&
	PhotonTools::PassSinglePhotonPreselPFISO(fixPh1st[iPair],fElectrons,fConversions,bsp,
						 fTracks,theVtx[iPair],rho2012,fPFCands,fApplyEleVeto,
						 fInvertElectronVeto);      

      if (pass1)
	pass2 = fixPh2nd[iPair]->Pt() > trailptcut &&
	  PhotonTools::PassSinglePhotonPreselPFISO(fixPh2nd[iPair],fElectrons,fConversions,bsp,
						   fTracks,theVtx[iPair],rho2012,fPFCands,fApplyEleVeto,
						   fInvertElectronVeto);           
      
      pass1 = pass1 && PhotonTools::PassCiCPFIsoSelection(fixPh1st[iPair], theVtx[iPair], fPFCands,
							  fPV, rho2012, leadptcut);
      if (pass1)
        pass2 = pass2 && PhotonTools::PassCiCPFIsoSelection(fixPh2nd[iPair], theVtx[iPair], fPFCands,
							    fPV, rho2012, trailptcut);
      break;      
    default:
      pass1 = true;
      pass2 = true;
    }
    
    //match to good electrons if requested
    //     if (fInvertElectronVeto) {
    //       pass1 &= !PhotonTools::PassElectronVeto(fixPh1st[iPair],fGoodElectrons);
    //       pass2 &= !PhotonTools::PassElectronVeto(fixPh2nd[iPair],fGoodElectrons);
    //     }
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

  // not used at all....
  //PhotonTools::CiCBaseLineCats catPh1 = PhotonTools::kCiCNoCat;
  //PhotonTools::CiCBaseLineCats catPh2 = PhotonTools::kCiCNoCat;

  double maxSumEt = 0.;
  for (unsigned int iPair=0; iPair<passPairs.size(); ++iPair) {
    double sumEt  = fixPh1st[passPairs[iPair]]->Et()
      +             fixPh2nd[passPairs[iPair]]->Et();
    if (sumEt > maxSumEt) {
      maxSumEt = sumEt;
      phHard   = fixPh1st[passPairs[iPair]];
      phSoft   = fixPh2nd[passPairs[iPair]];
      //catPh1   = cat1st  [passPairs[iPair]];
      //catPh2   = cat2nd  [passPairs[iPair]];
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

  // we also store the chosen Vtx, so later modules can use it

  Vertex* chosenVtx = NULL;
  if ( vtx ) {
    chosenVtx = new Vertex( *vtx );
    ChosenVtx->AddOwned( chosenVtx );
  }
  
  // sort according to pt
  GoodPhotons->Sort();
  
  // delete auxiliary photon collection...
  delete preselPh;
  //delete[] theVtx;

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
  if (fUseSingleLegConversions) ReqEventObject(fPFConversionName,     fPFConversions,  true);
  ReqEventObject(fBeamspotName,       fBeamspot,     true);
  ReqEventObject(fPFCandName,         fPFCands,      true);
  ReqEventObject(fJetsName,         fJets,      false);  
  ReqEventObject(fPFMetName,         fPFMet,      true);    
  if (!fIsData) {
    //ReqBranch(fPileUpName,            fPileUp);
    ReqBranch(fMCParticleName,        fMCParticles);
  }

  // determine photon selection type
  if      (fPhotonSelType.CompareTo("CiCSelection") == 0)
    fPhSelType =       kCiCPhSelection;
  else if      (fPhotonSelType.CompareTo("CiCPFSelection") == 0)
    fPhSelType =       kCiCPFPhSelection;  
  else if (fPhotonSelType.CompareTo("MVASelection") == 0) //MVA
    fPhSelType =       kMVAPhSelection;
  else if (fPhotonSelType.CompareTo("MITSelection") == 0)
    fPhSelType =       kMITPhSelection;
  else if (fPhotonSelType.CompareTo("MITPFSelection") == 0)
    fPhSelType =       kMITPFPhSelection;  
  else if (fPhotonSelType.CompareTo("NoSelection") == 0)
    fPhSelType =       kNoPhSelection;
  else {
    std::cerr<<" Photon Seclection "<<fPhotonSelType<<" not implemented."<<std::endl;
    return;    
  }

  if      (fVertexSelType.CompareTo("CiCSelection") == 0) {
    fVtxSelType =       kCiCVtxSelection;
    fVtxTools.InitP(1);
  }
  else if (fVertexSelType.CompareTo("MITSelection") == 0)
    fVtxSelType =       kMITVtxSelection;
  else if (fVertexSelType.CompareTo("CiCMVASelection") == 0) {
    fVtxSelType =       kCiCMVAVtxSelection;
    fVtxTools.InitP(1);
  }
  else if (fVertexSelType.CompareTo("CiCMVA2012Selection") == 0) {
    fVtxSelType =       kCiCMVAVtxSelection;
    fVtxTools.InitP(2);
  }  
  else if (fVertexSelType.CompareTo("MetSigSelection") == 0)
    fVtxSelType =       kMetSigVtxSelection;  
  else if (fVertexSelType.CompareTo("ZeroVtxSelection") == 0)
    fVtxSelType =       kStdVtxSelection;
  else {
    std::cerr<<" Vertex Seclection "<<fVertexSelType<<" not implemented."<<std::endl;
    return;
  }
  
  if      (fIdMVAType.CompareTo("2011IdMVA") == 0)
    fIdType =       MVATools::k2011IdMVA;
  else if (fIdMVAType.CompareTo("2012IdMVA_globe") == 0)
    fIdType =       MVATools::k2012IdMVA_globe;
  else if (fIdMVAType.CompareTo("2012IdMVA") == 0)
    fIdType =       MVATools::k2012IdMVA;
  else if (fIdMVAType.CompareTo("2011IdMVA_HZg") == 0)
    fIdType =       MVATools::k2011IdMVA_HZg;
  else if (fIdMVAType.CompareTo("None") == 0)
    fIdType =       MVATools::kNone;
  else {
    std::cerr<<" Id MVA "<<fIdMVAType<<" not implemented."<<std::endl;
    return;
  }
  

  if      (fShowerShapeType.CompareTo("None")            == 0)
    fSSType =       PhotonTools::kNoShowerShapeScaling;
  else if (fShowerShapeType.CompareTo("2011ShowerShape") == 0)
    fSSType =       PhotonTools::k2011ShowerShape;
  else if (fShowerShapeType.CompareTo("2012ShowerShape") == 0)
    fSSType =       PhotonTools::k2012ShowerShape;
  else {
    std::cerr<<"shower shape scale "<<fShowerShapeType<<" not implemented."<<std::endl;
    return;
  }
  
  if (fIsData)
    fPhFixFile = gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFixGRPV22.dat");
  else
    fPhFixFile = gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFixSTART42V13.dat");

  printf("initialize photon pair selector\n");

  // ----------------------------------------------------------------
  // anopther replace block (fab):
  //   delegate the choice of the weight files to the MVATools...
  fTool.InitializeMVA(fIdType);

//   switch (fIdType) {
//   case k2011IdMVA:
//     fTool.InitializeMVA(fVariableType_2011,fEndcapWeights_2011,fBarrelWeights_2011);
//     break;
    
//   case k2012IdMVA_globe:
//     fTool.InitializeMVA(fVariableType_2012_globe,fEndcapWeights_2012_globe,fBarrelWeights_2012_globe); 
//     break;
//   }
  // ----------------------------------------------------------------
  
  if (fVtxSelType==kMetSigVtxSelection) {
    //fMVAMet.Initialize();
//     fMVAMet.Initialize(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_lowpt.weights.xml"))),
//                         TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_highpt.weights.xml"))),
//                         TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/Utils/python/JetIdParams_cfi.py")),
//                         TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmet_42.root"))),
//                         TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetphi_42.root"))),
//                         TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu1_42.root"))),
//                         TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu2_42.root")))
//                        );

    fMVAMet.Initialize(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_lowpt.weights.xml"))),
                        TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_highpt.weights.xml"))),
                        TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/Utils/python/JetIdParams_cfi.py")),
                        TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmet_52.root"))),
                        TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetphi_52.root"))),
                        TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu1cov_52.root"))),
                        TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu2cov_52.root")))
                        );

  }

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

