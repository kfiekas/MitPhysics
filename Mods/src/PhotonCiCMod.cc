// $Id: PhotonCiCMod.cc,v 1.1 2011/06/01 18:11:52 fabstoec Exp $

#include "MitPhysics/Mods/interface/PhotonCiCMod.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"

using namespace mithep;

ClassImp(mithep::PhotonCiCMod)

//--------------------------------------------------------------------------------------------------
PhotonCiCMod::PhotonCiCMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPhotonBranchName  (Names::gkPhotonBrn),
  fGoodPhotonsName   (ModNames::gkGoodPhotonsName),
  fTrackBranchName   (Names::gkTrackBrn),
  fPileUpDenName     (Names::gkPileupEnergyDensityBrn),
  fElectronName      ("Electrons"),
  fPhotonPtMin(20.0),
  fApplySpikeRemoval(kFALSE),
  fAbsEtaMax(999.99),
  fPhotons(0),
  fTracks(0),
  fPileUpDen(0),
  fElectrons(0),
  
  fPVName("PrimaryVertexes"),
  fPV(0),
  fPVFromBranch(kTRUE)

{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void PhotonCiCMod::Process()
{
  // Process entries of the tree. 

  LoadEventObject(fPhotonBranchName,   fPhotons);

  PhotonOArr *GoodPhotons = new PhotonOArr;
  GoodPhotons->SetName(fGoodPhotonsName);
  
  
  Double_t _tRho = -1.;
  if (fPhotons->GetEntries()>0) {
    LoadEventObject(fTrackBranchName,    fTracks);
    LoadEventObject(fPileUpDenName,      fPileUpDen);
    LoadEventObject(fElectronName,       fElectrons);
    LoadEventObject(fPVName,             fPV);    

    if(fPileUpDen->GetEntries() > 0)
      _tRho = (Double_t) fPileUpDen->At(0)->Rho();    
  
  } else {
    AddObjThisEvt(GoodPhotons);
    return;
  }
    
  

  std::vector<const Photon*> phCat1;
  std::vector<const Photon*> phCat2;
  std::vector<const Photon*> phCat3;
  std::vector<const Photon*> phCat4;


  const BaseVertex* theVtx = NULL;
  if(fPV->GetEntries() > 0)
    theVtx = fPV->At(0);
  else {
    AddObjThisEvt(GoodPhotons);
    return;
  }

  float cic4_allcuts_temp_sublead[] = { 
    3.77459,     2.18305,     1.76811,     1.30029,
    11.5519,     3.48306,     3.87394,     1.89038,
    3.47103,      2.1822,     2.26931,     1.43769,
    0.0105631,  0.00974116,   0.0282572,   0.0265096,
    0.0834648,   0.0821447,   0.0648449,   0.0476437,
    0.94,    0.355935,    0.94,    0.316358,
    98.9979,     1.94497,     96.2292,     96.2294,
    1.5,         1.5,         1.5,         1.5 };   // THIS IS ????

  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {    
    const Photon *ph = fPhotons->At(i);        

    if (ph->Pt() <= fPhotonPtMin) continue;
    if (ph->AbsEta() >= fAbsEtaMax) continue;
    
    if (  ph->SCluster()->AbsEta()>=2.5 || (ph->SCluster()->AbsEta()>=1.4442 && ph->SCluster()->AbsEta()<=1.566)  ) 
      continue;
    
    Bool_t isbarrel = ph->SCluster()->AbsEta()<1.5;
    
    Bool_t passSpikeRemovalFilter = kTRUE;
    
    if (ph->SCluster() && ph->SCluster()->Seed()) {
      if(ph->SCluster()->Seed()->Energy() > 5.0 && 
         ph->SCluster()->Seed()->EMax() / ph->SCluster()->Seed()->E3x3() > 0.95 )
        passSpikeRemovalFilter = kFALSE;
    }
    
    if (fApplySpikeRemoval && !passSpikeRemovalFilter) continue;
    
    // compute all relevant observables first
    double ecalIso = ph->EcalRecHitIsoDr03();
    double hcalIso = ph->HcalTowerSumEtDr03();
    double trackIso1 = IsolationTools::CiCTrackIsolation(ph, theVtx, 0.3, 0.02, 0.0, 0.0, 0.1, 1.0, fTracks, NULL);
    double trackIso2 = IsolationTools::CiCTrackIsolation(ph, theVtx, 0.3, 0.02, 0.0, 0.0, 0.1, 1.0, fTracks, fPV);
    double trackIso3 = IsolationTools::CiCTrackIsolation(ph, theVtx, 0.4, 0.02, 0.0, 0.0, 0.1, 1.0, fTracks, NULL);

    double combIso1 = ecalIso+hcalIso+trackIso1 - 0.17*_tRho;
    double combIso2 = ecalIso+hcalIso+trackIso2 - 0.52*_tRho;
    
    double tIso1 = (combIso1)*50./ph->Et();
    double tIso2 = (combIso2)*50./ph->Et();
    double tIso3 = (trackIso3)*50./ph->Et();    
    
    double covIEtaIEta  =ph->CoviEtaiEta();
    double HoE = ph->HadOverEm();
    
    double R9 = ph->R9();
    
    double dRTrack = PhotonTools::ElectronVetoCiC(ph, fElectrons);

    // stupid hard coded, will update later...
    if( tIso1                       < cic4_allcuts_temp_sublead[0+0*4]  && 
	tIso2                       < cic4_allcuts_temp_sublead[0+1*4]  && 
	tIso3                       < cic4_allcuts_temp_sublead[0+2*4]  && 
	covIEtaIEta                 < cic4_allcuts_temp_sublead[0+3*4]  && 
	HoE                         < cic4_allcuts_temp_sublead[0+4*4]  && 
	R9                          > cic4_allcuts_temp_sublead[0+5*4]  && 
	dRTrack*100.                > cic4_allcuts_temp_sublead[0+6*4]  && 
	isbarrel ) phCat1.push_back(ph);

    if( tIso1                       < cic4_allcuts_temp_sublead[1+0*4]  && 
	tIso2                       < cic4_allcuts_temp_sublead[1+1*4]  && 
	tIso3                       < cic4_allcuts_temp_sublead[1+2*4]  && 
	covIEtaIEta                 < cic4_allcuts_temp_sublead[1+3*4]  && 
	HoE                         < cic4_allcuts_temp_sublead[1+4*4]  && 
	R9                          > cic4_allcuts_temp_sublead[1+5*4]  && 
	dRTrack*100.                > cic4_allcuts_temp_sublead[1+6*4]  && 
	isbarrel ) phCat2.push_back(ph);

    if( tIso1                       < cic4_allcuts_temp_sublead[2+0*4]  && 
	tIso2                       < cic4_allcuts_temp_sublead[2+1*4]  && 
	tIso3                       < cic4_allcuts_temp_sublead[2+2*4]  && 
	covIEtaIEta                 < cic4_allcuts_temp_sublead[2+3*4]  && 
	HoE                         < cic4_allcuts_temp_sublead[2+4*4]  && 
	R9                          > cic4_allcuts_temp_sublead[2+5*4]  && 
	dRTrack*100.                > cic4_allcuts_temp_sublead[2+6*4]     )  phCat3.push_back(ph);

    if( tIso1                       < cic4_allcuts_temp_sublead[3+0*4]  && 
	tIso2                       < cic4_allcuts_temp_sublead[3+1*4]  && 
	tIso3                       < cic4_allcuts_temp_sublead[3+2*4]  && 
	covIEtaIEta                 < cic4_allcuts_temp_sublead[3+3*4]  && 
	HoE                         < cic4_allcuts_temp_sublead[3+4*4]  && 
	R9                          > cic4_allcuts_temp_sublead[3+5*4]  && 
	dRTrack*100.                > cic4_allcuts_temp_sublead[3+6*4]     )  phCat4.push_back(ph);

  }

  // check which category has two vgalid passes
  bool foundCat=false;
  if(phCat1.size() > 1) {
    // pick the first two guys, assuming they are sorted...
    GoodPhotons->Add(phCat1[0]);
    GoodPhotons->Add(phCat1[1]);
  } else if(phCat2.size() > 1) {
    GoodPhotons->Add(phCat2[0]);
    GoodPhotons->Add(phCat2[1]);
  } else {
    if(phCat3.size() > 1) {
      // one must be a Endcap photons... (stupid, but that's how it is in CiC)
      if(phCat3[0]->SCluster()->AbsEta() > 1.5 || phCat3[1]->SCluster()->AbsEta() > 1.5) {
	GoodPhotons->Add(phCat3[0]);
	GoodPhotons->Add(phCat3[1]);
	foundCat=true;
      }
      if (!foundCat) {
	if(phCat4.size() > 1) {
	  if(phCat4[0]->SCluster()->AbsEta() > 1.5 || phCat4[1]->SCluster()->AbsEta() > 1.5) {
	    GoodPhotons->Add(phCat4[0]);
	    GoodPhotons->Add(phCat4[1]);
	  }
	}
      }
    }
  }

  // sort according to pt
  GoodPhotons->Sort();

  // add to event for other modules to use
  AddObjThisEvt(GoodPhotons);  
}

//--------------------------------------------------------------------------------------------------
void PhotonCiCMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the photon collection branch.

  ReqEventObject(fPhotonBranchName,   fPhotons,   kTRUE);
  ReqEventObject(fTrackBranchName,    fTracks,    kTRUE);
  ReqEventObject(fElectronName,       fElectrons, kTRUE);  
  ReqEventObject(fPileUpDenName,      fPileUpDen, kTRUE);
  ReqEventObject(fPVName,             fPV,        fPVFromBranch);

}
