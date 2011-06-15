// $Id: PhotonCiCMod.cc,v 1.2 2011/06/02 14:18:07 bendavid Exp $

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

  int numVertices = fPV->GetEntries();
  if(numVertices == 0) {
    AddObjThisEvt(GoodPhotons);
    return;
  }    
  
  PhotonOArr* preselPh  = new PhotonOArr;
  PhotonOArr* failingPh = new PhotonOArr;
  
  const Vertex** vtxRanked = new const Vertex*[numVertices];
  
  // 1. we do the pre-selection; but keep the non-passing photons in a secont container...
  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {    
    const Photon *ph = fPhotons->At(i);
    if(ph->AbsEta() > 2.5) {
      failingPh->Add(ph);
      continue;
    }
    if(ph->Et() < 30.) {
      failingPh->Add(ph);
      continue;
    }
    if(ph->HadOverEm() > 0.15) {
      failingPh->Add(ph);
      continue;
    }
    if(ph->AbsEta() < 1.5) {
      if(ph->CoviEtaiEta() > 0.013) {
	failingPh->Add(ph);
	continue;
      }
      if(ph->EcalRecHitIsoDr03() > 10.) {
	failingPh->Add(ph);
	continue;
      }
    } else {
      if(ph->CoviEtaiEta() > 0.03) {
	failingPh->Add(ph);
	continue;
      }
      if(ph->EcalRecHitIsoDr03() > 10.) {
	failingPh->Add(ph);
	continue;
      }
    }
    preselPh->Add(ph);
  }

  // sort both by pt
  preselPh->Sort();
  failingPh->Sort();
  
  // if we don't have two photons in total, nothing to do...
  if( (preselPh->GetEntries() + failingPh->GetEntries() ) < 2) {
    AddObjThisEvt(GoodPhotons);
    delete preselPh;
    delete failingPh;
    return;
  }

  // now we assign the two hardest photons to be the ones to select the Vertex...
  double ptggV[2];
  ptggV[0] = 0.;
  ptggV[1] = 0.;  
  for(unsigned int iPh =0; iPh < preselPh->GetEntries() && iPh < 2 ; ++iPh) {
    ptggV[0] += preselPh->At(iPh)->Px();
    ptggV[1] += preselPh->At(iPh)->Px();
  }
  for(unsigned int iPh =preselPh->GetEntries(); iPh < (preselPh->GetEntries()+failingPh->GetEntries()) && iPh < 2 ; ++iPh) {
    ptggV[0] += failingPh->At(iPh)->Px();
    ptggV[1] += failingPh->At(iPh)->Px();
  }  
  double ptgg = TMath::Sqrt(ptggV[0]*ptggV[0]+ptggV[1]+ptggV[1]);
  
  // loop over all vertices and assigne the ranks
  int* ptbal_rank  = new int[fPV->GetEntries()];
  int* ptasym_rank = new int[fPV->GetEntries()];
  int* total_rank  = new int[fPV->GetEntries()];
  double* ptbal = new double[fPV->GetEntries()];
  double* ptasym = new double[fPV->GetEntries()];

  for(int iVtx = 0; iVtx < numVertices; ++iVtx) {
    const Vertex* _tVtx = fPV->At(iVtx);
    ptbal [iVtx] = 0.0;
    ptasym[iVtx] = 0.0;
    ptbal_rank [iVtx] = 1;
    ptasym_rank[iVtx] = 1;
    double ptvtx = 0.;
    for(unsigned int iTrk = 0; iTrk < _tVtx->NTracks(); ++iTrk) {
      const Track* _tTrk = _tVtx->Trk(iTrk);
      ptvtx += _tTrk->Pt();
      ptbal[iVtx]+= (_tTrk->Px()*ptggV[0])+(_tTrk->Py()*ptggV[1]);
    }
    ptbal [iVtx] = -ptbal[iVtx]/ptgg;
    ptasym[iVtx] = (ptvtx - ptgg)*(ptvtx+ptgg);
    for(int cVtx =0; cVtx < iVtx; ++cVtx) {
      if(ptbal [iVtx] > ptbal [cVtx])
	ptbal_rank [cVtx]++;
      else
	ptbal_rank [iVtx]++;
      if(ptasym [iVtx] > ptasym [cVtx])
	ptasym_rank [cVtx]++;
      else
	ptasym_rank [iVtx]++;
    }
  }

  // compute the total rank
  for(int iVtx = 0; iVtx < numVertices; ++iVtx) {
    ptbal_rank [iVtx] = ptbal_rank [iVtx]*ptasym_rank [iVtx]*(iVtx+1);  
    total_rank [iVtx] = 0;
    for(int cVtx =0; cVtx < iVtx; ++cVtx) {
      if(ptbal_rank [iVtx] > ptbal_rank [cVtx])
	total_rank[iVtx]++;
      else
	total_rank[cVtx]++;
    }
  }

  // sort the Vtx by total rank
  for(int iVtx = 0; iVtx < numVertices; ++iVtx)
    vtxRanked[total_rank[iVtx]] = fPV->At(iVtx);
  
  // pick the higest ranked vertex
  const BaseVertex* theVtx = vtxRanked[0];
  
  // delete all utility arrays
  delete ptbal_rank  ;
  delete ptasym_rank ;
  delete total_rank  ;
  delete ptbal       ;
  delete ptasym      ;

  // fix all the photon kinematics using the new Vtx
  PhotonOArr *fixedPhotons = new PhotonOArr;
  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {    
    const Photon *ph = fPhotons->At(i);
    // first copy this guy
    Photon* fPh = new Photon(*ph);
    const SuperCluster* sc = ph->SCluster();
    // compute the new Photon momentum (assuming zero mass)
    ThreeVectorC tPh_3mom  = sc->Point() - theVtx->Position();
    tPh_3mom = tPh_3mom/tPh_3mom.R();
    fPh->SetMom(sc->Energy()*tPh_3mom.X(), sc->Energy()*tPh_3mom.Y(), sc->Energy()*tPh_3mom.Z(), sc->Energy());
    fixedPhotons->Add(fPh);
  }
  
  // these values are taken from the H2GGlobe code... (actually from Marco/s mail)
  float cic4_allcuts_temp_sublead[] = { 
    3.8,         2.2,         1.77,        1.29,
    11.7,        3.4,         3.9,         1.84,
    3.5,         2.2,         2.3,         1.45,
    0.0106,      0.0097,      0.028,       0.027,
    0.082,       0.062,       0.065,       0.048,
    0.94,        0.36,        0.94,        0.32,
    1.,          0.062,       0.97,        0.97,
    1.5,         1.5,         1.5,         1.5 };  // the last line is POixelmatchVeto and un-used
  
  for (UInt_t i=0; i<fixedPhotons->GetEntries(); ++i) {    
    const Photon *ph = fixedPhotons->At(i);        

    // cut on Et instead of Pt???
    if ( ph->Pt() <= fPhotonPtMin   ) continue;
    if ( ph->AbsEta() >= fAbsEtaMax ) continue;
    
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
    
    double tIso1 = (combIso1) *50./ph->Et();
    double tIso2 = (combIso2) *50./ph->Et();
    double tIso3 = (trackIso3)*50./ph->Et();    
    
    double covIEtaIEta  =ph->CoviEtaiEta();
    double HoE = ph->HadOverEm();
    
    double R9 = ph->R9();
    
    double dRTrack = PhotonTools::ElectronVetoCiC(ph, fElectrons);

    // check which category it is ...
    int _tCat = 1;
    if ( !isbarrel ) _tCat = 3;
    if ( R9 < 0.94 ) _tCat++;
    
    if( tIso1                       < cic4_allcuts_temp_sublead[_tCat-1+0*4]  && 
	tIso2                       < cic4_allcuts_temp_sublead[_tCat-1+1*4]  && 
	tIso3                       < cic4_allcuts_temp_sublead[_tCat-1+2*4]  && 
	covIEtaIEta                 < cic4_allcuts_temp_sublead[_tCat-1+3*4]  && 
	HoE                         < cic4_allcuts_temp_sublead[_tCat-1+4*4]  && 
	R9                          > cic4_allcuts_temp_sublead[_tCat-1+5*4]  && 
	dRTrack*100.                > cic4_allcuts_temp_sublead[_tCat-1+6*4]  && 
	isbarrel ) GoodPhotons->Add(ph);
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

