// $Id: IsolationTools.cc,v 1.24 2011/12/21 13:03:23 sixie Exp $

#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

ClassImp(mithep::IsolationTools)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::TrackIsolation(const Track *p, Double_t extRadius, Double_t intRadius,
                                        Double_t ptLow, Double_t maxVtxZDist, 
                                        const Collection<Track> *tracks) 
{
  //Computes the Track Isolation: Summed Transverse Momentum of all tracks inside an 
  //annulus around the electron seed track.  

  Double_t ptSum =0.;  
  for (UInt_t i=0; i<tracks->GetEntries();i++) {   
    Double_t tmpPt = tracks->At(i)->Pt();
    Double_t deltaZ = fabs(p->Z0() - tracks->At(i)->Z0());

    //ignore the track if it is below the pt threshold
    if (tmpPt < ptLow) 
      continue;    
    //ingore the track if it is too far away in Z
    if (deltaZ > maxVtxZDist) 
      continue;
           
    Double_t dr = MathUtils::DeltaR(p->Phi(),p->Eta(),tracks->At(i)->Phi(), tracks->At(i)->Eta());
    //add the track pt if it is inside the annulus
    if ( dr < extRadius && 
	 dr >= intRadius ) {
      ptSum += tmpPt;
    }
  }
  return ptSum;  
}

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::EcalIsolation(const SuperCluster *sc, Double_t coneSize, Double_t etLow, 
                                       const Collection<BasicCluster> *basicClusters) 
{
  //Computes the Ecal Isolation: Summed Transverse Energy of all Basic Clusters inside a  
  //cone around the electron, excluding those that are inside the electron super cluster.

  Double_t ecalIsol=0.;
  const BasicCluster *basicCluster= 0;
  for (UInt_t i=0; i<basicClusters->GetEntries();i++) {    
    basicCluster = basicClusters->At(i);    
    Double_t basicClusterEnergy    = basicCluster->Energy();
    Double_t basicClusterEta  = basicCluster->Eta();
    Double_t basicClusterEt   = basicClusterEnergy*sin(2*atan(exp(basicClusterEta)));           

    if (basicClusterEt > etLow) {            
      bool inSuperCluster = false;	  
      
      // loop over the basic clusters of the supercluster
      // to make sure that the basic cluster is not inside
      // the super cluster. We exclude those.
      for (UInt_t j=0; j<sc->ClusterSize(); j++) {
        const BasicCluster *tempBasicClusterInSuperCluster = sc->Cluster(j);	
        if (tempBasicClusterInSuperCluster == basicCluster) {
          inSuperCluster = true;	    
        }
      }
      
      if (!inSuperCluster) {	    
        Double_t dr = MathUtils::DeltaR(sc->Phi(), sc->Eta(),
                                      basicCluster->Phi(),basicCluster->Eta());
        if(dr < coneSize) {
          ecalIsol += basicClusterEt;
        }
      }
    }
  } 
  return ecalIsol;
}

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::CaloTowerHadIsolation(const ThreeVector *p, Double_t extRadius, 
                                               Double_t intRadius, Double_t etLow, 
                                               const Collection<CaloTower> *caloTowers) 
{
  //Computes the CaloTower Had Et Isolation: Summed Hadronic Transverse Energy of all Calo Towers 
  //inside an annulus around the electron super cluster position.

  Double_t sumEt = 0;
  for (UInt_t i=0; i<caloTowers->GetEntries();i++) {    
    Double_t caloTowerEt = caloTowers->At(i)->HadEt();
    Double_t dr = MathUtils::DeltaR(caloTowers->At(i)->Phi(), caloTowers->At(i)->Eta(),
                                  p->Phi(), p->Eta());
    if (dr < extRadius && dr > intRadius && caloTowerEt > etLow) {
      sumEt += caloTowerEt;
    }
  }
  return sumEt;
}

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::CaloTowerEmIsolation(const ThreeVector *p, Double_t extRadius, 
                                              Double_t intRadius, Double_t etLow, 
                                              const Collection<CaloTower> *caloTowers) 
{
  //Computes the CaloTower Em Et Isolation: Summed Hadronic Transverse Energy of all Calo Towers 
  //inside an annulus around the electron super cluster position.

  Double_t sumEt = 0;
  for (UInt_t i=0; i<caloTowers->GetEntries();i++) {    
    Double_t caloTowerEt = caloTowers->At(i)->EmEt();
    Double_t dr = MathUtils::DeltaR(caloTowers->At(i)->Phi(), caloTowers->At(i)->Eta(),
                                   p->Phi(), p->Eta());
    if (dr < extRadius && dr > intRadius && caloTowerEt > etLow) {
      sumEt += caloTowerEt;
    }
  }
  return sumEt;
}

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::PFMuonIsolation(const Muon *p, const PFCandidateCol *PFCands, 
                                      	 const Vertex *vertex, Double_t  delta_z, Double_t ptMin,
				     	 Double_t extRadius, Double_t intRadiusGamma, Double_t intRadius)
{
  //Computes the PF Isolation: Summed Transverse Momentum of all PF candidates inside an 
  //annulus around the particle seed track.  

  Double_t zLepton = 0.0;
  if(p->BestTrk()) zLepton = p->BestTrk()->DzCorrected(*vertex);

  Double_t ptSum =0.;  
  for (UInt_t i=0; i<PFCands->GetEntries();i++) {   
    const PFCandidate *pf = PFCands->At(i);

    // exclude muon
    if(pf->TrackerTrk() && p->TrackerTrk() &&
       pf->TrackerTrk() == p->TrackerTrk()) continue;

    Double_t dr = MathUtils::DeltaR(p->Mom(), pf->Mom());
    
    // pt cut applied to neutrals
    if(!pf->HasTrk() && pf->Pt() <= ptMin) continue;

    // ignore the pf candidate if it is too far away in Z
    if(pf->HasTrk()) {
      Double_t deltaZ = TMath::Abs(pf->BestTrk()->DzCorrected(*vertex) - zLepton);
      if (deltaZ >= delta_z) continue;
    }
      
    // inner cone veto for gammas
    if (pf->PFType() == PFCandidate::eGamma && dr < intRadiusGamma) continue;
    
    // inner cone veto for tracks
    if (dr < intRadius) continue;

    // add the pf pt if it is inside the extRadius 
    if (dr < extRadius) ptSum += pf->Pt();
  
  }
  return ptSum;
}

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::PFMuonIsolation(const Muon *p, const PFCandidateCol *PFCands,
                                      	 const MuonCol *goodMuons, const ElectronCol *goodElectrons, 
                                      	 const Vertex *vertex, Double_t  delta_z, Double_t ptMin,
				      	 Double_t extRadius, Double_t intRadiusGamma, Double_t intRadius)
{
  //Computes the PF Isolation: Summed Transverse Momentum of all PF candidates inside an 
  //annulus around the particle seed track.  

  Double_t zLepton = 0.0;
  if(p->BestTrk()) zLepton = p->BestTrk()->DzCorrected(*vertex);

  Double_t ptSum =0.;  
  for (UInt_t i=0; i<PFCands->GetEntries();i++) {   
    const PFCandidate *pf = PFCands->At(i);

    // exclude muon
    if(pf->TrackerTrk() && p->TrackerTrk() &&
       pf->TrackerTrk() == p->TrackerTrk()) continue;

    Double_t dr = MathUtils::DeltaR(p->Mom(), pf->Mom());
    
    // pt cut applied to neutrals
    if(!pf->HasTrk() && pf->Pt() <= ptMin) continue;

    // ignore the pf candidate if it is too far away in Z
    if(pf->HasTrk()) {
      Double_t deltaZ = TMath::Abs(pf->BestTrk()->DzCorrected(*vertex) - zLepton);
      if (deltaZ >= delta_z) continue;
    }
      
    // inner cone veto for gammas
    if (pf->PFType() == PFCandidate::eGamma && dr < intRadiusGamma) continue;
    
    // inner cone veto for tracks
    if (dr < intRadius) continue;

    // add the pf pt if it is inside the extRadius and outside the intRadius
    if (dr < extRadius ) {
      Bool_t isLepton = kFALSE;
      if(goodMuons){
        for (UInt_t nl=0; nl<goodMuons->GetEntries();nl++) {
	  const Muon *m = goodMuons->At(nl);
          if(pf->TrackerTrk() && m->TrackerTrk() &&
	     pf->TrackerTrk() == m->TrackerTrk()) {
	    isLepton = kTRUE;
	    break;
	  }
	}
      }
      if(goodElectrons && isLepton == kFALSE){
        for (UInt_t nl=0; nl<goodElectrons->GetEntries();nl++) {
	  const Electron *e = goodElectrons->At(nl);
          if(pf->TrackerTrk() && e->TrackerTrk() &&
	     pf->TrackerTrk() == e->TrackerTrk()) {
	    isLepton = kTRUE;
	    break;
	  }
          if(pf->GsfTrk() && e->GsfTrk() &&
	     pf->GsfTrk() == e->GsfTrk()) {
	    isLepton = kTRUE;
	    break;
	  }
	}
      }
      if (isLepton == kTRUE) continue;
      ptSum += pf->Pt();            
    }  
  }
  return ptSum;
}

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::PFElectronIsolation(const Electron *p, const PFCandidateCol *PFCands, 
                                      	     const Vertex *vertex, Double_t delta_z, Double_t ptMin,
				     	     Double_t extRadius, Double_t intRadius) 
{

  //Computes the PF Isolation: Summed Transverse Momentum of all PF candidates inside an 
  //annulus around the particle seed track.  

  Double_t zLepton = 0.0;
  if(p->BestTrk()) zLepton = p->BestTrk()->DzCorrected(*vertex);

  Double_t ptSum =0.;  
  for (UInt_t i=0; i<PFCands->GetEntries();i++) {   
    const PFCandidate *pf = PFCands->At(i);
    
    // pt cut applied to neutrals
    if(!pf->HasTrk() && pf->Pt() <= ptMin) continue;

    if(pf->TrackerTrk() && p->TrackerTrk() &&
       pf->TrackerTrk() == p->TrackerTrk()) continue;

    if(pf->GsfTrk() && p->GsfTrk() &&
       pf->GsfTrk() == p->GsfTrk()) continue;

    // ignore the pf candidate if it is too far away in Z
    if(pf->BestTrk()) {
      Double_t deltaZ = TMath::Abs(pf->BestTrk()->DzCorrected(*vertex) - zLepton);
      if (deltaZ >= delta_z) continue;
    }

           
    Double_t dr = MathUtils::DeltaR(p->Mom(), pf->Mom());
    // add the pf pt if it is inside the extRadius and outside the intRadius
    if ( dr < extRadius && 
	 dr >= intRadius ) {

      //EtaStrip Veto for Gamma 
      if (pf->PFType() == PFCandidate::eGamma && fabs(p->Eta() - pf->Eta()) < 0.025) continue;

      //InnerCone (One Tower = dR < 0.07) Veto for non-gamma neutrals
      if (!pf->HasTrk() && pf->PFType() == PFCandidate::eNeutralHadron
          && MathUtils::DeltaR(p->Mom(), pf->Mom()) < 0.07 ) continue; 

      ptSum += pf->Pt();            

    }
  }
  return ptSum;
}


//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::PFElectronIsolation(const Electron *p, const PFCandidateCol *PFCands, 
                      	     		     const MuonCol *goodMuons, const ElectronCol *goodElectrons,
					     const Vertex *vertex, Double_t  delta_z, Double_t ptMin,
			     		     Double_t extRadius, Double_t intRadius)
{

  //Computes the PF Isolation: Summed Transverse Momentum of all PF candidates inside an 
  //annulus around the particle seed track.  

  Double_t zLepton = 0.0;
  if(p->BestTrk()) zLepton = p->BestTrk()->DzCorrected(*vertex);

  Double_t ptSum =0.;  
  for (UInt_t i=0; i<PFCands->GetEntries();i++) {   
    const PFCandidate *pf = PFCands->At(i);
    
    // pt cut applied to neutrals
    if(!pf->HasTrk() && pf->Pt() <= ptMin) continue;

    if(pf->TrackerTrk() && p->TrackerTrk() &&
       pf->TrackerTrk() == p->TrackerTrk()) continue;

    if(pf->GsfTrk() && p->GsfTrk() &&
       pf->GsfTrk() == p->GsfTrk()) continue;

    // ignore the pf candidate if it is too far away in Z
    if(pf->BestTrk()) {
      Double_t deltaZ = TMath::Abs(pf->BestTrk()->DzCorrected(*vertex) - zLepton);
      if (deltaZ >= delta_z) continue;
    }
           
    Double_t dr = MathUtils::DeltaR(p->Mom(), pf->Mom());
    // add the pf pt if it is inside the extRadius and outside the intRadius
    if ( dr < extRadius && 
	 dr >= intRadius ) {

      //EtaStrip Veto for Gamma 
      if (pf->PFType() == PFCandidate::eGamma && fabs(p->Eta() - pf->Eta()) < 0.025) continue;

      //InnerCone (One Tower = dR < 0.07) Veto for non-gamma neutrals
      if (!pf->HasTrk() && pf->PFType() == PFCandidate::eNeutralHadron
          && MathUtils::DeltaR(p->Mom(), pf->Mom()) < 0.07 ) continue; 

      Bool_t isLepton = kFALSE;
      if(goodMuons){
        for (UInt_t nl=0; nl<goodMuons->GetEntries();nl++) {
	  const Muon *m = goodMuons->At(nl);
          if(pf->TrackerTrk() && m->TrackerTrk() &&
	     pf->TrackerTrk() == m->TrackerTrk()) {
	    isLepton = kTRUE;
	    break;
	  }
	}
      }
      if(goodElectrons && isLepton == kFALSE){
        for (UInt_t nl=0; nl<goodElectrons->GetEntries();nl++) {
	  const Electron *e = goodElectrons->At(nl);
          if(pf->TrackerTrk() && e->TrackerTrk() &&
	     pf->TrackerTrk() == e->TrackerTrk()) {
	    isLepton = kTRUE;
	    break;
	  }
          if(pf->GsfTrk() && e->GsfTrk() &&
	     pf->GsfTrk() == e->GsfTrk()) {
	    isLepton = kTRUE;
	    break;
	  }
	}
      }

      if (isLepton == kTRUE) continue;
      ptSum += pf->Pt();            

    }
  }
  return ptSum;
}
//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::BetaM(const TrackCol *tracks, const Muon *p, const Vertex *vertex, 
                               Double_t ptMin, Double_t  delta_z, Double_t extRadius,
			       Double_t intRadius){

  if(!tracks) return 1.0;
  if(tracks->GetEntries() <= 0) return 1.0;

  double Pt_jets_X = 0. ;
  double Pt_jets_Y = 0. ;
  double Pt_jets_X_tot = 0. ;
  double Pt_jets_Y_tot = 0. ;

  for(int i=0;i<int(tracks->GetEntries());i++){
    const Track *pTrack = tracks->At(i);

    if(pTrack && p->TrackerTrk() &&
       pTrack == p->TrackerTrk()) continue;

    if(pTrack->Pt() <= ptMin) continue;

    Double_t dr = MathUtils::DeltaR(pTrack->Mom(),p->Mom());
    if ( dr < extRadius && dr >= intRadius ) {
      Pt_jets_X_tot += pTrack->Px();
      Pt_jets_Y_tot += pTrack->Py();  
      double pDz = TMath::Abs(pTrack->DzCorrected(*vertex));
      if(pDz < delta_z){
        Pt_jets_X += pTrack->Px();
        Pt_jets_Y += pTrack->Py();
      }
    }
  }

  if(sqrt(Pt_jets_X_tot*Pt_jets_X_tot + Pt_jets_Y_tot*Pt_jets_Y_tot) > 0)
    return sqrt(Pt_jets_X*Pt_jets_X + Pt_jets_Y*Pt_jets_Y) / sqrt(Pt_jets_X_tot*Pt_jets_X_tot + Pt_jets_Y_tot*Pt_jets_Y_tot);

  return 1.0;
}

//--------------------------------------------------------------------------------------------------
Double_t IsolationTools::BetaE(const TrackCol *tracks, const Electron *p, const Vertex *vertex, 
                               Double_t ptMin, Double_t  delta_z, Double_t extRadius,
			       Double_t intRadius){

  if(!tracks) return 1.0;
  if(tracks->GetEntries() <= 0) return 1.0;

  double Pt_jets_X = 0. ;
  double Pt_jets_Y = 0. ;
  double Pt_jets_X_tot = 0. ;
  double Pt_jets_Y_tot = 0. ;

  for(int i=0;i<int(tracks->GetEntries());i++){
    const Track *pTrack = tracks->At(i);

    if(pTrack && p->TrackerTrk() &&
       pTrack == p->TrackerTrk()) continue;

    if(pTrack && p->GsfTrk() &&
       pTrack == p->GsfTrk()) continue;

    if(pTrack->Pt() <= ptMin) continue;

    Double_t dr = MathUtils::DeltaR(pTrack->Mom(),p->Mom());
    if ( dr < extRadius && dr >= intRadius ) {
      Pt_jets_X_tot += pTrack->Px();
      Pt_jets_Y_tot += pTrack->Py();  
      double pDz = TMath::Abs(pTrack->DzCorrected(*vertex));
      if(pDz < delta_z){
        Pt_jets_X += pTrack->Px();
        Pt_jets_Y += pTrack->Py();
      }
    }
  }

  if(sqrt(Pt_jets_X_tot*Pt_jets_X_tot + Pt_jets_Y_tot*Pt_jets_Y_tot) > 0)
    return sqrt(Pt_jets_X*Pt_jets_X + Pt_jets_Y*Pt_jets_Y) / sqrt(Pt_jets_X_tot*Pt_jets_X_tot + Pt_jets_Y_tot*Pt_jets_Y_tot);

  return 1.0;
}


// method added by F.Stoeckli: computes the track isolation with NO constrint on the OV-track compatibility
Double_t IsolationTools::TrackIsolationNoPV(const mithep::Particle* p, const BaseVertex* bsp, 
					    Double_t extRadius, 
					    Double_t intRadius, 
					    Double_t ptLow, 
					    Double_t etaStrip,
					    Double_t maxD0,
					    mithep::TrackQuality::EQuality quality,
					    const mithep::Collection<mithep::Track> *tracks,
                                            UInt_t maxNExpectedHitsInner,
                                            const mithep::DecayParticleCol *conversions) {
  
  // loop over all tracks
  Double_t tPt = 0.;
  //std::cout<<"      *** TrackIso:"<<std::endl;
  for(UInt_t i=0; i<tracks->GetEntries(); ++i) {
    const Track* t = tracks->At(i);
    if(t->Pt()>1. && false) {
      std::cout<<"     "<<i<<"    pt = "<<t->Pt()<<"  ("<<ptLow<<")"<<std::endl;
      std::cout<<"                d0 = "<<fabs(t->D0Corrected( *bsp) )<<"  ("<<maxD0<<")"<<std::endl;
      //std::cout<<"              conv ? "<<PhotonTools::MatchedConversion(t,conversions,bsp)<<std::endl;
      std::cout<<"                dR = "<<MathUtils::DeltaR(t->Mom(),p->Mom())<<"  ("<<extRadius<<","<<intRadius<<")"<<std::endl;
      std::cout<<"             dEta  = "<<fabs(t->Eta()-p->Eta())<<"   ("<<etaStrip<<")"<<std::endl;
    }
    if ( t->Pt() < ptLow ) continue;
    if ( ! t->Quality().Quality(quality) ) continue;
    // only check for beamspot if available, otherwise ignore cut
    if ( bsp && fabs(t->D0Corrected( *bsp) ) > maxD0) continue;
    if (t->NExpectedHitsInner()>maxNExpectedHitsInner) continue;
    if (conversions && PhotonTools::MatchedConversion(t,conversions,bsp)) continue;
    Double_t dR   = MathUtils::DeltaR(t->Mom(),p->Mom());
    Double_t dEta = fabs(t->Eta()-p->Eta());
    if(dR < extRadius && dR > intRadius && dEta > etaStrip) tPt += t->Pt();
  }
  return tPt;
}


Double_t IsolationTools::CiCTrackIsolation(const mithep::Photon* p, 
					   const BaseVertex* theVtx, 
					   Double_t extRadius, 
					   Double_t intRadius, 
					   Double_t ptLow, 
					   Double_t etaStrip,
					   Double_t maxD0,
					   Double_t maxDZ,
					   const mithep::Collection<mithep::Track> *tracks,
					   unsigned int* worstVtxIndex,
					   const mithep::Collection<mithep::Vertex> *vtxs,
					   const mithep::Collection<mithep::Electron> *eles,
					   bool print,
					   double* ptmax, double* dRmax) {
  
  UInt_t numVtx = 1;
  const BaseVertex* iVtx = theVtx;
  if( vtxs ) { 
    numVtx = vtxs->GetEntries();
    if (numVtx > 0)
      iVtx = vtxs->At(0);
    else
      return 0.;
  }
  
  // NEW for Electron T&P: need to remove the electron Gsf Track (applied if eles != NULL)
  const Track* theGsfTrack = NULL;
  if ( eles ) {
    // find the electron that matches the Photon SC
    for(UInt_t j=0; j<eles->GetEntries(); ++j) {
      if ( eles->At(j)->SCluster() == p->SCluster() ) {
	if( eles->At(j)->HasTrackerTrk() )
	  theGsfTrack = eles->At(j)->TrackerTrk();
	break;
      }
    }
  }

  if(print) {
    std::cout<<" Testing photon with"<<std::endl;
    std::cout<<"             Et  = "<<p->Et()<<std::endl;
    std::cout<<"             Eta = "<<p->Eta()<<std::endl;
    std::cout<<"             Phi = "<<p->Phi()<<std::endl;
  }

  Double_t iIso = 0.;
  Double_t maxIso = 0.;
  
  if(worstVtxIndex)
    *worstVtxIndex=0;

  double t_ptmax = 0.;
  double t_dRmax = 0.;

  for(UInt_t i=0; i<numVtx; ++i) {

    if(i>0) iVtx = vtxs->At(i);


    if(print) {
      std::cout<<"   Vertex #"<<i<<std::endl;
      std::cout<<"       with X = "<<iVtx->X()<<std::endl;
      std::cout<<"       with Y = "<<iVtx->Y()<<std::endl;
      std::cout<<"       with Z = "<<iVtx->Z()<<std::endl;
    }

    Photon* phTemp = new Photon(*p);

    // RESET CALO_POS!
    phTemp->SetCaloPosXYZ(p->SCluster()->Point().X(),p->SCluster()->Point().Y(),p->SCluster()->Point().Z());

    // compute the ph momentum with respect to this Vtx
    //FourVectorM phMom = p->MomVtx(iVtx->Position());
    FourVectorM phMom = phTemp->MomVtx(iVtx->Position());

    delete phTemp;

    if(print) {
      std::cout<<"         photon has changed to:"<<std::endl;
      std::cout<<"             Et  = "<<phMom.Et()<<std::endl;
      std::cout<<"             eta = "<<phMom.Eta()<<std::endl;
      std::cout<<"             Phi = "<<phMom.Phi()<<std::endl;
    }
    
    iIso = 0.;
    for(UInt_t j=0; j<tracks->GetEntries(); ++j) {
      const Track* t = tracks->At(j);
      if( theGsfTrack && t == theGsfTrack ) continue;

      //Double_t dR   = MathUtils::DeltaR(t->Mom(),p->Mom());
      //Double_t dEta = TMath::Abs(t->Eta()-p->Eta());

      Double_t dR   = MathUtils::DeltaR(t->Mom(),phMom);
      Double_t dEta = TMath::Abs(t->Eta()-phMom.Eta());

      if(print && t->Pt()>1. && false) {
	  std::cout<<"              passing track #"<<j<<std::endl;
	  std::cout<<"                          pt = "<<t->Pt()<<std::endl;
	  std::cout<<"                         eta = "<<t->Eta()<<std::endl;
	  std::cout<<"                         phi = "<<t->Phi()<<std::endl;
	  std::cout<<"                          d0 = "<<fabs(t->D0Corrected( *iVtx ))<<std::endl;
	  std::cout<<"                          dZ = "<<fabs(t->DzCorrected( *iVtx ))<<std::endl;
	  std::cout<<"                          dR = "<<dR<<std::endl;
	  std::cout<<"                        dEta = "<<dEta<<std::endl;
	  std::cout<<"                          vx = "<<t->X0()<<std::endl;
	  std::cout<<"                          vy = "<<t->Y0()<<std::endl;
	  std::cout<<"                          vz = "<<t->Z0()<<std::endl;
      }

      if ( t->Pt() < ptLow ) continue;
      // only check for beamspot if available, otherwise ignore cut
      if ( fabs(t->D0Corrected( *iVtx )) > maxD0) continue;
      if ( fabs(t->DzCorrected( *iVtx )) > maxDZ) continue;


      if(dR < extRadius && dR > intRadius && dEta > etaStrip) {
	iIso += t->Pt();

	if(t->Pt() > t_ptmax) {
	  t_ptmax=t->Pt();
	  t_dRmax=dR;
	}

	if(print && t->Pt()>1.) {
	  std::cout<<"              passing track #"<<j<<std::endl;
	  std::cout<<"                          pt = "<<t->Pt()<<std::endl;
	  std::cout<<"                         eta = "<<t->Eta()<<std::endl;
	  std::cout<<"                         phi = "<<t->Phi()<<std::endl;
	  std::cout<<"                          d0 = "<<fabs(t->D0Corrected( *iVtx ))<<std::endl;
	  std::cout<<"                          dZ = "<<fabs(t->DzCorrected( *iVtx ))<<std::endl;
	  std::cout<<"                          dR = "<<dR<<std::endl;
	  std::cout<<"                        dEta = "<<dEta<<std::endl;
	  std::cout<<"                          vx = "<<t->X0()<<std::endl;
	  std::cout<<"                          vy = "<<t->Y0()<<std::endl;
	  std::cout<<"                          vz = "<<t->Z0()<<std::endl;
	  std::cout<<"                             new  tIso = "<<iIso<<std::endl;
	}
      }
    }
    if ( iIso > maxIso ) {
      maxIso = iIso;
      if(worstVtxIndex)
	*worstVtxIndex=i;
    }
  }

  if(ptmax) (*ptmax)=t_ptmax;
  if(dRmax) (*dRmax)=t_dRmax;

  if(print) {
    if(worstVtxIndex)
      std::cout<<"   max TrkIso is given by Vtx #"<<*worstVtxIndex<<" with an amount of tIso = "<<maxIso<<std::endl;
    else
      std::cout<<"   max TrkIso is given by Vtx #0 with an amount of tIso = "<<maxIso<<std::endl;
  } 
  return maxIso;
}
