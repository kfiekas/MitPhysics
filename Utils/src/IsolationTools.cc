// $Id: IsolationTools.cc,v 1.3 2009/03/03 18:11:19 bendavid Exp $

#include "MitPhysics/Utils/interface/IsolationTools.h"
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
