//--------------------------------------------------------------------------------------------------
// $Id: IsolationTools.h,v 1.6 2011/02/05 05:48:13 ceballos Exp $
//
// IsolationTools
//
// Isolation functions to compute various kinds of isolation.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_ISOLATIONTOOLS_H
#define MITPHYSICS_UTILS_ISOLATIONTOOLS_H

#include <TMath.h>
#include "MitAna/DataTree/interface/Track.h"
#include "MitAna/DataTree/interface/BasicCluster.h"
#include "MitAna/DataTree/interface/SuperCluster.h"
#include "MitAna/DataTree/interface/CaloTower.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"

namespace mithep
{
  class IsolationTools {
    public:
      static Double_t TrackIsolation(const mithep::Track *p, Double_t extRadius, 
                                     Double_t intRadius, Double_t ptLow, Double_t maxVtxZDist, 
                                     const mithep::Collection<mithep::Track> *tracks); 
      static Double_t EcalIsolation(const SuperCluster *sc, Double_t coneSize, Double_t etLow, 
                                    const mithep::Collection<mithep::BasicCluster> *basicClusters);
      static Double_t CaloTowerHadIsolation(const ThreeVector *p,  Double_t extRadius, 
                                            Double_t intRadius, Double_t etLow, 
                                            const mithep::Collection<mithep::CaloTower> 
                                            *caloTowers);
      static Double_t CaloTowerEmIsolation(const ThreeVector *p, Double_t extRadius, 
                                           Double_t intRadius, Double_t etLow, 
                                           const mithep::Collection<mithep::CaloTower> *caloTowers);
      static Double_t PFMuonIsolation(const Muon *p, const PFCandidateCol *PFCands, 
                      	   	      const VertexCol *vertices, Double_t  delta_z, double ptMin,
			   	      Double_t extRadius, Double_t intRadius, int isoType);
      static Double_t PFElectronIsolation(const Electron *p, const PFCandidateCol *PFCands, 
                      	     		  const VertexCol *vertices, Double_t  delta_z, double ptMin,
			     		  Double_t extRadius, Double_t intRadius, int isoType);
    ClassDef(IsolationTools, 0) // Isolation tools
  };
}
#endif
