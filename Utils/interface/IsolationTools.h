//--------------------------------------------------------------------------------------------------
// $Id: IsolationTools.h,v 1.17 2011/06/01 18:11:52 fabstoec Exp $
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
#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/BasicCluster.h"
#include "MitAna/DataTree/interface/SuperCluster.h"
#include "MitAna/DataTree/interface/CaloTower.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"

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
      static Double_t PFMuonIsolation(const Muon *p, const PFCandidateCol *PFCands, const Vertex *vertex, 
                      	   	      Double_t  delta_z = 0.1, Double_t ptMin = 1.0,
			   	      Double_t extRadius = 0.4, Double_t intRadius = 0.07);
      static Double_t PFMuonIsolation(const Muon *p, const PFCandidateCol *PFCands, const Vertex *vertex, 
				      const MuonCol *goodMuons, const ElectronCol *goodElectrons,
				      Double_t  delta_z = 0.1, Double_t ptMin = 1.0, Double_t extRadius = 0.4,
				      Double_t intRadius = 0.07, int isoType = 0, Double_t beta = 1.0);
      static Double_t PFElectronIsolation(const Electron *p, const PFCandidateCol *PFCands, 
                      	     		  const Vertex *vertex, Double_t  delta_z, Double_t ptMin,
			     		  Double_t extRadius, Double_t intRadius);
      static Double_t PFElectronIsolation(const Electron *p, const PFCandidateCol *PFCands, const Vertex *vertex, 
                                          const MuonCol *goodMuons, const ElectronCol *goodElectrons, 
                                          Double_t  delta_z = 0.1 , Double_t ptMin = 1.0,Double_t extRadius = 0.4,			     		  
                                          Double_t intRadius = 0.0, int isoType = 0, Double_t beta = 1.0);
      static Double_t BetaM(const TrackCol *tracks, const Muon *p, const Vertex *vertex, 
                            Double_t ptMin, Double_t  delta_z, Double_t extRadius,
			    Double_t intRadius);
      static Double_t BetaE(const TrackCol *tracks, const Electron *p, const Vertex *vertex, 
                            Double_t ptMin, Double_t  delta_z, Double_t extRadius,
			    Double_t intRadius);

      // method added by F.Stoeckli: computes the track isolation with NO constrint on the OV-track compatibility
      static Double_t TrackIsolationNoPV(const mithep::Particle*, const BaseVertex*, 
					 Double_t extRadius, 
					 Double_t intRadius, 
					 Double_t ptLow, 
					 Double_t etaStrip,
					 Double_t maxD0,
					 mithep::TrackQuality::EQuality,
					 const mithep::Collection<mithep::Track> *tracks,
                                         UInt_t maxNExpectedHitsInner = 999,
                                         const mithep::DecayParticleCol *conversions = 0);

      // methods for Hgg BaseLien Selection. These isoaltion are stupid, but what can we do.... ;(
      static Double_t CiCTrackIsolation(const mithep::Photon*, 
					const BaseVertex*, 
					Double_t extRadius, 
					Double_t intRadius, 
					Double_t ptLow, 
					Double_t etaStrip,
					Double_t maxD0,
					Double_t maxDZ,
					const mithep::Collection<mithep::Track> *tracks,
					unsigned int* worstVtxIdx = NULL,
					const mithep::Collection<mithep::Vertex> *vtxs = NULL,
					bool print=false);


    ClassDef(IsolationTools, 0) // Isolation tools
  };
}
#endif
