//--------------------------------------------------------------------------------------------------
// $Id $
//
// ElectronTools
//
// Helper Class for electron Identification decisions.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_ELECTRONTOOLS_H
#define MITPHYSICS_UTILS_ELECTRONTOOLS_H

#include "MitAna/DataTree/interface/Electron.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"

namespace mithep {
  class ElectronTools {
    public:
      ElectronTools();
  
     enum EElIdType {
        kIdUndef = 0,       //not defined
        kTight,             //"Tight"
        kLoose,             //"Loose"
        kLikelihood,        //"Likelihood"
        kNoId,              //"NoId"
        kZeeId,             //"ZeeId"
        kCustomIdLoose,     //"CustomLoose"
        kCustomIdTight,     //"CustomTight"
        kVBTFWorkingPointFakeableId,
        kVBTFWorkingPoint95Id,
        kVBTFWorkingPoint90Id,
        kVBTFWorkingPoint85Id,
        kVBTFWorkingPoint80Id,
        kVBTFWorkingPointLowPtId,
        kVBTFWorkingPoint70Id,
        kMVA_BDTG_V1,
        kMVA_BDTG_V2
      };

      enum EElIsoType {
        kIsoUndef = 0,      	       //"not defined"
        kTrackCalo,         	       //"TrackCalo"
        kTrackJura,         	       //"TrackJura"
        kTrackJuraCombined, 	       //"TrackJuraCombined"
        kTrackJuraSliding,  	       //"TrackJuraSliding"
        kTrackJuraSlidingNoCorrection, //"TrackJuraSlidingNoCorrection"
        kCombinedRelativeConeAreaCorrected, //"CombinedRelativeConeAreaCorrected"
        kNoIso,             	       //"NoIso"
        kPFIso,             	       //"PFIso"
        kPFIsoNoL,          	       //"PFIsoNoL"
        kZeeIso,            	       //"ZeeIso"
        kCustomIso,         	       //"Custom"
        kVBTFWorkingPoint95Iso,
        kVBTFWorkingPoint90Iso,
        kVBTFWorkingPoint85Iso,
        kVBTFWorkingPoint80Iso,
        kVBTFWorkingPoint70Iso
      };

      static Bool_t       PassChargeFilter(const Electron *el);
      static Bool_t       PassConversionFilter(const Electron *el, const DecayParticleCol *conversions, 
                                               const BaseVertex *vtx, UInt_t nWrongHitsMax=0, Double_t probMin=1e-6,
                                               Double_t lxyMin = 2.0, Bool_t matchCkf = kTRUE, Bool_t requireArbitratedMerged = kFALSE, Double_t trkptMin = -99.);
      static Bool_t       PassCustomID(const Electron *el, EElIdType idType);
      static Bool_t       PassCustomIso(const Electron *el, EElIsoType isoType,
                                        Bool_t useCombineIso = kTRUE);
      static Bool_t       PassD0Cut(const Electron *el, const VertexCol *vertices, Double_t fD0Cut, Int_t nVertex = 0);
      static Bool_t       PassD0Cut(const Electron *el, const BeamSpotCol *beamspots, Double_t fD0Cut);
      static Bool_t       PassDZCut(const Electron *el, const VertexCol *vertices, Double_t fDZCut, Int_t nVertex = 0);
      static Bool_t       PassSpikeRemovalFilter(const Electron *ele);
      static Bool_t       PassTriggerMatching(const Electron *ele, const TriggerObjectCol *trigobjs);
      static Int_t        Classify(const Electron *ele);
      static Int_t        PassTightId(const Electron *ele, const VertexCol *vertices, 
                                      const DecayParticleCol *conversions, const Int_t typeCuts,
                                      Double_t beta = 1.0);
      static bool         compute_cut(double x, double et, double cut_min, double cut_max, bool gtn=false);
      static Double_t     Likelihood(ElectronLikelihood *LH, const Electron *ele);

    ClassDef(ElectronTools, 0) // Muon tools
  };
}

#endif
