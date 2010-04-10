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
#include "MitCommon/MathTools/interface/MathUtils.h"

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
        kVBTFWorkingPoint95Id,
        kVBTFWorkingPoint90Id,
        kVBTFWorkingPoint80Id,
        kVBTFWorkingPoint70Id
      };

      enum EElIsoType {
        kIsoUndef = 0,      //not defined
        kTrackCalo,         //"TrackCalo"
        kTrackJura,         //"TrackJura"
        kTrackJuraCombined, //"TrackJuraCombined"
        kTrackJuraSliding,  //"TrackJuraSliding"
        kNoIso,             //"NoIso"
        kZeeIso,            //"ZeeIso"
        kCustomIso,          //"Custom"
        kVBTFWorkingPoint95Iso,
        kVBTFWorkingPoint90Iso,
        kVBTFWorkingPoint80Iso,
        kVBTFWorkingPoint70Iso
      };

      static Bool_t       PassChargeFilter(const Electron *el);
      static Bool_t       PassConversionFilter(const Electron *el, const DecayParticleCol *conversions, 
                                               Bool_t WrongHitsRequirement );
      static Bool_t       PassCustomID(const Electron *el, EElIdType idType);
      static Bool_t       PassCustomIso(const Electron *el, EElIsoType isoType);
      static Bool_t       PassD0Cut(const Electron *el, const VertexCol *vertices, 
                                    Double_t fD0Cut, Bool_t fReverseD0Cut);
      static Bool_t       PassD0Cut(const Electron *el, const BeamSpotCol *beamspots, 
                                    Double_t fD0Cut, Bool_t fReverseD0Cut);
      static Bool_t       PassSpikeRemovalFilter(const Electron *ele);
   

    ClassDef(ElectronTools, 0) // Muon tools
  };
}

#endif
