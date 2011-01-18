//-------------------
//
// Jet Tools
//
// S Markson
//
//-------------------

#ifndef MITPHYSICS_UTILS_JETTOOLS_H
#define MITPHYSICS_UTILS_JETTOOLS_H

#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/CaloTowerCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include <TVector3.h>
#include <TLorentzVector.h>

namespace mithep {
  class JetTools {
  public:
    JetTools();
    virtual ~JetTools();
    
    static Double_t NJettiness(const ParticleOArr *particles, const JetOArr *jets, double Q = 1, double Y = 0.0);
    static Double_t NJettiness(const PFCandidateOArr *particles, const JetOArr *jets, double Q = 1, double Y = 0.0);
    static Double_t NJettiness(const TrackOArr *tracks, const JetOArr *jets, double Q = 1, double Y = 0.0);
    static Double_t NJettiness(const JetOArr *jetsS, const JetOArr *jets, double Q = 1, double Y = 0.0);
    static Double_t NJettiness(const CaloTowerOArr *calos, const JetOArr *jets, double Q = 1, double Y = 0.0);
    
    static Double_t M_r(const ParticleOArr *particles);
    static Double_t Beta_r(const ParticleOArr *particles);
    static Double_t M_r_t(const ParticleOArr *particles, const Met *met);
    static Double_t Razor(const ParticleOArr *particles, const Met *met);
    static Double_t CosineOmega(const Particle *particles0, const Particle *particles1);
    static Double_t MtHiggs(const ParticleOArr *leptons, const Met *met, double metFraction[2], int nsel);
    static void     Alpha(Double_t AlphaVar[2], const TrackCol *tracks, Jet *jet, const VertexCol *vertices, Double_t  delta_z, Double_t delta_cone);
    
    ClassDef(JetTools, 0)
  };

}

#endif
