#ifndef MITPHYSICS_UTILS_JETTOOLS_H
#define MITPHYSICS_UTILS_JETTOOLS_H

#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
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
    static Double_t Beta(const TrackCol *tracks, Jet *jet, const Vertex *vertex, Double_t  delta_z, Double_t delta_cone);
    static Double_t Beta(const PFJet *jet, const Vertex *vertex, Double_t  delta_z);
    static Bool_t   PassBetaVertexAssociationCut(const PFJet *jet, const Vertex *referenceVertex, const VertexCol *vertices, Double_t delta_z);
    static Double_t Beta2(const PFJet *jet, const Vertex *vertex, Double_t  delta_z);
    static Bool_t   PassBeta2VertexAssociationCut(const PFJet *jet, const Vertex *referenceVertex, const VertexCol *vertices, Double_t delta_z);
    static Int_t    MaxBetaVertexIndex(const PFJet *jet, const VertexCol *vertices, Double_t  delta_z);
    static Int_t    MaxBeta2VertexIndex(const PFJet *jet, const VertexCol *vertices, Double_t  delta_z);    
    static Int_t    JetToPVAssociation(const PFJet *jet, const VertexCol *vertices, Double_t  delta_z);    

    static Double_t           impactParameter(const PFJet *iJet,const Vertex *iVertex,bool iDZ=false);
    static const PFCandidate* leadCand       (const PFJet *iJet,int iPFType,bool i2nd=false);
    static Double_t           dRMean         (const PFJet *iJet,int iPFType);
    static Bool_t             passPFLooseId  (const PFJet *iJet);
    ClassDef(JetTools, 0)
  };

}

#endif
