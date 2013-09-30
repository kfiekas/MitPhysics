#ifndef MITPHYSICS_UTILS_QGTAGGER_H
#define MITPHYSICS_UTILS_QGTAGGER_H

#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitPhysics/Utils/interface/QGLikelihoodCalculator.h"

#include <map>



namespace mithep {
  class QGTagger {
  public:
    QGTagger(bool useCHS);
    virtual ~QGTagger();
    
    // setters
    void SetRhoIso(Float_t rhoIso) {variables["rhoIso"]=rhoIso;}

    // compute the inputs for gluon-quark discrimination
    void CalculateVariables(const PFJet *jet, const VertexCol *vertices);
    
    // getters
    Float_t QGValue(); 

  private:
    QGLikelihoodCalculator *qgLikelihood;
    std::map<TString, Float_t> variables;


    ClassDef(QGTagger, 0)
  };

}

#endif
