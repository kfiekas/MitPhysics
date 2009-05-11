//--------------------------------------------------------------------------------------------------
// $Id: MatchingTools.h,v 1.3 2009/02/17 06:49:01 phedex Exp $
//
// MatchingTools
//
//
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_MATCHINGTOOLS_H
#define MITPHYSICS_UTILS_MATCHINGTOOLS_H

#include <TMath.h>
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataCont/interface/Collection.h"

namespace mithep
{
  class MatchingTools {
    public:
      template<class V1, class V2> 
      static const V2 *Closest(const V1 *v1, const Collection<V2> &col, Double_t maxR);
      template<class V1, class V2> 
      static const V2 *ClosestEtaPhiPt(const V1 *v1, const Collection<V2> &col, Double_t max); 
  };

  //------------------------------------------------------------------------------------------------
  template<class V1, class V2>
  const V2 *MatchingTools::Closest(const V1 *v1, const Collection<V2> &col, Double_t maxR)
  {
    // Return closest geometrical neighbor in eta-phi within maximum delta R of maxR.
    
    const V2 *res = 0;

    const UInt_t ents = col.GetEntries();
    if (ents>0) {
      Double_t dR = 1e30;
      for (UInt_t i = 0; i<ents; ++i) {
        const V2 *v2 = col.At(i);
        Double_t diff = MathUtils::DeltaR(*v1,*v2);
        if ((diff<maxR) && (diff<dR)) {
          res = v2;
          dR  = diff;
        }
      }
    }
    return res;
  }

  //------------------------------------------------------------------------------------------------
  template<class V1, class V2>
  const V2 *MatchingTools::ClosestEtaPhiPt(const V1 *v1, const Collection<V2> &col, Double_t max)
  {
    // Return closest geometrical neighbor in eta-phi-pt within maximum delta of max.
    
    const V2 *res = 0;

    const UInt_t ents = col.GetEntries();
    if (ents>0) {
      Double_t d = 1e30;
      for (UInt_t i = 0; i<ents; ++i) {
        const V2 *v2 = col.At(i);
        Double_t diff = MathUtils::DeltaR2(*v1,*v2);
        Double_t ptd = 1-v1->Pt()/v2->Pt();
        diff += ptd*ptd;
        diff = TMath::Sqrt(diff);
        if ((diff<max) && (diff<d)) {
          res = v2;
          d   = diff;
        }
      }
    }
    return res;
  }
}
#endif
