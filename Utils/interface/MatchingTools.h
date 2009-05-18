//--------------------------------------------------------------------------------------------------
// $Id: MatchingTools.h,v 1.1 2009/05/11 08:01:51 loizides Exp $
//
// MatchingTools
//
// This class implements a couple of simple geometrical matching functions.
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_MATCHINGTOOLS_H
#define MITPHYSICS_UTILS_MATCHINGTOOLS_H

#include <TObjArray.h>
#include <TMath.h>
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataCont/interface/Collection.h"

namespace mithep
{
  class MatchingTools {
    public:
      template<class V1, class V2> 
      static const V2 *Closest(const V1 *v1, const Collection<V2> &col, Double_t maxR, 
                               Bool_t self=kTRUE);
      template<class V1, class V2> 
      static const V2 *ClosestRPt(const V1 *v1, const Collection<V2> &col, Double_t max,
                                  Bool_t self=kTRUE);
      template<class V1, class V2> 
      static TObjArray *Closests(const V1 *v1, const Collection<V2> &col, Double_t maxR,
                                 Bool_t self=kTRUE);
  };

  //------------------------------------------------------------------------------------------------
  template<class V1, class V2>
  const V2 *MatchingTools::Closest(const V1 *v1, const Collection<V2> &col, Double_t maxR, 
                                   Bool_t self)
  {
    // Return closest geometrical neighbor in eta-phi within maximum delta R of maxR.
    // If self is kFALSE make sure that identical pointers are excluded.
    
    if (!v1)
      return 0;

    const V2 *res = 0;
    const UInt_t ents = col.GetEntries();
    if (ents>0) {
      Double_t dR = 1e30;
      for (UInt_t i = 0; i<ents; ++i) {
        const V2 *v2 = col.At(i);
        if (!v2) 
          continue;
        if (!self && (void*)v1==(void*)v2) 
          continue;
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
  const V2 *MatchingTools::ClosestRPt(const V1 *v1, const Collection<V2> &col, Double_t max,
                                      Bool_t self)
  {
    // Return closest geometrical neighbor in eta-phi-pt within maximum delta of max.
    // If self is kFALSE make sure that identical pointers are excluded.
    
    if (!v1)
      return 0;

    const V2 *res = 0;
    const UInt_t ents = col.GetEntries();
    if (ents>0) {
      Double_t d = 1e30;
      for (UInt_t i = 0; i<ents; ++i) {
        const V2 *v2 = col.At(i);
        if (!v2) 
          continue;
        if (!self && (void*)v1==(void*)v2) 
          continue;
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

  //------------------------------------------------------------------------------------------------
  template<class V1, class V2>
  TObjArray *MatchingTools::Closests(const V1 *v1, const Collection<V2> &col, Double_t maxR,
                                     Bool_t self)
  {
    // Return closest geometrical neighbors in eta-phi within maximum delta R of maxR.
    // If self is kFALSE make sure that identical pointers are excluded.

    if (!v1)
      return 0;

    TObjArray *res = new TObjArray;
    res->SetOwner(kFALSE);
    const UInt_t ents = col.GetEntries();
    for (UInt_t i = 0; i<ents; ++i) {
      const V2 *v2 = col.At(i);
      if (!v2)
        continue;
      if (!self && (void*)v1==(void*)v2) 
        continue;
      Double_t diff = MathUtils::DeltaR(*v1,*v2);
      if (diff<maxR)
        res->Add(const_cast<V2*>(v2));
    }
    return res;
  }
}
#endif
