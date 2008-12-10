//--------------------------------------------------------------------------------------------------
// $Id: GenericSelMod.h,v 1.1 2008/12/09 10:18:19 loizides Exp $
//
// GenericSelMod
// 
// This module allows trivial event selection based on counting the number of
// particles found in the given kinematic range.
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_SELMODS_GENERICSELMOD_H
#define MITPHYSICS_SELMODS_GENERICSELMOD_H

#include "MitAna/TreeMod/interface/BaseSelMod.h" 

namespace mithep 
{
  template<class T>
  class GenericSelMod : public BaseSelMod
  {
    public:
      GenericSelMod(const char *name="GenericSelMod", 
                    const char *title="Generic selection module");
      ~GenericSelMod() {}

      const char              *GetColName()              const { return fColName;     }
      Double_t                 GetEtaMin()               const { return fEtaMin;      }
      Double_t                 GetEtaMax()               const { return fEtaMax;      }
      const char              *GetInputName()            const { return GetColName(); }
      UInt_t                   GetMinCounts()            const { return fMinCounts;   }
      Double_t                 GetPtMin()                const { return fPtMin;       }
      Double_t                 GetPtMax()                const { return fPtMax;       }
      void                     SetColName(const char *n)       { fColName=n;          }
      void                     SetEtaMin(Double_t e)           { fEtaMin = e;         }
      void                     SetEtaMax(Double_t e)           { fEtaMax = e;         }
      void                     SetInputName(const char *n)     { SetColName(n);       }
      void                     SetPtMin(Double_t pt)           { fPtMin = pt;         }
      void                     SetPtMax(Double_t pt)           { fPtMax = pt;         }
      void                     SetMinCounts(UInt_t c)          { fMinCounts = c;      }

    protected:
      void                     Process();

      TString                  fColName;    //name of input collection
      Double_t                 fPtMin;      //minimum pt required  (def=0 GeV)
      Double_t                 fPtMax;      //maximum pt required  (def=5000 GeV)
      Double_t                 fEtaMin;     //minimum eta required (def=-10)
      Double_t                 fEtaMax;     //maximum eta required (def=+10) 
      UInt_t                   fMinCounts;  //minimum number of particles required to accept event
      const Collection<T>     *fCol;        //!pointer to collection 

      ClassDefT(GenericSelMod,1) // Generic selection module
  };
}

//--------------------------------------------------------------------------------------------------
template<class T>
mithep::GenericSelMod<T>::GenericSelMod(const char *name, const char *title) : 
  BaseSelMod(name,title),
  fColName("SetMe"),
  fPtMin(0),
  fPtMax(5000),
  fEtaMin(-10),
  fEtaMax(10),
  fMinCounts(1),
  fCol(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
template<class T>
void mithep::GenericSelMod<T>::Process()
{
  // Process entries of the tree: Just load the branch and fill the histograms.

  fCol = GetObjThisEvt<Collection<T> >(GetColName());
  if (!fCol) {
    this->SendError(kAbortModule, "Process", 
                    "Could not obtain collection with name %s!", GetColName());
    return;
  }

  UInt_t counter = 0;
  UInt_t ents=this->fCol->GetEntries();
  for(UInt_t i=0;i<ents;++i) {
     const T *p = this->fCol->At(i);
     Double_t pt = p->Pt();
     if (pt<this->fPtMin) 
       continue;
     if (pt>this->fPtMax)
       continue;
     Double_t eta = p->Eta();
     if (eta<this->fEtaMin)
       continue;
     if (eta>this->fEtaMax)
       continue;

     ++counter;
  }

  // skip event if not enough particles are found in kinematic region
  if (counter<GetMinCounts())
    this->SkipEvent();
}
#endif
