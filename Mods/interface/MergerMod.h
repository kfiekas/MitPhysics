//--------------------------------------------------------------------------------------------------
// $Id: MergeParticlesMod.h,v 1.2 2008/12/10 17:28:22 loizides Exp $
//
// MergerMod
//
// This templated module merges arbitrarily many collections of objects. It puts them into an 
// object array of the templated type. 
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_MERGERMOD_H
#define MITPHYSICS_MODS_MERGERMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"

namespace mithep 
{
  template <class T>
  class MergerMod : public BaseMod 
  {
    public:
      MergerMod(const char *name="MergerMod", 
                const char *title="Merger module");
      ~MergerMod() {}

      void                     AddInputName(const char *n);
      const char              *GetMergedName()              const { return fMergedName;     }
      const char              *GetOutputName()              const { return GetMergedName(); }
      void                     SetMergedName(const char *n)       { fMergedName=n;          }
      void                     SetOutputName(const char *n)       { SetMergedName(n);       }

    protected:
      void                     SlaveBegin();
      void                     Process();

      StringVec                fInList;        //list of to be merged objects (input)
      TString                  fMergedName;    //name of merged collection (output)

    ClassDefT(MergerMod, 1) // Merger module
  };
}
#endif

//--------------------------------------------------------------------------------------------------
template <class T>
mithep::MergerMod<T>::MergerMod<T>(const char *name, const char *title) : 
  BaseMod(name,title),
  fInList(0)
{
  // Constructor.
}

template <class T>
void mithep::MergerMod<T>::AddInputName(const char *n)
{
  // Add name to the list of input collections.

  if (!n)
    return;
  fInList.AddCopy(std::string(n));
}

//--------------------------------------------------------------------------------------------------
template <class T>
void mithep::MergerMod<T>::Process()
{
  // Merge the input objects and publish merged collection. 

  ObjArray<T> *colout = new mithep::ObjArray<T>(0);
  colout->SetName(GetMergedName());

  UInt_t n = fInList.GetEntries();
  for (UInt_t i=0; i<n; ++i) {
    const TObject *in = GetObjThisEvt<TObject>(fInList.At(i)->c_str());
    if (!in)
      continue;

    const BaseCollection *bcol = dynamic_cast<const BaseCollection*>(in);
    if (bcol)
      colout->Add(bcol);
    else {
      const T *optr = dynamic_cast<const T*>(in);
      if (optr)
        colout->Add(optr);
    }
  }

  // attempt sorting (if implemented by underlying objects)
  colout->Sort();

  // add to event for other modules to use
  AddObjThisEvt(colout);
}

//--------------------------------------------------------------------------------------------------
template <class T>
void mithep::MergerMod<T>::SlaveBegin()
{
  // Check if output name was set.

  if (!fMergedName.IsNull())
    return;

  SendError(kAbortModule, "SlaveBegin", "No output name given.");
}
