//--------------------------------------------------------------------------------------------------
// $Id: FakeObject.h,v 1.2 2009/07/20 19:05:04 loizides Exp $
//
// FakeObject
//
// Class for particles which have been faked from a fakeable object
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_FAKEMODS_FAKEOBJECT_H
#define MITPHYSICS_FAKEMODS_FAKEOBJECT_H
 
#include "MitAna/DataTree/interface/Particle.h"

namespace mithep 
{
  class FakeObject : public Particle
  {
    public:
      FakeObject() {}
      ~FakeObject() {}

      const Particle *FakeParticle()           const { return fParticle;         }
      Bool_t          FakeTag()                const { return fFakeTag;          }
      Bool_t          MCTag()                  const { return fMCTag;            }
      EObjType        ObjType()                const { return fFakeType;         }

      void            SetFakeType(EObjType type)     { fFakeType = type;         }
      void            SetFakeTag(Bool_t tag)         { fFakeTag  = tag;          }
      void            SetMCTag(Bool_t tag)           { fMCTag    = tag;          }
      void            SetParticle(const Particle *p) { fParticle = p; ClearMom(); ClearCharge(); }

    protected:
      Double_t        GetCharge()              const;
      Double_t        GetMass()                const { return fParticle->Mass(); }
      void            GetMom()                 const;

      Bool_t          fFakeTag;               //whether the fake object passed lepton ID criteria
      Bool_t          fMCTag;                 //whether the fake object was a fake in Monte Carlo
      EObjType        fFakeType;              //the type of object it faked
      const Particle *fParticle;              //pointer to the original particle
 
    ClassDef(FakeObject, 1) // Fake object class
  };
}

//--------------------------------------------------------------------------------------------------
inline Double_t mithep::FakeObject::GetCharge() const
{
  // Get charge of the fake object. Charge is calculated from the original particle.
  
  if (fParticle) 
    return fParticle->Charge();
  else 
    return 0;
}

//--------------------------------------------------------------------------------------------------
inline void mithep::FakeObject::GetMom() const
{
  // Get momentum of the electron. We use an explicitly stored three vector, with the pdg mass,
  // since the momentum vector may be computed non-trivially in cmssw
  if (fParticle) {
    fCachedMom.SetCoordinates(fParticle->Pt(),fParticle->Eta(),
                              fParticle->Phi(),fParticle->Mass());
  } else {
    fCachedMom.SetCoordinates(0,0,0,0);
  }
}

#endif
