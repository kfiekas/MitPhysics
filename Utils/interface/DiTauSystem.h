//--------------------------------------------------------------------------------------------------
// $Id: DiTauSystem.h,v 1.4 2008/12/05 09:35:01 ceballos Exp $
//
// DiTauSystem
//
// Class to calculate the mass to the di-tau system. It is assumed that the tau is boosted 
// and that the neutrinos have the same flight direction as the tau. 
// CMS note 2006/082
//
// Authors: G.Ceballos 
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_DITAUSYSTEM_H
#define MITPHYSICS_UTILS_DITAUSYSTEM_H

#include <Rtypes.h>

namespace mithep
{
  class ChargedParticle;
  class Met;

  class DiTauSystem {
    public:
      DiTauSystem(ChargedParticle *t1, ChargedParticle *t2, const Met *met);
      ~DiTauSystem() {}

      Double_t         RecoMass()        const { return fRecoMass;}
      Double_t         TransverseMass()  const { return fMT;      }
      Double_t         TransverseEll()   const { return fETll;    }
      Double_t         TransverseEnn()   const { return fETnn;    }
      Double_t         VisMass()         const { return fVisMass; }
      Double_t         XTau1()           const { return fXTau[0]; }
      Double_t         XTau2()           const { return fXTau[1]; }
  
    private:
      void             Init();
  
      ChargedParticle *fT1;       //first tau
      ChargedParticle *fT2;       //second tau
      const Met       *fMet;      //missing et
      Double_t         fXTau[2];  //visible fraction of the tau momenta
      Double_t         fRecoMass; //higgs mass
      Double_t         fVisMass;  //visible mass
      Double_t         fMT;       //transverse visible mass
      Double_t         fETll;     //transverse energy of tau products
      Double_t         fETnn;     //transverse missing energy
  };
}
#endif
