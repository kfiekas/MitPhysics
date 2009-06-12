//--------------------------------------------------------------------------------------------------
// $Id: DiTauSystem.h,v 1.8 2009/03/23 22:17:04 loizides Exp $
//
// DiTauSystem
//
// Class to calculate the mass to the di-tau system. It is assumed that the tau is boosted 
// and that the neutrinos have the same flight direction as the tau. See CMS note 2006/082.
//
// Authors: G.Ceballos 
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_DITAUSYSTEM_H
#define MITPHYSICS_UTILS_DITAUSYSTEM_H

#include <Rtypes.h>

namespace mithep
{
  class Particle;
  class Met;

  class DiTauSystem {
    public:
      DiTauSystem(const Particle *t1, const Particle *t2, const Met *met);

      Double_t         RecoMass()        const { return fRecoMass;}
      Double_t         TransverseMass()  const { return fMT;      }
      Double_t         TransverseEll()   const { return fETll;    }
      Double_t         TransverseEnn()   const { return fETnn;    }
      Double_t         VisMass()         const { return fVisMass; }
      Double_t         XTau1()           const { return fXTau[0]; }
      Double_t         XTau2()           const { return fXTau[1]; }
  
    private:
      void             Init();
  
      const Particle  *fT1;       //first tau
      const Particle  *fT2;       //second tau
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
