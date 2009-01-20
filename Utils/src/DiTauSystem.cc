// $Id: DiTauSystem.cc,v 1.5 2008/12/11 20:29:36 ceballos Exp $

#include "MitPhysics/Utils/interface/DiTauSystem.h"
#include "MitAna/DataTree/interface/CompositeParticle.h"
#include "MitAna/DataTree/interface/Met.h"
#include "MitAna/DataTree/interface/Particle.h"

using namespace mithep;

//--------------------------------------------------------------------------------------------------
DiTauSystem::DiTauSystem(const Particle *t1, const Particle *t2, const Met *met) :
  fT1(t1),
  fT2(t2),
  fMet(met),
  fRecoMass(0),
  fVisMass(0),
  fMT(0),
  fETll(0),
  fETnn(0)
{
  // Constructor.

  Init();
}

//--------------------------------------------------------------------------------------------------
void DiTauSystem::Init()
{
  // Calculate the kinematical variables.

  CompositeParticle tt;
  tt.AddDaughter(fT1);
  tt.AddDaughter(fT2);
  
  CompositeParticle higgs;
  higgs.AddDaughter(fT1);
  higgs.AddDaughter(fT2);
  higgs.AddDaughter(fMet);
  
  Double_t xvar[3];
  xvar[0] = higgs.Px()*fT2->Py()-higgs.Py()*fT2->Px();
  xvar[1] = higgs.Py()*fT1->Px()-higgs.Px()*fT1->Py();
  xvar[2] = fT1->Px()*fT2->Py()-fT1->Py()*fT2->Px();
  
  for (Int_t i=0; i<2; ++i)
    xvar[i]==0 ? fXTau[i]=0 : fXTau[i]=xvar[2]/xvar[i];
  
  fVisMass  = tt.Mass();
  if (fXTau[0] > 0 && fXTau[1] > 0)
    fRecoMass = fVisMass / TMath::Sqrt(fXTau[0]*fXTau[1]);
  else 
    fRecoMass = 0;
  
  Double_t visMassS = fVisMass*fVisMass;
  Double_t ptll     = tt.Pt();
  Double_t ptmis    = fMet->Pt();
  if (visMassS > 0) {
    fETll    = TMath::Sqrt(ptll*ptll   + visMassS);
    fETnn    = TMath::Sqrt(ptmis*ptmis + visMassS);
    fMT      = (fETll+fETnn)*(fETll+fETnn)-(ptll+ptmis)*(ptll+ptmis);
    (fMT > 0) ? fMT=TMath::Sqrt(fMT) : fMT=0;
  }
}
