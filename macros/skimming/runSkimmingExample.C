// $Id: runSkimmingExample.C,v 1.2 2009/03/11 10:14:52 loizides Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TRandom.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#endif

namespace mithep 
{
  class Sel : public BaseMod
  {
    public:
      Sel(Int_t r=10) : fRand(r) {}
    protected:
      void       Process() { if (gRandom->Rndm()>1./fRand) SkipEvent(); }
      Int_t      fRand;
    ClassDef(Sel, 1)
  };
}

//--------------------------------------------------------------------------------------------------
void runSkimmingExample(const char *files, const char *prefix="skimtest", UInt_t nev=0)
{
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  Sel *smod = new Sel(3);
  OutputMod *omod = new OutputMod;
  omod->SetFileName(prefix);
  omod->Keep("*");
  if (0) {
    omod->Drop("*");
    omod->Keep("Muons");
    omod->Keep("MCParticles");
    omod->Keep("Electrons");
    omod->Keep("*Tracks*");
    omod->Keep("*Clusters*");
  }
  smod->Add(omod);

  // set up analysis
  Analysis *ana = new Analysis;
  ana->SetSuperModule(smod);
  if (nev)
    ana->SetProcessNEvents(nev);
  ana->AddFile(files);
  
  // run the analysis after successful initialisation
  ana->Run(0);
}
