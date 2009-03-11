// $Id: runSkimmingExample.C,v 1.1 2008/12/10 17:31:35 loizides Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#endif

namespace mithep 
{
  class Sel : public BaseMod
  {
    protected:
      void Process() {SkipEvent(); }

    ClassDef(Sel, 1)
  };
}

//--------------------------------------------------------------------------------------------------
void runSkimmingExample(const char *files, UInt_t nev=0)
{
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  Sel *selmod = new Sel;
  OutputMod *omod = new OutputMod;
  omod->Keep("*");
//  omod->Drop("*");
//  omod->Keep("Muons");
//  omod->Keep("MCParticles");
  // omod->Keep("Electrons");
  // omod->Keep("*Tracks*");
  //omod->Keep("*Clusters*");
  //omod->Drop("*Conversion*");
  //omod->Drop("ProtonRefitTracks");
  selmod->Add(omod);

  // set up analysis
  Analysis *ana = new Analysis;
  ana->SetSuperModule(selmod);
  if (nev)
    ana->SetProcessNEvents(nev);
  ana->AddFile(files);
  
  // run the analysis after successful initialisation
  ana->Run(0);
}
