// $Id: skimtest.C,v 1.1 2008/12/04 14:14:59 loizides Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#include "MitPhysics/SelMods/interface/HwwEvtPreSelMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void runExampleSkim(const char *files, UInt_t nev=0)
{
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  HwwEvtPreSelMod *smod = new HwwEvtPreSelMod;
  smod->SetLeptonMinMaxPt(30);
  OutputMod *omod = new OutputMod;
  omod->Drop("*");
  omod->Keep("Muons");
  omod->Keep("MCParticles");
  omod->Keep("Electrons");
  omod->Keep("Muons");
  omod->Keep("*Tracks*");
  omod->Keep("*Clusters*");
  omod->Drop("*Conversion*");
  omod->Drop("ProtonRefitTracks");
  smod->Add(omod);

  // set up analysis
  Analysis *ana = new Analysis;
  ana->SetSuperModule(smod);
  if (nev)
    ana->SetProcessNEvents(nev);
  ana->AddFile(files);
  
  // run the analysis after successful initialisation
  ana->Run();
}
