// $Id: runSimpleExample.C,v 1.8 2008/11/25 14:31:19 loizides Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/PhysicsMod/interface/FullExampleMod.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#include "MitPhysics/SelMods/interface/HwwEvtPreSelMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void skimtest(const char *files = "")
{
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  HwwEvtPreSelMod *smod = new HwwEvtPreSelMod;
  smod->SetLeptonMinMaxPt(30);
  FullExampleMod *fmod = new FullExampleMod;
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
  smod->Add(fmod);
  fmod->Add(omod);

  // set up analysis
  Analysis *ana = new Analysis;
  ana->SetSuperModule(smod);
  ana->SetProcessNEvents(1000);
  ana->AddFile(files);
  if (gROOT->IsBatch()) 
    ana->SetOutputName("mit-example-hist.root");
  
  // run the analysis after successful initialisation
  ana->Run(!gROOT->IsBatch());
}
