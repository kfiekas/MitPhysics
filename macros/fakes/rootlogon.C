// $Id: rootlogon.C,v 1.1 2008/12/10 17:31:35 loizides Exp $
{
  gROOT->Macro("$CMSSW_BASE/src/MitAna/macros/setRootEnv.C+");

  loadLibraries("libMitPhysics*.so");
}
