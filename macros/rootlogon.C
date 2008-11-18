// $Id: rootlogon.C,v 1.1 2008/07/17 09:36:13 loizides Exp $
{
  gROOT->Macro("$CMSSW_BASE/src/MitAna/macros/setRootEnv.C+");

  loadLibraries("libMitPhysics*.so");
}
