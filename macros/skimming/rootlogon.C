// $Id: rootlogon.C,v 1.1 2008/11/18 21:32:51 loizides Exp $
{
  gROOT->Macro("$CMSSW_BASE/src/MitAna/macros/setRootEnv.C+");

  loadLibraries("libMitPhysics*.so");
}
