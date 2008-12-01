//--------------------------------------------------------------------------------------------------
// $Id: ModNames.h,v 1.1 2008/11/27 16:28:34 loizides Exp $
//
// Names
//
// This class defines the standard names for branches,
// collections and what else we will standardize.
//
// Authors: C.Loizides, C.Paus
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_INIT_MODNAMES_H
#define MITPHYSICS_INIT_MODNAMES_H
 
#include "MitAna/DataTree/interface/Names.h"

namespace mithep 
{
  class ModNames 
  {
    public:
      static const char *gkGoodElectronsName;
      static const char *gkGoodMuonsName;
      static const char *gkGoodPhotonsName;
      static const char *gkGoodTausName;
      static const char *gkGoodJetsName;
      static const char *gkCleanElectronsName;
      static const char *gkCleanMuonsName;
      static const char *gkCleanPhotonsName;
      static const char *gkCleanTausName;
      static const char *gkCleanJetsName;
      static const char *gkCleanFwdJetsName;
      static const char *gkCleanNoFwdJetsName;
      static const char *gkMCLeptonsName;
      static const char *gkMCAllLeptonsName;
      static const char *gkMCTausName;
      static const char *gkMCNeutrinosName;
      static const char *gkMCQuarksName;
      static const char *gkMCqqHsName;
      static const char *gkMCBosonsName;
      static const char *gkMCPhotonsName;
  };
}
#endif
