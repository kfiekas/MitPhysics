//--------------------------------------------------------------------------------------------------
// $Id: ModNames.h,v 1.8 2009/04/17 13:45:10 ceballos Exp $
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
      static const char *gkCleanCaloMetName;
      static const char *gkCleanElectronsName;
      static const char *gkCleanFwdJetsName;
      static const char *gkCleanJetsName;
      static const char *gkCleanMuonsName;
      static const char *gkCleanNoFwdJetsName;
      static const char *gkCleanPhotonsName;
      static const char *gkCleanTausName;
      static const char *gkGoodElectronsName;
      static const char *gkGoodJetsName;
      static const char *gkGoodMuonsName;
      static const char *gkGoodPhotonsName;
      static const char *gkGoodTausName;
      static const char *gkMCAllLeptonsName;
      static const char *gkMCBosonsName;
      static const char *gkMCISRPhotonsName;
      static const char *gkMCLeptonsName;
      static const char *gkMCMETName;
      static const char *gkMCNeutrinosName;
      static const char *gkMCPhotonsName;
      static const char *gkMCQuarksName;
      static const char *gkMCRadPhotonsName;
      static const char *gkMCTausName;
      static const char *gkMCqqHsName;
      static const char *gkMergedLeptonsName;
      static const char *gkPubJetsName;
  };
}
#endif
