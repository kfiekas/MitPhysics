// $Id: MitCommonMathToolsLinkDef.h,v 1.4 2009/07/20 03:12:22 loizides Exp $

#ifndef MITPHYSICS_UTILS_LINKDEF_H
#define MITPHYSICS_UTILS_LINKDEF_H

#include "MitPhysics/Utils/interface/DiTauSystem.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/MatchingTools.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::DiTauSystem;
#pragma link C++ class mithep::IsolationTools;
#pragma link C++ class mithep::MatchingTools;
#pragma link C++ class mithep::MuonTools;

#endif
