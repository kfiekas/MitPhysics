// $Id: MitPhysicsUtilsLinkDef.h,v 1.2 2010/03/27 13:38:44 sixie Exp $

#ifndef MITPHYSICS_UTILS_LINKDEF_H
#define MITPHYSICS_UTILS_LINKDEF_H

#include "MitPhysics/Utils/interface/DiTauSystem.h"
#include "MitPhysics/Utils/interface/GeneratorTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/MatchingTools.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::DiTauSystem;
#pragma link C++ class mithep::GeneratorTools;
#pragma link C++ class mithep::IsolationTools;
#pragma link C++ class mithep::MatchingTools;
#pragma link C++ class mithep::MuonTools;
#pragma link C++ class mithep::ElectronTools;
#pragma link C++ enum mithep::ElectronTools::EElIdType;
#pragma link C++ enum mithep::ElectronTools::EElIsoType;

#endif
