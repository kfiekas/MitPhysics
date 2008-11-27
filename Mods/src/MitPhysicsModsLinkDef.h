// $Id: MitPhysicsModsLinkDef.h,v 1.2 2008/11/19 15:44:50 loizides Exp $

#ifndef MITPHYSICS_MODS_LINKDEF_H
#define MITPHYSICS_MODS_LINKDEF_H
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::ElectronCleaningMod+;
#pragma link C++ class mithep::ElectronIDMod+;
#pragma link C++ enum mithep::ElectronIDMod::EElIdType;
#pragma link C++ enum mithep::ElectronIDMod::EElIsoType;
#pragma link C++ class mithep::GeneratorMod+;
#pragma link C++ class mithep::JetCleaningMod+;
#pragma link C++ class mithep::JetIDMod+;
#pragma link C++ class mithep::MuonIDMod+;
#pragma link C++ enum mithep::MuonIDMod::EMuIdType;
#pragma link C++ enum mithep::MuonIDMod::EMuIsoType;
#pragma link C++ enum mithep::MuonIDMod::EMuClassType;
#endif
