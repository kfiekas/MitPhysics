// $Id: MitPhysicsModsLinkDef.h,v 1.9 2009/04/08 10:13:37 ceballos Exp $

#ifndef MITPHYSICS_MODS_LINKDEF_H
#define MITPHYSICS_MODS_LINKDEF_H
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/MergerMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/HKFactorProducer.h"
#include "MitPhysics/Mods/interface/TauIDMod.h"
#include "MitPhysics/Mods/interface/TauCleaningMod.h"
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
#pragma link C++ class mithep::MergeLeptonsMod+;
#pragma link C++ class mithep::MergerMod<TObject>+;
#pragma link C++ class mithep::MergerMod<mithep::DataBase>+;
#pragma link C++ class mithep::MergerMod<mithep::DataObject>+;
#pragma link C++ class mithep::MergerMod<mithep::Particle>+;
#pragma link C++ class mithep::MuonIDMod+;
#pragma link C++ enum mithep::MuonIDMod::EMuClassType;
#pragma link C++ enum mithep::MuonIDMod::EMuIdType;
#pragma link C++ enum mithep::MuonIDMod::EMuIsoType;
#pragma link C++ class mithep::PhotonCleaningMod+;
#pragma link C++ class mithep::PhotonIDMod+;
#pragma link C++ class mithep::HKFactorProducer+;
#pragma link C++ class mithep::TauIDMod+;
#pragma link C++ class mithep::TauCleaningMod+;
#endif
