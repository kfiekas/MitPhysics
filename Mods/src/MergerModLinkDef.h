// $Id: MitPhysicsModsLinkDef.h,v 1.14 2009/04/30 08:09:32 loizides Exp $

#ifndef MITPHYSICS_MODS_MERGERMODLINKDEF_H
#define MITPHYSICS_MODS_MERGERMODLINKDEF_H
#include "MitPhysics/Mods/interface/MergerMod.h"
#include "MitAna/DataTree/interface/DataBase.h"
#include "MitAna/DataTree/interface/DataObject.h"
#include "MitAna/DataTree/interface/Particle.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::MergerMod<TObject>+;
#pragma link C++ class mithep::MergerMod<mithep::DataBase>+;
#pragma link C++ class mithep::MergerMod<mithep::DataObject>+;
#pragma link C++ class mithep::MergerMod<mithep::Particle>+;
#endif
