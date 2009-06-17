// $Id: MitPhysicsModsLinkDef.h,v 1.14 2009/04/30 08:09:32 loizides Exp $

#ifndef MITPHYSICS_MODS_MATCHMODLINKDEF_H
#define MITPHYSICS_MODS_MATCHMODLINKDEF_H

#include "MitAna/DataTree/interface/Electron.h"
#include "MitAna/DataTree/interface/Muon.h"
#include "MitAna/DataTree/interface/Particle.h"
#include "MitAna/DataTree/interface/TriggerObject.h"
#include "MitPhysics/Mods/interface/MatchMod.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::MatchMod<mithep::Particle>+;
#pragma link C++ class mithep::MatchMod<mithep::Particle,mithep::TriggerObject>+;
#pragma link C++ class mithep::MatchMod<mithep::Electron,mithep::TriggerObject>+;
#pragma link C++ class mithep::MatchMod<mithep::Muon,mithep::TriggerObject>+;
#endif
