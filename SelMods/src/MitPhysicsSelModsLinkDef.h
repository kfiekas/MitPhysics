// $Id: MitPhysicsSelModsLinkDef.h,v 1.1 2008/12/10 19:51:05 loizides Exp $

#ifndef MITPHYSICS_SELMODS_LINKDEF_H
#define MITPHYSICS_SELMODS_LINKDEF_H
#include "MitPhysics/SelMods/interface/GenericSelMod.h"
#endif
 
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::GenericSelMod<mithep::Particle>+;
#pragma link C++ class mithep::GenericSelMod<mithep::CompositeParticle>+;
#pragma link C++ class mithep::GenericSelMod<mithep::Conversion>+;
#pragma link C++ class mithep::GenericSelMod<mithep::Electron>+;
#pragma link C++ class mithep::GenericSelMod<mithep::GenJet>+;
#pragma link C++ class mithep::GenericSelMod<mithep::Jet>+;
#pragma link C++ class mithep::GenericSelMod<mithep::MCParticle>+;
#pragma link C++ class mithep::GenericSelMod<mithep::Met>+;
#pragma link C++ class mithep::GenericSelMod<mithep::Muon>+;
#pragma link C++ class mithep::GenericSelMod<mithep::Photon>+;
#pragma link C++ class mithep::GenericSelMod<mithep::Track>+;
#pragma link C++ class mithep::GenericSelMod<mithep::TriggerObject>+;
#endif
