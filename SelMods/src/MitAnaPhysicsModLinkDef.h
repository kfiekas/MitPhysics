// $Id: MitAnaPhysicsModLinkDef.h,v 1.1 2008/10/14 06:13:54 loizides Exp $

#ifndef MITPHYSICS_SELMODS_LINKDEF_H
#define MITPHYSICS_SELMODS_LINKDEF_H
#include "MitPhysics/SelMods/interface/FwdJetEvtSelMod.h"
#include "MitPhysics/SelMods/interface/HwwEvtPreSelMod.h"
#include "MitPhysics/SelMods/interface/HwwEvtSelMod.h"
#include "MitPhysics/SelMods/interface/ZXEvtSelMod.h"
#include "MitPhysics/SelMods/interface/ttEvtSelMod.h"
#endif
 
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::FwdJetEvtSelMod+;
#pragma link C++ class mithep::HwwEvtPreSelMod+;
#pragma link C++ class mithep::HwwEvtSelMod+;
#pragma link C++ class mithep::ZXEvtSelMod+;
#pragma link C++ class mithep::ttEvtSelMod+;
#endif
