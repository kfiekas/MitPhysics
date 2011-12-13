// $Id: MVASystematicsMod.cc,v 1.1 2011/11/18 19:10:43 fabstoec Exp $

#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <TRandom3.h>
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitPhysics/Mods/interface/MVASystematicsMod.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"

#include <iostream>
#include <sstream>

using namespace mithep;

ClassImp(mithep::MVASystematicsMod)

//--------------------------------------------------------------------------------------------------
MVASystematicsMod::MVASystematicsMod(const char *name, const char *title) :
  BaseMod             (name,title),

  // define all the Branches to load
  fMCParticleName    (Names::gkMCPartBrn),

  // ----------------------------------------
  // collections....
  fMCParticles       (0),
  // -------------------------------------------
  fIsData            (false),

  fTupleName("hMVAtuple")
  
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void MVASystematicsMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do anything here.
}

//--------------------------------------------------------------------------------------------------
void MVASystematicsMod::Process()
{  
  IncNEventsProcessed();
  if( !fIsData )
    LoadBranch(fMCParticleName);

  Float_t _pth    = -100.;
  Float_t _y      = -100.;
  Float_t _genmass = -100.;
  if( !fIsData ) FindHiggsPtAndY(_pth, _y, _genmass);  

  hMVAtuple->Fill(_pth, _y, _genmass);


  return;
}

//--------------------------------------------------------------------------------------------------
void MVASystematicsMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches or objects created by earlier modules.

  if (!fIsData) {
    ReqBranch(fMCParticleName,        fMCParticles);
  }

  hMVAtuple = new TNtuple(fTupleName.Data(),fTupleName.Data(),"hpt:hy:hm");

  AddOutput(hMVAtuple);
}

//--------------------------------------------------------------------------------------------------
void MVASystematicsMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis. For this
  // module, we dont do anything here.
}

//--------------------------------------------------------------------------------------------------
void MVASystematicsMod::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.
}

// ----------------------------------------------------------------------------------------
// some helpfer functions....
void MVASystematicsMod::FindHiggsPtAndY(Float_t& pt, Float_t& Y, Float_t& mass) {

  pt = -999.;
  Y  = -999;
  mass = -999.;

  // loop over all GEN particles and look for status 1 photons
  for(UInt_t i=0; i<fMCParticles->GetEntries(); ++i) {
    const MCParticle* p = fMCParticles->At(i);
    if( p->Is(MCParticle::kH) ) {
      pt=p->Pt();
      Y = p->Rapidity();
      mass = p->Mass();
      break;
    }
  }
  
  return;
 }
