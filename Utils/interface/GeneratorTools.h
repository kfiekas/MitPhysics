//--------------------------------------------------------------------------------------------------
// $Id: GeneratorTools.h,v 1.3 2009/02/17 06:49:01 phedex Exp $
//
// GeneratorTools
//
// Various Generator level tools
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_GENERATORTOOLS_H
#define MITPHYSICS_UTILS_GENERATORTOOLS_H

#include "TString.h"
#include <TMath.h>
#include "MitAna/DataTree/interface/MCParticle.h"
#include "MitAna/DataTree/interface/Electron.h"

namespace mithep
{
  class GeneratorTools {
    public:
      
      static void PrintHepMCTable(const mithep::Collection<mithep::MCParticle> *particles, 
                           Bool_t showDaughters, int suppressEntriesAfterThisIndex);

      static void PrintNearbyParticles(const mithep::Collection<mithep::MCParticle> *particles, 
                                       Double_t eta, Double_t phi, Double_t deltaR);

      static TString ConvertPdgIdToName(Int_t pdgId);

      static const mithep::MCParticle* MatchElectronToSimParticle(
        const mithep::Collection<mithep::MCParticle> *particles,
        const mithep::Track *eleTrack, Bool_t isGsfTrack, 
        Int_t printDebugLevel, Bool_t printHepMCTable);
      
      static const mithep::MCParticle* FindElectronFakeAncestor(
        const mithep::MCParticle *matchedSimChargedParticle);

      static Int_t CategorizeFakeElectron(const mithep::MCParticle *ele);

      static const mithep::MCParticle* MatchMuonToSimParticle(
        const mithep::Collection<mithep::MCParticle> *particles,const mithep::Track *muonTrack, 
        Bool_t isTrackerTrack, Int_t printDebugLevel, Bool_t printHepMCTable);

      static Int_t CategorizeFakeMuon(const mithep::MCParticle *mu);

  };
}
#endif
