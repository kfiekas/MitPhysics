//--------------------------------------------------------------------------------------------------
// $Id: GeneratorMod.h,v 1.19 2009/02/13 14:33:31 loizides Exp $
//
// GeneratorMod
//
// This module collects interesting generator information and publishes collections
// for subsequent modules.
//
// Authors: G.Ceballos
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_GENERATORMOD_H
#define MITPHYSICS_MODS_GENERATORMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class GeneratorMod : public BaseMod
  {
    public:
      GeneratorMod(const char *name="GeneratorMod", 
                   const char *title="Generator information module");

      const char          *GetMCPartName()	 const { return fMCPartName; }	
      const char          *GetMCLeptonsName()    const { return fMCLeptonsName; }	
      const char          *GetMCAllLeptonsName() const { return fMCAllLeptonsName; }	
      const char          *GetMCTausName()	 const { return fMCTausName; }	
      const char          *GetMCNeutrinosName()  const { return fMCNeutrinosName; }	
      const char          *GetMCQuarksName()     const { return fMCQuarksName; }	
      const char          *GetMCqqHsName()	 const { return fMCqqHsName; }	
      const char          *GetMCBosonsName()     const { return fMCBosonsName; }	
      const char          *GetMCPhotonsName()    const { return fMCPhotonsName; }	
      void                 SetMCPartName(const char *s)	       { fMCPartName       = s; }	
      void                 SetMCLeptonsName(const char * s)    { fMCLeptonsName    = s; }	
      void                 SetMCAllLeptonsName(const char * s) { fMCAllLeptonsName = s; }	
      void                 SetMCTausName(const char *s)	       { fMCTausName       = s; }	
      void                 SetMCNeutrinosName(const char *s)   { fMCNeutrinosName  = s; }	
      void                 SetMCQuarksName(const char *s)      { fMCQuarksName     = s; }	
      void                 SetMCqqHsName(const char *s)	       { fMCqqHsName       = s; }	
      void                 SetMCBosonsName(const char *s)      { fMCBosonsName     = s; }	
      void                 SetMCPhotonsName(const char *s)     { fMCPhotonsName    = s; }	
      void                 SetPtLeptonMin(Double_t x)          { fPtLeptonMin      = x; }     
      void                 SetEtaLeptonMax(Double_t x)         { fEtaLeptonMax     = x; }     
      void                 SetPtPhotonMin(Double_t x)          { fPtPhotonMin      = x; }     
      void                 SetEtaPhotonMax(Double_t x)         { fEtaPhotonMax     = x; }	  
      void                 SetPdgIdCut(UInt_t d)	       { fPdgIdCut         = d; }	 
      void                 SetMassMinCut(Double_t x)	       { fMassMinCut       = x; }	 
      void                 SetMassMaxCut(Double_t x)	       { fMassMaxCut	   = x; }	 

    protected:
      void                 Process();
      void                 SlaveBegin();

      Bool_t               fFillHist;           //=true then fill histos (def=0)
      TString              fMCPartName;         //name of MCParticle branch
      TString              fMCLeptonsName;      //name of lepton coll (from W)
      TString              fMCAllLeptonsName;   //name of lepton coll (all)
      TString              fMCTausName;         //name of tau coll (hadronic decays)
      TString              fMCNeutrinosName;    //name of neutrinos coll
      TString              fMCQuarksName;       //name of quarks coll
      TString              fMCqqHsName;         //name of qqH coll
      TString              fMCBosonsName;       //name of bosons coll
      TString              fMCPhotonsName;      //name of photons coll
      Double_t             fPtLeptonMin;        //pt min for leptons
      Double_t             fEtaLeptonMax;       //eta max for leptons
      Double_t             fPtPhotonMin;        //pt min for photons
      Double_t             fEtaPhotonMax;       //eta max for photons
      UInt_t               fPdgIdCut;           //pdg id for particle used to select on mass (0=off)
      Double_t             fMassMinCut;	        //mass min for given PdgId particle 
      Double_t             fMassMaxCut;	        //mass max for given PdgId particle
      const MCParticleCol *fParticles;	        //!MCParticle branch

      TH1D                *hDGenLeptons[40];    //!histos for W leptons
      TH1D                *hDGenAllLeptons[20]; //!histos for all leptons
      TH1D                *hDGenTaus[20];       //!histos for taus
      TH1D                *hDGenNeutrinos[20];  //!histos for neutrinos
      TH1D                *hDGenQuarks[20];     //!histos for quarks
      TH1D                *hDGenWBF[20];        //!histos for WBF
      TH1D                *hDGenBosons[20];     //!histos for bosons
      TH1D                *hDGenPhotons[20];    //!histos for photons
      TH1D                *hDVMass[20];         //!histos for auxiliar work

    ClassDef(GeneratorMod, 1) // Module to gather generator information
  };
}
#endif
