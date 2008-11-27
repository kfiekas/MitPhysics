//--------------------------------------------------------------------------------------------------
// $Id: GeneratorMod.h,v 1.6 2008/11/24 14:13:20 loizides Exp $
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
      ~GeneratorMod() {}

      Bool_t         GetFillHist()         const { return fFillHist; }
      const char    *GetMCPartName()	   const { return fMCPartName; }	
      const char    *GSetMCLeptonsName()   const { return fMCLeptonsName; }	
      const char    *GetMCAllLeptonsName() const { return fMCAllLeptonsName; }	
      const char    *GetMCTausName()	   const { return fMCTausName; }	
      const char    *GetMCNeutrinosName()  const { return fMCNeutrinosName; }	
      const char    *GetMCQuarksName()     const { return fMCQuarksName; }	
      const char    *GetMCqqHsName()	   const { return fMCqqHsName; }	
      const char    *GetMCBosonsName()     const { return fMCBosonsName; }	
      void           SetFillHist(Bool_t b)	         { fFillHist	     = b; }
      void           SetMCPartName(const char *s)	 { fMCPartName       = s; }	
      void           SetMCLeptonsName(const char * s)    { fMCLeptonsName    = s; }	
      void           SetMCAllLeptonsName(const char * s) { fMCAllLeptonsName = s; }	
      void           SetMCTausName(const char *s)	 { fMCTausName       = s; }	
      void           SetMCNeutrinosName(const char *s)   { fMCNeutrinosName  = s; }	
      void           SetMCQuarksName(const char *s)      { fMCQuarksName     = s; }	
      void           SetMCqqHsName(const char *s)	 { fMCqqHsName       = s; }	
      void           SetMCBosonsName(const char *s)      { fMCBosonsName     = s; }	

    protected:
      Bool_t         fFillHist; 		//=true then fill histos (def=0)
      TString        fMCPartName;		//name of MCParticle branch
      TString        fMCLeptonsName ;		//name of lepton coll (from W)
      TString        fMCAllLeptonsName ;	//name of lepton coll (all)
      TString        fMCTausName;		//name of tau coll (hadronic decays)
      TString        fMCNeutrinosName;  	//name of neutrinos coll
      TString        fMCQuarksName;		//name of quarks coll
      TString        fMCqqHsName;		//name of qqH coll
      TString        fMCBosonsName;		//name of bosons coll
      MCParticleCol *fParticles;		//!MCParticle branch
      TH1D          *hDGenLeptons[20];          //!histos for W leptons
      TH1D          *hDGenAllLeptons[20];       //!histos for all leptons
      TH1D          *hDGenTaus[20];             //!histos for taus
      TH1D          *hDGenNeutrinos[20];        //!histos for neutrinos
      TH1D          *hDGenQuarks[20];           //!histos for quarks
      TH1D          *hDGenWBF[20];              //!histos for WBF
      TH1D          *hDGenBosons[20];           //!histos for bosons

      void           Process();
      void           SlaveBegin();
    
      ClassDef(GeneratorMod,1) // Module to gather generator information
  };
}
#endif
