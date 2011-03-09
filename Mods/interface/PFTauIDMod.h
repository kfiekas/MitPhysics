//--------------------------------------------------------------------------------------------------
// $Id: PFTauIDMod.h,v 1.2 2011/02/17 14:08:49 ceballos Exp $
//
// PFTauIDMod
//
// This module applies tau identification criteria and exports a pointer to a collection
// of "good Taus" according to the specified identification scheme.
//
// Authors: G.Ceballos
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PFTauIDMod_H
#define MITPHYSICS_MODS_PFTauIDMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PFTauFwd.h"

namespace mithep 
{
  class PFTauIDMod : public BaseMod
  {
    public:
      PFTauIDMod(const char *name="PFTauIDMod", 
               const char *title="Tau identification module");

      const char         *GetPFTausName()          	    const { return fPFTausName; 		  }
      const char         *GetGoodPFTausName()        	    const { return fGoodPFTausName;		  }
      Double_t            GetPtMin()	     	   	    const { return fPtMin;			  }
      Double_t            GetPtLeadChargedHadronPFCandMin() const { return fPtLeadChargedHadronPFCandMin; }
      UInt_t              GetIsoChargedHadronPtSumMax()     const { return fIsoChargedHadronPtSumMax;     }
      UInt_t              GetIsoGammaEtSumMax()	   	    const { return fIsoGammaEtSumMax;		  }
      Double_t            GetSignalMassMin()   	   	    const { return fSignalMassMin;		  }
      Double_t            GetSignalMassMax()       	    const { return fSignalMassMax;		  }
      Bool_t              GetIsHPSSel()       	            const { return fIsHPSSel;		          }
      void                SetPFTausName(const char *n)               { fPFTausName	    	     = n; }
      void                SetPtMin(Double_t x)                       { fPtMin 	           	     = x; }
      void                SetPtLeadChargedHadronPFCandMin(Double_t x){ fPtLeadChargedHadronPFCandMin = x; }
      void                SetIsoChargedHadronPtSumMax(Double_t x)    { fIsoChargedHadronPtSumMax     = x; }
      void                SetIsoGammaEtSumMax(Double_t x)   	     { fIsoGammaEtSumMax    	     = x; }
      void                SetSignalMassMin(Double_t x)	    	     { fSignalMassMin	    	     = x; }
      void                SetSignalMassMax(Double_t x)	    	     { fSignalMassMax	    	     = x; }
      void                SetIsHPSSel(Bool_t b)	    	             { fIsHPSSel                     = b; }
      void                SetHPSIso(const char *i)                   { fHPSIso                       = i; }

    protected:
      void                Process();
      void                SlaveBegin();

      TString             fPFTausName;               	 //name of tau collection (input)
      TString             fGoodPFTausName;           	 //name of exported good PFTau collection (output)
      Double_t            fPtMin;                    	 //min tau pt cut
      Double_t            fPtLeadChargedHadronPFCandMin; //min LeadChargedHadronPFCand pt cut
      Double_t            fIsoChargedHadronPtSumMax; 	 //maximum of Pt iso tracks
      Double_t            fIsoGammaEtSumMax;         	 //maximum of Pt iso neutrals
      Double_t            fSignalMassMin;            	 //minimum of tau mass
      Double_t            fSignalMassMax;            	 //maximum of tau mass
      Bool_t              fIsHPSSel;            	 //apply HPS Tau selection
      TString             fHPSIso;                       //isolation tightness: "tight", "medium", "loose"
      const PFTauCol     *fPFTaus;                       //!tau branch
    
    ClassDef(PFTauIDMod, 1) // Tau identification module
  };
}
#endif
