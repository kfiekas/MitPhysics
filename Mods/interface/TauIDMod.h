//--------------------------------------------------------------------------------------------------
// $Id: TauIDMod.h,v 1.3 2009/04/09 08:45:48 loizides Exp $
//
// TauIDMod
//
// This module applies tau identification criteria and exports a pointer to a collection
// of "good Taus" according to the specified identification scheme.
//
// Authors: G.Ceballos
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_TauIDMOD_H
#define MITPHYSICS_MODS_TauIDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Collections.h"

namespace mithep 
{
  class TauIDMod : public BaseMod
  {
    public:
      TauIDMod(const char *name="TauIDMod", 
               const char *title="Tau identification module");

      const char         *GetCaloTausName()        const { return fCaloTausName;        }   
      Double_t            GetEnergyFractionEmMax() const { return fEnergyFractionEmMax; }
      const char         *GetGoodTausName()        const { return fGoodTausName;        }   
      Double_t            GetHCalEtOverPtMin()     const { return fHCalEtOverPtMin;     }
      Double_t            GetIsoTrackPtSumMax()    const { return fIsoTrackPtSumMax;    }
      Double_t            GetJetPtMin()            const { return fJetPtMin;            }
      Double_t            GetLeadTrackPtMin()      const { return fLeadTrackPtMin;      }
      UInt_t              GetNIsoTracksMax()	   const { return fNIsoTracksMax;	}
      UInt_t              GetNSignalTracksMax()    const { return fNSignalTracksMax;	}
      Double_t            GetPtMin()               const { return fPtMin;               }
      Double_t            GetSignalTracksMassMax() const { return fSignalTracksMassMax; }
      void                SetCaloTausName(const char *n)    { fCaloTausName	   = n; } 
      void                SetEnergyFractionEmMax(Double_t x){ fEnergyFractionEmMax = x; }
      void                SetGoodTausName(const char *n)    { fGoodTausName	   = n; } 
      void                SetHCalEtOverPtMin(Double_t x)    { fHCalEtOverPtMin     = x; }
      void                SetIsoTrackPtSumMax(Double_t x)   { fIsoTrackPtSumMax    = x; }
      void                SetJetPtMin(Double_t x)           { fJetPtMin 	   = x; }
      void                SetLeadTrackPtMin(Double_t x)     { fLeadTrackPtMin      = x; }
      void                SetNIsoTracksMax(Int_t d)         { fNIsoTracksMax	   = d; }
      void                SetNSignalTracksMax(Int_t d)      { fNSignalTracksMax    = d; }
      void                SetPtMin(Double_t x)              { fPtMin 	           = x; }
      void                SetSignalTracksMassMax(Double_t x){ fSignalTracksMassMax = x; }

    protected:
      void                Process();
      void                SlaveBegin();

      TString             fCaloTausName;        //name of tau collection (input)
      TString             fGoodTausName;        //name of exported "good Tau" collection (output)
      Double_t            fPtMin;               //min pt cut
      Double_t            fJetPtMin;            //min jet pt cut
      Double_t            fLeadTrackPtMin;      //min leading track pt cut
      UInt_t              fNSignalTracksMax;	//maximum of signal tracks
      UInt_t              fNIsoTracksMax;	//maximum of iso tracks
      Double_t            fSignalTracksMassMax; //maximum of mass for signal tracks
      Double_t            fIsoTrackPtSumMax;    //maximum of Pt iso tracks
      Double_t            fEnergyFractionEmMax; //maximum of EnergyFractionEm
      Double_t            fHCalEtOverPtMin;     //minimum of HCalEt over Pt for the leading track

      const CaloTauCol   *fCaloTaus;            //!tau branch
    
    ClassDef(TauIDMod, 1) // Tau identification module
  };
}
#endif
