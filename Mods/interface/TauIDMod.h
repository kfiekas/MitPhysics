//--------------------------------------------------------------------------------------------------
// $Id: TauIDMod.h,v 1.1 2009/04/08 10:11:44 ceballos Exp $
//
// TauIDMod
//
// This module applies Tau identification criteria and exports a pointer to a collection
// of "good Taus" according to the specified identification scheme.
//
// Authors: S.Xie, C.Loizides
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
      const char         *GetGoodTausName()        const { return fGoodTausName;        }   
      Double_t            GetPtMin()               const { return fTauPtMin;            }
      Double_t            GetJetPtMin()            const { return fTauJetPtMin;         }
      UInt_t              GetNSignalTracksMax()    const { return fNSignalTracksMax;	}
      UInt_t              GetNIsoTracksMax()	   const { return fNIsoTracksMax;	}
      Double_t            GetSignalTracksMassMax() const { return fSignalTracksMassMax; }
      Double_t            GetIsoTrackPtSumMax()    const { return fIsoTrackPtSumMax;    }
      Double_t            GetEnergyFractionEmMax() const { return fEnergyFractionEmMax; }

      void                SetCaloTausName(const char *n)    { fCaloTausName	   = n; } 
      void                SetGoodTausName(const char *n)    { fGoodTausName	   = n; } 
      void                SetPtMin(Double_t x)              { fTauPtMin 	   = x; }
      void                SetJetPtMin(Double_t x)           { fTauJetPtMin 	   = x; }
      void                SetNSignalTracksMax(Int_t d)      { fNSignalTracksMax    = d; }
      void                SetNIsoTracksMax(Int_t d)         { fNIsoTracksMax	   = d; }
      void                SetSignalTracksMassMax(Double_t x){ fSignalTracksMassMax = x; }
      void                SetIsoTrackPtSumMax(Double_t x)   { fIsoTrackPtSumMax    = x; }
      void                SetEnergyFractionEmMax(Double_t x){ fEnergyFractionEmMax = x; }

    protected:
      void                Process();
      void                SlaveBegin();

      TString             fCaloTausName;        //name of Tau collection (input)
      const CaloTauCol   *fCaloTaus;            //!Tau branch
      TString             fGoodTausName;        //name of exported "good Tau" collection
      Double_t            fTauPtMin;            //min pt cut
      Double_t            fTauJetPtMin;         //min jet pt cut
      UInt_t              fNSignalTracksMax;	//maximum of signal tracks
      UInt_t              fNIsoTracksMax;	//maximum of iso tracks
      Double_t            fSignalTracksMassMax; //maximum of mass for signal tracks
      Double_t            fIsoTrackPtSumMax;    //maximum of Pt iso tracks
      Double_t            fEnergyFractionEmMax; //maximum of EnergyFractionEm
    
    ClassDef(TauIDMod, 1) // Tau identification module
  };
}
#endif
