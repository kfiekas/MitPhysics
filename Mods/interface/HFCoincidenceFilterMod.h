//--------------------------------------------------------------------------------------------------
// $Id $
//
// HFCoincidenceFilterMod
//
// This module selects events with HF Calotowers on both sides of the collision center
// with configurable cuts on the energy threshold of the calotowers and the number of towers found
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITMODS_MODS_HFCOINCIDENCEFILTERMOD_H
#define MITMODS_MODS_HFCOINCIDENCEFILTERMOD_H

#include <string>
#include <TString.h>
#include <TH1F.h>
#include "MitAna/DataTree/interface/VertexFwd.h" 
#include "MitAna/DataTree/interface/CaloTowerFwd.h" 
#include "MitAna/TreeMod/interface/BaseMod.h" 

namespace mithep 
{
  class HFCoincidenceFilterMod : public BaseMod {
    public:
      
      enum ECuts {
        eNTracks,
        eZ,
        eRho
      };
      
      HFCoincidenceFilterMod(const char *name="HFCoincidenceFilterMod", const char *title="Good PV Filter Module");
      ~HFCoincidenceFilterMod();

      Int_t                       GetNEvents()      const { return fNEvents;       }
      Int_t                       GetNAccepted()    const { return fNAcceped;      }
      Int_t                       GetNFailed()      const { return fNFailed;       }
      void                        SetAbortIfNotAccepted(Bool_t b)   { fAbort         = b; }
      void                        SetMinNCoincidentCaloTowers(UInt_t n)     { fMinNCoincidentCaloTowers = n; }
      void                        SetCaloTowerEnergyThreshold(Double_t x)   { fCaloTowerEnergyThreshold = x; }

    protected:
      void                        BeginRun();
      virtual void                OnAccepted()  {/*could be implemented in derived classes*/}
      virtual void                OnFailed()    {/*could be implemented in derived classes*/}
      void                        Process();
      void                        SlaveBegin();
      void                        SlaveTerminate();

      Bool_t                      fAbort;         //=true then abort (sub-)modules if not accepted
      UInt_t                      fMinNCoincidentCaloTowers; //minimum number of tracks for the vertex
      Double_t                    fCaloTowerEnergyThreshold; //maximum abs(z) of the vertex
      TString                     fCaloTowersName;  //Name of CaloTower collection
      Int_t                       fNEvents;         //!number of processed events
      Int_t                       fNAcceped;        //!number of accepted events
      Int_t                       fNFailed;         //!number of failed events
      const CaloTowerCol         *fCaloTowers;      //!CaloTower collection

    ClassDef(HFCoincidenceFilterMod, 1) // L1 TAM module
  };
}
#endif
