//--------------------------------------------------------------------------------------------------
// $Id: PhotonCiCMod.h,v 1.3 2011/06/29 18:28:05 fabstoec Exp $
//
// PhotonCiCMod
//
// This module applies photon identification criteria and exports a pointer to a collection
// of "good photons" according to the specified identification scheme.
//
// Authors: S.Xie, C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PHOTONCICMOD_H
#define MITPHYSICS_MODS_PHOTONCICMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PhotonFwd.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"

class TNtuple;
class TRandom3;

namespace mithep 
{
  class PhotonCiCMod : public BaseMod
  {
    public:
      PhotonCiCMod(const char *name="PhotonCiCMod", 
                  const char *title="Photon identification module");

      ~PhotonCiCMod();

      Bool_t              GetApplySpikeRemoval()      const { return fApplySpikeRemoval;   }
      const char         *GetGoodName()               const { return GetGoodPhotonsName(); }   
      const char         *GetGoodPhotonsName()        const { return fGoodPhotonsName;     }   
      const char         *GetInputName()              const { return fPhotonBranchName;    }   
      const char         *GetOutputName()             const { return GetGoodPhotonsName(); }   
      Double_t            GetPtMin()                  const { return fPhotonPtMin;         }
      Double_t            GetAbsEtaMax()	      const { return fAbsEtaMax;	   }
      void                SetApplySpikeRemoval(Bool_t b)    { fApplySpikeRemoval  = b;     }
      void                SetGoodName(const char *n)        { SetGoodPhotonsName(n);       }   
      void                SetGoodPhotonsName(const char *n) { fGoodPhotonsName = n;        }   
      void                SetInputName(const char *n)       { fPhotonBranchName= n;        }   
      void                SetTrackName(const char *n)       { fTrackBranchName = n;        }   
      void                SetOutputName(const char *n)      { SetGoodPhotonsName(n);       }    
      void                SetPtMin(Double_t pt)             { fPhotonPtMin     = pt;       }
      void                SetAbsEtaMax(Double_t x)          { fAbsEtaMax       = x;	   }

      void                SetPVName  (TString s) {fPVName   = s; fPVFromBranch   = false;};

      void                SetIsData (Bool_t b) { fIsData = b;};

      double findHiggsPt();

    protected:
      void                Process();
      void                SlaveBegin();
      unsigned int        findBestVertex(Photon* ph1, Photon* ph2, const BaseVertex* bsp, bool print=false);

      TString             fPhotonBranchName;     //name of photon collection (input)
      TString             fGoodPhotonsName;      //name of exported "good photon" collection
      TString             fTrackBranchName;      // name of the track collection (only needed for PU corrected isolation)
      TString             fPileUpDenName;        //name of the PU density collection      
      TString             fElectronName;
      Double_t            fPhotonPtMin;          //min pt cut
      Bool_t              fApplySpikeRemoval;    //whether apply spike removal      
      Double_t            fAbsEtaMax;  	         //max Abs Eta
      const PhotonCol    *fPhotons;              //!photon branch
      const TrackCol     *fTracks;               //!track branch
      const PileupEnergyDensityCol *fPileUpDen;  //!rho branch
      const ElectronCol  *fElectrons;            //!electron branch
      
      TString fPVName;
      const VertexCol *fPV;
      Bool_t fPVFromBranch;

      const DecayParticleCol *fConversions;
      TString fConversionName;

      const BeamSpotCol *fBeamspot;
     
      Double_t fDataEnCorr_EB_hR9;
      Double_t fDataEnCorr_EB_lR9;
      Double_t fDataEnCorr_EE_hR9;
      Double_t fDataEnCorr_EE_lR9;

      Double_t fMCSmear_EB_hR9;
      Double_t fMCSmear_EB_lR9;
      Double_t fMCSmear_EE_hR9;
      Double_t fMCSmear_EE_lR9;

      Bool_t fIsData;

      TNtuple* hCiCTuple;

      TRandom3* rng;

      TString fMCParticleName;
      const MCParticleCol* fMCParticles;

      TString                  fPileUpName;
      const PileupInfoCol     *fPileUp;

    ClassDef(PhotonCiCMod, 1) // Photon identification module
  };
}
#endif
