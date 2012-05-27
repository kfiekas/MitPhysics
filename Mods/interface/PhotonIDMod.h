//--------------------------------------------------------------------------------------------------
// $Id: PhotonIDMod.h,v 1.25 2012/01/13 15:27:29 fabstoec Exp $
//
// PhotonIDMod
//
// This module applies photon identification criteria and exports a pointer to a collection
// of "good photons" according to the specified identification scheme.
//
// Authors: S.Xie, C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PHOTONIDMOD_H
#define MITPHYSICS_MODS_PHOTONIDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PhotonFwd.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"

#include "MitPhysics/Utils/interface/MVATools.h"

namespace mithep 
{
  class PhotonIDMod : public BaseMod
  {
    public:
      PhotonIDMod(const char *name="PhotonIDMod", 
                  const char *title="Photon identification module");

      Bool_t              GetApplySpikeRemoval()      const { return fApplySpikeRemoval;   }
      Bool_t              GetApplyPixelSeed()         const { return fApplyPixelSeed;      }
      const char         *GetGoodName()               const { return GetGoodPhotonsName(); }   
      const char         *GetGoodPhotonsName()        const { return fGoodPhotonsName;     }   
      Double_t            GetHadOverEmMax()           const { return fHadOverEmMax;        }
      const char         *GetIDType()                 const { return fPhotonIDType;        }
      const char         *GetInputName()              const { return fPhotonBranchName;    }   
      const char         *GetIsoType()                const { return fPhotonIsoType;       }
      const char         *GetOutputName()             const { return GetGoodPhotonsName(); }   
      Double_t            GetPtMin()                  const { return fPhotonPtMin;         }
      Bool_t              GetApplyFiduciality()       const { return fFiduciality;         }
      Double_t            GetEtaWidthEB()	      const { return fEtaWidthEB;	   }
      Double_t            GetEtaWidthEE()	      const { return fEtaWidthEE;	   }
      Double_t            GetAbsEtaMax()	      const { return fAbsEtaMax;	   }
      void                SetApplySpikeRemoval(Bool_t b)    { fApplySpikeRemoval  = b;     }
      void                SetApplyPixelSeed(Bool_t b)       { fApplyPixelSeed  = b;        }
      void                SetApplyElectronVeto(Bool_t b)    { fApplyElectronVeto = b;      }
      void                SetInvertElectronVeto(Bool_t b)   { fInvertElectronVeto = b;     }      
      void                SetApplyElectronVetoConvRecovery(Bool_t b) { fApplyElectronVetoConvRecovery = b; }
      void                SetApplyConversionId(Bool_t b)    { fApplyConversionId = b;      }
      void                SetApplyTriggerMatching(Bool_t b)      { fApplyTriggerMatching = b;  }      
      void                SetGoodName(const char *n)        { SetGoodPhotonsName(n);       }   
      void                SetGoodPhotonsName(const char *n) { fGoodPhotonsName = n;        }   
      void                SetHadOverEmMax(Double_t hoe)     { fHadOverEmMax    = hoe;      }
      void                SetIDType(const char *type)       { fPhotonIDType    = type;     }
      void                SetInputName(const char *n)       { fPhotonBranchName= n;        }   
      void                SetTrackName(const char *n)       { fTrackBranchName = n;        }   
      void                SetBeamspotName(const char *n)    { fBeamspotBranchName = n;     }   
      void                SetIsoType(const char *type)      { fPhotonIsoType   = type;     }
      void                SetOutputName(const char *n)      { SetGoodPhotonsName(n);       }    
      void                SetPtMin(Double_t pt)             { fPhotonPtMin     = pt;       }
      void                SetR9Min(Double_t x)              { fPhotonR9Min     = x;        }
      void                SetEtaWidthEB(Double_t x)	    { fEtaWidthEB      = x;	   }
      void                SetEtaWidthEE(Double_t x)         { fEtaWidthEE      = x;	   }
      void                SetAbsEtaMax(Double_t x)          { fAbsEtaMax       = x;	   }
      void                SetApplyR9Min(Bool_t b)           { fApplyR9Min      = b;        }
      void                SetApplyFiduciality(Bool_t b)     { fFiduciality = b;            }      
      void                SetEffAreas(Double_t ecal, Double_t hcal, Double_t track) { 
	fEffAreaEcalEE = ecal; fEffAreaHcalEE = hcal; fEffAreaTrackEE = track;
	fEffAreaEcalEB = ecal; fEffAreaHcalEB = hcal; fEffAreaTrackEB = track;
      }
      void                SetEffAreasEEEB(Double_t ecalEE, Double_t hcalEE, Double_t trackEE,
					  Double_t ecalEB, Double_t hcalEB, Double_t trackEB) { 
	fEffAreaEcalEE = ecalEE; fEffAreaHcalEE = hcalEE; fEffAreaTrackEE = trackEE;
	fEffAreaEcalEB = ecalEB; fEffAreaHcalEB = hcalEB; fEffAreaTrackEB = trackEB;
      }
      void                SetTriggerObjectsName(const char *n)   { fTrigObjectsName = n;       }
    

    void                SetPhotonsFromBranch(bool b)           { fPhotonsFromBranch = b;           }

    void                SetPVName(const char *n)          { fPVName = n;                 }
    void                SetPVFromBranch(bool b)           { fPVFromBranch = b;           }
    void                SetIsData (Bool_t b) { fIsData = b;};
    void                SetGoodElectronsFromBranch(Bool_t b) { fGoodElectronsFromBranch = b; }
    void                SetGoodElectronName(TString name) { fGoodElectronName = name; }

    void                SetBdtCutBarrel(double a) {fbdtCutBarrel = a; }
    void                SetBdtCutEndcap(double a) {fbdtCutEndcap = a; }

    void                SetDoMCR9Scaling(Bool_t b)        { fDoMCR9Scaling = b; }
    void                SetMCR9Scale(Double_t ebscale, Double_t eescale) { fMCR9ScaleEB = ebscale; fMCR9ScaleEE = eescale; }
    void                SetDoMCSigIEtaIEtaScaling(Bool_t b)        { fDoMCSigIEtaIEtaScaling = b; }
    void                SetDoMCWidthScaling(Bool_t b)        { fDoMCWidthScaling = b; }
    void                SetDoMCErrScaling(Bool_t b)        { fDoMCErrScaling = b; }
    void                SetMCErrScale(Double_t ebscale, Double_t eescale) { fMCErrScaleEB = ebscale; fMCErrScaleEE = eescale; }

      enum EPhIdType {
        kIdUndef = 0,       //not defined
        kTight,             //"Tight"
        kLoose,             //"Loose"
        kLooseEM,           //"LooseEM"
	kBaseLineCiC,        //"2011" Hgg BaseLine CiC
	kBaseLineCiCPF,        //"2012" Hgg BaseLine CiC
	kMITMVAId,          // MingMing MVA ID
	kMITPhSelection,    //MIT loose preselection (for mva)
	kMITPFPhSelection,    //MIT loose preselection (for mva)
        kCustomId           //"Custom"
      };

      enum EPhIsoType {
        kIsoUndef = 0,      //not defined        
        kNoIso,             //"NoIso"
        kCombinedIso,       //"CombinedIso"
        kCustomIso,         //"Custom"
	kMITPUCorrected     //PileUp Corrected Hgg Isolation
      };

    protected:
      void                Process();
      void                SlaveBegin();

      TString             fPhotonBranchName;     //name of photon collection (input)
      TString             fGoodPhotonsName;      //name of exported "good photon" collection
      TString             fTrackBranchName;      // name of the track collection (only needed for PU corrected isolation)
      TString             fBeamspotBranchName;   //name of the Beamspot collection (only needed for PU corrected isolation)
      TString             fPileUpDenName;        //name of the PU density collection      
      TString             fConversionName;       //name of conversion branch
      TString             fElectronName;
      TString             fGoodElectronName;
      TString             fTrigObjectsName;        //name of trigger object collection
      TString             fPVName;
      TString             fMCParticleName;
      TString             fPileUpName;
      TString             fPFCandsName;
      TString             fPhotonIDType;         //type of photon identification we impose
      TString             fPhotonIsoType;        //type of photon isolation we impose
      Double_t            fPhotonPtMin;          //min pt cut
      Double_t            fHadOverEmMax;         //maximum of hadronic/em energy
      Bool_t              fApplySpikeRemoval;    //whether apply spike removal      
      Bool_t              fApplyPixelSeed;       //=true then apply pixel seed constraint
      Bool_t              fApplyElectronVeto;    //=true then apply electron veto (with no conversion recovery)
      Bool_t              fInvertElectronVeto;    //=true then invert electron veto (for cic selection only atm)      
      Bool_t              fApplyElectronVetoConvRecovery; //=true then apply electron veto with conversion recovery
      Bool_t              fApplyConversionId;    //=true then apply conversion id cuts
      Bool_t              fApplyTriggerMatching;   //match to hlt photon (default=0)      
      Double_t            fPhotonR9Min;          //min R9 value
      EPhIdType           fPhIdType;             //!identification scheme
      EPhIsoType          fPhIsoType;            //!isolation scheme
      Bool_t              fFiduciality;          //=true then apply fiducual requirement
      Double_t            fEtaWidthEB;  	 //max Eta Width in ECAL Barrel
      Double_t            fEtaWidthEE;  	 //max Eta Width in ECAL End Cap
      Double_t            fAbsEtaMax;  	         //max Abs Eta
      Bool_t              fApplyR9Min;           //apply R9 min
      Double_t            fEffAreaEcalEE;
      Double_t            fEffAreaHcalEE;
      Double_t            fEffAreaTrackEE;
      Double_t            fEffAreaEcalEB;
      Double_t            fEffAreaHcalEB;
      Double_t            fEffAreaTrackEB;
      const PhotonCol    *fPhotons;              //!photon branch
      const TrackCol     *fTracks;               //!track branch
      const BeamSpotCol  *fBeamspots;            //!beamspot branch    
      const PileupEnergyDensityCol *fPileUpDen;  //!rho branch
      const DecayParticleCol *fConversions;      //!conversion branch
      const ElectronCol  *fElectrons;            //!electron branch
      const ElectronCol  *fGoodElectrons;        //!electron branch
      const VertexCol*    fPV;                   //!
      const MCParticleCol          *fMCParticles;//!
      const PileupInfoCol          *fPileUp;     //!  
      const PFCandidateCol          *fPFCands;     //!  

      Double_t fbdtCutBarrel;
      Double_t fbdtCutEndcap;
      int                         fVariableType;
      TString                     fEndcapWeights;
      TString                     fBarrelWeights;
      MVATools                    fTool;

      Bool_t fDoMCR9Scaling;
      Double_t fMCR9ScaleEB;
      Double_t fMCR9ScaleEE;
      
      Bool_t fDoMCSigIEtaIEtaScaling;
      Bool_t fDoMCWidthScaling;
      
      Bool_t fDoMCErrScaling;
      Double_t fMCErrScaleEB;
      Double_t fMCErrScaleEE;    

      Bool_t              fPhotonsFromBranch;
      Bool_t              fPVFromBranch;
      Bool_t              fGoodElectronsFromBranch;
      Bool_t              fIsData;



    ClassDef(PhotonIDMod, 1) // Photon identification module
  };
}
#endif
