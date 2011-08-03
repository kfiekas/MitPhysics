//--------------------------------------------------------------------------------------------------
// $Id: PhotonTreeWriter.h,v 1.6 2011/07/27 17:50:52 fabstoec Exp $
//
// PhotonTreeWriter
//
// Authors: J. Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PHOTONTREEWRITER_H
#define MITPHYSICS_MODS_PHOTONTREEWRITER_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PhotonFwd.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"

#include "MitPhysics/Utils/interface/PhotonTools.h"

class TNtuple;
class TRandom3;

namespace mithep 
{
  
  class PhotonTreeWriterPhoton
  {
    public:  
      void SetVars(const Photon *p, const DecayParticle *c = 0, const MCParticle *m = 0);
      Float_t Ecor()    const { return ecor;    };
      Float_t Ecorerr() const { return ecorerr; };
      
    private:  
      Float_t e;
      Float_t pt;
      Float_t eta;
      Float_t phi;
      Float_t r9;
      Float_t e3x3;
      Float_t e5x5;
      Float_t sce;
      Float_t scrawe;
      Float_t scpse;
      Float_t sceta;
      Float_t scphi;
      UInt_t scnclusters;
      UInt_t scnhits;
      Float_t hovere;
      Float_t sigietaieta;
      Bool_t isbarrel;
      Bool_t isr9reco;
      Bool_t isr9cat;
      UChar_t phcat;
      
      //quantities from seed basic cluster
      Float_t sigiphiphi;
      Float_t covietaiphi;
      Float_t emax;
      Float_t e2nd;
      Float_t etop;
      Float_t ebottom;
      Float_t eleft;
      Float_t eright;
      Float_t e1x3;
      Float_t e3x1;
      Float_t e1x5;
      Float_t e2x2;
      Float_t e4x4;
      Float_t e2x5max;
      Float_t e2x5top;
      Float_t e2x5bottom;
      Float_t e2x5left;
      Float_t e2x5right;      
      
      //energy correction quantities from PhotonFix
      Float_t ecor;
      Float_t ecorerr;
      Float_t etac;
      Float_t etas;
      Float_t etam;
      Float_t phic;
      Float_t phis;
      Float_t phim;  
      Float_t xz;
      Float_t xc;
      Float_t xs;
      Float_t xm;
      Float_t yz;
      Float_t yc;
      Float_t ys;
      Float_t ym;        
      
      //conversion quantities
      Bool_t hasconversion;
      Float_t convp;
      Float_t convpt;
      Float_t conveta;
      Float_t convphi;
      Float_t convdeta;
      Float_t convdphi;
      Float_t convvtxrho;
      Float_t convvtxz;
      Float_t convvtxphi;
      Float_t convleadpt;
      Float_t convtrailpt;
      Float_t convleadtrackpt;
      Char_t convleadtrackalgo;
      Char_t convleadtrackalgos;
      Char_t convleadtrackcharge;
      Float_t convtrailtrackpt;
      Char_t convtrailtrackalgo;
      Char_t convtrailtrackalgos;      
      Char_t trailtrackcharge;
      
      
      //generator level quantities
      Bool_t ispromptgen;
      Float_t gene;
      Float_t genpt;
      Float_t geneta;
      Float_t genphi;
      Float_t genz;
      

      
  };
  
  class PhotonTreeWriterDiphotonEvent
  {
    public:
      Float_t rho;
      Float_t genHiggspt;
      Float_t genHiggsZ;
      Float_t genmass;
      Float_t gencostheta;
      Float_t vtxZ;
      Int_t   nVtx;
      Int_t   numPU;
      Int_t   numPUminus;
      Int_t   numPUplus;
      Float_t mass;
      Float_t ptgg;
      Float_t costheta;
      UInt_t  evt;
      UInt_t  run;
      UInt_t  lumi;
      UChar_t evtcat;
      
      //corrected quantities from PhotonFix corrections
      Float_t masscor;
      Float_t masscorerr;
      
      PhotonTreeWriterPhoton photons[2];

    
  };
  
  class PhotonTreeWriter : public BaseMod
  {
  public:
    PhotonTreeWriter(const char *name ="PhotonTreeWriter", 
		       const char *title="Selecting PhotonPairs");
    
    ~PhotonTreeWriter();

    // setting all the input Names
    void                SetInputPhotonsName(const char *n){ fPhotonBranchName= n;        }
    void                SetPhotonsFromBranch(bool b)      { fPhotonsFromBranch = b;      }
    void                SetTrackName(const char *n)       { fTrackBranchName = n;        }
    void                SetElectronName(const char *n)    { fElectronName = n;           }
    void                SetConversionName(const char *n)  { fConversionName = n;         }
    void                SetPUDensityName(const char *n)   { fPileUpDenName = n;          }
    void                SetPVName(const char *n)          { fPVName = n;                 }
    void                SetPVFromBranch(bool b)           { fPVFromBranch = b;           }
    void                SetMCParticle(const char *n)      { fMCParticleName = n;         }
    void                SetPUInfoName(const char *n)      { fPileUpName = n;             }
    void                SetBeamspotName(const char *n)    { fBeamspotName = n;           }
    void                SetPFCandName(const char *n)      { fPFCandName = n;             }



    // set basic Cut variables (FOR PRE-SELECTION)

    // is Data Or Not?
    void                SetIsData (Bool_t b) { fIsData = b;};
    

    void                SetInvertElectronVeto(Bool_t b)   { fInvertElectronVeto = b;     }          

    void                SetTupleName(const char* c)     { fTupleName     = c; }
    void                SetGoodElectronsFromBranch(Bool_t b) { fGoodElectronsFromBranch = b; }
    void                SetGoodElectronName(TString name) { fGoodElectronName = name; }
    void                SetWriteDiphotonTree(Bool_t b)  { fWriteDiphotonTree = b; }
    void                SetWriteSingleTree(Bool_t b)    { fWriteSingleTree = b; }

  protected:
    void                Process();
    void                SlaveBegin();

    // private auxiliary methods...
    void                FindHiggsPtAndZ(Float_t& pt, Float_t& z, Float_t& mass);
    Float_t             GetEventCat(PhotonTools::CiCBaseLineCats cat1, PhotonTools::CiCBaseLineCats cat2);

    // Names for the input Collections
    TString             fPhotonBranchName;
    TString             fElectronName;
    TString             fGoodElectronName;
    TString             fConversionName;
    TString             fTrackBranchName;
    TString             fPileUpDenName;    
    TString             fPVName;
    TString             fBeamspotName;
    TString             fPFCandName;
    TString             fMCParticleName;
    TString             fPileUpName;
    
    // is it Data or MC?
    Bool_t              fIsData;
    
    // in case there's some PV pre-selection
    Bool_t              fPhotonsFromBranch;
    Bool_t              fPVFromBranch;
    Bool_t              fGoodElectronsFromBranch;

    const PhotonCol              *fPhotons;
    const ElectronCol            *fElectrons;
    const ElectronCol            *fGoodElectrons;    
    const DecayParticleCol       *fConversions;
    const TrackCol               *fTracks;
    const PileupEnergyDensityCol *fPileUpDen;
    const VertexCol              *fPV;
    const BeamSpotCol            *fBeamspot;
    const PFCandidateCol         *fPFCands;
    const MCParticleCol          *fMCParticles;
    const PileupInfoCol          *fPileUp;    
       
    // --------------------------------
    Bool_t              fInvertElectronVeto;    //=true then invert electron veto (for cic selection only atm)     
    Bool_t fWriteDiphotonTree;
    Bool_t fWriteSingleTree;


    // --------------------------------
    // validation Tuple
    TString fTupleName;
    PhotonTreeWriterDiphotonEvent* fDiphotonEvent;
    PhotonTreeWriterPhoton* fSinglePhoton;    
    TTree* hCiCTuple;
    TTree* hCiCTupleSingle;

    ClassDef(PhotonTreeWriter, 1) // Photon identification module
  };
}
#endif
