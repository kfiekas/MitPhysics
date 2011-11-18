//--------------------------------------------------------------------------------------------------
// $Id: PhotonTreeWriter.h,v 1.3 2011/10/07 09:56:36 bendavid Exp $
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
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitPhysics/Utils/interface/PhotonFix.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"

class TNtuple;
class TRandom3;

namespace mithep 
{
  
  class PhotonTreeWriterPhoton
  {
    public:  
      void SetVars(const Photon *p, const DecayParticle *c, const Electron *ele, const SuperCluster *pfsc, const MCParticle *m, PhotonFix &phfixph, PhotonFix &phfixele);
      Float_t Ecor()    const { return ecor;    };
      Float_t Ecorerr() const { return ecorerr; };
      Float_t Ecorele()    const { return ecorele;    };
      Float_t Ecoreleerr() const { return ecoreleerr; };      
      
    private:  
      Bool_t hasphoton;
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
      Float_t scetawidth;
      Float_t scphiwidth;
      Float_t hovere;
      Float_t sigietaieta;
      Bool_t isbarrel;
      Bool_t isr9reco;
      Bool_t isr9cat;
      Char_t phcat;
      Float_t eerr;
      Float_t eerrsmeared;
      Float_t esmearing;
      
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
      Float_t xseed;
      Float_t yseed;
      Float_t zseed;
      
      Float_t eseed;
      Float_t etaseed;
      Float_t phiseed;
      Int_t ietaseed;
      Int_t iphiseed;
      Int_t ixseed;
      Int_t iyseed;
      Float_t etacryseed;
      Float_t phicryseed;
      Float_t xcryseed;
      Float_t ycryseed;
      Float_t thetaaxisseed;
      Float_t phiaxisseed;
      
      //quantities from second basic cluster, if present
      Float_t ebc2;
      Float_t etabc2;
      Float_t phibc2;
      Int_t ietabc2;
      Int_t iphibc2;
      Int_t ixbc2;
      Int_t iybc2;
      Float_t etacrybc2;
      Float_t phicrybc2;
      Float_t xcrybc2;
      Float_t ycrybc2;
      Float_t thetaaxisbc2;
      Float_t phiaxisbc2;      
      
      //energy correction quantities from PhotonFix
      Float_t ecor;
      Float_t ecorerr;
      Float_t ecorele;
      Float_t ecoreleerr;
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
      
      //electron quantities
      Bool_t haselectron;
      Bool_t eleisecaldriven;      
      Bool_t eleistrackerdriven;     
      Float_t elee;
      Float_t elept;
      Float_t eleeta;
      Float_t elephi;
      Char_t elecharge;
      Float_t elefbrem;
      Float_t eledeta;
      Float_t eledphi;
      Float_t elep;
      Float_t elepin;
      Float_t elepout;
      
      //pf supercluster quantities
      Bool_t haspfsc;
      Float_t pfsce;
      Float_t pfscrawe;
      Float_t pfsceta;
      Float_t pfscphi;
      
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
      Float_t bsX;
      Float_t bsY;
      Float_t bsZ;
      Float_t vtxX;
      Float_t vtxY;      
      Float_t vtxZ;
      Int_t   nVtx;
      Int_t   numPU;
      Int_t   numPUminus;
      Int_t   numPUplus;
      Float_t mass;
      Float_t masserr;
      Float_t masserrsmeared;
      Float_t masserrwrongvtx;
      Float_t masserrsmearedwrongvtx;
      Float_t vtxprob;
      Float_t ptgg;
      Float_t costheta;
      Float_t massele;
      Float_t ptee;
      Float_t costhetaele;
      Float_t mt;
      Float_t cosphimet;
      Float_t mtele;
      Float_t cosphimetele;
      UInt_t  evt;
      UInt_t  run;
      UInt_t  lumi;
      UChar_t evtcat;
      UInt_t  nobj;
      Float_t pfmet;
      Float_t pfmetphi;
      Float_t pfmetx;
      Float_t pfmety;
      Bool_t  ismc;
      
      //corrected quantities from PhotonFix corrections
      Float_t masscor;
      Float_t masscorerr;
      Float_t masscorele;
      Float_t masscoreleerr;
      
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
    void                SetSuperClusterName(const char *n) { fSuperClusterName = n;      }
    void                SetPhFixDataFile(const char *n)   { fPhFixDataFile = n;          }


    // set basic Cut variables (FOR PRE-SELECTION)

    // is Data Or Not?
    void                SetIsData (Bool_t b) { fIsData = b;};
    

    void                SetInvertElectronVeto(Bool_t b)   { fInvertElectronVeto = b;     }          

    void                SetTupleName(const char* c)     { fTupleName     = c; }
    void                SetGoodElectronsFromBranch(Bool_t b) { fGoodElectronsFromBranch = b; }
    void                SetGoodElectronName(TString name) { fGoodElectronName = name; }
    void                SetWriteDiphotonTree(Bool_t b)  { fWriteDiphotonTree = b; }
    void                SetWriteSingleTree(Bool_t b)    { fWriteSingleTree = b; }
    void                SetLoopOnGoodElectrons(Bool_t b) { fLoopOnGoodElectrons = b; }
    void                SetExcludeSinglePrompt(Bool_t b) { fExcludeSinglePrompt = b; }
    void                SetExcludeDoublePrompt(Bool_t b) { fExcludeDoublePrompt = b; }

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
    TString             fSuperClusterName;
    TString             fPFMetName;
    
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
    const SuperClusterCol        *fSuperClusters;   
    const PFMetCol               *fPFMet;
    
    // --------------------------------
    Bool_t fLoopOnGoodElectrons; //primary loop over good electrons collection instead of photons
    Bool_t              fInvertElectronVeto;    //=true then invert electron veto (for cic selection only atm)     
    Bool_t fWriteDiphotonTree;
    Bool_t fWriteSingleTree;

    Bool_t fExcludeSinglePrompt;
    Bool_t fExcludeDoublePrompt;
    
    TString fPhFixDataFile;
    PhotonFix fPhfixph;
    PhotonFix fPhfixele;

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
