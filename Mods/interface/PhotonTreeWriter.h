//--------------------------------------------------------------------------------------------------
// $Id: PhotonTreeWriter.h,v 1.11 2012/05/02 16:57:20 fabstoec Exp $
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
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitPhysics/Utils/interface/PhotonFix.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"

class TNtuple;
class TRandom3;

namespace mithep 
{
  
  class PhotonTreeWriterPhoton
  {
    public:  
      void SetVars(const Photon *p, const DecayParticle *c, const Electron *ele, const SuperCluster *pfsc, const MCParticle *m, PhotonFix &phfixph, PhotonFix &phfixele, const TrackCol* trackCol,const VertexCol* vtxCol,Double_t _tRho, const ElectronCol* els=0, Bool_t applyElectronVeto=kTRUE);
      Float_t Ecor()    const { return ecor;    };
      Float_t Ecorerr() const { return ecorerr; };
      Float_t Ecorele()    const { return ecorele;    };
      Float_t Ecoreleerr() const { return ecoreleerr; };      
      
    private:  
      UChar_t hasphoton;
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
      UInt_t  scnclusters;
      UInt_t  scnhits;
      Float_t scetawidth;
      Float_t scphiwidth;
      Float_t hovere;
      Float_t sigietaieta;
      Bool_t  isbarrel;
      Bool_t  isr9reco;
      Bool_t  isr9cat;
      Char_t  phcat;
      Float_t eerr;
      Float_t eerrsmeared;
      Float_t esmearing;
      Float_t idmva;
      Float_t ecalisodr03;
      Float_t hcalisodr03;
      Float_t trkisohollowdr03;
      Float_t ecalisodr04;
      Float_t hcalisodr04;
      Float_t trkisohollowdr04;      
      Float_t trackiso1;
      Float_t trackiso2;
      Float_t combiso1;
      Float_t combiso2;

      //quantities from seed basic cluster      
      Float_t eseed;
      Float_t etaseed;
      Float_t phiseed;
      Int_t   ietaseed;
      Int_t   iphiseed;
      Int_t   ixseed;
      Int_t   iyseed;
      Float_t etacryseed;
      Float_t phicryseed;
      Float_t xcryseed;
      Float_t ycryseed;
      Float_t thetaaxisseed;
      Float_t phiaxisseed;
      Float_t sigietaietaseed;
      Float_t sigiphiphiseed;
      Float_t covietaiphiseed;
      Float_t e3x3seed;
      Float_t e5x5seed;      
      Float_t emaxseed;
      Float_t e2ndseed;
      Float_t etopseed;
      Float_t ebottomseed;
      Float_t eleftseed;
      Float_t erightseed;      
      Float_t e1x3seed;
      Float_t e3x1seed;
      Float_t e1x5seed;
      Float_t e2x2seed;
      Float_t e4x4seed;
      Float_t e2x5maxseed;
      Float_t e2x5topseed;
      Float_t e2x5bottomseed;
      Float_t e2x5leftseed;
      Float_t e2x5rightseed;      
      Float_t xseedseed;
      Float_t yseedseed;
      Float_t zseedseed;
      UInt_t  nhitsseed;
      
      //quantities from second basic cluster, if present
      Float_t ebc2;
      Float_t etabc2;
      Float_t phibc2;
      Int_t   ietabc2;
      Int_t   iphibc2;
      Int_t   ixbc2;
      Int_t   iybc2;
      Float_t etacrybc2;
      Float_t phicrybc2;
      Float_t xcrybc2;
      Float_t ycrybc2;
      Float_t thetaaxisbc2;
      Float_t phiaxisbc2;
      Float_t sigietaietabc2;
      Float_t sigiphiphibc2;
      Float_t covietaiphibc2;
      Float_t e3x3bc2;
      Float_t e5x5bc2;            
      Float_t emaxbc2;
      Float_t e2ndbc2;
      Float_t etopbc2;
      Float_t ebottombc2;
      Float_t eleftbc2;
      Float_t erightbc2;
      Float_t e1x3bc2;
      Float_t e3x1bc2;
      Float_t e1x5bc2;
      Float_t e2x2bc2;
      Float_t e4x4bc2;
      Float_t e2x5maxbc2;
      Float_t e2x5topbc2;
      Float_t e2x5bottombc2;
      Float_t e2x5leftbc2;
      Float_t e2x5rightbc2;      
      Float_t xbc2bc2;
      Float_t ybc2bc2;
      Float_t zbc2bc2;      
      UInt_t  nhitsbc2;
      
      //quantities from lowest energy basic cluster if present
      Float_t ebclast;
      Float_t etabclast;
      Float_t phibclast;
      Int_t   ietabclast;
      Int_t   iphibclast;
      Int_t   ixbclast;
      Int_t   iybclast;
      Float_t etacrybclast;
      Float_t phicrybclast;
      Float_t xcrybclast;
      Float_t ycrybclast;
      Float_t thetaaxisbclast;
      Float_t phiaxisbclast;
      Float_t sigietaietabclast;
      Float_t sigiphiphibclast;
      Float_t covietaiphibclast;
      Float_t e3x3bclast;
      Float_t e5x5bclast;             
      UInt_t  nhitsbclast;
      
      //quantities from second lowest energy basic cluster if present
      Float_t ebclast2;
      Float_t etabclast2;
      Float_t phibclast2;
      Int_t   ietabclast2;
      Int_t   iphibclast2;
      Int_t   ixbclast2;
      Int_t   iybclast2;
      Float_t etacrybclast2;
      Float_t phicrybclast2;
      Float_t xcrybclast2;
      Float_t ycrybclast2;
      Float_t thetaaxisbclast2;
      Float_t phiaxisbclast2;
      Float_t sigietaietabclast2;
      Float_t sigiphiphibclast2;
      Float_t covietaiphibclast2;
      Float_t e3x3bclast2;
      Float_t e5x5bclast2;
      UInt_t  nhitsbclast2;             
      
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
      UChar_t hasconversion;
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
      Char_t  convleadtrackalgo;
      Char_t  convleadtrackalgos;
      Char_t  convleadtrackcharge;
      Float_t convtrailtrackpt;
      Char_t  convtrailtrackalgo;
      Char_t  convtrailtrackalgos;      
      Char_t  trailtrackcharge;
      
      //electron quantities
      UChar_t haselectron;
      UChar_t eleisecaldriven;      
      UChar_t eleistrackerdriven;     
      Float_t elee;
      Float_t elept;
      Float_t eleeta;
      Float_t elephi;
      Char_t  elecharge;
      Float_t elefbrem;
      Float_t eledeta;
      Float_t eledphi;
      Float_t elep;
      Float_t elepin;
      Float_t elepout;
      
      //pf supercluster quantities
      UChar_t haspfsc;
      Float_t pfsce;
      Float_t pfscrawe;
      Float_t pfsceta;
      Float_t pfscphi;
      
      //generator level quantities
      UChar_t ispromptgen;
      Float_t gene;
      Float_t genpt;
      Float_t geneta;
      Float_t genphi;
      Float_t genz;
      Int_t   pdgid;
      Int_t   motherpdgid;
  };
  
  class PhotonTreeWriterDiphotonEvent
  {
    public:
      // ------------ BTAG STUFF -------------------
      Float_t btag;
      Float_t btagJetPt;
      Float_t btagJetEta;
      // ----------- LEPTON TAG STUFF -------------
      Int_t leptonTag;
      // ---------- MUON STUFF --------------------
      Float_t muonPt;
      Float_t muonEta;
      Float_t muDR1;
      Float_t muDR2;
      Float_t muIso1;
      Float_t muIso2;
      Float_t muIso3;
      Float_t muIso4;
      Float_t muD0;
      Float_t muDZ;
      Int_t muNhits;
      Float_t muChi2;
      Int_t muNpixhits;
      Int_t muNegs;
      Int_t muNMatch;      
      // ----------- ELECTRON STUFF --------------
      Float_t elePt;
      Float_t eleEta;
      Float_t eleSCEta;     
      Float_t eleIso1;
      Float_t eleIso2;
      Float_t eleIso3;
      Float_t eleIso4;
      Float_t eleDist;
      Float_t eleDcot;
      Float_t eleCoviee;
      Float_t eleDphiin;
      Float_t eleDetain;
      Float_t eleDR1;
      Float_t eleDR2;
      Float_t eleMass1;
      Float_t eleMass2;
      Int_t eleNinnerHits;     
      // -----------------------------------------
      Float_t rho;
      Float_t genHiggspt;
      Float_t genHiggsZ;
      Float_t genmass;
      Float_t gencostheta;
      Float_t bsX;
      Float_t bsY;
      Float_t bsZ;
      Float_t bsSigmaZ; 
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
      Float_t deltamvtx;
      Float_t ptgg;
      Float_t etagg;
      Float_t phigg;
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
      UChar_t  ismc;
      
      //corrected quantities from PhotonFix corrections
      Float_t masscor;
      Float_t masscorerr;
      Float_t masscorele;
      Float_t masscoreleerr;
      
      //jet quantities
      Float_t jet1pt;
      Float_t jet1eta;
      Float_t jet1phi;
      Float_t jet1mass;
      Float_t jet2pt;
      Float_t jet2eta;
      Float_t jet2phi;
      Float_t jet2mass;
      Float_t jetcentralpt;
      Float_t jetcentraleta;
      Float_t jetcentralphi;
      Float_t jetcentralmass;
      Float_t dijetpt;
      Float_t dijeteta;
      Float_t dijetphi;
      Float_t dijetmass;
      Float_t jetetaplus;
      Float_t jetetaminus;
      
      Float_t zeppenfeld;
      Float_t dphidijetgg;
      
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
    void                SetPFJetName(const char *n)       { fPFJetName = n;              }
    void                SetPFJetsFromBranch(Bool_t b)     { fPFJetsFromBranch = b;       }
    void                SetEnableJets(Bool_t b)           { fEnableJets = b;             }
    void                SetApplyLeptonTag(Bool_t b)       { fApplyLeptonTag = b;         }
    void                SetApplyBTag(Bool_t b)            { fApplyBTag = b;              }
    void                SetPhFixDataFile(const char *n)   { fPhFixDataFile = n;          }


    // set basic Cut variables (FOR PRE-SELECTION)

    // is Data Or Not?
    void                SetIsData (Bool_t b)                 { fIsData = b; };
    

    void                SetApplyElectronVeto(Bool_t b)   { fApplyElectronVeto = b;     }          

    void                SetTupleName(const char* c)          { fTupleName = c; }
    void                SetGoodElectronsFromBranch(Bool_t b) { fGoodElectronsFromBranch = b; }
    void                SetGoodElectronName(TString name)    { fGoodElectronName = name; }
    void                SetWriteDiphotonTree(Bool_t b)       { fWriteDiphotonTree = b; }
    void                SetWriteSingleTree(Bool_t b)         { fWriteSingleTree = b; }
    void                SetLoopOnGoodElectrons(Bool_t b)     { fLoopOnGoodElectrons = b; }
    void                SetExcludeSinglePrompt(Bool_t b)     { fExcludeSinglePrompt = b; }
    void                SetExcludeDoublePrompt(Bool_t b)     { fExcludeDoublePrompt = b; }

    void                SetLeptonTagElectronsName(TString name) { fLeptonTagElectronsName = name; }
    void                SetLeptonTagMuonsName    (TString name) { fLeptonTagMuonsName     = name; }

  protected:
    void                Process();
    void                SlaveBegin();
    // Private auxiliary methods...
    void                FindHiggsPtAndZ(Float_t& pt, Float_t& z, Float_t& mass);
    Float_t             GetEventCat    (PhotonTools::CiCBaseLineCats cat1,
					PhotonTools::CiCBaseLineCats cat2);

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
    TString             fPFJetName;

    TString             fLeptonTagElectronsName;
    TString             fLeptonTagMuonsName;

    
    // is it Data or MC?
    Bool_t              fIsData;
    
    // in case there's some PV pre-selection
    Bool_t              fPhotonsFromBranch;
    Bool_t              fPVFromBranch;
    Bool_t              fGoodElectronsFromBranch;
    Bool_t              fPFJetsFromBranch;

    const PhotonCol               *fPhotons;
    const ElectronCol             *fElectrons;
    const ElectronCol             *fGoodElectrons;    
    const DecayParticleCol        *fConversions;
    const TrackCol                *fTracks;
    const PileupEnergyDensityCol  *fPileUpDen;
    const VertexCol               *fPV;
    const BeamSpotCol             *fBeamspot;
    const PFCandidateCol          *fPFCands;
    const MCParticleCol           *fMCParticles;
    const PileupInfoCol           *fPileUp;    
    const SuperClusterCol         *fSuperClusters;   
    const PFMetCol                *fPFMet;
    const JetCol                  *fPFJets;
    
    const ElectronCol             *fLeptonTagElectrons;
    const MuonCol                 *fLeptonTagMuons;

    // --------------------------------
    Bool_t                         fLoopOnGoodElectrons; //loop over good elecs instead of photons
    Bool_t                         fApplyElectronVeto;   //invert elec veto (for cic sel. only atm)
    Bool_t                         fWriteDiphotonTree;
    Bool_t                         fWriteSingleTree;

    Bool_t                         fExcludeSinglePrompt;
    Bool_t                         fExcludeDoublePrompt;
    
    Bool_t                         fEnableJets;

    Bool_t                         fApplyLeptonTag;
    Bool_t                         fApplyBTag;

    TString                        fPhFixDataFile;
    PhotonFix                      fPhfixph;
    PhotonFix                      fPhfixele;

    // --------------------------------
    // validation Tuple
    TString                        fTupleName;
    PhotonTreeWriterDiphotonEvent* fDiphotonEvent;
    PhotonTreeWriterPhoton*        fSinglePhoton;    
    TTree*                         hCiCTuple;
    TTree*                         hCiCTupleSingle;

    ClassDef(PhotonTreeWriter, 1) // Photon identification module
  };
}
#endif
