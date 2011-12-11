//--------------------------------------------------------------------------------------------------
// $Id $
//
// MVATools
//
// Helper Class for photon Identification decisions.
//
// Authors: M.Yang 2011/10/12
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_MVATOOLS_H
#define MITPHYSICS_UTILS_MVATOOLS_H

#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/Electron.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/BaseVertex.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/SuperCluster.h"
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"


class TRandom3;
namespace TMVA {//MVA
  class Reader;
}

namespace mithep {
  class MVATools {
    public:
    MVATools();
      
    
    //--------------------------
    //MVA
    //--------------------------
    
    void InitializeMVA(int VariableType, TString EndcapWeights,TString BarrelWeights);
    Bool_t PassMVASelection(const Photon* p,const Vertex* vtx,const TrackCol* trackCol,const VertexCol* vtxCol,Double_t _tRho,Float_t bdtCutBarrel, Float_t bdtCutEndcap, const ElectronCol* els=0, Bool_t applyElectronVeto=kTRUE);
    Int_t PassElectronVetoInt(const Photon* p, const ElectronCol* els);
      Float_t GetMVAbdtValue(const Photon* p,const Vertex* vtx,const TrackCol* trackCol,const VertexCol* vtxCol,Double_t _tRho, const ElectronCol* els=0, Bool_t applyElectronVeto=kTRUE);
      
      TMVA::Reader *fReaderEndcap;
      TMVA::Reader *fReaderBarrel;    
      
      //MVA Variables
      float HoE;
      float covIEtaIEta;
      float tIso1abs;
      float tIso3abs;
      float tIso2abs;
      float R9;
    
      float absIsoEcal;
      float absIsoHcal;
      float RelEMax;
      float RelETop;
      float RelEBottom;
      float RelELeft;
      float RelERight;
      float RelE2x5Max;
      float RelE2x5Top;
      float RelE2x5Bottom;
      float RelE2x5Left;
      float RelE2x5Right;
      float RelE5x5;

      float EtaWidth;
      float PhiWidth;
      float CoviEtaiPhi;
      float CoviPhiiPhi;

      float NVertexes;
      float RelPreshowerEnergy;

      //variable for v2 and v1
      float RelIsoEcal;
      float RelIsoHcal;
      float tIso1;
      float tIso3;
      float tIso2;
      
      //variables used to compute mva variables
      
      Bool_t PassElecVeto;  
      
      double ecalIso3;
      double ecalIso4;
      double hcalIso4; 
      
      unsigned int wVtxInd;
      
      double trackIso1;
      
      // track iso only
      double trackIso3;
      
      // track iso worst vtx
      double trackIso2; 
      
      double combIso1;
      double combIso2;
      
      double RawEnergy;
      
      double dRTrack;
      
      //spectator variables
      double Pt_MVA;
      double ScEta_MVA;
      
      Bool_t isbarrel;
      
      // check which category it is ...
      int _tCat;
      
      //MVA  
      Bool_t PassMVA;
      TMVA::Reader *reader;
      Float_t bdt;
      Int_t PassElecVetoInt;
      
      ClassDef(MVATools, 0) 
	};
}

#endif
