// ======================================================================================================
// $Id: ZGTools.h,v 1.0 2012/08/20 16:11:01 ksingh
// 
// ZGTools (Helper functions for Higgs to Z + gamma analysis)
//
// Authors: K. Hahn, K. Singh
//
// angles following conventions in Gainer, et. al.  
// arXiv:1112.1405v2 [hep-ph], 6 May 2012
// see also : Gainer, et. al. 
// arVix:1108.2274v1 [hep-ph], 10 Aug 2011
//
// ======================================================================================================
#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
#include "MitAna/DataTree/interface/PileupEnergyDensity.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitPhysics/Utils/interface/ZGTools.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitPhysics/Utils/interface/MVATools.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/StableParticle.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/Muon.h"
#include "TDataMember.h"
#include "TFile.h"
#include <TNtuple.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <math.h>
#include <string>
#include <iostream>
#include <vector>
#include <assert.h>



using namespace std;

ClassImp(mithep::ZGTools)

using namespace mithep;

ZGTools::ZGTools()
{
  // Constructor.
}

// --------------------------------------------------------------------------------------------------
// set angles from a struct of lab-frame 4vecs 
//
ZGAngles ZGTools::getZGAngles( ZGLabVectors &l, bool debug ) 
// --------------------------------------------------------------------------------------------------
{
    ZGAngles a;
    if( debug ) { 
    cout << endl << endl;
    cout << "veczg(x:y:z) :: " << l.veczg.X() <<":"<< l.veczg.Y() <<":"<< l.veczg.Z() << endl;
    cout << "\tpt: " << l.veczg.Pt() << "\teta: " << l.veczg.Eta() << "\tphi:" << l.veczg.Phi() << endl;
    cout << "vecz(x:y:z) :: "  << l.vecz.X() <<":"<< l.vecz.Y() <<":"<< l.vecz.Z() << endl;
    cout << "\tpt: " << l.vecz.Pt() << "\teta: " << l.vecz.Eta() << "\tphi:" << l.vecz.Phi() << endl;
    cout << "veclp(x:y:z) :: "  << l.veclp.X() <<":"<< l.veclp.Y() <<":"<< l.veclp.Z() << endl;
    cout << "\tpt: " << l.veclp.Pt() << "\teta: " << l.veclp.Eta() << "\tphi:" << l.veclp.Phi() << endl;
    cout << "veclm(x:y:z) :: "  << l.veclm.X() <<":"<< l.veclm.Y() <<":"<< l.veclm.Z() << endl;
    cout << "\tpt: " << l.veclm.Pt() << "\teta: " << l.veclm.Eta() << "\tphi:" << l.veclm.Phi() << endl;
    cout << "vecg(x:y:z) :: "  << l.vecg.X() <<":"<< l.vecg.Y() <<":"<< l.vecg.Z() << endl;
    cout << "\tpt: " << l.vecg.Pt() << "\teta: " << l.vecg.Eta() << "\tphi:" << l.vecg.Phi() << endl;
  }



  //
  // define the frames & coordinate systems
  // --------------------------------------------------------------------------------------------------

  TVector3 Xframe  = l.veczg.BoostVector();
  TVector3 Z1frame = l.vecz.BoostVector();

  if( debug ) { 
    cout << "Xboost(x:y:z) :: " << Xframe.X() <<":"<< Xframe.Y() <<":"<< Xframe.Z() << endl;
    cout << "Zboost(x:y:z) :: " << Z1frame.X() <<":"<< Z1frame.Y() <<":"<< Z1frame.Z() << endl;
  }

  // "partons" 
  TLorentzVector  kq( 0, 0, (l.veczg.E()+l.veczg.Pz())/2, (l.veczg.E()+l.veczg.Pz())/2 );
  TLorentzVector  kqbar( 0, 0, (l.veczg.Pz()-l.veczg.E())/2, (l.veczg.E()-l.veczg.Pz())/2 );
  TLorentzVector veckq_in_Xframe(kq);
  TLorentzVector veckqbar_in_Xframe(kqbar);
  veckq_in_Xframe.Boost(-1*Xframe);
  veckqbar_in_Xframe.Boost(-1*Xframe);
  if( debug ) { 
    cout << "kq_X(x:y:z) :: " << veckq_in_Xframe.X() <<":"<< veckq_in_Xframe.Y() <<":"<< veckq_in_Xframe.Z() << endl;
    cout << "kqb_X(x:y:z) :: " << veckqbar_in_Xframe.X() <<":"<< veckqbar_in_Xframe.Y() <<":"<< veckqbar_in_Xframe.Z() << endl;
  }

  // Z vectors
  if( debug ) { cout << "boosting Z,g ..." << endl;}
  TLorentzVector vecz_in_Xframe (l.vecz); 
  TLorentzVector vecg_in_Xframe (l.vecg); 
  TLorentzVector vecz_in_Z1frame (l.vecz); 
  vecz_in_Xframe.Boost(-1*Xframe);
  vecg_in_Xframe.Boost(-1*Xframe);
  vecz_in_Z1frame.Boost(-1*Z1frame);

  if( debug ) { cout << "rotating ..." << endl;}
  // coord system in the CM frame
  TVector3 uz_in_Xframe = vecz_in_Xframe.Vect().Unit();
  TVector3 uy_in_Xframe = (veckq_in_Xframe.Vect().Unit().Cross(uz_in_Xframe.Unit())).Unit();
  TVector3 ux_in_Xframe = (uy_in_Xframe.Unit().Cross(uz_in_Xframe.Unit())).Unit();
  TRotation rotation;
  rotation = rotation.RotateAxes( ux_in_Xframe,uy_in_Xframe,uz_in_Xframe ).Inverse();

  //
  // for going to the Z frames from the CM frame,
  // boost after transform
  //
  if( debug ) { cout << "xform z ..." << endl; }
  TLorentzVector vecz_in_Xframe_newcoords(vecz_in_Xframe);
  vecz_in_Xframe_newcoords.Transform(rotation);
  TVector3 Z1frame_from_Xframe_newcoords = vecz_in_Xframe_newcoords.BoostVector();
  if( debug ) { 
    cout << "vecz_in_Xframe_newcoords (x:y:z) :: " 
	 << vecz_in_Xframe_newcoords.X() <<":"
	 << vecz_in_Xframe_newcoords.Y() <<":"
	 << vecz_in_Xframe_newcoords.Z() << endl;
    
    cout << "Z1frame_from_Xframe_newcoords (x:y:z) :: " 
	 << Z1frame_from_Xframe_newcoords.X() <<":"
	 << Z1frame_from_Xframe_newcoords.Y() <<":"
	 << Z1frame_from_Xframe_newcoords.Z() << endl;
  }


  //  
  // theta(lm), phi(lm) in Z1 frame
  // --------------------------------------------------------------------------------------------------
  TLorentzVector veclm_in_Z1frame(l.veclm);
  TLorentzVector veclp_in_Z1frame(l.veclp);

  // first boost to CM, then redefine coords  
  if( debug ) { cout << "xform lm..." << endl; }
  veclm_in_Z1frame.Boost(-1*Xframe);
  veclm_in_Z1frame.Transform(rotation);
  veclp_in_Z1frame.Boost(-1*Xframe);
  veclp_in_Z1frame.Transform(rotation);
  
  // then boost to Z1
  if ( debug ) { cout << "boost lm ..." << endl;}
  veclm_in_Z1frame.Boost(-1*Z1frame_from_Xframe_newcoords);
  veclp_in_Z1frame.Boost(-1*Z1frame_from_Xframe_newcoords);

  // now get angles
  a.phi      = veclm_in_Z1frame.Phi();
  a.costheta_lm = veclm_in_Z1frame.CosTheta();
  a.costheta_lp = veclp_in_Z1frame.CosTheta();
  
  //
  // or, just go directly to Z frame, skip CM frame
  //
  /*
  // direct
  TLorentzVector veclm_in_Z1frame_direct(l.veclm);
  veclm_in_Z1frame_direct.Boost(-1*Z1frame);
  TLorentzVector veclp_in_Z1frame_direct(l.veclp);
  veclp_in_Z1frame_direct.Boost(-1*Z1frame);
  a.costheta_lm = veclm_in_Z1frame_direct.Vect().Unit().Dot((-1*l.vecz.Vect()).Unit());
  a.costheta_lp = veclp_in_Z1frame_direct.Vect().Unit().Dot((-1*l.vecz.Vect()).Unit());
  */
    
  if( debug ) { 
    TVector3 unitz(0,0,1);
    cout << "compare :: ct: " << a.costheta_lm
	 << "\tv.ct: " << veclm_in_Z1frame.Vect().CosTheta()
	 << "\tc(t): " << TMath::Cos(veclm_in_Z1frame.Theta())
	 << "\tc(v.t): " << TMath::Cos(veclm_in_Z1frame.Theta())
	 << "\tdot: " << veclm_in_Z1frame.Vect().Unit().Dot(unitz)
	 << endl;
    cout << "costheta lp vs lm : " << veclp_in_Z1frame.CosTheta() << " vs " << veclm_in_Z1frame.CosTheta() << endl;
  }


  //  
  // Theta(Z1,Zg-lab) in X frame 
  // --------------------------------------------------------------------------------------------------
  TLorentzVector veczg_in_Xframe(l.veczg);
  //  veczg_in_Xframe.Boost(-1*Xframe);
  veczg_in_Xframe.Transform(rotation);

  if( debug ) { 
    cout << "veczg_X (x:y:z) ::  " << veczg_in_Xframe.X() 
	 <<":"<< veczg_in_Xframe.Y()  
	 <<":"<< veczg_in_Xframe.Z()  
	 << endl;
    cout << "vecz_X (x:y:z) ::  " << vecz_in_Xframe.X() 
	 <<":"<< vecz_in_Xframe.Y()  
	 <<":"<< vecz_in_Xframe.Z()  
	 << endl;
    cout << "vecg_X (x:y:z) ::  " << vecg_in_Xframe.X() 
	 <<":"<< vecg_in_Xframe.Y()  
	 <<":"<< vecg_in_Xframe.Z()  
	 << endl;

    // check
    TLorentzVector ztmp(l.vecz);
    ztmp.Boost(-1*Xframe);
    ztmp.Transform(rotation);
    
    TLorentzVector gtmp(l.vecg);
    gtmp.Boost(-1*Xframe);
    gtmp.Transform(rotation);
    
    cout << "ztmp.CosTheta(): " << ztmp.CosTheta() << endl;
    cout << "gtmp.CosTheta(): " << gtmp.CosTheta() << endl;
  }
  

  TLorentzVector veczg_in_Xframe_newcoords(l.veczg);
  //  veczg_in_Xframe_newcoords.Boost(-1*Xframe);
  veczg_in_Xframe_newcoords.Transform(rotation);
  a.cosTheta = (-1*veczg_in_Xframe_newcoords.Vect()).CosTheta();

  if( debug ) { cout << "phi: " << a.phi << "\tct: " << a.costheta_lm << "\tcT:" << a.cosTheta << endl;}

  a.ptg   = l.vecg.Pt();
  a.etag  = l.vecg.Eta();
  a.ptl1  = l.veclp.Pt();
  a.etal1 = l.veclp.Eta();
  a.ptl2  = l.veclm.Pt();
  a.etal2 = l.veclm.Eta();
  a.ptz   = (l.veclp+l.veclm).Pt();
  a.etaz  = (l.veclp+l.veclm).Eta();
  a.mzg   = l.veczg.M();
  a.mz    = l.vecz.M();

 return a;
}


bool ZGTools::electronCutBasedIDLoose(   const mithep::Electron *ele,
                                         const mithep::Vertex * vtx,
                                         const mithep::DecayParticleCol *conversions,
					 const float year)
//----------------------------------------------------------------------------------------
{

  bool pass=true;

  // KH, replace w/ EG definition
  //  float ooemoop = fabs(1 - ele->ESuperClusterOverP())/(ele->SCluster()->Et()*TMath::CosH(ele->SCluster()->Eta()));
  float ooemoop       = (1.0/ele->EcalEnergy() - ele->ESuperClusterOverP()/ele->EcalEnergy());

  pass =  EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::LOOSE,
				      (fabs(ele->Eta())<1.566),//             (ele->IsEB()||ele->IsEBEtaGap()), 
				      ele->Pt(),
				      ele->SCluster()->Eta(),
				      ele->DeltaEtaSuperClusterTrackAtVtx(),
				      ele->DeltaPhiSuperClusterTrackAtVtx(),
				      ele->CoviEtaiEta(),
				      ele->HadronicOverEm(),
				      ooemoop,
				      ele->BestTrk()->D0Corrected(*vtx),
				      ele->BestTrk()->DzCorrected(*vtx),
				      0.,//const float iso_ch, 
				      0.,//const float iso_em, 
				      0.,//const float iso_nh, 
				      !mithep::ElectronTools::PassConversionFilter(ele, conversions, vtx, 0, 1e-6, 2.0, kTRUE, kFALSE),
				      ele->CorrectedNExpectedHitsInner(),
				      0.
				      );
  /*
  pass =  ZGTools::PassWP(2,
                 (fabs(ele->Eta())<1.566),//             (ele->IsEB()||ele->IsEBEtaGap()), 
                  ele->Pt(),
                 ele->SCluster()->Eta(),
                 ele->DeltaEtaSuperClusterTrackAtVtx(),
                 ele->DeltaPhiSuperClusterTrackAtVtx(),
                 ele->CoviEtaiEta(),
                 ele->HadronicOverEm(),
                 ooemoop,
                 ele->BestTrk()->D0Corrected(*vtx),
                 ele->BestTrk()->DzCorrected(*vtx),
                 0.,//const float iso_ch, 
                 0.,//const float iso_em, 
                 0.,//const float iso_nh, 
                 !mithep::ElectronTools::PassConversionFilter(ele, conversions, vtx, 0, 1e-6, 2.0, kTRUE, kFALSE),
                 ele->CorrectedNExpectedHitsInner(),
                 0.,
		//const double rho,
		 year
                 );
  */
  
  if( (fabs(ele->SCluster()->Eta()) <= 1.566 && fabs(ele->SCluster()->Eta()) >= 1.4442) )
      pass = false;

  return pass;

}


bool ZGTools::PassWP(int workingPoint, const bool isEB, const float pt, const float eta,
    const float dEtaIn, const float dPhiIn, const float sigmaIEtaIEta, const float hoe,
    const float ooemoop, const float d0vtx, const float dzvtx, const float iso_ch, const float iso_em, const float iso_nh,
    const bool vtxFitConversion, const unsigned int mHits, const double rho, const float y)
{
    float cut_dEtaIn[2]         = {999.9, 999.9};
    float cut_dPhiIn[2]         = {999.9, 999.9};
    float cut_sigmaIEtaIEta[2]  = {999.9, 999.9};
    float cut_hoe[2]            = {999.9, 999.9};
    float cut_ooemoop[2]        = {999.9, 999.9};
    float cut_d0vtx[2]          = {999.9, 999.9};
    float cut_dzvtx[2]          = {999.9, 999.9};
    float cut_iso[2]            = {999.9, 999.9};
    bool cut_vtxFit[2]          = {false, false};
    unsigned int cut_mHits[2]   = {999, 999};

    if (workingPoint == 1) {//VETO
        cut_dEtaIn[0]        = 0.007; cut_dEtaIn[1]        = 0.010;
        cut_dPhiIn[0]        = 0.800; cut_dPhiIn[1]        = 0.700;
        cut_sigmaIEtaIEta[0] = 0.010; cut_sigmaIEtaIEta[1] = 0.030;
        cut_hoe[0]           = 0.150; cut_hoe[1]           = 999.9;
        cut_ooemoop[0]       = 999.9; cut_ooemoop[1]       = 999.9;
        cut_d0vtx[0]         = 0.040; cut_d0vtx[1]         = 0.040;
        cut_dzvtx[0]         = 0.200; cut_dzvtx[1]         = 0.200;
        cut_vtxFit[0]        = false; cut_vtxFit[1]        = false;
        cut_mHits[0]         = 999  ; cut_mHits[1]         = 999;
        cut_iso[0]           = 0.150; cut_iso[1]           = 0.150;
    }
    else if (workingPoint == 2) {//LOOSE
        cut_dEtaIn[0]        = 0.007; cut_dEtaIn[1]        = 0.009;
        cut_dPhiIn[0]        = 0.150; cut_dPhiIn[1]        = 0.100;
        cut_sigmaIEtaIEta[0] = 0.010; cut_sigmaIEtaIEta[1] = 0.030;
        cut_hoe[0]           = 0.120; cut_hoe[1]           = 0.100;
        cut_ooemoop[0]       = 0.050; cut_ooemoop[1]       = 0.050;
        cut_d0vtx[0]         = 0.020; cut_d0vtx[1]         = 0.020;
        cut_dzvtx[0]         = 0.200; cut_dzvtx[1]         = 0.200;
        cut_vtxFit[0]        = true ; cut_vtxFit[1]        = true;
        cut_mHits[0]         = 1    ; cut_mHits[1]         = 1;
        if (pt >= 20.0) {
            cut_iso[0] = 0.150; cut_iso[1] = 0.150;
        }
        else {
            cut_iso[0] = 0.150; cut_iso[1] = 0.100;
        }
    }
    else if (workingPoint == 3) {//MEDIUM
        cut_dEtaIn[0]        = 0.004; cut_dEtaIn[1]        = 0.007;
        cut_dPhiIn[0]        = 0.060; cut_dPhiIn[1]        = 0.030;
        cut_sigmaIEtaIEta[0] = 0.010; cut_sigmaIEtaIEta[1] = 0.030;
        cut_hoe[0]           = 0.120; cut_hoe[1]           = 0.100;
        cut_ooemoop[0]       = 0.050; cut_ooemoop[1]       = 0.050;
        cut_d0vtx[0]         = 0.020; cut_d0vtx[1]         = 0.020;
        cut_dzvtx[0]         = 0.100; cut_dzvtx[1]         = 0.100;
        cut_vtxFit[0]        = true ; cut_vtxFit[1]        = true;
        cut_mHits[0]         = 1    ; cut_mHits[1]         = 1;
        if (pt >= 20.0) {
            cut_iso[0] = 0.150; cut_iso[1] = 0.150;
        }
        else {
            cut_iso[0] = 0.150; cut_iso[1] = 0.100;
        }
    }
    else if (workingPoint == 4) {//TIGHT
        cut_dEtaIn[0]        = 0.004; cut_dEtaIn[1]        = 0.005;
        cut_dPhiIn[0]        = 0.030; cut_dPhiIn[1]        = 0.020;
        cut_sigmaIEtaIEta[0] = 0.010; cut_sigmaIEtaIEta[1] = 0.030;
        cut_hoe[0]           = 0.120; cut_hoe[1]           = 0.100;
        cut_ooemoop[0]       = 0.050; cut_ooemoop[1]       = 0.050;
        cut_d0vtx[0]         = 0.020; cut_d0vtx[1]         = 0.020;
        cut_dzvtx[0]         = 0.100; cut_dzvtx[1]         = 0.100;
        cut_vtxFit[0]        = true ; cut_vtxFit[1]        = true;
        cut_mHits[0]         = 0    ; cut_mHits[1]         = 0;
        if (pt >= 20.0) {
            cut_iso[0] = 0.100; cut_iso[1] = 0.100;
        }
        else {
            cut_iso[0] = 0.100; cut_iso[1] = 0.070;
        }
    }
    else {
        std::cout << "[EgammaCutBasedEleId::TestWP] Undefined working point" << std::endl;
    }

    // choose cut if barrel or endcap
    unsigned int idx = isEB ? 0 : 1;

    // effective area for isolation
    mithep::ElectronTools eleT;
    mithep::ElectronTools::EElectronEffectiveAreaTarget EffectiveAreaVersion;
if (y == 2011){
    EffectiveAreaVersion  = mithep::ElectronTools::kEleEAData2011;}
if (y == 2012){
    EffectiveAreaVersion  = mithep::ElectronTools::kEleEAData2012;}

///
    float AEff = eleT.ElectronEffectiveArea(eleT.kEleGammaAndNeutralHadronIso04,eta,EffectiveAreaVersion);
 
   // apply to neutrals
    double rhoPrime = std::max(rho, 0.0);
    double iso_n = std::max(iso_nh + iso_em - rhoPrime * AEff, 0.0);

    // compute final isolation
    double iso = (iso_n + iso_ch) / pt;

    bool DETAIN = false;
    bool DPHIIN = false;
    bool SIGMAIETAIETA = false;
    bool HOE = false;
    bool OOEMOOP = false;
    bool D0VTX = false;
    bool DZVTX = false;
    bool VTXFIT = false;
    bool MHITS = false;
    bool ISO = false;
    // test cuts
    if (fabs(dEtaIn) < cut_dEtaIn[idx])             DETAIN = true;
    if (fabs(dPhiIn) < cut_dPhiIn[idx])             DPHIIN =true;
    if (sigmaIEtaIEta < cut_sigmaIEtaIEta[idx])     SIGMAIETAIETA = true;
    if (hoe < cut_hoe[idx])                         HOE = true;
    if (fabs(ooemoop) < cut_ooemoop[idx])           OOEMOOP = true;
    if (fabs(d0vtx) < cut_d0vtx[idx])               D0VTX = true;
    if (fabs(dzvtx) < cut_dzvtx[idx])               DZVTX = true;
    if (!cut_vtxFit[idx] || !vtxFitConversion)      VTXFIT = true;
    if (mHits <= cut_mHits[idx])                    MHITS = true;
    if (iso < cut_iso[idx])                         ISO = true;

    if (DETAIN && DPHIIN && SIGMAIETAIETA && HOE && OOEMOOP && D0VTX && DZVTX && VTXFIT && MHITS && ISO) return true;
    return false;

}
//--------------------------------------------------------------------------------------------------
bool ZGTools::electronPFIso04(const mithep::Electron * ele,
                      const mithep::Vertex * vtx,
                      const mithep::PFCandidateCol * fPFCandidates,
                      const mithep::PileupEnergyDensityCol * fPUEnergyDensity,
                      mithep::ElectronTools::EElectronEffectiveAreaTarget EffectiveAreaVersion,
		      vector<bool> PFnoPUflag,
		      float y)
//--------------------------------------------------------------------------------------------------
{

  //
  // final iso
  //
  Double_t fChargedIso = 0.0;
  Double_t fGammaIso = 0.0;
  Double_t fNeutralHadronIso = 0.0;


  //
  //Loop over PF Candidates
  //
  for(UInt_t k=0; k<fPFCandidates->GetEntries(); ++k) {

    const mithep::PFCandidate *pf = (mithep::PFCandidate*)((*fPFCandidates)[k]);

    //
    // veto FSR recovered photons
    //
    //bool vetoPhoton = false;
    //for( int p=0; p<photonsToVeto.size(); p++ ) {
    //  if( pf == photonsToVeto[p] ) {
    //    vetoPhoton = true;
    //    break;
    //  }
    //} if( vetoPhoton ) continue;

    //Double_t deta = (ele->Eta() - pf->Eta());
    //Double_t dphi = mithep::MathUtils::DeltaPhi(Double_t(ele->Phi()),Double_t(pf->Phi()));
    Double_t dr = mithep::MathUtils::DeltaR(ele->Phi(),ele->Eta(), pf->Phi(), pf->Eta());

    if (dr > 0.4) continue;
    //if (dr > 0.3) continue;
    if( !(PFnoPUflag[k]) ) continue; // my PF no PU hack

    //if(ctrl.debug) {
    //  cout << "pf :: type: " << pf->PFType() << "\tpt: " << pf->Pt() << "\tdR: " << dr;
    //  if( pf->HasTrackerTrk() ) cout << "\tdZ: " << pf->TrackerTrk()->DzCorrected(*vtx)
    //                                << "\ttrk: " << pf->HasTrackerTrk()
    //                                 << "\tgsf: " << pf->HasGsfTrk();
    //  cout << endl;
    //}
    // sync : I don't think theyre doing this ...
    //
    //     if ( (pf->HasTrackerTrk() && (pf->TrackerTrk() == ele->TrackerTrk())) ||
    //   (pf->HasGsfTrk() && (pf->GsfTrk() == ele->GsfTrk()))) { 
    //       if( ctrl.debug ) cout << "\tskipping, matches to the electron ..."  << endl;
    //       continue;
    //     }
    //
    // Lepton Footprint Removal
    //
    //Bool_t IsLeptonFootprint = kFALSE;
    if (dr < 1.0) {
    //
    // Charged Iso
    //
      //    if (pf->Charge() != 0 && (pf->HasTrackerTrk()||pf->HasGsfTrk()) ) {
      if (pf->Charge() != 0 && (pf->HasTrackerTrk()) ) {
    //  if( ctrl.debug) cout << "charged:: pt: " << pf->Pt()
    //                       << "\ttype: " << pf->PFType()
    //                       << "\ttrk: " << pf->TrackerTrk()
    //                       << "\tdeta: " << deta
    //                       << "\tdphi: " << dphi
    //                       << "\tdz: " << pf->TrackerTrk()->DzCorrected(*vtx)
    //                       << "\tdzy: " << pf->TrackerTrk()->D0Corrected(*vtx)
    //                       << "\tdR: " << dr
    //                       << "\teleSCeta: " << ele->SCluster()->Eta()
    //                       << endl;
      // check against ggA
      // Veto any PFmuon, or PFEle
      if (abs(pf->PFType()) == PFCandidate::eElectron || abs(pf->PFType()) == PFCandidate::eMuon) {
    //    if( ctrl.debug ) cout << "\t skipping, pf is and ele or mu .." <<endl;
        continue;
      }


      // for Zg, PV association from POG twiki
      /*
      ThreeVector sv = pf->SourceVertex();
      ThreeVector pv = vtx->Position();
      cout << "chIso sourceVtx: " <<  sv.X() <<":" << sv.Y() << ":" << sv.Z() << endl;
      cout << "chIso pv: " <<  pv.X() <<":" << pv.Y() << ":" << pv.Z() << endl;
      cout << "chIso pvpf: " <<  pvpf.X() <<":" << pvpf.Y() << ":" << pvpf.Z() << endl;
      if( sv != pvpf ) continue;
      */

      /*
      // PV association, take 2 
      const Track * pftrk = pf->TrackerTrk();
      bool found_trk=false;
      for(int t=0; t<vtx->NTracks(); t++ ) { 
        if( vtx->Trk(t) == pftrk ) { 
          found_trk = true;
          break;
        }
      }
      if( !found_trk ) { 
        cout << "\t\tno PV assoc, skipping ..." << endl;
        continue;
      }
      */

      /*
      // Zg, in E/G code
      if( fabs(pf->TrackerTrk()->DzCorrected(*vtx)) > 1. || 
          fabs(pf->TrackerTrk()->D0Corrected(*vtx)) >0.1  ) {
        cout << "\t\tpf fails dz,dxy ..." << endl;
          continue;
      }
      */

      // Footprint Veto
      if (fabs(ele->SCluster()->Eta()) > 1.479 && dr < 0.015) {
//        cout << "\t\tpf fails footprint dR ..." << endl;
        continue;
      }


      fChargedIso += pf->Pt();
    }

    //
    // Gamma Iso 
    //
    else if (abs(pf->PFType()) == PFCandidate::eGamma) {

      bool match=false;

      if( ele->PFSCluster() == pf->SCluster() ) {
        match=true;
//        cout << "gammaiso :: matching SC top level ... " << endl;
      }
//      if( ctrl.debug) cout << "gamma:: " << pf->Pt()
//                           << "\tdeta: " << deta
//                           << "\tdphi: " << dphi
//                           << "\tdr: " << dr
//                           << "\tmatch: " << match
//                           << "\tmva: " << pf->MvaGamma()
//                           << "\telemishits: " << ele->CorrectedNExpectedHitsInner()
//                           << "\tpfHasSC: " << pf->HasSCluster()
//                           << "\tused: ";


      // misshits for Zg, see http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc?revision=1.22&view=markup

      if( ele->CorrectedNExpectedHitsInner() > 0 &&
          pf->MvaGamma() > 0.99 ) { // pf->HasSCluster() always 0????
        if( match ) {
//          cout << 0  << endl;
          continue;
        }
      }


      if (fabs(ele->SCluster()->Eta()) > 1.479) {
        if (mithep::MathUtils::DeltaR(ele->Phi(),ele->Eta(), pf->Phi(), pf->Eta()) < 0.08) {
//          cout << 0  << endl;
          continue;
        }
      }


      // KH, add to sync 
      //      if( pf->Pt() > 0.5 ) 
        fGammaIso += pf->Pt();
//        cout << 1  << endl;
    }

    //
    // Neutral Iso 
    //
    else if( abs(pf->PFType()) == PFCandidate::eNeutralHadron ) {
//      if( ctrl.debug) cout << "neutral:: " << pf->Pt() << " "
//                           << dr << endl;
      // KH, add to sync 
      //      if( pf->Pt() > 0.5 ) 
        fNeutralHadronIso += pf->Pt();
    }


    }
  }

  double rho=0;
if (y == 2011){
  if( (EffectiveAreaVersion == mithep::ElectronTools::kEleEAFall11MC) ||
      (EffectiveAreaVersion == mithep::ElectronTools::kEleEAData2011) ) {
    if (!(std::isnan(fPUEnergyDensity->At(0)->RhoKt6PFJetsForIso25()) ||
          std::isinf(fPUEnergyDensity->At(0)->RhoKt6PFJetsForIso25())))
      rho = fPUEnergyDensity->At(0)->RhoKt6PFJetsForIso25();
    // !!!!!!!!!!!!! TMP HACK FOR SYNC !!!!!!!!!!!!!!!!!!!!!
    EffectiveAreaVersion  = mithep::ElectronTools::kEleEAData2011;
    // !!!!!!!!!!!!! TMP HACK FOR SYNC !!!!!!!!!!!!!!!!!!!!!
  } else {
    if (!(std::isnan(fPUEnergyDensity->At(0)->RhoKt6PFJets()) ||
          std::isinf(fPUEnergyDensity->At(0)->RhoKt6PFJets())))
      rho = fPUEnergyDensity->At(0)->RhoKt6PFJets();
    // !!!!!!!!!!!!! TMP HACK FOR SYNC !!!!!!!!!!!!!!!!!!!!!
    EffectiveAreaVersion  = mithep::ElectronTools::kEleEAData2012;
    // !!!!!!!!!!!!! TMP HACK FOR SYNC !!!!!!!!!!!!!!!!!!!!!
  }
}
if (y == 2012){
      rho = fPUEnergyDensity->At(0)->RhoKt6PFJets();
    EffectiveAreaVersion  = mithep::ElectronTools::kEleEAData2012;

}
//  if(ctrl.debug) cout << "rho: " << rho << endl;
  mithep::ElectronTools eleT;
  double pfIso = fChargedIso + fmax(0.0,(fGammaIso + fNeutralHadronIso
                                        -rho*eleT.ElectronEffectiveArea(eleT.kEleGammaAndNeutralHadronIso04,
                                                                   ele->Eta(),EffectiveAreaVersion)));

//  gChargedIso = fChargedIso;
//  gGammaIso = fGammaIso;
//  gNeutralIso = fNeutralHadronIso;

//  if( ctrl.debug ) {
//    cout << "PFiso: " << pfIso
//         << "\tfChargedIso: " << fChargedIso
//         << "\tfGammaIso: " << fGammaIso
//         << "\tfNeutralHadronIso: " << fNeutralHadronIso
//         << endl;
//  }

  if( (pfIso/ele->Pt()) < 0.4 ) return true;
  return false;
}

bool ZGTools::photonCutBasedLoose2012ID(const mithep::Photon *pho,
                                          const mithep::Vertex * vtx,
                                          const mithep::PFCandidateCol * fPFCandidates,
                                          const mithep::PileupEnergyDensityCol * fPUEnergyDensity,
                                          const mithep::ElectronCol * fElectrons,
                                          const mithep::DecayParticleCol * fConversions,
                                          const mithep::VertexCol * fVertexes )
//--------------------------------------------------------------------------------------------------
{
  bool pass=true;

  //
  // kinematics 
  //
  if( pho->Pt() < 10. )                            pass=false;
  if(fabs(pho->SCluster()->Eta())>2.5 ||
     ( fabs(pho->SCluster()->Eta())>1.4442 &&
       fabs(pho->SCluster()->Eta())<1.566 )  )     pass=false;
  if( !pass ) {
    return false;
  }

  //
  // ID
  //
  if( !PhotonTools::PassElectronVetoConvRecovery(pho, fElectrons, fConversions, vtx ) ) pass = false;
  if( pho->HadOverEmTow()  > 0.05 )                   pass=false;
  bool isEB = (pho->SCluster()->AbsEta()<1.5);
  if( pho->CoviEtaiEta() > (isEB ? 0.012 : 0.034) ) pass=false;

  return pass;
}

bool ZGTools::photonCutBasedLoose2012Isolation(const mithep::Photon * ph,
                                                 const mithep::VertexCol * vtxArr,
                                                 //const mithep::Vertex * vtx,
                                                 const mithep::PFCandidateCol * fPFCandidates,
                                                 const mithep::PileupEnergyDensityCol * fPUEnergyDensity,
						vector<bool> PFnoPUflag)
//-------------------------------------------------------------------------------------------------------
{

  bool pass=true;

  //  photonPFIso03( ctrl, ph, vtx, fPFCandidates, fPUEnergyDensity); 
  vector<double> isolations;
  isolations = photonPFIso03( ph, vtxArr, fPFCandidates, fPUEnergyDensity, PFnoPUflag);


  double rho = fPUEnergyDensity->At(0)->RhoKt6PFJets();
  //  double rho = fPUEnergyDensity->At(0)->RhoKt6PFJetsForIso25();  
  double fChargedIso = max(isolations[0] - rho*photonEffectiveEra_Ch(ph->SCluster()->Eta()), 0.);
  double fGammaIso   = max(isolations[1] - rho*photonEffectiveEra_Ga(ph->SCluster()->Eta()), 0.);
  double fNeutralIso = max(isolations[2] - rho*photonEffectiveEra_Nh(ph->SCluster()->Eta()), 0.);
  //cout << rho << endl;
  bool isEB = (ph->SCluster()->AbsEta()<1.5);
  if(isEB) {
    if( fChargedIso > 2.6 )                    {// cout << "eb, fail ch ..." << endl; 
pass=false;}
    if( fNeutralIso > (3.5 + 0.04*ph->Pt())  ) {// cout << "eb, fail nh ..." << (3.5 + 0.04*ph->Pt())<< endl; 
pass=false;}
    if( fGammaIso > (1.3 + 0.005*ph->Pt())  )  {// cout << "eb, fail ga ..." << (1.3 + 0.005*ph->Pt()) << endl; 
pass=false;}
  } else {
    if( fChargedIso > 2.3 )                    {// cout << "ee, fail ch ..." << endl; 
pass=false;}
    if( fNeutralIso > (2.9 + 0.04*ph->Pt())  ) {// cout << "ee, fail nh ..." << (1.5 + 0.04*ph->Pt()) << endl; 
pass=false;}
    //if( fGammaIso > (1.0 + 0.005*ph->Pt())  )  { cout << "ee, fail ga ..." << (1.0 + 0.005*ph->Pt()) << endl; pass=false;}
  }

  return pass;
}


bool ZGTools::photonCutBasedMedium2012ID(const mithep::Photon *pho,
                                          const mithep::Vertex * vtx,
                                          const mithep::PFCandidateCol * fPFCandidates,
                                          const mithep::PileupEnergyDensityCol * fPUEnergyDensity,
                                          const mithep::ElectronCol * fElectrons,
                                          const mithep::DecayParticleCol * fConversions,
                                          const mithep::VertexCol * fVertexes )
//--------------------------------------------------------------------------------------------------
{
  bool pass=true;

  //
  // kinematics 
  //
  if( pho->Pt() < 10. )                            pass=false;
  if(fabs(pho->SCluster()->Eta())>2.5 ||
     ( fabs(pho->SCluster()->Eta())>1.4442 &&
       fabs(pho->SCluster()->Eta())<1.566 )  )     pass=false;
  if( !pass ) {
    return false;
  }

  //
  // ID
  //
  if( !PhotonTools::PassElectronVetoConvRecovery(pho, fElectrons, fConversions, vtx ) ) pass = false;
  if( pho->HadOverEmTow()  > 0.05 )                   pass=false;
  bool isEB = (pho->SCluster()->AbsEta()<1.5);
  if( pho->CoviEtaiEta() > (isEB ? 0.011 : 0.033) ) pass=false;

  return pass;
}

bool ZGTools::photonCutBasedMedium2012Isolation(const mithep::Photon * ph,
                                                 const mithep::VertexCol * vtxArr,
                                                 //const mithep::Vertex * vtx,
                                                 const mithep::PFCandidateCol * fPFCandidates,
                                                 const mithep::PileupEnergyDensityCol * fPUEnergyDensity,
						vector<bool> PFnoPUflag,
						float year)
//-------------------------------------------------------------------------------------------------------
{

  bool pass=true;

  //  photonPFIso03( ctrl, ph, vtx, fPFCandidates, fPUEnergyDensity); 
  vector<double> isolations;
  isolations = photonPFIso03( ph, vtxArr, fPFCandidates, fPUEnergyDensity, PFnoPUflag);
  //cout << "isolations: " << isolations[0] << " " << isolations[1] << " " << isolations[2] << "\n";

  double rho = 0;

  if (year == 2011) {
    if (!(std::isnan(fPUEnergyDensity->At(0)->RhoKt6PFJetsForIso25()) ||
          std::isinf(fPUEnergyDensity->At(0)->RhoKt6PFJetsForIso25())))
      rho = fPUEnergyDensity->At(0)->RhoKt6PFJetsForIso25();
  }
  if (year == 2012) {
    rho = fPUEnergyDensity->At(0)->RhoKt6PFJets();
  }

  double fChargedIso = max(isolations[0] - rho*photonEffectiveEra_Ch(ph->SCluster()->Eta()), 0.);
  double fGammaIso   = max(isolations[1] - rho*photonEffectiveEra_Ga(ph->SCluster()->Eta()), 0.);
  double fNeutralIso = max(isolations[2] - rho*photonEffectiveEra_Nh(ph->SCluster()->Eta()), 0.);
  //cout << rho << endl;
  bool isEB = (ph->SCluster()->AbsEta()<1.5);

  //cout << "debug : ph iso: " << fChargedIso << " " << fGammaIso << " " << fNeutralIso << "\n";
  //cout << "ph iso corr: " << rho*photonEffectiveEra_Ch(ph->SCluster()->Eta()) << " " << rho*photonEffectiveEra_Ga(ph->SCluster()->Eta()) << " " << rho*photonEffectiveEra_Nh(ph->SCluster()->Eta()) << "\n";

  if(isEB) {
    if( fChargedIso > 1.5 )                    { //cout << "eb, fail ch ..." << endl; 
      pass=false;}
    if( fNeutralIso > (1.0 + 0.04*ph->Pt())  ) { //cout << "eb, fail nh ..." << (3.5 + 0.04*ph->Pt())<< endl; 
      pass=false;}
    if( fGammaIso > (0.7 + 0.005*ph->Pt())  )  { //cout << "eb, fail ga ..." << (1.3 + 0.005*ph->Pt()) << endl; 
      pass=false;}
  } else {
    if( fChargedIso > 1.2 )                    { //cout << "ee, fail ch ..." << endl; 
      pass=false;}
    if( fNeutralIso > (1.5 + 0.04*ph->Pt())  ) { //cout << "ee, fail nh ..." << (1.5 + 0.04*ph->Pt()) << endl; 
      pass=false;}
    if( fGammaIso > (1.0 + 0.005*ph->Pt())  )  { //cout << "ee, fail ga ..." << (1.0 + 0.005*ph->Pt()) << endl; 
      pass=false;}
  }

  return pass;
}



float ZGTools::photonEffectiveEra_Ga( float eta ) {
  if( fabs(eta)<=1.0 )          return 0.148;
  else if(fabs(eta)>1.0 &&
          fabs(eta)<=1.479 )    return 0.130;
  else if(fabs(eta)>1.479 &&
          fabs(eta)<=2.0 )      return 0.112;
  else if(fabs(eta)>2.0 &&
          fabs(eta)<=2.2 )      return 0.216;
  else if(fabs(eta)>2.2 &&
          fabs(eta)<=2.3 )      return 0.262;
  else if(fabs(eta)>2.3 &&
          fabs(eta)<=2.4 )      return 0.260;
  else                          return 0.266;
} 

  
float ZGTools::photonEffectiveEra_Nh( float eta ) {
  if( fabs(eta)<=1.0 )          return 0.030;
  else if(fabs(eta)>1.0 &&
          fabs(eta)<=1.479 )    return 0.057;
  else if(fabs(eta)>1.479 &&  
          fabs(eta)<=2.0 )      return 0.039;
  else if(fabs(eta)>2.0 &&
          fabs(eta)<=2.2 )      return 0.015;
  else if(fabs(eta)>2.2 &&
          fabs(eta)<=2.3 )      return 0.024;
  else if(fabs(eta)>2.3 &&
          fabs(eta)<=2.4 )      return 0.039;
  else                          return 0.072;
}

float ZGTools::photonEffectiveEra_Ch( float eta ) {
  if( fabs(eta)<=1.0 )          return 0.012;
  else if(fabs(eta)>1.0 &&
          fabs(eta)<=1.479 )    return 0.010;
  else if(fabs(eta)>1.479 &&
          fabs(eta)<=2.0 )      return 0.014;
  else if(fabs(eta)>2.0 &&
          fabs(eta)<=2.2 )      return 0.012;
  else if(fabs(eta)>2.2 &&
          fabs(eta)<=2.3 )      return 0.016;
  else if(fabs(eta)>2.3 &&
          fabs(eta)<=2.4 )      return 0.020;
  else                          return 0.012;

}

//--------------------------------------------------------------------------------------------------
vector<double> ZGTools::photonPFIso03(const mithep::Photon * ph,
                    const mithep::VertexCol * vtxArr,
                      const mithep::PFCandidateCol * fPFCandidates,
                      const mithep::PileupEnergyDensityCol * fPUEnergyDensity,
		      vector<bool> PFnoPUflag	)
//--------------------------------------------------------------------------------------------------
{

  bool isEB = (ph->SCluster()->AbsEta()<1.5);

  //
  // final iso 
  //
  Double_t fChargedIso  = 0.0;
  Double_t fGammaIso  = 0.0;
  Double_t fNeutralHadronIso  = 0.0;

  //
  //Loop over PF Candidates
  //
  for(UInt_t k=0; k<fPFCandidates->GetEntries(); ++k) {

    if( !(PFnoPUflag[k]) ) continue; // my PF no PU hack
    const mithep::PFCandidate *pf = (mithep::PFCandidate*)((*fPFCandidates)[k]);

    /*
    Double_t deta = (ph->SCluster()->Eta() - pf->Eta());
    Double_t dphi = mithep::MathUtils::DeltaPhi(Double_t(ph->SCluster()->Phi()),Double_t(pf->Phi()));
    Double_t dr = mithep::MathUtils::DeltaR(ph->SCluster()->Phi(),ph->SCluster()->Eta(), pf->Phi(), pf->Eta());

    Double_t deta2 = (ph->Eta() - pf->Eta());
    Double_t dphi2 = mithep::MathUtils::DeltaPhi(Double_t(ph->Phi()),Double_t(pf->Phi()));
    Double_t dr2 = mithep::MathUtils::DeltaR(ph->Phi(),ph->Eta(), pf->Phi(), pf->Eta());
    */

    const mithep::Vertex * vtx = (mithep::Vertex *)(vtxArr->At(0));

    double fEta, fPhi;
    if( pf->PFType() == PFCandidate::eGamma ||  pf->PFType() == PFCandidate::eNeutralHadron ) {
      TVector3 direction( ph->SCluster()->Point().X() - pf->SourceVertex().X(),
                          ph->SCluster()->Point().Y() - pf->SourceVertex().Y(),
                          ph->SCluster()->Point().Z() - pf->SourceVertex().Z() );

      fEta = direction.Eta();
      fPhi = direction.Phi();
    } else {
      TVector3 direction( ph->SCluster()->Point().X() - vtx->X(),
                          ph->SCluster()->Point().Y() - vtx->Y(),
                          ph->SCluster()->Point().Z() - vtx->Z() );

      fEta = direction.Eta();
      fPhi = direction.Phi();
    }

    double deta = (fEta - pf->Eta());
    double dphi = mithep::MathUtils::DeltaPhi(fPhi,Double_t(pf->Phi()));
    double dr   = mithep::MathUtils::DeltaR(fPhi,fEta, pf->Phi(), pf->Eta());


    if (dr > 0.3) continue;
    if (abs(pf->PFType()) == PFCandidate::eElectron || abs(pf->PFType()) == PFCandidate::eMuon) continue;

    //
    // Charged Iso 
    //
    if (pf->Charge() != 0 ) {
      if( fabs(pf->TrackerTrk()->D0Corrected(*vtx))>0.1 ||
          fabs(pf->TrackerTrk()->DzCorrected(*vtx))>0.2   ) continue;
      if( dr < 0.02 ) continue; // from the PF iso page???
      fChargedIso += pf->Pt();
      //cout << "debug: charged " << pf->Pt() << "\n";
    }

    //
    // Gamma Iso 
    //
    else if (abs(pf->PFType()) == PFCandidate::eGamma  ) {

      if( isEB && ( fabs(deta)<0.015 && fabs(dphi)<1.)) continue; // && fabs(dphi2)<1.))
      //if( !isEB &&  dr<0.07 ) continue;
      if( !isEB &&  dr<0.00864*fabs(sinh(ph->SCluster()->Eta()))*4 ) continue;

      //cout << "debug: gamma iso: " << dr << " " << pf->Pt() << " : " << deta << " " << dphi << " " << 0.00864*fabs(sinh(ph->SCluster()->Eta()))*4 << " " << isEB << "\n";
      fGammaIso += pf->Pt();
    }

    //
    // Other Neutrals
    //
    else {
      // not used for Zg ...
      //if( pf->Pt() > 0.5 ) {
        fNeutralHadronIso += pf->Pt();   
	//cout << "debug: neutral hadron " << pf->Pt() << "\n";

        //}
    }
  }


  //double rho = fPUEnergyDensity->At(0)->RhoKt6PFJets();

//   gChargedIso = max(fChargedIso - rho*photonEffectiveEra_Ch(ph->SCluster()->Eta()), 0.);
//   gGammaIso   = max(fGammaIso - rho*photonEffectiveEra_Ga(ph->SCluster()->Eta()), 0.);
//   gNeutralIso = max(fNeutralHadronIso - rho*photonEffectiveEra_Nh(ph->SCluster()->Eta()), 0.);
  
vector<double> iso;
iso.push_back(fChargedIso);
iso.push_back(fGammaIso);
iso.push_back(fNeutralHadronIso);

//cout << "iso: " << iso[0] << " " << iso[1] << " " << iso[2] << " \n";

return iso;
}

bool ZGTools::muonIDPOGTightSelection(int year,
                                        const mithep::Muon *mu,
                                        const mithep::Vertex * vtx,
                                        const mithep::PFCandidateCol * pfCandidates )
//--------------------------------------------------------------------------------------------------
{
  bool pass = true;

  //
  // 2012 tight ID : https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Tight_Muon
  //
  if( year == 2012 ) {
    if(!(mu->IsGlobalMuon()) )                          pass = false;
    if(!(mu->IsPFMuon()) )                              pass = false;
    if(mu->GlobalTrk()->RChi2() > 10.)                  pass = false;
    if(mu->NValidHits() <= 0 )                          pass = false;
    if(mu->NSegments() <= 1)                             pass = false;
    if(fabs(mu->TrackerTrk()->D0Corrected(*vtx)) > 0.2) pass = false;
    if(fabs(mu->TrackerTrk()->DzCorrected(*vtx)) > 0.5) pass = false;
    if(mu->BestTrk()->NPixelHits() <= 0)                pass = false;
    if(mu->NTrkLayersHit() <= 5)                        pass = false;
  }
  //
  // 2011 tight ID : https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Tight_Muon_selection
  //
  else if (year == 2011 ) {
    //cout << "muon ID era is 2011..." << endl;
    if(!(mu->IsGlobalMuon()) )                          pass = false;
    if(mu->GlobalTrk()->RChi2() > 10.)                  pass = false;
    if(mu->NValidHits() <= 0 )                          pass = false;
    if(mu->NSegments() <= 1)                            pass = false;
    if(fabs(mu->TrackerTrk()->D0Corrected(*vtx)) > 0.2) pass = false;
    if(mu->BestTrk()->NPixelHits() <= 0)                pass = false;
    if(mu->NTrkLayersHit() <= 8)                        pass = false;
  } else {
    pass = false;
  }

  return pass;
}


bool ZGTools::muonPFIso04(const mithep::Muon * mu,
                    const mithep::Vertex * vtx,
                    const mithep::PFCandidateCol * fPFCandidates,
                    const mithep::PileupEnergyDensityCol * fPUEnergyDensity,
                    mithep::MuonTools::EMuonEffectiveAreaTarget EffectiveAreaVersion,
		    vector<bool> PFnoPUflag,
		    Float_t y)
//--------------------------------------------------------------------------------------------------
{

  //
  // final iso 
  //
  Double_t fChargedIso  = 0.0;
  Double_t fGammaIso  = 0.0;
  Double_t fNeutralHadronIso  = 0.0;

  //
  //Loop over PF Candidates
  //
 // std::cout << fPFCandidates->GetEntries() << " " << PFnoPUflag.size() << std::endl;
  for(UInt_t k=0; k<fPFCandidates->GetEntries(); ++k) {
    
    if( !(PFnoPUflag[k]) ) continue; // my PF no PU hack
    const mithep::PFCandidate *pf = (mithep::PFCandidate*)((*fPFCandidates)[k]);


    //Double_t deta = (mu->Eta() - pf->Eta());
    //Double_t dphi = mithep::MathUtils::DeltaPhi(Double_t(mu->Phi()),Double_t(pf->Phi()));
    Double_t dr = mithep::MathUtils::DeltaR(mu->Phi(),mu->Eta(), pf->Phi(), pf->Eta());


    if (dr > 0.4) continue;
    if (pf->HasTrackerTrk() && (pf->TrackerTrk() == mu->TrackerTrk()) ) continue;

    //
    // Charged Iso 
    //
    if (pf->Charge() != 0 && (pf->HasTrackerTrk()||pf->HasGsfTrk()) ) {

      //if( dr < 0.01 ) continue; // only for muon iso mva?
      if (abs(pf->PFType()) == PFCandidate::eElectron || abs(pf->PFType()) == PFCandidate::eMuon) continue;
      fChargedIso += pf->Pt();
    }

    //
    // Gamma Iso 
    //
    else if (abs(pf->PFType()) == PFCandidate::eGamma) {
      //      if( ctrl.era==2011  ) { 
      //      fGammaIso += pf->Pt();
      //      cout << "\tused"; 
      //      }
      //      else if( pf->Pt() > 0.5 && dr > 0.01 ) {
      if( pf->Pt() > 0.5 && dr > 0.01 ) {
      fGammaIso += pf->Pt();
      }
    }

    //
    // Other Neutrals
    //
    else {
      //       if( ctrl.era == 2011 ) { 
      //        fNeutralHadronIso += pf->Pt();
      //       }
      //      else if( pf->Pt() > 0.5 ) {
      if( pf->Pt() > 0.5 && dr > 0.01) {
        fNeutralHadronIso += pf->Pt();
      }
    }
  }

  double rho=0;
  if (y == 2011){
    if( (EffectiveAreaVersion == mithep::MuonTools::kMuEAFall11MC) ||
	(EffectiveAreaVersion == mithep::MuonTools::kMuEAData2011) ) {
      if (!(std::isnan(fPUEnergyDensity->At(0)->RhoKt6PFJetsForIso25()) ||
	    std::isinf(fPUEnergyDensity->At(0)->RhoKt6PFJetsForIso25())))
	rho = fPUEnergyDensity->At(0)->RhoKt6PFJetsForIso25();
      //rho = fPUEnergyDensity->At(0)->Rho();
      // !!!!!!!!!!!!! TMP HACK FOR SYNC !!!!!!!!!!!!!!!!!!!!!
      EffectiveAreaVersion  = mithep::MuonTools::kMuEAData2011;
      // !!!!!!!!!!!!! TMP HACK FOR SYNC !!!!!!!!!!!!!!!!!!!!!
    } else {
      if (!(std::isnan(fPUEnergyDensity->At(0)->RhoKt6PFJetsCentralNeutral()) ||
	    std::isinf(fPUEnergyDensity->At(0)->RhoKt6PFJetsCentralNeutral())))
	rho = fPUEnergyDensity->At(0)->RhoKt6PFJetsCentralNeutral();
      // !!!!!!!!!!!!! TMP HACK FOR SYNC !!!!!!!!!!!!!!!!!!!!!
      EffectiveAreaVersion  = mithep::MuonTools::kMuEAData2012;
      // !!!!!!!!!!!!! TMP HACK FOR SYNC !!!!!!!!!!!!!!!!!!!!!
    }
  }
  if (y == 2012){
    rho = fPUEnergyDensity->At(0)->RhoKt6PFJetsCentralNeutral();
    EffectiveAreaVersion  = mithep::MuonTools::kMuEAData2012;
    
    
  }
  //  TLorentzVector  tmpvec;
  //  tmpvec.SetPtEtaPhiM(mu->Pt(),mu->Eta(),mu->Phi(),mu->Mass());
  //  for( int p=0; p<photonsToVeto.size(); p++ ) {
  //    const mithep::PFCandidate * pf  = photonsToVeto[p];
  //    TLorentzVector pfvec;
  //    pfvec.SetPtEtaPhiM(pf->Pt(),pf->Eta(),pf->Phi(),0.);
  //    tmpvec += pfvec;
  //  }
  mithep::MuonTools muT;
  double pfIso = fChargedIso + fmax(0.0,(fGammaIso + fNeutralHadronIso
					 -rho*muT.MuonEffectiveArea(muT.kMuGammaAndNeutralHadronIso04,
								    //tmpvec.Eta(),EffectiveAreaVersion)));
								    mu->Eta(),EffectiveAreaVersion)));
  bool pass = false;
  if( (pfIso/mu->Pt()) < 0.12 ) pass = true;
  return pass;
}

//--------------------------------------------------------------------------------------------------
// KH implementation of phosphor corrections.  NB: r9 dependence is not yet provieded
// so it's a no-op ATM
//
std::pair<float,float> ZGTools::getPhosphorScaleRes(unsigned year, 
						    bool isMC,
						    double scEta, 
						    double pt, 
						    double r9) 
//--------------------------------------------------------------------------------------------------
{
  
  
  bool isEB = (fabs(scEta)<=1.4442); 
  int ptbin=-1;
  if     ( pt >10. && pt<= 12.  ) ptbin=0;
  else if( pt >12. && pt<= 15.  ) ptbin=1;
  else if( pt >15. && pt<= 20.  ) ptbin=2;
  else if( pt >20. && pt<= 999. ) ptbin=3;
  
  assert(ptbin>=0 );

  if( year == 2011 ) {
    if( isMC ) { 
      if( isEB ) { 
	if( ptbin==0 )      return std::pair<float,float> (1.43,5.71);
	else if( ptbin==1 ) return std::pair<float,float> (1.29,5.16);
	else if( ptbin==2 ) return std::pair<float,float> (0.86,4.15);
	else                return std::pair<float,float> (0.42,2.68);
      } else { // isEE 
	if( ptbin==0 )      return std::pair<float,float> (2.90,7.64);
	else if( ptbin==1 ) return std::pair<float,float> (1.91,6.12);
	else if( ptbin==2 ) return std::pair<float,float> (1.83,5.15);
	else                return std::pair<float,float> (0.87,3.44);
      } // isEE
    } else { // isData
      if( isEB ) { 
	if( ptbin==0 )      return std::pair<float,float> (-0.13,6.03);
	else if( ptbin==1 ) return std::pair<float,float> (0.41,5.16);
	else if( ptbin==2 ) return std::pair<float,float> (0.01,4.10);
	else                return std::pair<float,float> (0.42,2.88);
      } else { // isEE 
	if( ptbin==0 )      return std::pair<float,float> (1.95,7.60);
	else if( ptbin==1 ) return std::pair<float,float> (0.86,6.79);
	else if( ptbin==2 ) return std::pair<float,float> (0.65,6.05);
	else                return std::pair<float,float> (0.28,4.34);
      } // isEE
    } // isData


  } else if ( year == 2012 ) { 
    if( isMC ) { 
      if( isEB ) { 
	if( ptbin==0 )      return std::pair<float,float> (0.85,4.00);
	else if( ptbin==1 ) return std::pair<float,float> (0.36,4.32);
	else if( ptbin==2 ) return std::pair<float,float> (0.48,3.36);
	else                return std::pair<float,float> (0.30,1.92);
      } else { // isEE
	if( ptbin==0 )      return std::pair<float,float> (1.57,6.68);
	else if( ptbin==1 ) return std::pair<float,float> (1.79,5.59);
	else if( ptbin==2 ) return std::pair<float,float> (0.96,4.61);
	else                return std::pair<float,float> (0.60,3.02);
      }
    } else { //isData
      if( isEB ) { 
	if( ptbin==0 )      return std::pair<float,float> (-0.18,5.27);
	else if( ptbin==1 ) return std::pair<float,float> (-0.22,5.07);
	else if( ptbin==2 ) return std::pair<float,float> (0.01,4.84);
	else                return std::pair<float,float> (0.47,2.49);
      } else { // isEE
	if( ptbin==0 )      return std::pair<float,float> (-1.24,10.94);
	else if( ptbin==1 ) return std::pair<float,float> (-1.47,8.08);
	else if( ptbin==2 ) return std::pair<float,float> (-0.69,5.37);
	else                return std::pair<float,float> (1.12,5.45);
      } // isEE
    } // isData
  } // 2012

  return std::pair<float,float> (0.,0.);

}
