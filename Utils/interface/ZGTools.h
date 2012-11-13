#ifndef MITPHYSICS_MODS_ZGTOOLS_H
#define MITPHYSICS_MODS_ZGTOOLS_H

#include "TLorentzVector.h"
#include "MitAna/DataTree/interface/PileupEnergyDensity.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
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
#include <utility> // for std::pair 
namespace mithep
{

class ZGAngles { 
 public:
  float costheta_lm, costheta_lp, phi, cosTheta, cosThetaG;
  float ptg, ptz, ptl1, ptl2, etag, etaz, etal1, etal2, mzg, mz;
};



class ZGLabVectors { 
 public:
  TLorentzVector veckq;
  TLorentzVector veckqbar;

  TLorentzVector veczg;
  TLorentzVector vecz;
  TLorentzVector vecg;
  TLorentzVector veclp;
  TLorentzVector veclm;
};

class ZGTools
	{
	public:
	ZGTools();
	
	static ZGAngles getZGAngles( ZGLabVectors &l, bool );

	static bool electronCutBasedIDLoose(   	const mithep::Electron *ele,
                                         	const mithep::Vertex * vtx,
                                         	const mithep::DecayParticleCol *conversions,
						const float year);

	static bool PassWP(	int workingPoint, const bool isEB, const float pt, const float eta,
    				const float dEtaIn, const float dPhiIn, const float sigmaIEtaIEta, const float hoe,
    				const float ooemoop, const float d0vtx, const float dzvtx, const float iso_ch, const float iso_em, const float iso_nh,
    				const bool vtxFitConversion, const unsigned int mHits, const double rho, const float y);

	static bool electronPFIso04(	const mithep::Electron * ele,
                      			const mithep::Vertex * vtx,
                      			const mithep::PFCandidateCol * fPFCandidates,
                      			const mithep::PileupEnergyDensityCol * fPUEnergyDensity,
                      			mithep::ElectronTools::EElectronEffectiveAreaTarget EffectiveAreaVersion,
                      			vector<bool> PFnoPUflag,
					float y);

	static bool photonCutBasedLoose2012ID(	const mithep::Photon *pho,
                                          	const mithep::Vertex * vtx,
                                          	const mithep::PFCandidateCol * fPFCandidates,
                                          	const mithep::PileupEnergyDensityCol * fPUEnergyDensity,
                                          	const mithep::ElectronCol * fElectrons,
                                          	const mithep::DecayParticleCol * fConversions,
                                          	const mithep::VertexCol * fVertexes );


	static bool photonCutBasedLoose2012Isolation(	const mithep::Photon * ph,
                                                 	const mithep::VertexCol * vtxArr,
                                                 	const mithep::PFCandidateCol * fPFCandidates,
                                                 	const mithep::PileupEnergyDensityCol * fPUEnergyDensity,
							vector<bool> PFnoPUflag);

	static bool photonCutBasedMedium2012ID(	const mithep::Photon *pho,
                                          	const mithep::Vertex * vtx,
                                          	const mithep::PFCandidateCol * fPFCandidates,
                                          	const mithep::PileupEnergyDensityCol * fPUEnergyDensity,
                                          	const mithep::ElectronCol * fElectrons,
                                          	const mithep::DecayParticleCol * fConversions,
                                          	const mithep::VertexCol * fVertexes );


	static bool photonCutBasedMedium2012Isolation(	const mithep::Photon * ph,
                                                 	const mithep::VertexCol * vtxArr,
                                                 	const mithep::PFCandidateCol * fPFCandidates,
                                                 	const mithep::PileupEnergyDensityCol * fPUEnergyDensity,
							vector<bool> PFnoPUflag, 
							float year );

	static float photonEffectiveEra_Ga( float eta ); 

	static float photonEffectiveEra_Nh( float eta );

	static float photonEffectiveEra_Ch( float eta );

	static vector<double> photonPFIso03(	const mithep::Photon * ph,
                    				const mithep::VertexCol * vtxArr,
                      				const mithep::PFCandidateCol * fPFCandidates,
                      				const mithep::PileupEnergyDensityCol * fPUEnergyDensity,
						vector<bool> PFnoPUflag);

	static bool muonIDPOGTightSelection(	int year,
                                        	const mithep::Muon *mu,
                                        	const mithep::Vertex * vtx,
                                        	const mithep::PFCandidateCol * pfCandidates );

	static bool muonPFIso04(const mithep::Muon * mu,
                    		const mithep::Vertex * vtx,
                    		const mithep::PFCandidateCol * fPFCandidates,
                    		const mithep::PileupEnergyDensityCol * fPUEnergyDensity,
                    		mithep::MuonTools::EMuonEffectiveAreaTarget EffectiveAreaVersion,
				vector<bool> PFnoPUflag,
				float y);

	static std::pair<float,float> getPhosphorScaleRes(unsigned year, 
							  bool isMC,
							  double scEta, 
							  double pt, 
							  double r9); 

ClassDef(ZGTools,0)
};

}
#endif
