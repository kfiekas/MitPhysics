#include "MitPhysics/Utils/interface/RecoilTools.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include <TFile.h>

ClassImp(mithep::RecoilTools)

using namespace mithep;

RecoilTools::RecoilTools(TString iJetMVAFile) { 
  fJetIDMVA = new JetIDMVA();
  fJetIDMVA->Initialize( "JetIDMVA",iJetMVAFile,JetIDMVA::kBaseline);
}
//--------------------------------------------------------------------------------------------------
RecoilTools::~RecoilTools() { 
  delete fJetIDMVA;
}
//--------------------------------------------------------------------------------------------------
bool RecoilTools::filter(const PFJet *iJet,Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2) { 
  double pDEta1 = iJet->Eta() - iEta1;
  double pDPhi1 = fabs(iJet->Phi() - iPhi1); if(pDPhi1 > 2.*TMath::Pi()-pDPhi1) pDPhi1 = 2.*TMath::Pi()-pDPhi1;
  double pDR1   = sqrt(pDEta1*pDEta1 + pDPhi1*pDPhi1);
  if(pDR1 < 0.5) return false;
  double pDEta2 = iJet->Eta() - iEta2;
  double pDPhi2 = fabs(iJet->Phi() - iPhi2); if(pDPhi2 > 2.*TMath::Pi()-pDPhi2) pDPhi2 = 2.*TMath::Pi()-pDPhi2;
  double pDR2   = sqrt(pDEta2*pDEta2 + pDPhi2*pDPhi2);
  if(pDR2 < 0.5) return false;
  return true;
}

//--------------------------------------------------------------------------------------------------
Met RecoilTools::pfRecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,
			  const PFMet *iMet) { 
  Met lPFMet(iMet->Px(),iMet->Py()); 
  lPFMet.SetSumEt(iMet->SumEt());
  lPFMet.SetMex  (lPFMet.Mex()+iVisPt*cos(iVisPhi));  
  lPFMet.SetMey  (lPFMet.Mey()+iVisPt*sin(iVisPhi));
  lPFMet.SetSumEt(lPFMet.SumEt()-iVisSumEt);
  return lPFMet;
}
//--------------------------------------------------------------------------------------------------
Met RecoilTools::trackMet(const PFCandidateCol *iCands,const Vertex *iVertex,Double_t iDZCut) { 
  double trkMetx  = 0;
  double trkMety  = 0;
  double trkSumEt = 0; 
  for(UInt_t i=0; i<iCands->GetEntries(); ++i) {
    const PFCandidate *pfcand = iCands->At(i);
    if( (pfcand->HasTrackerTrk() && (fabs(pfcand->TrackerTrk()->DzCorrected(*iVertex))< iDZCut)) ||
        (pfcand->HasGsfTrk()     && (fabs(pfcand->GsfTrk()->DzCorrected(*iVertex))    < iDZCut)) ) {
      trkMetx  -= pfcand->Px();
      trkMety  -= pfcand->Py();
      trkSumEt += pfcand->Pt();
    }
  }
  Met lMet(trkMetx,trkMety);
  lMet.SetSumEt(trkSumEt);
  return lMet;
}

//--------------------------------------------------------------------------------------------------
//Compute the recoil => here this requires the vector sum of the visible components
//VisPt    => Vector sum pT  of the visible non-recoiling components
//VisPhi   => Vector sum Phi of the visible non-recoiling components
//visSumEt => Vector sum Et  =f the visible non-recoiling components
Met RecoilTools::trackRecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,
			     const PFCandidateCol *iCands,const Vertex *iVertex,double iDZCut) { 
  Met lTrkMet = trackMet(iCands,iVertex,iDZCut);
  lTrkMet.SetMex  (lTrkMet.Mex()+iVisPt*cos(iVisPhi));  
  lTrkMet.SetMey  (lTrkMet.Mey()+iVisPt*sin(iVisPhi));
  lTrkMet.SetSumEt(lTrkMet.SumEt()-iVisSumEt);
  return lTrkMet;
}
//--------------------------------------------------------------------------------------------------
void RecoilTools::addNeut(const PFJet *iJet,FourVectorM &iVec,Double_t &iSumEt, 
			  FactorizedJetCorrector *iJetCorrector,const PileupEnergyDensityCol *iPUEnergyDensity,
			  int iSign) { 
  FourVectorM lVec(0,0,0,0);
  double lPt = fJetIDMVA->correctedPt(iJet,iJetCorrector,iPUEnergyDensity);
  lPt *= (iJet->NeutralEmEnergy()/iJet->E() + iJet->NeutralHadronEnergy()/iJet->E());
  lVec.SetPt(lPt); lVec.SetEta(iJet->Eta()); lVec.SetPhi(iJet->Phi()); lVec.SetM(iJet->Mass());
  if(iSign > 0) iVec -= lVec;
  if(iSign < 0) iVec += lVec;
  iSumEt += lPt;
}

//--------------------------------------------------------------------------------------------------
//Corrected Jets
void RecoilTools::addNeut(const PFJet *iJet,FourVectorM &iVec,Double_t &iSumEt,int iSign) { 
  FourVectorM lVec(0,0,0,0);
  double lPt = iJet->Pt();
  lPt *= (iJet->NeutralEmEnergy()/iJet->E() + iJet->NeutralHadronEnergy()/iJet->E());
  lVec.SetPt(lPt); lVec.SetEta(iJet->Eta()); lVec.SetPhi(iJet->Phi()); lVec.SetM(iJet->Mass());
  if(iSign > 0) iVec -= lVec;
  if(iSign < 0) iVec += lVec;
  iSumEt += lPt;
}

//--------------------------------------------------------------------------------------------------
Met RecoilTools::NoPUMet( const PFJetCol       *iJets,FactorizedJetCorrector *iJetCorrector,
			  const PileupEnergyDensityCol *iPileupEnergyDensity,
			  const PFCandidateCol *iCands,const Vertex *iVertex,
			  Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2,Double_t iDZCut) { 

  FourVectorM lVec        (0,0,0,0); double lSumEt          = 0; 
  for(UInt_t i0 = 0; i0 < iCands->GetEntries(); i0++) { 
    const PFCandidate *pPF = iCands->At(i0);
    const Track* pTrack = pPF->TrackerTrk();
    if(pPF->GsfTrk()) pTrack = pPF->GsfTrk();
    if(pTrack        ==  0                           ) continue;
    if(	 !((pPF->HasTrackerTrk() && (fabs(pPF->TrackerTrk()->DzCorrected(*iVertex))<iDZCut)) ||
	   (pPF->HasGsfTrk()     && (fabs(pPF->GsfTrk()->DzCorrected(*iVertex))    <iDZCut)))) continue; 
    lVec     -= pPF->Mom();
    lSumEt   += pPF->Pt();
  }
  for(UInt_t i0 = 0; i0 < iJets->GetEntries(); i0++) {
    const PFJet *pJet = iJets->At(i0);
    if(!fJetIDMVA->pass(pJet,iVertex,iJetCorrector,iPileupEnergyDensity)) continue;
    if(!filter(pJet,iPhi1,iEta1,iPhi2,iEta2))                             continue; //Quick cleaning==> if not done already
    addNeut(pJet,lVec,lSumEt,iJetCorrector,iPileupEnergyDensity);
  }
  Met lMet(lVec.Px(),lVec.Py());
  lMet.SetSumEt(lSumEt);
  return lMet;
}
//--------------------------------------------------------------------------------------------------
//Corrected Jets
Met RecoilTools::NoPUMet( const PFJetCol       *iJets,const PFCandidateCol *iCands,const Vertex *iVertex,
			  Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2,Double_t iDZCut) { 

  FourVectorM lVec        (0,0,0,0); double lSumEt          = 0; 
  for(UInt_t i0 = 0; i0 < iCands->GetEntries(); i0++) { 
    const PFCandidate *pPF = iCands->At(i0);
    const Track* pTrack = pPF->TrackerTrk();
    if(pPF->GsfTrk()) pTrack = pPF->GsfTrk();
    if(pTrack        ==  0                           ) continue;
    if(	 !((pPF->HasTrackerTrk() && (fabs(pPF->TrackerTrk()->DzCorrected(*iVertex))<iDZCut)) ||
	   (pPF->HasGsfTrk()     && (fabs(pPF->GsfTrk()->DzCorrected(*iVertex))    <iDZCut)))) continue; 
    lVec     -= pPF->Mom();
    lSumEt   += pPF->Pt();
  }
  for(UInt_t i0 = 0; i0 < iJets->GetEntries(); i0++) {
    const PFJet *pJet = iJets->At(i0);
    if(!fJetIDMVA->pass(pJet,iVertex))                                    continue;
    if(!filter(pJet,iPhi1,iEta1,iPhi2,iEta2))                             continue; //Quick cleaning==> if not done already
    addNeut(pJet,lVec,lSumEt);
  }
  Met lMet(lVec.Px(),lVec.Py());
  lMet.SetSumEt(lSumEt);
  return lMet;
}
//--------------------------------------------------------------------------------------------------
Met RecoilTools::NoPURecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,
			    const PFJetCol       *iJets,FactorizedJetCorrector *iJetCorrector,
			    const PileupEnergyDensityCol *iPileupEnergyDensity,
			    const PFCandidateCol *iCands,const Vertex *iVertex,
			    Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2,
			    Double_t iDZCut) { 
  
  Met lNoPUMet = NoPUMet(iJets,iJetCorrector,iPileupEnergyDensity,iCands,iVertex,iPhi1,iEta1,iPhi2,iEta2,iDZCut);
  lNoPUMet.SetMex  (lNoPUMet.Mex()+iVisPt*cos(iVisPhi));  
  lNoPUMet.SetMey  (lNoPUMet.Mey()+iVisPt*sin(iVisPhi));
  lNoPUMet.SetSumEt(lNoPUMet.SumEt()-iVisSumEt);
  return lNoPUMet;
}
//--------------------------------------------------------------------------------------------------
//Corrected Jets
Met RecoilTools::NoPURecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,
			    const PFJetCol       *iJets,const PFCandidateCol *iCands,const Vertex *iVertex,
			    Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2,
			    Double_t iDZCut) { 
  
  Met lNoPUMet = NoPUMet(iJets,iCands,iVertex,iPhi1,iEta1,iPhi2,iEta2,iDZCut);
  lNoPUMet.SetMex  (lNoPUMet.Mex()+iVisPt*cos(iVisPhi));  
  lNoPUMet.SetMey  (lNoPUMet.Mey()+iVisPt*sin(iVisPhi));
  lNoPUMet.SetSumEt(lNoPUMet.SumEt()-iVisSumEt);
  return lNoPUMet;
}
//--------------------------------------------------------------------------------------------------
Met RecoilTools::PUCMet( const PFJetCol       *iJets,FactorizedJetCorrector *iJetCorrector,
			 const PileupEnergyDensityCol *iPileupEnergyDensity,
			 const PFCandidateCol *iCands,const Vertex *iVertex,
			 Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2,
			 Double_t iDZCut) { 

  FourVectorM lVec        (0,0,0,0); double lSumEt          = 0; 
  for(UInt_t i0 = 0; i0 < iCands->GetEntries(); i0++) { 
    const PFCandidate *pPF = iCands->At(i0);
    const Track* pTrack = pPF->TrackerTrk();
    if(pPF->GsfTrk()) pTrack = pPF->GsfTrk();
    if(pTrack == 0                                   && 
       (pPF->PFType() == PFCandidate::eGamma         || 
	pPF->PFType() == PFCandidate::eEGammaHF      || 
	pPF->PFType() == PFCandidate::eNeutralHadron || 
	pPF->PFType() == PFCandidate::eHadronHF      ))
      {lVec -= pPF->Mom(); lSumEt += pPF->Pt();}
    if(pTrack        ==  0                           ) continue;
    if(	 !((pPF->HasTrackerTrk() && (fabs(pPF->TrackerTrk()->DzCorrected(*iVertex))<iDZCut)) ||
	   (pPF->HasGsfTrk()     && (fabs(pPF->GsfTrk()->DzCorrected(*iVertex))    <iDZCut)))) continue; 
    lVec     -= pPF->Mom();
    lSumEt   += pPF->Pt();
  }
  for(UInt_t i0 = 0; i0 < iJets->GetEntries(); i0++) {
    const PFJet *pJet = iJets->At(i0);
    if(fJetIDMVA->correctedPt(pJet,iJetCorrector,iPileupEnergyDensity) < fJetIDMVA->fJetPtMin
       && pJet->TrackCountingHighEffBJetTagsDisc() == -100)              continue; //This line is a bug in the Met training//
    if(!JetTools::passPFLooseId(pJet))                                   continue;
    if(fJetIDMVA->pass(pJet,iVertex,iJetCorrector,iPileupEnergyDensity)) continue;
    if(!filter(pJet,iPhi1,iEta1,iPhi2,iEta2))                            continue; //Quick cleaning==> if not done already
    addNeut(pJet,lVec,lSumEt,iJetCorrector,iPileupEnergyDensity,-1);
  }
  Met lMet(lVec.Px(),lVec.Py());
  lMet.SetSumEt(lSumEt);
  return lMet;
}
//--------------------------------------------------------------------------------------------------
//Corrected jets
Met RecoilTools::PUCMet( const PFJetCol       *iJets,const PFCandidateCol *iCands,const Vertex *iVertex,
			 Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2,
			 Double_t iDZCut) { 

  FourVectorM lVec        (0,0,0,0); double lSumEt          = 0; 
  for(UInt_t i0 = 0; i0 < iCands->GetEntries(); i0++) { 
    const PFCandidate *pPF = iCands->At(i0);
    const Track* pTrack = pPF->TrackerTrk();
    if(pPF->GsfTrk()) pTrack = pPF->GsfTrk();
    if(pTrack == 0                                   && 
       (pPF->PFType() == PFCandidate::eGamma         || 
	pPF->PFType() == PFCandidate::eEGammaHF      || 
	pPF->PFType() == PFCandidate::eNeutralHadron || 
	pPF->PFType() == PFCandidate::eHadronHF      ))
      {lVec -= pPF->Mom(); lSumEt += pPF->Pt();}
    if(pTrack        ==  0                           ) continue;
    if(	 !((pPF->HasTrackerTrk() && (fabs(pPF->TrackerTrk()->DzCorrected(*iVertex))<iDZCut)) ||
	   (pPF->HasGsfTrk()     && (fabs(pPF->GsfTrk()->DzCorrected(*iVertex))    <iDZCut)))) continue; 
    lVec     -= pPF->Mom();
    lSumEt   += pPF->Pt();
  }
  for(UInt_t i0 = 0; i0 < iJets->GetEntries(); i0++) {
    const PFJet *pJet = iJets->At(i0);
    if(pJet->Pt() < fJetIDMVA->fJetPtMin
       && pJet->TrackCountingHighEffBJetTagsDisc() == -100)              continue; //This line is a bug in the Met training//
    if(!JetTools::passPFLooseId(pJet))                                   continue;
    if(fJetIDMVA->pass(pJet,iVertex))                                    continue;
    if(!filter(pJet,iPhi1,iEta1,iPhi2,iEta2))                            continue; //Quick cleaning==> if not done already
    addNeut(pJet,lVec,lSumEt,-1);
  }
  Met lMet(lVec.Px(),lVec.Py());
  lMet.SetSumEt(lSumEt);
  return lMet;
}
//----> This MET is a bug need to fix it
//--------------------------------------------------------------------------------------------------
Met RecoilTools::PUCRecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,
			   const PFJetCol       *iJets,FactorizedJetCorrector *iJetCorrector,
			   const PileupEnergyDensityCol *iPileupEnergyDensity,
			   const PFCandidateCol *iCands,const Vertex *iVertex,
			   Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2,
			   Double_t iDZCut) { 
  Met lPUCMet = PUCMet(iJets,iJetCorrector,iPileupEnergyDensity,iCands,iVertex,iPhi1,iEta1,iPhi2,iEta2,iDZCut);
  lPUCMet.SetMex  (lPUCMet.Mex()+iVisPt*cos(iVisPhi));  
  lPUCMet.SetMey  (lPUCMet.Mey()+iVisPt*sin(iVisPhi));
  lPUCMet.SetSumEt(lPUCMet.SumEt()-iVisSumEt);
  return lPUCMet;
}
//--------------------------------------------------------------------------------------------------
//Corrected Jets
Met RecoilTools::PUCRecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,
			   const PFJetCol       *iJets,
			   const PFCandidateCol *iCands,const Vertex *iVertex,
			   Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2,
			   Double_t iDZCut) { 
  Met lPUCMet = PUCMet(iJets,iCands,iVertex,iPhi1,iEta1,iPhi2,iEta2,iDZCut);
  lPUCMet.SetMex  (lPUCMet.Mex()+iVisPt*cos(iVisPhi));  
  lPUCMet.SetMey  (lPUCMet.Mey()+iVisPt*sin(iVisPhi));
  lPUCMet.SetSumEt(lPUCMet.SumEt()-iVisSumEt);
  return lPUCMet;
}
//--------------------------------------------------------------------------------------------------
Met RecoilTools::PUMet( const PFJetCol       *iJets,FactorizedJetCorrector *iJetCorrector,
			const PileupEnergyDensityCol *iPileupEnergyDensity,
			const PFCandidateCol *iCands,const Vertex *iVertex,
			Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2,
			Double_t iDZCut) { 

  FourVectorM lVec        (0,0,0,0); double lSumEt          = 0; 
  for(UInt_t i0 = 0; i0 < iCands->GetEntries(); i0++) { 
    const PFCandidate *pPF = iCands->At(i0);
    const Track* pTrack = pPF->TrackerTrk();
    if(pPF->GsfTrk()) pTrack = pPF->GsfTrk();
    if(pTrack        ==  0                           ) continue;
    if(	  ((pPF->HasTrackerTrk() && (fabs(pPF->TrackerTrk()->DzCorrected(*iVertex))<iDZCut)) ||
	   (pPF->HasGsfTrk()     && (fabs(pPF->GsfTrk()->DzCorrected(*iVertex))    <iDZCut)))) continue; 
    lVec     -= pPF->Mom();
    lSumEt   += pPF->Pt();
  }
  for(UInt_t i0 = 0; i0 < iJets->GetEntries(); i0++) {
    const PFJet *pJet = iJets->At(i0);
    if(fJetIDMVA->correctedPt(pJet,iJetCorrector,iPileupEnergyDensity) < fJetIDMVA->fJetPtMin
       && pJet->TrackCountingHighEffBJetTagsDisc() == -100)              continue; //This line is a bug in the Met training//
    if(!JetTools::passPFLooseId(pJet))                                   continue;
    if(!filter(pJet,iPhi1,iEta1,iPhi2,iEta2))                            continue; //Quick cleaning
    if(fJetIDMVA->pass(pJet,iVertex,iJetCorrector,iPileupEnergyDensity)) continue;
    addNeut(pJet,lVec,lSumEt,iJetCorrector,iPileupEnergyDensity);
  }
  Met lMet(lVec.Px(),lVec.Py());
  lMet.SetSumEt(lSumEt);
  return lMet;
}
//--------------------------------------------------------------------------------------------------
//Corrected Jets
Met RecoilTools::PUMet( const PFJetCol       *iJets,
			const PFCandidateCol *iCands,const Vertex *iVertex,
			Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2,
			Double_t iDZCut) { 

  FourVectorM lVec        (0,0,0,0); double lSumEt          = 0; 
  for(UInt_t i0 = 0; i0 < iCands->GetEntries(); i0++) { 
    const PFCandidate *pPF = iCands->At(i0);
    const Track* pTrack = pPF->TrackerTrk();
    if(pPF->GsfTrk()) pTrack = pPF->GsfTrk();
    if(pTrack        ==  0                           ) continue;
    if(	  ((pPF->HasTrackerTrk() && (fabs(pPF->TrackerTrk()->DzCorrected(*iVertex))<iDZCut)) ||
	   (pPF->HasGsfTrk()     && (fabs(pPF->GsfTrk()->DzCorrected(*iVertex))    <iDZCut)))) continue; 
    lVec     -= pPF->Mom();
    lSumEt   += pPF->Pt();
  }
  for(UInt_t i0 = 0; i0 < iJets->GetEntries(); i0++) {
    const PFJet *pJet = iJets->At(i0);
    if(pJet->Pt() < fJetIDMVA->fJetPtMin
       && pJet->TrackCountingHighEffBJetTagsDisc() == -100)              continue; //This line is a bug in the Met training//
    if(!JetTools::passPFLooseId(pJet))                                   continue;
    if(!filter(pJet,iPhi1,iEta1,iPhi2,iEta2))                            continue; //Quick cleaning
    if(fJetIDMVA->pass(pJet,iVertex))                                    continue;
    addNeut(pJet,lVec,lSumEt);
  }
  Met lMet(lVec.Px(),lVec.Py());
  lMet.SetSumEt(lSumEt);
  return lMet;
}
//--------------------------------------------------------------------------------------------------
Met RecoilTools::PURecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,
			  const PFJetCol       *iJets,FactorizedJetCorrector *iJetCorrector,
			  const PileupEnergyDensityCol *iPileupEnergyDensity,
			  const PFCandidateCol *iCands,const Vertex *iVertex,
			  Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2,
			  Double_t iDZCut) { 
  Met lPUMet = PUMet(iJets,iJetCorrector,iPileupEnergyDensity,iCands,iVertex,iPhi1,iEta1,iPhi2,iEta2,iDZCut);
  lPUMet.SetMex  (lPUMet.Mex()+iVisPt*cos(iVisPhi));  
  lPUMet.SetMey  (lPUMet.Mey()+iVisPt*sin(iVisPhi));
  lPUMet.SetSumEt(lPUMet.SumEt()-iVisSumEt);
  return lPUMet;
}
//--------------------------------------------------------------------------------------------------
//Corrected Jets
Met RecoilTools::PURecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,
			  const PFJetCol       *iJets,
			  const PFCandidateCol *iCands,const Vertex *iVertex,
			  Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2,
			  Double_t iDZCut) { 
  Met lPUMet = PUMet(iJets,iCands,iVertex,iPhi1,iEta1,iPhi2,iEta2,iDZCut);
  lPUMet.SetMex  (lPUMet.Mex()+iVisPt*cos(iVisPhi));  
  lPUMet.SetMey  (lPUMet.Mey()+iVisPt*sin(iVisPhi));
  lPUMet.SetSumEt(lPUMet.SumEt()-iVisSumEt);
  return lPUMet;
}
