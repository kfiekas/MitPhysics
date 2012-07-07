// $Id: MetTools.cc,v 1.12 2012/04/28 19:10:01 ceballos Exp $

#include "MitPhysics/Utils/interface/MetTools.h"
#include <TFile.h>

ClassImp(mithep::MetTools)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
MetTools::MetTools(const MuonCol *fMuons, const PFCandidateCol *fPFCandidates, const Vertex *fVertex, 
		   float deltaZCut, float ptCut, float etaCut) {

  float trackNumeratorX  =0, trackNumeratorY  =0;
  float neutralNumeratorX=0, neutralNumeratorY=0;
  float CHSNumeratorX  =0, CHSNumeratorY  =0, NHSNumeratorX  =0, NHSNumeratorY  =0;

  // muons Pt
  for (UInt_t m = 0; m < fMuons->GetEntries(); ++m) {
    trackNumeratorX -= fMuons->At(m)->Px();
    trackNumeratorY -= fMuons->At(m)->Py();
  }

  // PF candidates pT
  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {

    // CHS computation
    if (fPFCandidates->At(i)->HasTrk()){
	CHSNumeratorX -= fPFCandidates->At(i)->Px();
	CHSNumeratorY -= fPFCandidates->At(i)->Py();
    }
    // NHS computation
    if (fPFCandidates->At(i)->HasTrk() ||
       (fPFCandidates->At(i)->AbsEta() < 3.0 && fPFCandidates->At(i)->Pt() > 4.0)){
	NHSNumeratorX -= fPFCandidates->At(i)->Px();
	NHSNumeratorY -= fPFCandidates->At(i)->Py();
    }

    // charged
    if (fPFCandidates->At(i)->HasTrackerTrk()){

      bool isMuonTrack = false;
      for (UInt_t m = 0; m < fMuons->GetEntries(); ++m) {
	if (fMuons->At(m)->TrackerTrk() == fPFCandidates->At(i)->TrackerTrk()) {
	  isMuonTrack = true;
	  break;
	}
      }      
      if (isMuonTrack) continue;

      if (fabs(fPFCandidates->At(i)->TrackerTrk()->DzCorrected(*fVertex)) < deltaZCut) {
	trackNumeratorX -= fPFCandidates->At(i)->Px();
	trackNumeratorY -= fPFCandidates->At(i)->Py();
        fRecoil += fPFCandidates->At(i)->Mom();
        fChargedRecoil += fPFCandidates->At(i)->Mom();
      }
    }

    // neutral    
    if (fPFCandidates->At(i)->PFType()== PFCandidate::eNeutralHadron || fPFCandidates->At(i)->PFType()== PFCandidate::eGamma) {
      if (fPFCandidates->At(i)->Pt() > ptCut and fabs(fPFCandidates->At(i)->Eta()) < etaCut ) {
	neutralNumeratorX -= fPFCandidates->At(i)->Px();
	neutralNumeratorY -= fPFCandidates->At(i)->Py();
        fRecoil += fPFCandidates->At(i)->Mom();
      }
    }
  }

  fCorrectedMet = mithep::Met(trackNumeratorX+neutralNumeratorX, trackNumeratorY+neutralNumeratorY);
  fCorrectedTrackMet = mithep::Met(trackNumeratorX, trackNumeratorY);
  fCHSMet = mithep::Met(CHSNumeratorX, CHSNumeratorY);
  fNHSMet = mithep::Met(NHSNumeratorX, NHSNumeratorY);
}

MetTools::MetTools(const ElectronCol *fElectrons, const PFCandidateCol *fPFCandidates, const Vertex *fVertex, 
		   float deltaZCut, float ptCut, float etaCut) {

  float trackNumeratorX  =0, trackNumeratorY  =0;
  float neutralNumeratorX=0, neutralNumeratorY=0;
  float CHSNumeratorX  =0, CHSNumeratorY  =0, NHSNumeratorX  =0, NHSNumeratorY  =0;

  // electrons Pt
  for (UInt_t m = 0; m < fElectrons->GetEntries(); ++m) {
    trackNumeratorX -= fElectrons->At(m)->Px();
    trackNumeratorY -= fElectrons->At(m)->Py();
  }

  // PF candidates pT
  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {

    // CHS computation
    if (fPFCandidates->At(i)->HasTrk()){
	CHSNumeratorX -= fPFCandidates->At(i)->Px();
	CHSNumeratorY -= fPFCandidates->At(i)->Py();
    }
    // NHS computation
    if (fPFCandidates->At(i)->HasTrk() ||
       (fPFCandidates->At(i)->AbsEta() < 3.0 && fPFCandidates->At(i)->Pt() > 4.0)){
	NHSNumeratorX -= fPFCandidates->At(i)->Px();
	NHSNumeratorY -= fPFCandidates->At(i)->Py();
    }

    // charged
    if (fPFCandidates->At(i)->HasTrackerTrk() || fPFCandidates->At(i)->HasGsfTrk()){
      bool isElectronTrack = false;
      for (UInt_t m = 0; m < fElectrons->GetEntries(); ++m) {
	if ( (fElectrons->At(m)->TrackerTrk() == fPFCandidates->At(i)->TrackerTrk()) or
	     (fElectrons->At(m)->HasGsfTrk() and fElectrons->At(m)->GsfTrk() == fPFCandidates->At(i)->GsfTrk()) ) {
	  isElectronTrack = true;
	  break;
	}
      }      
      if (isElectronTrack) continue;

      if ((fPFCandidates->At(i)->HasTrackerTrk() && fabs(fPFCandidates->At(i)->TrackerTrk()->DzCorrected(*fVertex)) < deltaZCut) ||
          (fPFCandidates->At(i)->HasGsfTrk()     && fabs(fPFCandidates->At(i)->GsfTrk()->DzCorrected(*fVertex)    ) < deltaZCut)) {
	trackNumeratorX -= fPFCandidates->At(i)->Px();
	trackNumeratorY -= fPFCandidates->At(i)->Py();
        fRecoil += fPFCandidates->At(i)->Mom();
        fChargedRecoil += fPFCandidates->At(i)->Mom();
      }
    }

    // neutral    
    if (fPFCandidates->At(i)->PFType()== PFCandidate::eNeutralHadron || fPFCandidates->At(i)->PFType()== PFCandidate::eGamma) {
      if (fPFCandidates->At(i)->Pt() > ptCut and fabs(fPFCandidates->At(i)->Eta()) < etaCut ) {
	neutralNumeratorX -= fPFCandidates->At(i)->Px();
	neutralNumeratorY -= fPFCandidates->At(i)->Py();
        fRecoil += fPFCandidates->At(i)->Mom();
      }
    }
  }

  fCorrectedMet = mithep::Met(trackNumeratorX+neutralNumeratorX, trackNumeratorY+neutralNumeratorY);
  fCorrectedTrackMet = mithep::Met(trackNumeratorX, trackNumeratorY);
  fCHSMet = mithep::Met(CHSNumeratorX, CHSNumeratorY);
  fNHSMet = mithep::Met(NHSNumeratorX, NHSNumeratorY);
}


//--------------------------------------------------------------------------------------------------
MetTools::MetTools(const MuonCol *fMuons, const PFCandidateCol *fPFCandidates, const PFJetCol *pfJets, const Vertex *fVertex, 
		   float deltaZCut, float ptCut, float etaCut) {

  float trackNumeratorX  =0, trackNumeratorY  =0;
  float neutralNumeratorX=0, neutralNumeratorY=0;
  float CHSNumeratorX  =0, CHSNumeratorY  =0, NHSNumeratorX  =0, NHSNumeratorY  =0;

  // muons Pt
  for (UInt_t m = 0; m < fMuons->GetEntries(); ++m) {
    trackNumeratorX -= fMuons->At(m)->Px();
    trackNumeratorY -= fMuons->At(m)->Py();
  }
  
  // jets Pt
  for (UInt_t j = 0; j < pfJets->GetEntries(); ++j) {
    trackNumeratorX -= pfJets->At(j)->Px();
    trackNumeratorY -= pfJets->At(j)->Py();
  }
  
  // PF candidates pT
  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {
    
    // CHS computation
    if (fPFCandidates->At(i)->HasTrk()){
	CHSNumeratorX -= fPFCandidates->At(i)->Px();
	CHSNumeratorY -= fPFCandidates->At(i)->Py();
    }
    // NHS computation
    if (fPFCandidates->At(i)->HasTrk() ||
       (fPFCandidates->At(i)->AbsEta() < 3.0 && fPFCandidates->At(i)->Pt() > 4.0)){
	NHSNumeratorX -= fPFCandidates->At(i)->Px();
	NHSNumeratorY -= fPFCandidates->At(i)->Py();
    }

    // jets
    bool inTheJet = false;
    for (UInt_t j = 0; j < pfJets->GetEntries(); ++j) {
      for (UInt_t c=0;c<pfJets->At(j)->NPFCands();++c){
	if (pfJets->At(j)->PFCand(c) == fPFCandidates->At(i)) {
	  inTheJet=true;
	  break;
	}
      }
      if (inTheJet) break; 
    }
    if (inTheJet) continue;
    
    // charged
    if (fPFCandidates->At(i)->HasTrackerTrk()){

      bool isMuonTrack = false;
      for (UInt_t m = 0; m < fMuons->GetEntries(); ++m) {
	if (fMuons->At(m)->TrackerTrk() == fPFCandidates->At(i)->TrackerTrk()) {
	  isMuonTrack = true;
	  break;
	}
      }      
      if (isMuonTrack) continue;
      
      if (fabs(fPFCandidates->At(i)->TrackerTrk()->DzCorrected(*fVertex)) < deltaZCut) {
	trackNumeratorX -= fPFCandidates->At(i)->Px();
	trackNumeratorY -= fPFCandidates->At(i)->Py();
        fRecoil += fPFCandidates->At(i)->Mom();
        fChargedRecoil += fPFCandidates->At(i)->Mom();
      }
    }

    // neutral    
    if (fPFCandidates->At(i)->PFType()== PFCandidate::eNeutralHadron || fPFCandidates->At(i)->PFType()== PFCandidate::eGamma) {
      if (fPFCandidates->At(i)->Pt() > ptCut and fabs(fPFCandidates->At(i)->Eta()) < etaCut ) {
	neutralNumeratorX -= fPFCandidates->At(i)->Px();
	neutralNumeratorY -= fPFCandidates->At(i)->Py();
        fRecoil += fPFCandidates->At(i)->Mom();
      }
    }
  }

  fCorrectedMet = mithep::Met(trackNumeratorX+neutralNumeratorX, trackNumeratorY+neutralNumeratorY);
  fCorrectedTrackMet = mithep::Met(trackNumeratorX, trackNumeratorY);
  fCHSMet = mithep::Met(CHSNumeratorX, CHSNumeratorY);
  fNHSMet = mithep::Met(NHSNumeratorX, NHSNumeratorY);
}


MetTools::MetTools(const ElectronCol *fElectrons, const PFCandidateCol *fPFCandidates, const PFJetCol *pfJets, const Vertex *fVertex, 
		   float deltaZCut, float ptCut, float etaCut) {

  float trackNumeratorX  =0, trackNumeratorY  =0;
  float neutralNumeratorX=0, neutralNumeratorY=0;
  float CHSNumeratorX  =0, CHSNumeratorY  =0, NHSNumeratorX  =0, NHSNumeratorY  =0;

  // electrons Pt
  for (UInt_t m = 0; m < fElectrons->GetEntries(); ++m) {
    trackNumeratorX -= fElectrons->At(m)->Px();
    trackNumeratorY -= fElectrons->At(m)->Py();
  }

  // jets Pt
  for (UInt_t j = 0; j < pfJets->GetEntries(); ++j) {
    trackNumeratorX -= pfJets->At(j)->Px();
    trackNumeratorY -= pfJets->At(j)->Py();
  }

  // PF candidates pT
  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {

    // CHS computation
    if (fPFCandidates->At(i)->HasTrk()){
	CHSNumeratorX -= fPFCandidates->At(i)->Px();
	CHSNumeratorY -= fPFCandidates->At(i)->Py();
    }
    // NHS computation
    if (fPFCandidates->At(i)->HasTrk() ||
       (fPFCandidates->At(i)->AbsEta() < 3.0 && fPFCandidates->At(i)->Pt() > 4.0)){
	NHSNumeratorX -= fPFCandidates->At(i)->Px();
	NHSNumeratorY -= fPFCandidates->At(i)->Py();
    }

    // jets
    bool inTheJet = false;
    for (UInt_t j = 0; j < pfJets->GetEntries(); ++j) {
      for (UInt_t c=0;c<pfJets->At(j)->NPFCands();++c){
	if (pfJets->At(j)->PFCand(c) == fPFCandidates->At(i)) {
	  inTheJet=true;
	  break;
	}
      }
      if (inTheJet) break; 
    }
    if (inTheJet) continue;

    // charged
    if (fPFCandidates->At(i)->HasTrackerTrk() || fPFCandidates->At(i)->HasGsfTrk()){
      bool isElectronTrack = false;
      for (UInt_t m = 0; m < fElectrons->GetEntries(); ++m) {
	if ( (fElectrons->At(m)->HasTrackerTrk() and fElectrons->At(m)->TrackerTrk() == fPFCandidates->At(i)->TrackerTrk()) or
	     (fElectrons->At(m)->HasGsfTrk() and fElectrons->At(m)->GsfTrk() == fPFCandidates->At(i)->GsfTrk()) ) {
	  isElectronTrack = true;
	  break;
	}
      }      
      if (isElectronTrack) continue;

      if ((fPFCandidates->At(i)->HasTrackerTrk() && fabs(fPFCandidates->At(i)->TrackerTrk()->DzCorrected(*fVertex)) < deltaZCut) ||
          (fPFCandidates->At(i)->HasGsfTrk()     && fabs(fPFCandidates->At(i)->GsfTrk()->DzCorrected(*fVertex)    ) < deltaZCut)) {
	trackNumeratorX -= fPFCandidates->At(i)->Px();
	trackNumeratorY -= fPFCandidates->At(i)->Py();
        fRecoil += fPFCandidates->At(i)->Mom();
        fChargedRecoil += fPFCandidates->At(i)->Mom();
      }
    }

    // neutral    
    if (fPFCandidates->At(i)->PFType()== PFCandidate::eNeutralHadron || fPFCandidates->At(i)->PFType()== PFCandidate::eGamma) {
      if (fPFCandidates->At(i)->Pt() > ptCut and fabs(fPFCandidates->At(i)->Eta()) < etaCut ) {
	neutralNumeratorX -= fPFCandidates->At(i)->Px();
	neutralNumeratorY -= fPFCandidates->At(i)->Py();
        fRecoil += fPFCandidates->At(i)->Mom();
      }
    }
  }

  fCorrectedMet = mithep::Met(trackNumeratorX+neutralNumeratorX, trackNumeratorY+neutralNumeratorY);
  fCorrectedTrackMet = mithep::Met(trackNumeratorX, trackNumeratorY);
  fCHSMet = mithep::Met(CHSNumeratorX, CHSNumeratorY);
  fNHSMet = mithep::Met(NHSNumeratorX, NHSNumeratorY);
}


void MetTools::AddToCorrectedTrackMet( const Particle *p, bool debug ) {
  float MetX = fCorrectedTrackMet.Mex();
  float MetY = fCorrectedTrackMet.Mey();

  if (debug) std::cout << "AddToCorrectedTrackMet:\n";
  if (debug) std::cout << "Before: " << MetX << " " << MetY << std::endl;

  MetX -= p->Px();
  MetY -= p->Py();

  if (debug) std::cout << "Add : " << p->Px() << " " << p->Py() << std::endl;
  if (debug) std::cout << "After : " << MetX << " " << MetY << std::endl;

  fCorrectedTrackMet.SetMex(MetX);
  fCorrectedTrackMet.SetMey(MetY);

  return;
}

void MetTools::AddToCorrectedMet( const Particle *p) {
  float MetX=fCorrectedMet.Mex();
  float MetY=fCorrectedMet.Mey();
  
  MetX -= p->Px();
  MetY -= p->Py();

  fCorrectedMet.SetMex(MetX);
  fCorrectedMet.SetMey(MetY);

  return;
}

void MetTools::AddToRecoil( const Particle *p) {

  if (p->Charge() != 0) {
    fChargedRecoil += p->Mom();
  }
  fRecoil += p->Mom();

  return;
}



void MetTools::RemoveParticleInIsoConeFromTrackMet( const Particle *p, const PFCandidateCol *fPFCandidates, const Vertex *fVertex, float deltaZCut , float deltaR, bool debug ) {
  float MetX = fCorrectedTrackMet.Mex();
  float MetY = fCorrectedTrackMet.Mey();
  
  if (debug) std::cout << "RemoveParticleInIsoConeFromTrackMet:\n";
  if (debug) std::cout << "Before: " << MetX << " " << MetY << std::endl;

  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {
    //charged
    if (fPFCandidates->At(i)->HasTrackerTrk() || fPFCandidates->At(i)->HasGsfTrk()) {
      //passes dZ cut
      if ((fPFCandidates->At(i)->HasTrackerTrk() && fabs(fPFCandidates->At(i)->TrackerTrk()->DzCorrected(*fVertex)) < deltaZCut) ||
          (fPFCandidates->At(i)->HasGsfTrk()     && fabs(fPFCandidates->At(i)->GsfTrk()->DzCorrected(*fVertex)    ) < deltaZCut)) {        
        //inside cone
        if (MathUtils::DeltaR(fPFCandidates->At(i)->Mom(), p->Mom()) < deltaR ) {
          MetX += fPFCandidates->At(i)->Px();
          MetY += fPFCandidates->At(i)->Py();
          if (debug) std::cout << "Subtract : " << fPFCandidates->At(i)->Px() << " " << fPFCandidates->At(i)->Py() << " : " << std::endl;
        }
      }
    }
  }

  if (debug) std::cout << "After : " << MetX << " " << MetY << std::endl;

  fCorrectedTrackMet.SetMex(MetX);
  fCorrectedTrackMet.SetMey(MetY);

  return;
}


void MetTools::RemoveParticleInIsoConeFromCorrectedMet( const Particle *p, const PFCandidateCol *fPFCandidates, const Vertex *fVertex, float deltaZCut, float ptCut, float etaCut , float deltaR) {
  float MetX = fCorrectedMet.Mex();
  float MetY = fCorrectedMet.Mey();
  
  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {

    //inside cone
    if (MathUtils::DeltaR(fPFCandidates->At(i)->Mom(), p->Mom()) < deltaR ) {

      //charged
      if (fPFCandidates->At(i)->HasTrackerTrk() || fPFCandidates->At(i)->HasGsfTrk()) {
        //passes dZ cut
        if ((fPFCandidates->At(i)->HasTrackerTrk() && fabs(fPFCandidates->At(i)->TrackerTrk()->DzCorrected(*fVertex)) < deltaZCut) ||
            (fPFCandidates->At(i)->HasGsfTrk()     && fabs(fPFCandidates->At(i)->GsfTrk()->DzCorrected(*fVertex)    ) < deltaZCut)) {                  
          MetX += p->Px();
          MetY += p->Py();
        }
      }

      //neutrals
      if (fPFCandidates->At(i)->PFType()== PFCandidate::eNeutralHadron || fPFCandidates->At(i)->PFType()== PFCandidate::eGamma) {
        if (fPFCandidates->At(i)->Pt() > ptCut and fabs(fPFCandidates->At(i)->Eta()) < etaCut ) {
          MetX += fPFCandidates->At(i)->Px();
          MetY += fPFCandidates->At(i)->Py();
        }
      }
    }

  }

  fCorrectedMet.SetMex(MetX);
  fCorrectedMet.SetMey(MetY);

  return;
}


void MetTools::RemoveParticleInIsoConeFromRecoil( const Particle *p, const PFCandidateCol *fPFCandidates, const Vertex *fVertex, float deltaZCut, float ptCut, float etaCut , float deltaR) {
  
  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {

    //inside cone
    if (MathUtils::DeltaR(fPFCandidates->At(i)->Mom(), p->Mom()) < deltaR ) {

      //charged
      if (fPFCandidates->At(i)->HasTrackerTrk() || fPFCandidates->At(i)->HasGsfTrk()) {
        //passes dZ cut
        if ((fPFCandidates->At(i)->HasTrackerTrk() && fabs(fPFCandidates->At(i)->TrackerTrk()->DzCorrected(*fVertex)) < deltaZCut) ||
            (fPFCandidates->At(i)->HasGsfTrk()     && fabs(fPFCandidates->At(i)->GsfTrk()->DzCorrected(*fVertex)    ) < deltaZCut)) {                  
          fChargedRecoil -= fPFCandidates->At(i)->Mom();
          fRecoil -= fPFCandidates->At(i)->Mom();
        }
      }

      //neutrals
      if (fPFCandidates->At(i)->PFType()== PFCandidate::eNeutralHadron || fPFCandidates->At(i)->PFType()== PFCandidate::eGamma) {
        if (fPFCandidates->At(i)->Pt() > ptCut and fabs(fPFCandidates->At(i)->Eta()) < etaCut ) {          
          fRecoil -= fPFCandidates->At(i)->Mom();
        }
      }
    }
  }

  return;
}


MetTools::MetTools(const MuonCol *fMuons, const ElectronCol *fElectrons, const PFCandidateCol *fPFCandidates, 
                   const Vertex *fVertex, float deltaZCut, float ptCut, float etaCut, float intRadius, 
		   const GenericParticle *genP) {

  float trackNumeratorX  =0, trackNumeratorY  =0;
  float neutralNumeratorX=0, neutralNumeratorY=0;
  float CHSNumeratorX  =0, CHSNumeratorY  =0, NHSNumeratorX  =0, NHSNumeratorY  =0;

  float trackNumeratorGenPX  =0, trackNumeratorGenPY  =0;
  // added to the track quantities
  if (genP) {
    trackNumeratorGenPX -= genP->Px();
    trackNumeratorGenPY -= genP->Py();
  }

  // muons Pt
  for (UInt_t m = 0; m < fMuons->GetEntries(); ++m) {
    trackNumeratorX -= fMuons->At(m)->Px();
    trackNumeratorY -= fMuons->At(m)->Py();
  }

  // electrons Pt
  for (UInt_t m = 0; m < fElectrons->GetEntries(); ++m) {
    trackNumeratorX -= fElectrons->At(m)->Px();
    trackNumeratorY -= fElectrons->At(m)->Py();
  }

  // PF candidates pT
  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {

    // CHS computation
    if (fPFCandidates->At(i)->HasTrk()){
	CHSNumeratorX -= fPFCandidates->At(i)->Px();
	CHSNumeratorY -= fPFCandidates->At(i)->Py();
    }
    // NHS computation
    if (fPFCandidates->At(i)->HasTrk() ||
       (fPFCandidates->At(i)->AbsEta() < 3.0 && fPFCandidates->At(i)->Pt() > 4.0)){
	NHSNumeratorX -= fPFCandidates->At(i)->Px();
	NHSNumeratorY -= fPFCandidates->At(i)->Py();
    }

    // charged
    if (fPFCandidates->At(i)->HasTrackerTrk() || fPFCandidates->At(i)->HasGsfTrk()){
      bool isMuonTrack = false;
      for (UInt_t m = 0; m < fMuons->GetEntries(); ++m) {
	if (fMuons->At(m)->TrackerTrk() == fPFCandidates->At(i)->TrackerTrk()) {
	  isMuonTrack = true;
	  break;
	}
	if (intRadius > 0.0 && MathUtils::DeltaR(fPFCandidates->At(i)->Mom(), fMuons->At(m)->Mom()) < intRadius) {
	  isMuonTrack = true;
	  break;
	}
      }    
      if (isMuonTrack) continue;
      
      bool isElectronTrack = false;
      for (UInt_t m = 0; m < fElectrons->GetEntries(); ++m) {
	if ( (fElectrons->At(m)->HasTrackerTrk() and fElectrons->At(m)->TrackerTrk() == fPFCandidates->At(i)->TrackerTrk()) or
	     (fElectrons->At(m)->HasGsfTrk() and fElectrons->At(m)->GsfTrk() == fPFCandidates->At(i)->GsfTrk()) ) {
	  isElectronTrack = true;
	  break;
	}
	if (intRadius > 0.0 && MathUtils::DeltaR(fPFCandidates->At(i)->Mom(), fElectrons->At(m)->Mom()) < intRadius) {
	  isElectronTrack = true;
	  break;
	}
      }
      if (isElectronTrack) continue;

      if ((fPFCandidates->At(i)->HasTrackerTrk() && fabs(fPFCandidates->At(i)->TrackerTrk()->DzCorrected(*fVertex)) < deltaZCut) ||
          (fPFCandidates->At(i)->HasGsfTrk()     && fabs(fPFCandidates->At(i)->GsfTrk()->DzCorrected(*fVertex)    ) < deltaZCut)) {
	trackNumeratorX -= fPFCandidates->At(i)->Px();
	trackNumeratorY -= fPFCandidates->At(i)->Py();
        fRecoil += fPFCandidates->At(i)->Mom();
        fChargedRecoil += fPFCandidates->At(i)->Mom();
      }
    }

    // neutral    
    if (fPFCandidates->At(i)->PFType()== PFCandidate::eNeutralHadron || fPFCandidates->At(i)->PFType()== PFCandidate::eGamma) {
      if (fPFCandidates->At(i)->Pt() > ptCut and fabs(fPFCandidates->At(i)->Eta()) < etaCut ) {
	neutralNumeratorX -= fPFCandidates->At(i)->Px();
	neutralNumeratorY -= fPFCandidates->At(i)->Py();
        fRecoil += fPFCandidates->At(i)->Mom();
      }
    }
  }
  fCorrectedMet = mithep::Met(trackNumeratorX+neutralNumeratorX, trackNumeratorY+neutralNumeratorY);
  fCorrectedTrackMet = mithep::Met(trackNumeratorX+trackNumeratorGenPX, trackNumeratorY+trackNumeratorGenPY);
  fCHSMet = mithep::Met(CHSNumeratorX+trackNumeratorGenPX, CHSNumeratorY+trackNumeratorGenPY);
  fNHSMet = mithep::Met(NHSNumeratorX, NHSNumeratorY);
}


MetTools::MetTools(const MuonCol *fMuons, const ElectronCol *fElectrons, const PFCandidateCol *fPFCandidates, 
		   const PFJetCol *pfJets,
                   const Vertex *fVertex, float deltaZCut, float ptCut, float etaCut, float intRadius, 
		   const GenericParticle *genP) {

  float trackNumeratorX  =0, trackNumeratorY  =0;
  float neutralNumeratorX=0, neutralNumeratorY=0;
  float CHSNumeratorX  =0, CHSNumeratorY  =0, NHSNumeratorX  =0, NHSNumeratorY  =0;

  float trackNumeratorGenPX  =0, trackNumeratorGenPY  =0;
  // added to the track quantities
  if (genP) {
    trackNumeratorGenPX -= genP->Px();
    trackNumeratorGenPY -= genP->Py();
  }

  // muons Pt
  for (UInt_t m = 0; m < fMuons->GetEntries(); ++m) {
    trackNumeratorX -= fMuons->At(m)->Px();
    trackNumeratorY -= fMuons->At(m)->Py();
  }

  // electrons Pt
  for (UInt_t m = 0; m < fElectrons->GetEntries(); ++m) {
    trackNumeratorX -= fElectrons->At(m)->Px();
    trackNumeratorY -= fElectrons->At(m)->Py();
  }

  // jets Pt
  for (UInt_t j = 0; j < pfJets->GetEntries(); ++j) {
    trackNumeratorX -= pfJets->At(j)->Px();
    trackNumeratorY -= pfJets->At(j)->Py();
  }

  // PF candidates pT
  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {

    // CHS computation
    if (fPFCandidates->At(i)->HasTrk()){
	CHSNumeratorX -= fPFCandidates->At(i)->Px();
	CHSNumeratorY -= fPFCandidates->At(i)->Py();
    }
    // NHS computation
    if (fPFCandidates->At(i)->HasTrk() ||
       (fPFCandidates->At(i)->AbsEta() < 3.0 && fPFCandidates->At(i)->Pt() > 4.0)){
	NHSNumeratorX -= fPFCandidates->At(i)->Px();
	NHSNumeratorY -= fPFCandidates->At(i)->Py();
    }

    // jets
    bool inTheJet = false;
    for (UInt_t j = 0; j < pfJets->GetEntries(); ++j) {
      for (UInt_t c=0;c<pfJets->At(j)->NPFCands();++c){
	if (pfJets->At(j)->PFCand(c) == fPFCandidates->At(i)) {
	  inTheJet=true;
	  break;
	}
      }
      if (inTheJet) break; 
    }
    if (inTheJet) continue;

    // charged
    if (fPFCandidates->At(i)->HasTrackerTrk() || fPFCandidates->At(i)->HasGsfTrk()){
      bool isMuonTrack = false;
      for (UInt_t m = 0; m < fMuons->GetEntries(); ++m) {
	if (fMuons->At(m)->TrackerTrk() == fPFCandidates->At(i)->TrackerTrk()) {
	  isMuonTrack = true;
	  break;
	}
	if (intRadius > 0.0 && MathUtils::DeltaR(fPFCandidates->At(i)->Mom(), fMuons->At(m)->Mom()) < intRadius) {
	  isMuonTrack = true;
	  break;
	}
      }      
      if (isMuonTrack) continue;
      
      bool isElectronTrack = false;
      for (UInt_t m = 0; m < fElectrons->GetEntries(); ++m) {
	if ( (fElectrons->At(m)->HasTrackerTrk() and fElectrons->At(m)->TrackerTrk() == fPFCandidates->At(i)->TrackerTrk()) or
	     (fElectrons->At(m)->HasGsfTrk() and fElectrons->At(m)->GsfTrk() == fPFCandidates->At(i)->GsfTrk()) ) {
	  isElectronTrack = true;
	  break;
	}
	if (intRadius > 0.0 && MathUtils::DeltaR(fPFCandidates->At(i)->Mom(), fElectrons->At(m)->Mom()) < intRadius) {
	  isElectronTrack = true;
	  break;
	}
      }      
      if (isElectronTrack) continue;

      if ((fPFCandidates->At(i)->HasTrackerTrk() && fabs(fPFCandidates->At(i)->TrackerTrk()->DzCorrected(*fVertex)) < deltaZCut) ||
          (fPFCandidates->At(i)->HasGsfTrk()     && fabs(fPFCandidates->At(i)->GsfTrk()->DzCorrected(*fVertex)    ) < deltaZCut)) {
	trackNumeratorX -= fPFCandidates->At(i)->Px();
	trackNumeratorY -= fPFCandidates->At(i)->Py();
        fRecoil += fPFCandidates->At(i)->Mom();
        fChargedRecoil += fPFCandidates->At(i)->Mom();
      }
    }

    // neutral    
    if (fPFCandidates->At(i)->PFType()== PFCandidate::eNeutralHadron || fPFCandidates->At(i)->PFType()== PFCandidate::eGamma) {
      if (fPFCandidates->At(i)->Pt() > ptCut and fabs(fPFCandidates->At(i)->Eta()) < etaCut ) {
	neutralNumeratorX -= fPFCandidates->At(i)->Px();
	neutralNumeratorY -= fPFCandidates->At(i)->Py();
        fRecoil += fPFCandidates->At(i)->Mom();
      }
    }
  }
  fCorrectedMet = mithep::Met(trackNumeratorX+neutralNumeratorX, trackNumeratorY+neutralNumeratorY);
  fCorrectedTrackMet = mithep::Met(trackNumeratorX+trackNumeratorGenPX, trackNumeratorY+trackNumeratorGenPY);
  fCHSMet = mithep::Met(CHSNumeratorX+trackNumeratorGenPX, CHSNumeratorY+trackNumeratorGenPY);
  fNHSMet = mithep::Met(NHSNumeratorX, NHSNumeratorY);
}


Met MetTools::GetMinimumMet(const Met *UncorrectedMet) {

  return (fCorrectedMet.Pt() < UncorrectedMet->Pt()) ?  fCorrectedMet : *UncorrectedMet; 
}

Met MetTools::GetMinimumTrackMet(const Met *UncorrectedMet) {

  return (fCorrectedTrackMet.Pt() < UncorrectedMet->Pt()) ?  fCorrectedTrackMet : *UncorrectedMet; 
}

