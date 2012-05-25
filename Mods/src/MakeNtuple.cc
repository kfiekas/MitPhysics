//$Id:MakeNtuple.cc,v 1.1 2011/07/12 Mingming Yang 
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TRandom3.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/GeneratorTools.h"//tools to get the generator level information
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitPhysics/Mods/interface/MakeNtuple.h"

using namespace mithep;
using namespace TMath;

ClassImp(mithep::MakeNtuple)
  
//---------------------------------------------------------------------------
  MakeNtuple::MakeNtuple(const char *name, const char *title):
    BaseMod (name,title),
    //------------------------------- 
   fTrigObjsName       (Names::gkHltObjBrn),
    //-------------------------------
    fPhotonName         (Names::gkPhotonBrn),
    fElectronName       ("Electrons"),
    fTrackBranchName    (Names::gkTrackBrn),
    fPVName             (Names::gkPVBeamSpotBrn),
    fPileUpDenName      (Names::gkPileupEnergyDensityBrn),
    fPhotonsFromBranch  (kTRUE),
    //-------------------------------
    fIsData             (kTRUE),
    //----scale removal--------------
    fOverlapCut         (-1.0),
    //----MC Info--------------------
    fMcEventInfoName    (Names::gkMCEvtInfoBrn),
    fMcParticleName     (Names::gkMCPartBrn),
    fPileupInfoName     (Names::gkPileupInfoBrn),
    //----PF-------------------------
    fPFCandName        (Names::gkPFCandidatesBrn),
    //-------------------------------
    fPhotons            (0),
    fElectrons          (0),
    fMcParticles        (0),  
    fMcEventInfo        (0),
    fPileupInfos        (0),
    fPFCands            (0),
    //-------------------------------
    fPhotonPtMin       (25.0),
    fApplyElectronVeto (kTRUE),
    finvertElectronVeto(kFALSE),
    //-------------------------------
    fBeamSpotName(Names::gkBeamSpotBrn),
    fConversionName(Names::gkMvfConversionBrn),
    fBeamSpot(0),
    fConversions(0),
    fIsSignal(kTRUE)
{
  //constructor
}
//----------------------------------------------------------------------------
void MakeNtuple::Process()
{
  // count the events we have processed
  IncNEventsProcessed();
  
  // access the MC Information if needed
  if (!fIsData) {
    LoadBranch(fMcEventInfoName);
    LoadBranch(fMcParticleName);
    LoadBranch(fPileupInfoName);
    fprocessid = fMcEventInfo->ProcessId();
    for (UInt_t i=0; i<fPileupInfos->GetEntries(); ++i) {
      const PileupInfo *puinfo = fPileupInfos->At(i);
      if (puinfo->GetBunchCrossing()==0)  fNVertexesGenPile= puinfo->GetPU_NumInteractions();
    } 
    if (fOverlapCut > 0. && fMcEventInfo->Scale() > fOverlapCut) {
      MDB(kModules,1)
	printf(" Reject event with scale %f  above cut: %f\n",fMcEventInfo->Scale(),fOverlapCut);
      return;
    }
    if( fMcParticles  != NULL && fMcParticles->GetEntries()>0){
      printf("fMcParticles exists!\n");
    }
  }
  
  //load branch
  LoadEventObject(fPhotonName,fPhotons);
  LoadEventObject(fElectronName,fElectrons);
  LoadEventObject(fTrackBranchName,fTracks);
  LoadEventObject(fPVName,fPV);
  LoadEventObject(fPileUpDenName,fPileUpDen);
  LoadEventObject(fPFCandName,         fPFCands);
  LoadBranch(fBeamSpotName);
  LoadBranch(fConversionName);

  Float_t _pth     = -100.;
  Float_t _decayZ  = -100.;
  Float_t _genmass = -100.; 
  
  //select photon  
  if(!fIsData){
    FindHiggsPtAndZ(_pth, _decayZ, _genmass);
    Float_t EventNum=GetEventHeader()->EvtNum();
    if(fIsSignal){
      int PhotonColSize=fPhotons->GetEntries();
      if(PhotonColSize>=2){
	Float_t ptmax1=0;
	Float_t ptmax2=0;
	int imax1=0;
	int imax2=0;
	for(int i=0;i<PhotonColSize;i++){
	  if((fPhotons->At(i)->Pt())>ptmax1){
	    ptmax1=fPhotons->At(i)->Pt();
	    imax1=i;
	  }
	}
	for (int i=0;i<PhotonColSize;i++) {
	  if(((fPhotons->At(i)->Pt())>ptmax2)&&((fPhotons->At(i)->Pt()<ptmax1))){
	    ptmax2=fPhotons->At(i)->Pt();
	    imax2=i;}
	}
	const Photon *p1= fPhotons->At(imax1);
	const Photon *p2= fPhotons->At(imax2);
	Float_t GenPhotonID1=-999;
	Float_t GenPhotonMotherID1=-999;
	Float_t GenPhotonID2=-999;
	Float_t GenPhotonMotherID2=-999;
	int match1=0;
        int match2=0;
	match1=MatchRecPhotonsToGenPhotonsReal(p1);
        if(match1==1){
	  GenPhotonID1=fMatchedGenPhotonID;
	  GenPhotonMotherID1=fMatchedGenPhotonMotherID;
	}
        match2=MatchRecPhotonsToGenPhotonsReal(p2);
        if(match2==1){
	  GenPhotonID2=fMatchedGenPhotonID;
	  GenPhotonMotherID2=fMatchedGenPhotonMotherID;
	}
     	int p1real=0;
	int p2real=0;
        if( PhotonTools::MatchMC(p1,fMcParticles,kFALSE)){p1real=1;}
	if( PhotonTools::MatchMC(p2,fMcParticles,kFALSE)){p2real=1;}
	
        if(p1real==1 && p2real==1 && ptmax1>25 && ptmax2>25){
	  Double_t VtxProb = 1.0;
	  const Vertex* SelVtx = fVtxTools.findVtxBasicRanking(p1,p2,fBeamSpot->At(0),fPV,fConversions,kTRUE,VtxProb);
          //printf("VtxProb:%f\n",VtxProb);
          VtxProb = (Float_t)VtxProb;
	  FillPhotonTree(p1,p2,SelVtx,fPFCands,GenPhotonID1,GenPhotonMotherID1,VtxProb,EventNum,_decayZ);
	  FillPhotonTree(p2,p1,SelVtx,fPFCands,GenPhotonID2,GenPhotonMotherID2,VtxProb,EventNum,_decayZ);
  	}
      }
    }
    
    if(!fIsSignal){
      int PhotonColSize=fPhotons->GetEntries();
      if(PhotonColSize>=2){
	Float_t ptmax1=0;
	Float_t ptmax2=0;
	int imax1=0;
	int imax2=0;
	for(int i=0;i<PhotonColSize;i++){
	  if((fPhotons->At(i)->Pt())>ptmax1){
	    ptmax1=fPhotons->At(i)->Pt();
	    imax1=i;
	  }
	}
	for (int i=0;i<PhotonColSize;i++) {
	  if(((fPhotons->At(i)->Pt())>ptmax2)&&((fPhotons->At(i)->Pt()<ptmax1))){
	    ptmax2=fPhotons->At(i)->Pt();
	    imax2=i;}
	}
	const Photon *p1= fPhotons->At(imax1);
	const Photon *p2= fPhotons->At(imax2);
	Float_t GenPhotonID1=-999;
	Float_t GenPhotonMotherID1=-999;
	Float_t GenPhotonID2=-999;
	Float_t GenPhotonMotherID2=-999;
	int match1=0;
        int match2=0;
	match1=MatchRecPhotonsToGenPhotonsReal(p1);
        if(match1==1){
	  GenPhotonID1=fMatchedGenPhotonID;
	  GenPhotonMotherID1=fMatchedGenPhotonMotherID;
	}
        match2=MatchRecPhotonsToGenPhotonsReal(p2);
        if(match2==1){
	  GenPhotonID2=fMatchedGenPhotonID;
	  GenPhotonMotherID2=fMatchedGenPhotonMotherID;
	}
        //printf("GenPhotonID1:%f  ",GenPhotonID1);
        //printf("GenPhotonMotherID1:%f  ",GenPhotonMotherID1);
        //printf("GenPhotonID2:%f  ",GenPhotonID2);
        //printf("GenPhotonMotherID2:%f  ",GenPhotonMotherID2);
	int p1real=0;
	int p2real=0;
        if( PhotonTools::MatchMC(p1,fMcParticles,kFALSE)){p1real=1;}
	if( PhotonTools::MatchMC(p2,fMcParticles,kFALSE)){p2real=1;}
        int p1p2=p1real+p2real;
        if(p1p2==1 && ptmax1>25 && ptmax2>25){
	  Double_t VtxProb = 1.0;
	  const Vertex* SelVtx = fVtxTools.findVtxBasicRanking(p1,p2,fBeamSpot->At(0),fPV,fConversions,kTRUE,VtxProb);
	  // printf("VtxProb:%f\n",VtxProb);
          VtxProb = (Float_t)VtxProb;
          if(p1real==0){
	    FillPhotonTree(p1,p2,SelVtx,fPFCands,GenPhotonID1,GenPhotonMotherID1,VtxProb,EventNum,_decayZ);
	  }
          if(p2real==0){
	    FillPhotonTree(p2,p1,SelVtx,fPFCands,GenPhotonID2,GenPhotonMotherID2,VtxProb,EventNum,_decayZ);
	  }
  	}
      }
    }  
  }

  if(fIsData){
    Float_t EventNum=GetEventHeader()->EvtNum();
    int PhotonColSize=fPhotons->GetEntries();
    if(PhotonColSize>=2){
      Float_t ptmax1=0;
      Float_t ptmax2=0;
      int imax1=0;
      int imax2=0;
      for(int i=0;i<PhotonColSize;i++){
	if((fPhotons->At(i)->Pt())>ptmax1){
	  ptmax1=fPhotons->At(i)->Pt();
	  imax1=i;
	}
      }
      for (int i=0;i<PhotonColSize;i++) {
	if(((fPhotons->At(i)->Pt())>ptmax2)&&((fPhotons->At(i)->Pt()<ptmax1))){
	  ptmax2=fPhotons->At(i)->Pt();
	  imax2=i;}
      }
      const Photon *p1= fPhotons->At(imax1);
      const Photon *p2= fPhotons->At(imax2);
      Float_t GenPhotonID1=-999;
      Float_t GenPhotonMotherID1=-999;
      Float_t GenPhotonID2=-999;
      Float_t GenPhotonMotherID2=-999;
      int match1=0;
      int match2=0;
      if(ptmax1>25 && ptmax2>25){
	Double_t VtxProb = 1.0;
	const Vertex* SelVtx = fVtxTools.findVtxBasicRanking(p1,p2,fBeamSpot->At(0),fPV,fConversions,kTRUE,VtxProb);
	VtxProb = (Float_t)VtxProb;
	FillPhotonTree(p1,p2,SelVtx,fPFCands,GenPhotonID1,GenPhotonMotherID1,VtxProb,EventNum,_decayZ);
	FillPhotonTree(p2,p1,SelVtx,fPFCands,GenPhotonID2,GenPhotonMotherID2,VtxProb,EventNum,_decayZ);
      }
    }
    
  }
}
//-----------------------------------------------------------------------------
void MakeNtuple::SlaveBegin()
{
  //require brunch
  ReqEventObject(fPhotonName,fPhotons,fPhotonsFromBranch);
  ReqEventObject(fElectronName,       fElectrons, kTRUE);
  ReqEventObject(fTrackBranchName,    fTracks,    kTRUE);
  ReqEventObject(fPVName,             fPV,        kTRUE);
  ReqEventObject(fPileUpDenName,      fPileUpDen, kTRUE);
  ReqEventObject(fPileupInfoName,     fPileupInfos, kTRUE);
  ReqEventObject(fPFCandName,         fPFCands,    true);
  ReqBranch(fBeamSpotName,fBeamSpot);
  ReqBranch(fConversionName,fConversions);
  
  // for MC only to adjust potential overlaps from generation
  if (! fIsData) {
    ReqBranch(fMcEventInfoName,fMcEventInfo);
    printf(" Monte Carlo Information block %p\n",(void*) fMcEventInfo);
    printf(" --> this is no data. Access the McEventInfo.\n\n");
    
    ReqBranch(fMcParticleName,fMcParticles);
    
  }
  else
    printf(" --> this is data. Drop the McEventInfo.\n\n");
  
  //book ntuple
  hPhotonNtuple = new TNtuple("hPhotonNtuple","hPhotonNtuple","fprocessid:MatchGenPhotonID:MatchGenPhotonMotherID:Category:PassOfficial:PassElectronVeto:PassEleVetoConvRecovery:dRTrack:Pt:Eta:ScEta:ScPhi:ScEta_accompany:Et:EcalIsoDr03:HcalIsoDr03:TrkIsoHollowDr03:HoE:R9:covIEtaIEta:sigIEtaIEta:tIso1:tIso2:tIso3:tIso1abs:tIso2abs:tIso3abs:RelIsoEcal:RelIsoHcal:RelEMax:RelETop:RelEBottom:RelELeft:RelERight:RelE2x5Max:RelE2x5Top:RelE2x5Bottom:RelE2x5Left:RelE2x5Right:RelE5x5:RelE1x3:RelE3x1:RelE1x5:RelE2x2:RelE3x2:RelE4X4:EtaWidth:PhiWidth:RelPreshowerEnergy:CoviEtaiPhi:CoviPhiiPhi:NVertexes:_tRho:_NewRho_42:_RhoKt6PFJets:VtxProb:fNVertexesGenPile:EventNum:passpre:RawEnergy:AbsIsoEcal:AbsIsoHcal:GammaIso_DR0045To0p01:GammaIso_DR007To0p01:GammaIso_DR0To0p001:GammaIso_DR001To0p002:GammaIso_DR002To0p005:GammaIso_DR005To0p01:GammaIso_DR01To0p02:GammaIso_DR02To0p03:GammaIso_DR03To0p04:GammaIso_DR04To0p05:GammaIso_DR05To0p06:NeutralHadronIso_DR0To0p001:NeutralHadronIso_DR001To0p002:NeutralHadronIso_DR002To0p005:NeutralHadronIso_DR005To0p01:NeutralHadronIso_DR01To0p02:NeutralHadronIso_DR02To0p03:NeutralHadronIso_DR03To0p04:NeutralHadronIso_DR04To0p05:NeutralHadronIso_DR05To0p06:ChargedIso_selvtx_DR0To0p001:ChargedIso_selvtx_DR001To0p002:ChargedIso_selvtx_DR002To0p005:ChargedIso_selvtx_DR005To0p01:ChargedIso_selvtx_DR01To0p02:ChargedIso_selvtx_DR02To0p03:ChargedIso_selvtx_DR03To0p04:ChargedIso_selvtx_DR04To0p05:ChargedIso_selvtx_DR05To0p06:ChargedIso_selvtx_DR002To0p01:ChargedIso_worstvtx_DR0To0p001:ChargedIso_worstvtx_DR001To0p002:ChargedIso_worstvtx_DR002To0p005:ChargedIso_worstvtx_DR005To0p01:ChargedIso_worstvtx_DR01To0p02:ChargedIso_worstvtx_DR02To0p03:ChargedIso_worstvtx_DR03To0p04:ChargedIso_worstvtx_DR04To0p05:ChargedIso_worstvtx_DR05To0p06:ChargedIso_worstvtx_DR002To0p01:ChargedCount_selvtx_DR0To0p001:ChargedCount_selvtx_DR001To0p002:ChargedCount_selvtx_DR002To0p005:ChargedCount_selvtx_DR005To0p01:ChargedCount_selvtx_DR01To0p02:ChargedCount_selvtx_DR02To0p03:ChargedCount_selvtx_DR03To0p04:ChargedCount_selvtx_DR04To0p05:ChargedCount_selvtx_DR05To0p06:ChargedCount_selvtx_DR002To0p01:ChargedCount_worstvtx_DR0To0p001:ChargedCount_worstvtx_DR001To0p002:ChargedCount_worstvtx_DR002To0p005:ChargedCount_worstvtx_DR005To0p01:ChargedCount_worstvtx_DR01To0p02:ChargedCount_worstvtx_DR02To0p03:ChargedCount_worstvtx_DR03To0p04:ChargedCount_worstvtx_DR04To0p05:ChargedCount_worstvtx_DR05To0p06:ChargedCount_worstvtx_DR002To0p01:DiphotonMass:RelEmaxOverE33:RelE22:CovEtaPhi:PsEffWidthSigmaXX:PsEffWidthSigmaYY:vtxZ:_decayZ:s4ratio:lambdaratio");

  AddOutput(hPhotonNtuple);
  
  fVtxTools.InitP();
}

//-----------------------------------------------------------------------------------------------------------
int MakeNtuple::MatchRecPhotonsToGenPhotonsReal(const Photon *photonRec)//ming: match the reconstructed photon to the generated photon
{
  int MatchReal=0;
  
  Float_t MinDeltaR=0.3;
  int MinDeltaRj=-1;
  Float_t photonRec_pt=photonRec->Pt();
  
  // loop through all generated photons and try to find a match for the reconstructed photon being considerated
  for(UInt_t j=0; j<fMcParticles->GetEntries(); j++){
    const MCParticle *photonGen =fMcParticles->At(j);
    const MCParticle *Mother=photonGen->DistinctMother();
    if(photonGen->AbsPdgId()==22 && photonGen->IsGenerated()){//photon PdgId=22
      Float_t photonGen_pt=photonGen->Pt();
      Float_t diffPtRatio=fabs(photonGen_pt-photonRec_pt)/photonRec_pt;
      Float_t deltaR=MathUtils::DeltaR(*photonRec,*photonGen);
      if (deltaR<MinDeltaR && diffPtRatio<0.5) {
	MinDeltaR=deltaR;
	MinDeltaRj=j;
      }
    }
  }
  
  if(MinDeltaRj!=-1){
    const MCParticle *MatchedGenPhoton=fMcParticles->At(MinDeltaRj);
    const MCParticle *MatchedGenPhotonMother=MatchedGenPhoton->DistinctMother();
    fMatchedGenPhotonID=MatchedGenPhoton->PdgId();
    fMatchedGenPhotonMotherID=MatchedGenPhotonMother->PdgId();
    if(fabs(fMatchedGenPhotonMotherID)==25 ||(fabs(fMatchedGenPhotonMotherID)>0 && fabs(fMatchedGenPhotonMotherID)<10) || fabs(fMatchedGenPhotonMotherID)==21){
      MatchReal=1;
    }
  }
  return MatchReal;
}

void MakeNtuple::FillPhotonTree(const Photon *p,const Photon *p_accompany,const Vertex *SelVtx,const PFCandidateCol *fPFCands,Float_t GenPhotonID,Float_t GenPhotonMotherID,Float_t VtxProb, Float_t EventNum,Float_t _decayZ) {

  // these values are taken from the H2GGlobe code... (actually from Marco/s mail)
  float cic4_allcuts_temp_sublead[] = { 
    3.8,         2.2,         1.77,        1.29,
    11.7,        3.4,         3.9,         1.84,
    3.5,         2.2,         2.3,         1.45,
    0.0106,      0.0097,      0.028,       0.027,
    0.082,       0.062,       0.065,       0.048,
    0.94,        0.36,        0.94,        0.32,
    1.,          0.062,       0.97,        0.97,
    1.5,         1.5,         1.5,         1.5 };  // the last line is PixelmatchVeto and un-used
 
  Float_t DiphotonMass=0; 
  Float_t vtxZ=SelVtx->Z(); 
   
  Float_t MatchGenPhotonID=GenPhotonID;
  Float_t MatchGenPhotonMotherID=GenPhotonMotherID;
  Float_t Category=0;
  Float_t PassOfficial=0;
  Float_t PassElectronVeto=0;
  Float_t PassEleVetoConvRecovery=0;
  Float_t dRTrack=0;
  Float_t Pt=0;
  Float_t Eta=0;
  Float_t ScEta=0;
  Float_t ScPhi=0;
  Float_t ScEta_accompany=0;
  Float_t Et=0;
  Float_t EcalIsoDr03=0;
  Float_t HcalIsoDr03=0;
  Float_t TrkIsoHollowDr03=0;
  Float_t HoE=0;
  Float_t R9=0;
  Float_t covIEtaIEta=0;
  Float_t sigIEtaIEta=0;
  Float_t tIso1=0;
  Float_t tIso2=0;
  Float_t tIso3=0;
  Float_t tIso1abs=0;
  Float_t tIso2abs=0;
  Float_t tIso3abs=0;
  Float_t RelIsoEcal=0;
  Float_t RelIsoHcal=0;
  Float_t RelEMax=0;
  Float_t RelETop=0;
  Float_t RelEBottom=0;
  Float_t RelELeft=0;
  Float_t RelERight=0;
  Float_t RelE2x5Max=0;
  Float_t RelE2x5Top=0;
  Float_t RelE2x5Bottom=0;
  Float_t RelE2x5Left=0;
  Float_t RelE2x5Right=0;
  Float_t RelE5x5=0;
  Float_t RelE1x3=0;
  Float_t RelE3x1=0;
  Float_t RelE1x5=0;
  Float_t RelE2x2=0;
  Float_t RelE3x2=0;
  Float_t RelE4X4=0;
  Float_t EtaWidth=0;
  Float_t PhiWidth=0;
  Float_t RelPreshowerEnergy=0;
  Float_t CoviEtaiPhi=0;
  Float_t CoviPhiiPhi=0;
  Float_t NVertexes=0;
  Float_t _tRho=0;
  Float_t _NewRho_42=0;
  Float_t _RhoKt6PFJets=0;
  Float_t AbsIsoEcal=0;
  Float_t AbsIsoHcal=0;

   //get the variables used to compute MVA variables
   Float_t ecalIso3 = p->EcalRecHitIsoDr03();
   Float_t ecalIso4 = p->EcalRecHitIsoDr04();
   Float_t hcalIso4 = p->HcalTowerSumEtDr04();
   unsigned int wVtxInd = 0;
   Float_t trackIso1 = IsolationTools::CiCTrackIsolation(p,SelVtx, 0.3, 0.02, 0.0, 0.0, 0.1, 1.0, fTracks);
   Float_t trackIso2 = IsolationTools::CiCTrackIsolation(p,SelVtx, 0.4, 0.02, 0.0, 0.0, 0.1, 1.0, fTracks, &wVtxInd,fPV);
   Float_t trackIso3 = IsolationTools::CiCTrackIsolation(p,SelVtx, 0.3, 0.02, 0.0, 0.0, 0.1, 1.0, fTracks);
   _tRho = (Float_t) fPileUpDen->At(0)->RhoRandomLowEta();
   _NewRho_42 = (Float_t) fPileUpDen->At(0)->Rho();
   _RhoKt6PFJets = (Float_t) fPileUpDen->At(0)->RhoKt6PFJets();
   Float_t combIso1 = ecalIso3+hcalIso4+trackIso1 - 0.17*_tRho;
   Float_t combIso2 = ecalIso4+hcalIso4+trackIso2 - 0.52*_tRho;
   Float_t RawEnergy = p->SCluster()->RawEnergy();

   //mass
   DiphotonMass=(p->Mom()+p_accompany->Mom()).M();

   //preselection variables
   Et=p->Et();
   EcalIsoDr03=p->EcalRecHitIsoDr03();  
   HcalIsoDr03=p->HcalTowerSumEtDr03();
   TrkIsoHollowDr03=p->HollowConeTrkIsoDr03();
   
   //compute MVA variables
   HoE = p->HadOverEm();
   covIEtaIEta = p->CoviEtaiEta();
   sigIEtaIEta = sqrt(p->SCluster()->Seed()->CoviEtaiEta()); 
   dRTrack = PhotonTools::ElectronVetoCiC(p,fElectrons); 
   tIso1 = (combIso1) *50./p->Et();
   tIso3 = (trackIso3)*50./p->Et();
   tIso2 = (combIso2) *50./(p->MomVtx(fPV->At(wVtxInd)->Position()).Pt());
   tIso1abs = (combIso1);
   tIso3abs = (trackIso3);
   tIso2abs = (combIso2);
   R9 = p->R9();

   //newly added MVA variables 1
   RelIsoEcal=(ecalIso3-0.17*_tRho)/p->Et();
   RelIsoHcal=(hcalIso4-0.17*_tRho)/p->Et();

   AbsIsoEcal=ecalIso3;
   AbsIsoHcal=hcalIso4;

   RelEMax=p->SCluster()->Seed()->EMax()/RawEnergy;
   RelETop=p->SCluster()->Seed()->ETop()/RawEnergy;
   RelEBottom=p->SCluster()->Seed()->EBottom()/RawEnergy;
   RelELeft=p->SCluster()->Seed()->ELeft()/RawEnergy;
   RelERight=p->SCluster()->Seed()->ERight()/RawEnergy;
   RelE2x5Max=p->SCluster()->Seed()->E2x5Max()/RawEnergy;
   RelE2x5Top=p->SCluster()->Seed()->E2x5Top()/RawEnergy;
   RelE2x5Bottom=p->SCluster()->Seed()->E2x5Bottom()/RawEnergy;
   RelE2x5Left=p->SCluster()->Seed()->E2x5Left()/RawEnergy;
   RelE2x5Right=p->SCluster()->Seed()->E2x5Right()/RawEnergy;
   RelE5x5=p->SCluster()->Seed()->E5x5()/RawEnergy;

   //newly added MVA variables 2
   EtaWidth=p->SCluster()->EtaWidth();
   PhiWidth=p->SCluster()->PhiWidth();
   RelPreshowerEnergy=p->SCluster()->PreshowerEnergy()/RawEnergy;

   CoviEtaiPhi=p->SCluster()->Seed()->CoviEtaiPhi(); 
   CoviPhiiPhi=p->SCluster()->Seed()->CoviPhiiPhi();

   RelE1x3=p->SCluster()->Seed()->E1x3()/RawEnergy;
   RelE3x1=p->SCluster()->Seed()->E3x1()/RawEnergy;
   RelE1x5=p->SCluster()->Seed()->E1x5()/RawEnergy;
   RelE2x2=p->SCluster()->Seed()->E2x2()/RawEnergy;
   RelE3x2=p->SCluster()->Seed()->E3x2()/RawEnergy;
   RelE4X4=p->SCluster()->Seed()->E4x4()/RawEnergy;

   //spectator variables 1
   Pt=p->Pt();
   Eta=p->Eta();
   ScEta=p->SCluster()->Eta();
   ScPhi=p->SCluster()->Phi();
   ScEta_accompany=p_accompany->SCluster()->Eta();;
   NVertexes=fPV->GetEntries();//not for mva yet

   //Official selection variables
   if(fabs(ScEta)<1.4442 && p->R9()>=0.94){
     Category=1;  
     if( PhotonTools::PassCiCSelection(p, SelVtx, fTracks, fElectrons, fPV, _tRho, fPhotonPtMin, fApplyElectronVeto) ){
       PassOfficial=1;
     }
   }
   if(fabs(ScEta)<1.4442 && p->R9()<0.94){
     Category=2;  
     if( PhotonTools::PassCiCSelection(p, SelVtx, fTracks, fElectrons, fPV, _tRho, fPhotonPtMin, fApplyElectronVeto) ){
       PassOfficial=1;
     }
   }
   if(fabs(ScEta)>1.566 && fabs(ScEta)<2.5 && p->R9()>=0.94){
     Category=3; 
     if( PhotonTools::PassCiCSelection(p, SelVtx, fTracks, fElectrons, fPV, _tRho, fPhotonPtMin, fApplyElectronVeto) ){
       PassOfficial=1;
     }
     
   }
   if(fabs(ScEta)>1.566 && fabs(ScEta)<2.5 && p->R9()<0.94){
     Category=4; 
     if( PhotonTools::PassCiCSelection(p, SelVtx, fTracks, fElectrons, fPV, _tRho, fPhotonPtMin, fApplyElectronVeto) ){
       PassOfficial=1;
     }
   }
   
   //cic ele veto
   Bool_t isbarrel = (fabs(ScEta)<1.4442);
   int _tCat = 1;
   if ( !isbarrel ) _tCat = 3;
   if ( R9 < 0.94 ) _tCat++;
   if(dRTrack > cic4_allcuts_temp_sublead[_tCat-1+6*4]){
     PassElectronVeto=1;
   }

   //mit ele veto
   if(PhotonTools::PassElectronVetoConvRecovery(p,fElectrons,fConversions,fBeamSpot->At(0))){
     PassEleVetoConvRecovery=1;
   }

   float passpre=0;

   if(PhotonTools::PassSinglePhotonPresel(p,fElectrons,fConversions,fBeamSpot->At(0), fTracks,SelVtx,_tRho,fApplyElectronVeto,finvertElectronVeto)){
     passpre=1; 
   }

   //PF ISO
   //photon
   Float_t GammaIso_DR0To0p001=0;
   Float_t GammaIso_DR001To0p002=0;
   Float_t GammaIso_DR002To0p005=0;
   Float_t GammaIso_DR005To0p01=0;
   Float_t GammaIso_DR01To0p02=0;
   Float_t GammaIso_DR02To0p03=0; 
   Float_t GammaIso_DR03To0p04=0;
   Float_t GammaIso_DR04To0p05=0; 
   Float_t GammaIso_DR05To0p06=0; 

   Float_t GammaIso_DR0045To0p01=0;
   Float_t GammaIso_DR007To0p01=0; 

   //Neutral Hadron
   Float_t NeutralHadronIso_DR0To0p001=0;
   Float_t NeutralHadronIso_DR001To0p002=0;
   Float_t NeutralHadronIso_DR002To0p005=0;
   Float_t NeutralHadronIso_DR005To0p01=0;
   Float_t NeutralHadronIso_DR01To0p02=0;
   Float_t NeutralHadronIso_DR02To0p03=0; 
   Float_t NeutralHadronIso_DR03To0p04=0;
   Float_t NeutralHadronIso_DR04To0p05=0; 
   Float_t NeutralHadronIso_DR05To0p06=0;
   //Charged iso
   Float_t ChargedIso_selvtx_DR0To0p001=0;
   Float_t ChargedIso_selvtx_DR001To0p002=0;
   Float_t ChargedIso_selvtx_DR002To0p005=0;
   Float_t ChargedIso_selvtx_DR005To0p01=0;
   Float_t ChargedIso_selvtx_DR01To0p02=0;
   Float_t ChargedIso_selvtx_DR02To0p03=0; 
   Float_t ChargedIso_selvtx_DR03To0p04=0;
   Float_t ChargedIso_selvtx_DR04To0p05=0; 
   Float_t ChargedIso_selvtx_DR05To0p06=0;
   Float_t ChargedIso_selvtx_DR002To0p01=0;

   Float_t ChargedIso_worstvtx_DR0To0p001=0;
   Float_t ChargedIso_worstvtx_DR001To0p002=0;
   Float_t ChargedIso_worstvtx_DR002To0p005=0;
   Float_t ChargedIso_worstvtx_DR005To0p01=0;
   Float_t ChargedIso_worstvtx_DR01To0p02=0;
   Float_t ChargedIso_worstvtx_DR02To0p03=0; 
   Float_t ChargedIso_worstvtx_DR03To0p04=0;
   Float_t ChargedIso_worstvtx_DR04To0p05=0; 
   Float_t ChargedIso_worstvtx_DR05To0p06=0;
   Float_t ChargedIso_worstvtx_DR002To0p01=0;

   //Charged count
   Float_t ChargedCount_selvtx_DR0To0p001=0;
   Float_t ChargedCount_selvtx_DR001To0p002=0;
   Float_t ChargedCount_selvtx_DR002To0p005=0;
   Float_t ChargedCount_selvtx_DR005To0p01=0;
   Float_t ChargedCount_selvtx_DR01To0p02=0;
   Float_t ChargedCount_selvtx_DR02To0p03=0; 
   Float_t ChargedCount_selvtx_DR03To0p04=0;
   Float_t ChargedCount_selvtx_DR04To0p05=0; 
   Float_t ChargedCount_selvtx_DR05To0p06=0;
   Float_t ChargedCount_selvtx_DR002To0p01=0;
   
   Float_t ChargedCount_worstvtx_DR0To0p001=0;
   Float_t ChargedCount_worstvtx_DR001To0p002=0;
   Float_t ChargedCount_worstvtx_DR002To0p005=0;
   Float_t ChargedCount_worstvtx_DR005To0p01=0;
   Float_t ChargedCount_worstvtx_DR01To0p02=0;
   Float_t ChargedCount_worstvtx_DR02To0p03=0; 
   Float_t ChargedCount_worstvtx_DR03To0p04=0;
   Float_t ChargedCount_worstvtx_DR04To0p05=0; 
   Float_t ChargedCount_worstvtx_DR05To0p06=0;
   Float_t ChargedCount_worstvtx_DR002To0p01=0;

   //new shower shape
   Float_t RelEmaxOverE33=0; 
   Float_t RelE22=0; 
   Float_t CovEtaPhi=0; 
   Float_t PsEffWidthSigmaXX=0; 
   Float_t PsEffWidthSigmaYY=0; 

   Float_t s4ratio=0;
   Float_t lambdaratio=0;

   for(UInt_t i=0; i<fPFCands->GetEntries(); i++){

     const PFCandidate *pf= fPFCands->At(i);
     Float_t dr = mithep::MathUtils::DeltaR(p->Phi(),p->Eta(), pf->Phi(), pf->Eta());
     Float_t deta = fabs(p->Eta()-pf->Eta());

     if(deta>0.015 && (fabs(ScEta)<1.4442 && pf->Pt()>0.08) || (fabs(ScEta)>1.566 && fabs(ScEta)<2.5 && pf->Pt()>0.1) ){
       if(pf->PFType()==PFCandidate::eGamma){
	 if(dr>0 && dr<=0.01) GammaIso_DR0To0p001+=pf->Pt();
	 if(dr>0.01 && dr<=0.02) GammaIso_DR001To0p002+=pf->Pt();
	 if(dr>0.02 && dr<=0.05) GammaIso_DR002To0p005+=pf->Pt();
	 if(dr>0.05 && dr<=0.1) GammaIso_DR005To0p01+=pf->Pt();
	 if(dr>0.1 && dr<=0.2) GammaIso_DR01To0p02+=pf->Pt();
	 if(dr>0.2 && dr<=0.3) GammaIso_DR02To0p03+=pf->Pt();
	 if(dr>0.3 && dr<=0.4) GammaIso_DR03To0p04+=pf->Pt();
	 if(dr>0.4 && dr<=0.5) GammaIso_DR04To0p05+=pf->Pt();
	 if(dr>0.5 && dr<=0.6) GammaIso_DR05To0p06+=pf->Pt();
	 
	 if(dr>0.045 && dr<=0.1) GammaIso_DR0045To0p01+=pf->Pt();
	 if(dr>0.07 && dr<=0.1) GammaIso_DR007To0p01+=pf->Pt(); 
       }
     }
     
     if(pf->PFType()==PFCandidate::eNeutralHadron){
       if(dr>0 && dr<=0.01) NeutralHadronIso_DR0To0p001+=pf->Pt();
       if(dr>0.01 && dr<=0.02) NeutralHadronIso_DR001To0p002+=pf->Pt();
       if(dr>0.02 && dr<=0.05) NeutralHadronIso_DR002To0p005+=pf->Pt(); 
       if(dr>0.05 && dr<=0.1) NeutralHadronIso_DR005To0p01+=pf->Pt();
       if(dr>0.1 && dr<=0.2) NeutralHadronIso_DR01To0p02+=pf->Pt();
       if(dr>0.2 && dr<=0.3) NeutralHadronIso_DR02To0p03+=pf->Pt();
       if(dr>0.3 && dr<=0.4) NeutralHadronIso_DR03To0p04+=pf->Pt();
       if(dr>0.4 && dr<=0.5) NeutralHadronIso_DR04To0p05+=pf->Pt();
       if(dr>0.5 && dr<=0.6) NeutralHadronIso_DR05To0p06+=pf->Pt();
     }
   }

   //iso
   ChargedIso_selvtx_DR0To0p001=IsolationTools::PFChargedIsolation(p,SelVtx, 0.01, 0,1, 0.0, 0.1, 0.2,fPFCands);
   ChargedIso_selvtx_DR001To0p002=IsolationTools::PFChargedIsolation(p,SelVtx, 0.02, 0.01, 1, 0.0, 0.1, 0.2,fPFCands);
   ChargedIso_selvtx_DR002To0p005=IsolationTools::PFChargedIsolation(p,SelVtx, 0.05, 0.02, 1, 0.0, 0.1, 0.2,fPFCands);
   ChargedIso_selvtx_DR005To0p01=IsolationTools::PFChargedIsolation(p,SelVtx, 0.1, 0.05, 1, 0.0, 0.1, 0.2,fPFCands);
   ChargedIso_selvtx_DR01To0p02=IsolationTools::PFChargedIsolation(p,SelVtx, 0.2, 0.1, 1, 0.0, 0.1, 0.2,fPFCands);
   ChargedIso_selvtx_DR02To0p03=IsolationTools::PFChargedIsolation(p,SelVtx, 0.3, 0.2, 1, 0.0, 0.1, 0.2,fPFCands);
   ChargedIso_selvtx_DR03To0p04=IsolationTools::PFChargedIsolation(p,SelVtx, 0.4, 0.3, 1, 0.0, 0.1, 0.2,fPFCands);
   ChargedIso_selvtx_DR04To0p05=IsolationTools::PFChargedIsolation(p,SelVtx, 0.5, 0.4, 1, 0.0, 0.1, 0.2,fPFCands);
   ChargedIso_selvtx_DR05To0p06=IsolationTools::PFChargedIsolation(p,SelVtx, 0.6, 0.5, 1, 0.0, 0.1, 0.2,fPFCands);
   ChargedIso_selvtx_DR002To0p01=IsolationTools::PFChargedIsolation(p,SelVtx, 0.1, 0.02, 1, 0.0, 0.1, 0.2,fPFCands);
   
   ChargedIso_worstvtx_DR0To0p001=IsolationTools::PFChargedIsolation(p,SelVtx, 0.01, 0,1, 0.0, 0.1, 0.2,fPFCands,&wVtxInd,fPV);
   ChargedIso_worstvtx_DR001To0p002=IsolationTools::PFChargedIsolation(p,SelVtx, 0.02, 0.01, 1, 0.0, 0.1, 0.2,fPFCands,&wVtxInd,fPV);
   ChargedIso_worstvtx_DR002To0p005=IsolationTools::PFChargedIsolation(p,SelVtx, 0.05, 0.02, 1, 0.0, 0.1, 0.2,fPFCands,&wVtxInd,fPV);
   ChargedIso_worstvtx_DR005To0p01=IsolationTools::PFChargedIsolation(p,SelVtx, 0.1, 0.05, 1, 0.0, 0.1, 0.2,fPFCands, &wVtxInd,fPV);
   ChargedIso_worstvtx_DR01To0p02=IsolationTools::PFChargedIsolation(p,SelVtx, 0.2, 0.1, 1, 0.0, 0.1, 0.2,fPFCands, &wVtxInd,fPV);
   ChargedIso_worstvtx_DR02To0p03=IsolationTools::PFChargedIsolation(p,SelVtx, 0.3, 0.2, 1, 0.0, 0.1, 0.2,fPFCands, &wVtxInd,fPV);
   ChargedIso_worstvtx_DR03To0p04=IsolationTools::PFChargedIsolation(p,SelVtx, 0.4, 0.3, 1, 0.0, 0.1, 0.2,fPFCands, &wVtxInd,fPV);
   ChargedIso_worstvtx_DR04To0p05=IsolationTools::PFChargedIsolation(p,SelVtx, 0.5, 0.4, 1, 0.0, 0.1, 0.2,fPFCands, &wVtxInd,fPV);
   ChargedIso_worstvtx_DR05To0p06=IsolationTools::PFChargedIsolation(p,SelVtx, 0.6, 0.5, 1, 0.0, 0.1, 0.2,fPFCands, &wVtxInd,fPV);
   ChargedIso_worstvtx_DR002To0p01=IsolationTools::PFChargedIsolation(p,SelVtx, 0.1, 0.02, 1, 0.0, 0.1, 0.2,fPFCands, &wVtxInd,fPV);

   //count
   ChargedCount_selvtx_DR0To0p001=IsolationTools::PFChargedCount(p,SelVtx, 0.01, 0,1.6, 0.0, 0.1, 0.2,fPFCands);
   ChargedCount_selvtx_DR001To0p002=IsolationTools::PFChargedCount(p,SelVtx, 0.02, 0.01, 1.6, 0.0, 0.1, 0.2,fPFCands);
   ChargedCount_selvtx_DR002To0p005=IsolationTools::PFChargedCount(p,SelVtx, 0.05, 0.02, 1.6, 0.0, 0.1, 0.2,fPFCands);
   ChargedCount_selvtx_DR005To0p01=IsolationTools::PFChargedCount(p,SelVtx, 0.1, 0.05, 1.6, 0.0, 0.1, 0.2,fPFCands);
   ChargedCount_selvtx_DR01To0p02=IsolationTools::PFChargedCount(p,SelVtx, 0.2, 0.1, 1.6, 0.0, 0.1, 0.2,fPFCands);
   ChargedCount_selvtx_DR02To0p03=IsolationTools::PFChargedCount(p,SelVtx, 0.3, 0.2, 1.6, 0.0, 0.1, 0.2,fPFCands);
   ChargedCount_selvtx_DR03To0p04=IsolationTools::PFChargedCount(p,SelVtx, 0.4, 0.3, 1.6, 0.0, 0.1, 0.2,fPFCands);
   ChargedCount_selvtx_DR04To0p05=IsolationTools::PFChargedCount(p,SelVtx, 0.5, 0.4, 1.6, 0.0, 0.1, 0.2,fPFCands);
   ChargedCount_selvtx_DR05To0p06=IsolationTools::PFChargedCount(p,SelVtx, 0.6, 0.5, 1.6, 0.0, 0.1, 0.2,fPFCands);
   ChargedCount_selvtx_DR002To0p01=IsolationTools::PFChargedCount(p,SelVtx, 0.1, 0.02, 1.6, 0.0, 0.1, 0.2,fPFCands);
   
   ChargedCount_worstvtx_DR0To0p001=IsolationTools::PFChargedCount(p,SelVtx, 0.01, 0,1.6, 0.0, 0.1, 0.2,fPFCands,&wVtxInd,fPV);
   ChargedCount_worstvtx_DR001To0p002=IsolationTools::PFChargedCount(p,SelVtx, 0.02, 0.01, 1.6, 0.0, 0.1, 0.2,fPFCands,&wVtxInd,fPV);
   ChargedCount_worstvtx_DR002To0p005=IsolationTools::PFChargedCount(p,SelVtx, 0.05, 0.02, 1.6, 0.0, 0.1, 0.2,fPFCands,&wVtxInd,fPV);
   ChargedCount_worstvtx_DR005To0p01=IsolationTools::PFChargedCount(p,SelVtx, 0.1, 0.05, 1.6, 0.0, 0.1, 0.2,fPFCands, &wVtxInd,fPV);
   ChargedCount_worstvtx_DR01To0p02=IsolationTools::PFChargedCount(p,SelVtx, 0.2, 0.1, 1.6, 0.0, 0.1, 0.2,fPFCands, &wVtxInd,fPV);
   ChargedCount_worstvtx_DR02To0p03=IsolationTools::PFChargedCount(p,SelVtx, 0.3, 0.2, 1.6, 0.0, 0.1, 0.2,fPFCands, &wVtxInd,fPV);
   ChargedCount_worstvtx_DR03To0p04=IsolationTools::PFChargedCount(p,SelVtx, 0.4, 0.3, 1.6, 0.0, 0.1, 0.2,fPFCands, &wVtxInd,fPV);
   ChargedCount_worstvtx_DR04To0p05=IsolationTools::PFChargedCount(p,SelVtx, 0.5, 0.4, 1.6, 0.0, 0.1, 0.2,fPFCands, &wVtxInd,fPV);
   ChargedCount_worstvtx_DR05To0p06=IsolationTools::PFChargedCount(p,SelVtx, 0.6, 0.5, 1.6, 0.0, 0.1, 0.2,fPFCands, &wVtxInd,fPV);
   ChargedCount_worstvtx_DR002To0p01=IsolationTools::PFChargedCount(p,SelVtx, 0.1, 0.02, 1.6, 0.0, 0.1, 0.2,fPFCands, &wVtxInd,fPV);
  

   RelEmaxOverE33=p->SCluster()->Seed()->EMax()/p->SCluster()->Seed()->E3x3(); 
   RelE22=p->SCluster()->Seed()->E2x2()/p->SCluster()->RawEnergy(); 
   CovEtaPhi=p->SCluster()->Seed()->CovEtaPhi(); 
   PsEffWidthSigmaXX=p->SCluster()->PsEffWidthSigmaXX(); 
   PsEffWidthSigmaYY=p->SCluster()->PsEffWidthSigmaYY(); 
   
   s4ratio=p->SCluster()->Seed()->E2x2()/p->SCluster()->Seed()->E5x5(); 

   Float_t CovEtaEta=p->SCluster()->Seed()->CovEtaEta();
   Float_t CovPhiPhi=p->SCluster()->Seed()->CovPhiPhi();

   lambdaratio=(CovEtaEta+CovPhiPhi-sqrt((CovEtaEta-CovPhiPhi)*(CovEtaEta-CovPhiPhi) +4*CovEtaPhi*CovEtaPhi))/(CovEtaEta+CovPhiPhi+sqrt((CovEtaEta-CovPhiPhi)*(CovEtaEta-CovPhiPhi)+4*CovEtaPhi*CovEtaPhi));

   //fill the ntuple
   Float_t fillPhotonMatchRealGenTree[]={fprocessid,MatchGenPhotonID,MatchGenPhotonMotherID,Category,PassOfficial,PassElectronVeto,PassEleVetoConvRecovery,dRTrack,Pt,Eta,ScEta,ScPhi,ScEta_accompany,Et,EcalIsoDr03,HcalIsoDr03,TrkIsoHollowDr03,HoE,R9,covIEtaIEta,sigIEtaIEta,tIso1,tIso2,tIso3,tIso1abs,tIso2abs,tIso3abs,RelIsoEcal,RelIsoHcal,RelEMax,RelETop,RelEBottom,RelELeft,RelERight,RelE2x5Max,RelE2x5Top,RelE2x5Bottom,RelE2x5Left,RelE2x5Right,RelE5x5,RelE1x3,RelE3x1,RelE1x5,RelE2x2,RelE3x2,RelE4X4,EtaWidth,PhiWidth,RelPreshowerEnergy,CoviEtaiPhi,CoviPhiiPhi,NVertexes,_tRho,_NewRho_42,_RhoKt6PFJets,VtxProb,fNVertexesGenPile,EventNum,passpre,RawEnergy,AbsIsoEcal,AbsIsoHcal,GammaIso_DR0045To0p01,GammaIso_DR007To0p01,GammaIso_DR0To0p001,GammaIso_DR001To0p002,GammaIso_DR002To0p005,GammaIso_DR005To0p01,GammaIso_DR01To0p02,GammaIso_DR02To0p03,GammaIso_DR03To0p04,GammaIso_DR04To0p05,GammaIso_DR05To0p06,NeutralHadronIso_DR0To0p001,NeutralHadronIso_DR001To0p002,NeutralHadronIso_DR002To0p005,NeutralHadronIso_DR005To0p01,NeutralHadronIso_DR01To0p02,NeutralHadronIso_DR02To0p03,NeutralHadronIso_DR03To0p04,NeutralHadronIso_DR04To0p05,NeutralHadronIso_DR05To0p06,ChargedIso_selvtx_DR0To0p001,ChargedIso_selvtx_DR001To0p002,ChargedIso_selvtx_DR002To0p005,ChargedIso_selvtx_DR005To0p01,ChargedIso_selvtx_DR01To0p02,ChargedIso_selvtx_DR02To0p03,ChargedIso_selvtx_DR03To0p04,ChargedIso_selvtx_DR04To0p05,ChargedIso_selvtx_DR05To0p06,ChargedIso_selvtx_DR002To0p01,ChargedIso_worstvtx_DR0To0p001,ChargedIso_worstvtx_DR001To0p002,ChargedIso_worstvtx_DR002To0p005,ChargedIso_worstvtx_DR005To0p01,ChargedIso_worstvtx_DR01To0p02,ChargedIso_worstvtx_DR02To0p03,ChargedIso_worstvtx_DR03To0p04,ChargedIso_worstvtx_DR04To0p05,ChargedIso_worstvtx_DR05To0p06,ChargedIso_worstvtx_DR002To0p01,ChargedCount_selvtx_DR0To0p001,ChargedCount_selvtx_DR001To0p002,ChargedCount_selvtx_DR002To0p005,ChargedCount_selvtx_DR005To0p01,ChargedCount_selvtx_DR01To0p02,ChargedCount_selvtx_DR02To0p03,ChargedCount_selvtx_DR03To0p04,ChargedCount_selvtx_DR04To0p05,ChargedCount_selvtx_DR05To0p06,ChargedCount_selvtx_DR002To0p01,ChargedCount_worstvtx_DR0To0p001,ChargedCount_worstvtx_DR001To0p002,ChargedCount_worstvtx_DR002To0p005,ChargedCount_worstvtx_DR005To0p01,ChargedCount_worstvtx_DR01To0p02,ChargedCount_worstvtx_DR02To0p03,ChargedCount_worstvtx_DR03To0p04,ChargedCount_worstvtx_DR04To0p05,ChargedCount_worstvtx_DR05To0p06,ChargedCount_worstvtx_DR002To0p01,DiphotonMass,RelEmaxOverE33,RelE22,CovEtaPhi,PsEffWidthSigmaXX,PsEffWidthSigmaYY,vtxZ,_decayZ,s4ratio,lambdaratio};
	
   hPhotonNtuple->Fill(fillPhotonMatchRealGenTree);
}

void MakeNtuple::FindHiggsPtAndZ(Float_t& pt, Float_t& decayZ, Float_t& mass)
{
  pt     = -999.;
  decayZ = -999.;
  mass   = -999.;

  // loop over all GEN particles and look for status 1 photons
  for(UInt_t i=0; i<fMcParticles->GetEntries(); ++i) {
    const MCParticle* p = fMcParticles->At(i);
    if (p->Is(MCParticle::kH) || (!fApplyElectronVeto &&
				  (p->AbsPdgId()==23 || p->AbsPdgId()==24))) {
      pt     = p->Pt();
      decayZ = p->DecayVertex().Z();
      mass   = p->Mass();
      break;
    }
  }
  
  return;
}
