// $Id: GenFakeableObjsMod.cc,v 1.1 2009/06/30 10:47:17 loizides Exp $

#include "MitPhysics/FakeMods/interface/GenFakeableObjsMod.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"

using namespace mithep;

ClassImp(mithep::GenFakeableObjsMod)

//--------------------------------------------------------------------------------------------------
GenFakeableObjsMod::GenFakeableObjsMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fVetoTriggerJet(false),
  fVetoGenLeptons(true),
  fVetoCleanLeptons(false),
  fElectronFOType("GsfPlusSC"),
  fMuonFOType("IsoTrack"),
  fTriggerName("NotSpecified"),
  fTriggerObjectsName("NotSpecified"),
  fElectronBranchName(Names::gkElectronBrn),
  fMuonBranchName(Names::gkMuonBrn),
  fTrackBranchName(Names::gkTrackBrn),
  fGsfTrackBranchName(Names::gkGsfTrackBrn),
  fBarrelSuperClusterBranchName(Names::gkBarrelSuperClusterBrn),
  fEndcapSuperClusterBranchName(Names::gkEndcapSuperClusterBrn),
  fVertexBranchName(Names::gkPVBeamSpotBrn),
  fConversionBranchName(Names::gkMvfConversionBrn),
  fGoodJetsName(ModNames::gkGoodJetsName),
  fCleanElectronsName(ModNames::gkCleanElectronsName),
  fCleanMuonsName(ModNames::gkCleanMuonsName),
  fMCLeptonsName(ModNames::gkMCLeptonsName),
  fMCTausName(ModNames::gkMCTausName),
  fElFakeableObjsName(ModNames::gkElFakeableObjsName),
  fMuFakeableObjsName(ModNames::gkMuFakeableObjsName),
  fElFOType(kElFOUndef),
  fMuFOType(kMuFOUndef),
  fBarrelSuperClusters(0),
  fEndcapSuperClusters(0),
  fTracks(0),
  fGsfTracks(0),
  fVertices(0),
  fConversions(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void GenFakeableObjsMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  //------------------------------------------------------------------------------------------------
  // Request the branches (no significant time safed by not doing this)
  //------------------------------------------------------------------------------------------------
  ReqBranch(fElectronBranchName,              fElectrons);
  ReqBranch(fMuonBranchName,                  fMuons);
  ReqBranch(fTrackBranchName,                 fTracks);
  ReqBranch(fGsfTrackBranchName,              fGsfTracks);
  ReqBranch(fBarrelSuperClusterBranchName,    fBarrelSuperClusters);
  ReqBranch(fEndcapSuperClusterBranchName,    fEndcapSuperClusters);
  ReqBranch(fConversionBranchName,            fConversions);
  ReqBranch(fVertexBranchName, fVertices);

  if (fElectronFOType.CompareTo("GsfPlusSC") == 0) 
    fElFOType = kElFOGsfPlusSC;
  else if (fElectronFOType.CompareTo("Reco") == 0) 
    fElFOType = kElFOReco;
  else if (fElectronFOType.CompareTo("Loose") == 0) 
    fElFOType = kElFOLoose;
  else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified electron fakeable object %s is not defined.",
              fElectronFOType.Data());
    return;
  }

  if (fMuonFOType.CompareTo("IsoTrack") == 0) 
    fMuFOType = kMuFOIsoTrack;
  else if (fMuonFOType.CompareTo("Global") == 0) 
    fMuFOType = kMuFOGlobal;
  else if (fMuonFOType.CompareTo("TrackerMuon") == 0) 
    fMuFOType = kMuFOTrackerMuon;
  else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified muon fakeable object %s is not defined.",
              fMuonFOType.Data());
    return;
  }

}

//--------------------------------------------------------------------------------------------------
void GenFakeableObjsMod::Process()
{
  // Process entries of the tree.
  LoadBranch(fElectronBranchName);
  LoadBranch(fMuonBranchName);
  LoadBranch(fTrackBranchName);
  LoadBranch(fGsfTrackBranchName);
  LoadBranch(fBarrelSuperClusterBranchName);
  LoadBranch(fEndcapSuperClusterBranchName);
  LoadBranch(fVertexBranchName);
  LoadBranch(fConversionBranchName);

  //Load Trigger Objects
  const TriggerObjectCol *triggerObjects = GetHLTObjects(fTriggerObjectsName);
  if (!triggerObjects && fVetoTriggerJet) {
    cout << "Error: Could not load Trigger Object Collection with name " 
         << fTriggerObjectsName << endl;
  }

  // get input clean object collections
  mithep::ParticleOArr *CleanLeptons = dynamic_cast<mithep::ParticleOArr*>
    (FindObjThisEvt(ModNames::gkMergedLeptonsName));
  const ElectronCol *CleanElectrons = 0;
  if (!fCleanElectronsName.IsNull())
    CleanElectrons = GetObjThisEvt<ElectronCol>(fCleanElectronsName);
  const MuonCol *CleanMuons = 0;
  if (!fCleanMuonsName.IsNull())
    CleanMuons = GetObjThisEvt<MuonCol>(fCleanMuonsName);
  const JetCol *GoodJets = 0;
  if (!fGoodJetsName.IsNull())
    GoodJets = GetObjThisEvt<JetCol>(fGoodJetsName);

  //get input MC collections
  const MCParticleCol *GenLeptons = 0;
  if (!fMCLeptonsName.IsNull())
    GenLeptons = GetObjThisEvt<MCParticleCol>(fMCLeptonsName);
  const MCParticleCol *GenTaus = 0;
  if (!fMCTausName.IsNull())
    GenTaus = GetObjThisEvt<MCParticleCol>(fMCTausName);
  ObjArray<MCParticle> *GenLeptonsAndTaus = new ObjArray<MCParticle>;
  for (UInt_t i=0; i<GenLeptons->GetEntries(); i++)
    GenLeptonsAndTaus->Add(GenLeptons->At(i));
  for (UInt_t i=0; i<GenTaus->GetEntries(); i++)
    GenLeptonsAndTaus->Add(GenTaus->At(i));

  //Combine Barrel and Endcap superclusters into the same ObjArray
  ObjArray <SuperCluster> *SuperClusters = new ObjArray <SuperCluster>;
  for (UInt_t i=0; i<fBarrelSuperClusters->GetEntries(); i++) {
    SuperClusters->Add(fBarrelSuperClusters->At(i));
  }
  for (UInt_t i=0; i<fEndcapSuperClusters->GetEntries(); i++) {
    SuperClusters->Add(fEndcapSuperClusters->At(i));
  }

  //collections for duplicate removed electrons
  ObjArray<Electron> *DuplicateRemovedElectrons = new ObjArray<Electron>;
  std::vector<const Electron*> tmpDuplicateRemovedElectrons;

  // create final output collection
  ElectronArr *ElFakeableObjs = new ElectronArr;
  MuonArr *MuFakeableObjs = new MuonArr;


  //***********************************************************************************************
  //First do duplicate electron removal
  //***********************************************************************************************
  for (UInt_t i=0; i<fElectrons->GetEntries(); ++i) {  
    const Electron *e = fElectrons->At(i);   
    Bool_t isElectronOverlap = kFALSE;

    for (UInt_t j=0; j<tmpDuplicateRemovedElectrons.size(); ++j) {
      Double_t deltaR = MathUtils::DeltaR(tmpDuplicateRemovedElectrons[j]->Mom(), e->Mom());
      if (e->SCluster() == tmpDuplicateRemovedElectrons[j]->SCluster() ||
          e->GsfTrk() == tmpDuplicateRemovedElectrons[j]->GsfTrk() ||
          deltaR < 0.1) {
        isElectronOverlap = kTRUE;
      }
    
    
      if (isElectronOverlap) {
        if (TMath::Abs(tmpDuplicateRemovedElectrons[j]->ESuperClusterOverP() - 1) > 
            TMath::Abs(e->ESuperClusterOverP() - 1)) {	   
          tmpDuplicateRemovedElectrons[j] = e;
        }      
        break;
      }
    }
    
    if (!isElectronOverlap) {
      tmpDuplicateRemovedElectrons.push_back(e);
    }
  }
  for (UInt_t i=0; i<tmpDuplicateRemovedElectrons.size(); ++i) {  
    DuplicateRemovedElectrons->Add(tmpDuplicateRemovedElectrons[i]);
  }

  //***********************************************************************************************
  //Fakeable Objects for Electron Fakes
  //Supercluster matched to nearest isolated track.  
  //***********************************************************************************************
  if (fElFOType == kElFOGsfPlusSC) {

    std::vector<const Electron*> GsfTrackSCDenominators;

    //loop over all super clusters
    for (UInt_t i=0; i<SuperClusters->GetEntries(); i++) {
      const SuperCluster *cluster = SuperClusters->At(i);
      
      //find best matching track based on DR to the cluster 
      const Track *EOverPMatchedTrk = NULL;
      double BestEOverP = 5000.0;
      for (UInt_t j=0; j<fGsfTracks->GetEntries(); j++) {
        const Track *trk = fGsfTracks->At(j);
        
        //Use best E/P matching within dR of 0.3
        double dR = MathUtils::DeltaR(cluster->Phi(), cluster->Eta(), trk->Phi(), trk->Eta());
        Double_t EOverP = cluster->Energy() / trk->P();
        if( fabs(1-EOverP) < fabs(1-BestEOverP)  && dR < 0.3 ) { 
          BestEOverP = EOverP;
          EOverPMatchedTrk = trk;
        }
      }
      
      //****************************************************************************************
      //Use Best E/P Matching and require the track is within 0.3 of the super cluster position
      //****************************************************************************************
      if( EOverPMatchedTrk ) {
        
        //calculate track isolation around the matched track
        Double_t matchiso=0;
        matchiso = IsolationTools::TrackIsolation(EOverPMatchedTrk, 0.3, 0.015, 1.0, 0.2,
                                                  fTracks);
        
        //Veto denominators matching to real electrons      
        Bool_t IsGenLepton = false;
        for (UInt_t l=0; l<GenLeptonsAndTaus->GetEntries(); l++) {
          if (MathUtils::DeltaR(EOverPMatchedTrk->Phi(), EOverPMatchedTrk->Eta(),
                                GenLeptonsAndTaus->At(l)->Phi(), 
                                GenLeptonsAndTaus->At(l)->Eta()) < 0.1) {
            IsGenLepton = true;
          }
        }

        //Veto denominators matching to clean leptons
        Bool_t IsCleanLepton = false;
        for (UInt_t l=0; l<CleanLeptons->GetEntries(); l++) {
          if (MathUtils::DeltaR(EOverPMatchedTrk->Phi(), EOverPMatchedTrk->Eta(),
                                CleanLeptons->At(l)->Phi(), 
                                CleanLeptons->At(l)->Eta()) < 0.1) {
            IsCleanLepton = true;
          }
        }

        //Veto on Leading jet
        Bool_t IsTriggerJet = false;
        if (fVetoTriggerJet) {
          for (UInt_t l=0; l<triggerObjects->GetEntries(); l++) {      
            Double_t deltaR = MathUtils::DeltaR(EOverPMatchedTrk->Phi(), EOverPMatchedTrk->Eta(),
                                                triggerObjects->At(l)->Phi(), triggerObjects->At(l)->Eta());
            if (triggerObjects->At(l)->TrigName() == fTriggerName.Data() 
                && triggerObjects->At(l)->Type() == TriggerObject::TriggerJet
                && deltaR < 0.3
              ) {
              IsTriggerJet = true;
              break;
            }
          }
        }

        //create new electron object for the denominator under consideration
        Bool_t denominatorSaved = false;
        Electron *denominator = new Electron();
        Double_t p = TMath::Sqrt(cluster->Energy()*cluster->Energy() 
                               - denominator->Mass()*denominator->Mass());
        denominator->SetPtEtaPhi(TMath::Abs(p*TMath::Cos(EOverPMatchedTrk->Lambda())), 
                                 EOverPMatchedTrk->Eta(),EOverPMatchedTrk->Phi());
        denominator->SetGsfTrk(EOverPMatchedTrk);
        denominator->SetSuperCluster(cluster);
     
        //****************************************************************************************
        // conversion filter
        //****************************************************************************************        
        Bool_t isGoodConversion = kFALSE;
        for (UInt_t ifc=0; ifc<fConversions->GetEntries(); ifc++) {          
          Bool_t ConversionMatchFound = kFALSE;
          for (UInt_t d=0; d<fConversions->At(ifc)->NDaughters(); d++) {
            const Track *trk = dynamic_cast<const ChargedParticle*>
              (fConversions->At(ifc)->Daughter(d))->Trk();
            if (EOverPMatchedTrk == trk) {
              ConversionMatchFound = kTRUE;
              break;
            }
          }
          
          // if match between the gsftrack and one of the conversion legs
          if (ConversionMatchFound == kTRUE){
            isGoodConversion =  (fConversions->At(ifc)->Prob() > 0.0005) &&
              (fConversions->At(ifc)->Lxy() > 0) &&
              (fConversions->At(ifc)->Lz() > 0);
            
            if (isGoodConversion == kTRUE) {
              for (UInt_t d=0; d<fConversions->At(ifc)->NDaughters(); d++) {
                const Track *trk = dynamic_cast<const ChargedParticle*>
                  (fConversions->At(ifc)->Daughter(d))->Trk();
                
                if (trk) {
                  // These requirements are not used for the GSF track (d == 1)
                  if (!(trk->NHits() > 8 && trk->Prob() > 0.005) && d == 0)
                    isGoodConversion = kFALSE;
                  
                  const StableData *sd = dynamic_cast<const StableData*>
                    (fConversions->At(ifc)->DaughterDat(d));
                  if (sd->NWrongHits() != 0)
                    isGoodConversion = kFALSE;
                  
                } else {
                  isGoodConversion = kFALSE;
                }
              }
            }
          }          
          if (isGoodConversion == kTRUE) break;          
        } // loop over all conversions 

        //****************************************************************************************
        // D0 Cut        
        //****************************************************************************************        
        double d0Min = 99999;
        for(UInt_t i0 = 0; i0 < fVertices->GetEntries(); i0++) {
          double pD0 = EOverPMatchedTrk->D0Corrected(*fVertices->At(i0));
          if(TMath::Abs(pD0) < TMath::Abs(d0Min)) d0Min = TMath::Abs(pD0);
        }

        //****************************************************************************************
        // Make denominator object cuts
        //****************************************************************************************
        if( matchiso <= 10.0 && denominator->Pt() > 10.0
            && !isGoodConversion
            && d0Min < 0.025
            && !(fVetoCleanLeptons && IsCleanLepton)
            && !(fVetoGenLeptons && IsGenLepton)
            && !(fVetoTriggerJet && IsTriggerJet)
          ) {
          
          //check whether we have duplicate denominators. If yes then choose best E/P one.
          Bool_t foundDuplicate = false;
          for (UInt_t d=0; d<GsfTrackSCDenominators.size();++d) {
            if (GsfTrackSCDenominators[d]->GsfTrk() == denominator->GsfTrk()) {
              if (fabs(denominator->SCluster()->Energy()/denominator->GsfTrk()->P() - 1) 
                  < fabs(GsfTrackSCDenominators[d]->SCluster()->Energy()/
                         GsfTrackSCDenominators[d]->GsfTrk()->P() - 1)) {
                //swap this one with previous one and delete the previous one
                const Electron *denominatorToBeDeleted = GsfTrackSCDenominators[d];
                GsfTrackSCDenominators[d] = denominator;
                denominatorSaved = true;
                foundDuplicate = true;
                delete denominatorToBeDeleted;
                break;
              }
            }
          }
          if (!foundDuplicate) {
            GsfTrackSCDenominators.push_back(denominator);
            denominatorSaved = true;
          }                              
        } //end if candidate passes denominator cuts
        //delete denominator candidate object
        if (!denominatorSaved) {
          delete denominator;
        }
      }//end if track -> SC match was found      
    } //loop over SC

    //Save denominators permanently for export
    for (UInt_t d=0; d<GsfTrackSCDenominators.size() ; ++d) {
      Electron *tmpElectron = ElFakeableObjs->AddNew();
      tmpElectron->SetPtEtaPhi(GsfTrackSCDenominators[d]->Pt(), 
                               GsfTrackSCDenominators[d]->Eta(),GsfTrackSCDenominators[d]->Phi());
      tmpElectron->SetGsfTrk(GsfTrackSCDenominators[d]->GsfTrk());
      tmpElectron->SetSuperCluster(GsfTrackSCDenominators[d]->SCluster());
      delete GsfTrackSCDenominators[d];
    }
  } else if (fElFOType == kElFOReco) {
    for (UInt_t i=0; i<DuplicateRemovedElectrons->GetEntries(); i++) {  
      const Electron *denominator = DuplicateRemovedElectrons->At(i);
      Double_t denominatorIso =  
        denominator->TrackIsolation() + denominator->EcalJurassicIsolation() - 1.5;

      //Veto denominators matching to real electrons      
      Bool_t IsGenLepton = false;
      for (UInt_t l=0; l<GenLeptonsAndTaus->GetEntries(); l++) {
        if (MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                              GenLeptonsAndTaus->At(l)->Phi(), 
                              GenLeptonsAndTaus->At(l)->Eta()) < 0.1) {
          IsGenLepton = true;
        }
      }

      //Veto denominators matching to clean leptons
      Bool_t IsCleanLepton = false;
      for (UInt_t l=0; l<CleanLeptons->GetEntries(); l++) {
        if (MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                              CleanLeptons->At(l)->Phi(), 
                              CleanLeptons->At(l)->Eta()) < 0.1) {
          IsCleanLepton = true;
          }
      }
      
      //Veto on Trigger jet
      Bool_t IsTriggerJet = false;
      if (fVetoTriggerJet) {
        for (UInt_t l=0; l<triggerObjects->GetEntries(); l++) {      
          Double_t deltaR = MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                                              triggerObjects->At(l)->Phi(), triggerObjects->At(l)->Eta());
          if (triggerObjects->At(l)->TrigName() == fTriggerName.Data() 
              && triggerObjects->At(l)->Type() == TriggerObject::TriggerJet
              && deltaR < 0.3
            ) {
            IsTriggerJet = true;
            break;
          }
        }
      }

      //*************************************************************
      // conversion filter
      //*************************************************************        
      Bool_t isGoodConversion = kFALSE;
      for (UInt_t ifc=0; ifc<fConversions->GetEntries(); ifc++) {          
        Bool_t ConversionMatchFound = kFALSE;
        for (UInt_t d=0; d<fConversions->At(ifc)->NDaughters(); d++) {
          const Track *trk = dynamic_cast<const ChargedParticle*>
            (fConversions->At(ifc)->Daughter(d))->Trk();
          if (denominator->GsfTrk() == trk) {
            ConversionMatchFound = kTRUE;
            break;
          }
        }
        
        // if match between the gsftrack and one of the conversion legs
        if (ConversionMatchFound == kTRUE){
          isGoodConversion =  (fConversions->At(ifc)->Prob() > 0.0005) &&
            (fConversions->At(ifc)->Lxy() > 0) &&
            (fConversions->At(ifc)->Lz() > 0);
          
          if (isGoodConversion == kTRUE) {
            for (UInt_t d=0; d<fConversions->At(ifc)->NDaughters(); d++) {
              const Track *trk = dynamic_cast<const ChargedParticle*>
                (fConversions->At(ifc)->Daughter(d))->Trk();
              
              if (trk) {
                // These requirements are not used for the GSF track (d == 1)
                if (!(trk->NHits() > 8 && trk->Prob() > 0.005) && d == 0)
                  isGoodConversion = kFALSE;
                
                const StableData *sd = dynamic_cast<const StableData*>
                  (fConversions->At(ifc)->DaughterDat(d));
                if (sd->NWrongHits() != 0)
                  isGoodConversion = kFALSE;
                
              } else {
                isGoodConversion = kFALSE;
              }
            }
          }
        }          
        if (isGoodConversion == kTRUE) break;          
      } // loop over all conversions 
      
      //****************************************************************************************
      // D0 Cut        
      //****************************************************************************************        
      double d0Min = 99999;
      for(UInt_t i0 = 0; i0 < fVertices->GetEntries(); i0++) {
        double pD0 = denominator->GsfTrk()->D0Corrected(*fVertices->At(i0));
        if(TMath::Abs(pD0) < TMath::Abs(d0Min)) d0Min = TMath::Abs(pD0);
      }

      //****************************************************************************************
      // Apply Denominator Cuts
      //****************************************************************************************        
      if (denominator->Pt() > 10.0 && 
          denominatorIso < 10.0
          && !isGoodConversion
          && d0Min < 0.025
          && !(fVetoGenLeptons && IsGenLepton)
          && !(fVetoCleanLeptons && IsCleanLepton)
          && !(fVetoTriggerJet && IsTriggerJet)
        ) {
        Electron *tmpElectron = ElFakeableObjs->AddNew();
        tmpElectron->SetPtEtaPhi(denominator->Pt(), denominator->Eta(),denominator->Phi());
        tmpElectron->SetGsfTrk(denominator->GsfTrk());
        tmpElectron->SetSuperCluster(denominator->SCluster());        
      }
    }
  } else if (fElFOType == kElFOLoose) {
    for (UInt_t i=0; i<DuplicateRemovedElectrons->GetEntries(); i++) {  
      const Electron *denominator = DuplicateRemovedElectrons->At(i);
      Double_t denominatorIso =  denominator->TrackIsolation() + denominator->EcalJurassicIsolation() - 1.5;

      //Veto denominators matching to real electrons      
      Bool_t IsGenLepton = false;
      for (UInt_t l=0; l<GenLeptonsAndTaus->GetEntries(); l++) {
        if (MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                              GenLeptonsAndTaus->At(l)->Phi(), 
                              GenLeptonsAndTaus->At(l)->Eta()) < 0.1) {
          IsGenLepton = true;
        }
      }

      //Veto denominators matching to clean leptons
      Bool_t IsCleanLepton = false;
      for (UInt_t l=0; l<CleanLeptons->GetEntries(); l++) {
        if (MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                              CleanLeptons->At(l)->Phi(), 
                              CleanLeptons->At(l)->Eta()) < 0.1) {
          IsCleanLepton = true;
          }
      }

      //Veto on Trigger jet
      Bool_t IsTriggerJet = false;
      if (fVetoTriggerJet) {
        for (UInt_t l=0; l<triggerObjects->GetEntries(); l++) {      
          Double_t deltaR = MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                                              triggerObjects->At(l)->Phi(), triggerObjects->At(l)->Eta());
          if (triggerObjects->At(l)->TrigName() == fTriggerName.Data() 
              && triggerObjects->At(l)->Type() == TriggerObject::TriggerJet
              && deltaR < 0.3
            ) {
            IsTriggerJet = true;
            break;
          }
        }
      }

      //****************************************************************************************
      // conversion filter
      //****************************************************************************************        
      Bool_t isGoodConversion = kFALSE;
      for (UInt_t ifc=0; ifc<fConversions->GetEntries(); ifc++) {          
        Bool_t ConversionMatchFound = kFALSE;
        for (UInt_t d=0; d<fConversions->At(ifc)->NDaughters(); d++) {
          const Track *trk = dynamic_cast<const ChargedParticle*>
            (fConversions->At(ifc)->Daughter(d))->Trk();
          if (denominator->GsfTrk() == trk) {
            ConversionMatchFound = kTRUE;
            break;
          }
        }
        
        // if match between the gsftrack and one of the conversion legs
        if (ConversionMatchFound == kTRUE){
          isGoodConversion =  (fConversions->At(ifc)->Prob() > 0.0005) &&
            (fConversions->At(ifc)->Lxy() > 0) &&
            (fConversions->At(ifc)->Lz() > 0);
          
          if (isGoodConversion == kTRUE) {
            for (UInt_t d=0; d<fConversions->At(ifc)->NDaughters(); d++) {
              const Track *trk = dynamic_cast<const ChargedParticle*>
                (fConversions->At(ifc)->Daughter(d))->Trk();
              
              if (trk) {
                // These requirements are not used for the GSF track (d == 1)
                if (!(trk->NHits() > 8 && trk->Prob() > 0.005) && d == 0)
                  isGoodConversion = kFALSE;
                
                const StableData *sd = dynamic_cast<const StableData*>
                  (fConversions->At(ifc)->DaughterDat(d));
                if (sd->NWrongHits() != 0)
                  isGoodConversion = kFALSE;
                
              } else {
                isGoodConversion = kFALSE;
              }
            }
          }
        }          
        if (isGoodConversion == kTRUE) break;          
      } // loop over all conversions 
      
      //****************************************************************************************
      // D0 Cut        
      //****************************************************************************************        
      double d0Min = 99999;
      for(UInt_t i0 = 0; i0 < fVertices->GetEntries(); i0++) {
        double pD0 = denominator->GsfTrk()->D0Corrected(*fVertices->At(i0));
        if(TMath::Abs(pD0) < TMath::Abs(d0Min)) d0Min = TMath::Abs(pD0);
      }

      if (denominator->Pt() > 10.0 && denominatorIso < 10.0
          && !isGoodConversion
          && d0Min < 0.025
          && denominator->PassLooseID()
          && !(fVetoGenLeptons && IsGenLepton)
          && !(fVetoCleanLeptons && IsCleanLepton)
          && !(fVetoTriggerJet && IsTriggerJet)
        ) {
        Electron *tmpElectron = ElFakeableObjs->AddNew();
        tmpElectron->SetPtEtaPhi(denominator->Pt(), denominator->Eta(),denominator->Phi());
        tmpElectron->SetGsfTrk(denominator->GsfTrk());
        tmpElectron->SetSuperCluster(denominator->SCluster());
      }
    }
  }

  //***********************************************************************************************
  //Fakeable Objects for Muon Fakes
  //***********************************************************************************************
  if (fMuFOType == kMuFOIsoTrack) {
    for (UInt_t i=0; i<fTracks->GetEntries(); i++) {
      const Track *track = fTracks->At(i);
      Double_t trackIsolation = IsolationTools::TrackIsolation(track, 0.4, 0.015, 1.0, 
                                                               0.2, fTracks);
      //Determine if muon fakeable object matches a gen lepton
      Bool_t IsGenLepton = false;
      for (UInt_t l=0; l<GenLeptonsAndTaus->GetEntries(); l++) {
        if (MathUtils::DeltaR(track->Phi(), track->Eta(),
                              GenLeptonsAndTaus->At(l)->Phi(), 
                              GenLeptonsAndTaus->At(l)->Eta()) < 0.3) {
          IsGenLepton = true;
        }
      }
      
      //Veto denominators matching to clean leptons
      Bool_t IsCleanLepton = false;
      for (UInt_t l=0; l<CleanLeptons->GetEntries(); l++) {
        if (MathUtils::DeltaR(track->Phi(), track->Eta(),
                              CleanLeptons->At(l)->Phi(), 
                              CleanLeptons->At(l)->Eta()) < 0.1) {
          IsCleanLepton = true;
          }
      }

      //Veto on Trigger jet
      Bool_t IsTriggerJet = false;
      if (fVetoTriggerJet) {
        for (UInt_t l=0; l<triggerObjects->GetEntries(); l++) {      
          Double_t deltaR = MathUtils::DeltaR(track->Phi(), track->Eta(),
                                              triggerObjects->At(l)->Phi(), 
                                              triggerObjects->At(l)->Eta());
          if (triggerObjects->At(l)->TrigName() == fTriggerName.Data() 
              && triggerObjects->At(l)->Type() == TriggerObject::TriggerJet
              && deltaR < 0.3
            ) {
            IsTriggerJet = true;
            break;
          }
        }
      }
      
      //****************************************************************************************
      // D0 Cut        
      //****************************************************************************************
      double d0Min = 99999;
      for(UInt_t i0 = 0; i0 < fVertices->GetEntries(); i0++) {
        double pD0 = track->D0Corrected(*fVertices->At(i0));
        if(TMath::Abs(pD0) < TMath::Abs(d0Min)) d0Min = TMath::Abs(pD0);
      }


      //define denominator cuts
      if (track->Pt() > 10.0 && trackIsolation < 10.0
          && d0Min < 0.025
          && !(fVetoGenLeptons && IsGenLepton)  
          && !(fVetoCleanLeptons && IsCleanLepton)  
          && !(fVetoTriggerJet && IsTriggerJet)
        ) {
        //add to fakeable objects
        Muon* tmpMuon = MuFakeableObjs->AddNew();
        tmpMuon->SetTrackerTrk(track);
      }
    }
  } else if (fMuFOType == kMuFOGlobal) {
    for (UInt_t i=0; i<fMuons->GetEntries(); i++) {
      const Muon *denominator = fMuons->At(i);
      Double_t totalIsolation = denominator->IsoR03SumPt() + denominator->IsoR03EmEt() + 
        denominator->IsoR03HadEt();

      //Determine if muon fakeable object matches a gen lepton
      Bool_t IsGenLepton = false;
      for (UInt_t l=0; l<GenLeptonsAndTaus->GetEntries(); l++) {
        if (MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                              GenLeptonsAndTaus->At(l)->Phi(), 
                              GenLeptonsAndTaus->At(l)->Eta()) < 0.3) {
          IsGenLepton = true;
        }
      }
      
      //Veto denominators matching to clean leptons
      Bool_t IsCleanLepton = false;
      for (UInt_t l=0; l<CleanLeptons->GetEntries(); l++) {
        if (MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                              CleanLeptons->At(l)->Phi(), 
                              CleanLeptons->At(l)->Eta()) < 0.1) {
          IsCleanLepton = true;
          }
      }

      //Veto on Trigger jet
      Bool_t IsTriggerJet = false;
      if (fVetoTriggerJet) {
        for (UInt_t l=0; l<triggerObjects->GetEntries(); l++) {      
          Double_t deltaR = MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                                              triggerObjects->At(l)->Phi(), triggerObjects->At(l)->Eta());
          if (triggerObjects->At(l)->TrigName() == fTriggerName.Data() 
              && triggerObjects->At(l)->Type() == TriggerObject::TriggerJet
              && deltaR < 0.3
            ) {
            IsTriggerJet = true;
            break;
          }
        }
      }

      //****************************************************************************************
      // D0 Cut        
      //****************************************************************************************        
      double d0Min = 99999;
      for(UInt_t i0 = 0; i0 < fVertices->GetEntries(); i0++) {
        double pD0 = denominator->Trk()->D0Corrected(*fVertices->At(i0));
        if(TMath::Abs(pD0) < TMath::Abs(d0Min)) d0Min = TMath::Abs(pD0);
      }
      
      if (denominator->Pt() > 10.0 && totalIsolation < 10.0 && denominator->HasGlobalTrk()
          && d0Min < 0.025
          && !(fVetoGenLeptons && IsGenLepton)
          && !(fVetoCleanLeptons && IsCleanLepton)
          && !(fVetoTriggerJet && IsTriggerJet)
        ) {
        Muon* tmpMuon = MuFakeableObjs->AddNew();
        tmpMuon->SetGlobalTrk(denominator->GlobalTrk());
        tmpMuon->SetTrackerTrk(denominator->TrackerTrk());
      }
    }
  } else if (fMuFOType == kMuFOTrackerMuon) {
    for (UInt_t i=0; i<fMuons->GetEntries(); i++) {
      const Muon *denominator = fMuons->At(i);
      Double_t totalIsolation = denominator->IsoR03SumPt() + denominator->IsoR03EmEt() + denominator->IsoR03HadEt();

      //Determine if muon fakeable object matches a gen lepton
      Bool_t IsGenLepton = false;
      for (UInt_t l=0; l<GenLeptonsAndTaus->GetEntries(); l++) {
        if (MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                              GenLeptonsAndTaus->At(l)->Phi(), 
                              GenLeptonsAndTaus->At(l)->Eta()) < 0.3) {
          IsGenLepton = true;
        }
      }

      //Veto denominators matching to clean leptons
      Bool_t IsCleanLepton = false;
      for (UInt_t l=0; l<CleanLeptons->GetEntries(); l++) {
        if (MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                              CleanLeptons->At(l)->Phi(), 
                              CleanLeptons->At(l)->Eta()) < 0.1) {
          IsCleanLepton = true;
          }
      }

      //Veto on Trigger jet
      Bool_t IsTriggerJet = false;
      if (fVetoTriggerJet) {
        for (UInt_t l=0; l<triggerObjects->GetEntries(); l++) {      
          Double_t deltaR = MathUtils::DeltaR(denominator->Phi(), denominator->Eta(),
                                              triggerObjects->At(l)->Phi(), triggerObjects->At(l)->Eta());
          if (triggerObjects->At(l)->TrigName() == fTriggerName.Data() 
              && triggerObjects->At(l)->Type() == TriggerObject::TriggerJet
              && deltaR < 0.3
            ) {
            IsTriggerJet = true;
            break;
          }
        }
      }

      //****************************************************************************************
      // D0 Cut        
      //****************************************************************************************        
      double d0Min = 99999;
      for(UInt_t i0 = 0; i0 < fVertices->GetEntries(); i0++) {
        double pD0 = denominator->Trk()->D0Corrected(*fVertices->At(i0));
        if(TMath::Abs(pD0) < TMath::Abs(d0Min)) d0Min = TMath::Abs(pD0);
      }

      if (denominator->Pt() > 10.0 && totalIsolation < 10.0 && denominator->HasTrackerTrk()
          && d0Min < 0.025
          && !(fVetoGenLeptons && IsGenLepton)
          && !(fVetoCleanLeptons && IsCleanLepton)
          && !(fVetoTriggerJet && IsTriggerJet)
        ) {
        Muon* tmpMuon = MuFakeableObjs->AddNew();
        tmpMuon->SetTrackerTrk(denominator->TrackerTrk());
      }
    }
  }

  //***********************************************************************************************
  //Export the fakeable object collections for other modules to use
  //***********************************************************************************************
  ElFakeableObjs->SetName(fElFakeableObjsName);
  AddObjThisEvt(ElFakeableObjs);
  MuFakeableObjs->SetName(fMuFakeableObjsName);
  AddObjThisEvt(MuFakeableObjs);

  delete GenLeptonsAndTaus;
  delete SuperClusters;
  delete DuplicateRemovedElectrons;
}
