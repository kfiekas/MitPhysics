// $Id: ElectronTools.cc,v 1.15 2010/10/09 20:15:33 bendavid Exp $

#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include <TFile.h>

ClassImp(mithep::ElectronTools)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
ElectronTools::ElectronTools()  
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassCustomID(const Electron *ele, EElIdType idType) {

  Double_t  fCuts[6][8];             //!custom id cuts

  Double_t tightcuts[6][8]={
    {0.086, 0.1, 0.052, 0.0, 0.050, 0.059, 0.061, 0.0},   //hovere
    {0.011, 0.011, 0.011, 0.0, 0.033, 0.029, 0.030, 0.0},     //sigmaetaeta
    {0.038, 0.024, 0.045, 0.0, 0.034, 0.017, 0.026, 0.0},   //deltaphiin
    {0.0081, 0.0029, 0.0051, 0.0, 0.0070, 0.0062, 0.0088, 0.0}, //deltaetain
    {0.0,    0.9,    0.0,    0.0, 0.0,    0.78,   0.0,    0.0},             //eoverp
    {0.8,0.2,0.9,0,0,0,0,0}};                              //extra cuts fbrem and E_Over_P 

  Double_t loosecuts[6][8]={
    {0.076, 0.033, 0.07, 0.0, 0.083,0.148, 0.033, 0.0},         //hovere
    {0.0101, 0.0095, 0.0097, 0.0, 0.03, 0.03, 0.03, 0.0},       //sigmaetaeta
    {0.053, 0.0189, 0.059, 0.099, 0.0278,0.0157, 0.042, 0.080}, //deltaphiin
    {0.0078, 0.00259, 0.0062, 0.0, 0.0078,0.0061, 0.0061, 0.0}, //deltaetain
    {0.3, 0.92, 0.211, 0.0, 0.42, 0.88, 0.68, 0.0},             //eoverp
    {0.8,0.2,0,0,0,0,0,0}};                                     //extra cuts fbrem and E_Over_P 

  Double_t VBTFWorkingPoint95[6][8] = {
    {0.15,  0.15,  0.15,  0.15,  0.07,   0.07,   0.07,   0.07  }, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03  }, //sigmaetaeta
    {0.8,   0.8,   0.8,   0.8,   0.7,    0.7,    0.7,    0.7   }, //deltaphiin
    {0.007, 0.007, 0.007, 0.007, 0.010,  0.010,  0.010,  0.010 }, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0   }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0     }  //extra cuts fbrem and E_Over_P 
  };            

  Double_t VBTFWorkingPoint90[6][8] = {
    {0.12,  0.12,  0.12,  0.12,  0.05,   0.05,   0.05,   0.05  }, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03  }, //sigmaetaeta
    {0.8,   0.8,   0.8,   0.8,   0.7,    0.7,    0.7,    0.7   }, //deltaphiin
    {0.007, 0.007, 0.007, 0.007, 0.009,  0.009,  0.009,  0.009 }, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0   }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0     }  //extra cuts fbrem and E_Over_P 
  };            

  Double_t VBTFWorkingPoint85[6][8] = {
    {0.04,  0.04,  0.04,  0.04,  0.025,  0.025,  0.025,  0.025 }, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03  }, //sigmaetaeta
    {0.06,  0.06,  0.06,  0.06,  0.04,   0.04,   0.04,   0.04  }, //deltaphiin
    {0.006, 0.006, 0.006, 0.006, 0.007,  0.007,  0.007,  0.007 }, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0   }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0     }  //extra cuts fbrem and E_Over_P 
  };            

  Double_t VBTFWorkingPoint80[6][8] = {
    {0.04,  0.04,  0.04,  0.04,  0.025,  0.025,  0.025,  0.025}, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03 }, //sigmaetaeta
    {0.06,  0.06,  0.06,  0.06,  0.03,   0.03,   0.03,   0.03 }, //deltaphiin
    {0.004, 0.004, 0.004, 0.004, 0.007,  0.007,  0.007,  0.007}, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0  }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0    }  //extra cuts fbrem and E_Over_P 
  };            

  Double_t VBTFWorkingPoint70[6][8] = {
    {0.025, 0.025, 0.025, 0.025, 0.012,  0.012,  0.012,  0.012}, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03 }, //sigmaetaeta
    {0.03,  0.03,  0.03,  0.03,  0.02,   0.02,   0.02,   0.02 }, //deltaphiin
    {0.004, 0.004, 0.004, 0.004, 0.005,  0.005,  0.005,  0.005}, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0  }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0    }  //extra cuts fbrem and E_Over_P 
  };            

  switch (idType) {
    case kCustomIdTight:    
      memcpy(fCuts,tightcuts,sizeof(fCuts));
      break;
    case kCustomIdLoose:
      memcpy(fCuts,loosecuts,sizeof(fCuts));
      break;
    case kVBTFWorkingPoint95Id:
      memcpy(fCuts,VBTFWorkingPoint95,sizeof(fCuts));
      break;
    case kVBTFWorkingPoint90Id:
      memcpy(fCuts,VBTFWorkingPoint90,sizeof(fCuts));
      break;
    case kVBTFWorkingPoint85Id:
      memcpy(fCuts,VBTFWorkingPoint85,sizeof(fCuts));
      break;
    case kVBTFWorkingPoint80Id:
      memcpy(fCuts,VBTFWorkingPoint80,sizeof(fCuts));
      break;
    case kVBTFWorkingPoint70Id:
      memcpy(fCuts,VBTFWorkingPoint70,sizeof(fCuts));
      break;
    default:
      memset(fCuts,0,sizeof(fCuts));
      break;
  }


  // Based on RecoEgamma/ElectronIdentification/src/CutBasedElectronID.cc.
  Double_t eOverP = ele->ESuperClusterOverP();
  Double_t fBrem  = ele->FBrem();

  if ( (fCuts[5][0]>0.0) && (eOverP < fCuts[5][0]) && (fCuts[5][1]>0.0) && (fBrem < fCuts[5][1]))
    return kFALSE;

  if ( (fCuts[5][2]>0.0) && (eOverP < fCuts[5][2]*(1-fBrem)))
    return kFALSE;

  Int_t cat = 2;
  if ((ele->IsEB() && fBrem<0.06) || (ele->IsEE() && fBrem<0.1)) 
    cat=1;
  else if (eOverP < 1.2 && eOverP > 0.8) 
    cat=0;
  
  if(ele->SCluster() == 0)
    return kFALSE;
  Double_t eSeedOverPin = ele->ESeedClusterOverPIn(); 
  Double_t hOverE       = ele->HadronicOverEm();
  Double_t sigmaee      = ele->CoviEtaiEta();
  Double_t deltaPhiIn   = TMath::Abs(ele->DeltaPhiSuperClusterTrackAtVtx());
  Double_t deltaEtaIn   = TMath::Abs(ele->DeltaEtaSuperClusterTrackAtVtx());

  Int_t eb = 1;
  if (ele->IsEB()) 
    eb = 0;
 
  if (hOverE>fCuts[0][cat+4*eb])
    return kFALSE;

  if (sigmaee>fCuts[1][cat+4*eb])
    return kFALSE;

  if (deltaPhiIn>fCuts[2][cat+4*eb])
    return kFALSE; 

  if(deltaEtaIn>fCuts[3][cat+4*eb])
    return kFALSE; 

  if(eSeedOverPin<fCuts[4][cat+4*eb])
    return kFALSE;

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassCustomIso(const Electron *ele, EElIsoType isoType,
                                    Bool_t useCombineIso) 
{
  Bool_t pass = kTRUE;
  Double_t fIsoCuts[4][2];          //!custom isolation cuts
  Double_t VBTFWorkingPoint95[4][2] = {
    {0.15 , 0.08   },   //TrkIso
    {2.00 , 0.06   },   //ECALIso
    {0.12 , 0.05   },   //HCALIso
    {0.15,  0.10   }   //Combined    
  };            

  Double_t VBTFWorkingPoint90[4][2] = {
    {0.12 , 0.05   },   //TrkIso
    {0.09 , 0.06   },   //ECALIso
    {0.10 , 0.03   },   //HCALIso
    {0.10,  0.07   }   //Combined    
  };            

  Double_t VBTFWorkingPoint85[4][2] = {
    {0.09 , 0.05   },   //TrkIso
    {0.08 , 0.05   },   //ECALIso
    {0.10 , 0.025  },   //HCALIso
    {0.09,  0.06   }   //Combined    
  };            

  Double_t VBTFWorkingPoint80[4][2] = {
    {0.09 , 0.04   },   //TrkIso
    {0.07 , 0.05   },   //ECALIso
    {0.10 , 0.025  },   //HCALIso
    {0.07,  0.06   }   //Combined    
  };            

  Double_t VBTFWorkingPoint70[4][2] = {
    {0.05 , 0.025  },   //TrkIso
    {0.06 , 0.025  },   //ECALIso
    {0.03 , 0.020  },   //HCALIso
    {0.04,  0.030  }   //Combined
  };            

  switch (isoType) {
    case kVBTFWorkingPoint95Iso:
      memcpy(fIsoCuts,VBTFWorkingPoint95,sizeof(fIsoCuts));
      break;
    case kVBTFWorkingPoint90Iso:
      memcpy(fIsoCuts,VBTFWorkingPoint90,sizeof(fIsoCuts));
      break;
    case kVBTFWorkingPoint85Iso:
      memcpy(fIsoCuts,VBTFWorkingPoint85,sizeof(fIsoCuts));
      break;
    case kVBTFWorkingPoint80Iso:
      memcpy(fIsoCuts,VBTFWorkingPoint80,sizeof(fIsoCuts));
      break;
    case kVBTFWorkingPoint70Iso:
      memcpy(fIsoCuts,VBTFWorkingPoint70,sizeof(fIsoCuts));
      break;
    default:
      memset(fIsoCuts,0,sizeof(fIsoCuts));
      break;
  }

  Double_t trkIso  = ele->TrackIsolationDr03() / ele->Pt();
  Double_t ecalIso = ele->EcalRecHitIsoDr03() / ele->Pt();
  Double_t hcalIso = ele->HcalTowerSumEtDr03() / ele->Pt();
  Double_t combinedIso = ele->TrackIsolationDr03() + ele->EcalRecHitIsoDr03() + ele->HcalTowerSumEtDr03();
  if(ele->IsEB()) combinedIso = ele->TrackIsolationDr03() + TMath::Max(ele->EcalRecHitIsoDr03() - 1.0, 0.0) + ele->HcalTowerSumEtDr03();
  combinedIso = combinedIso / ele->Pt();

  Int_t eb = 1;
  if (ele->IsEB()) 
    eb = 0;
 
  if(useCombineIso == kFALSE){
    if (trkIso>fIsoCuts[0][eb])
      pass = kFALSE;
    if (ecalIso>fIsoCuts[1][eb])
      pass = kFALSE;
    if (hcalIso>fIsoCuts[2][eb])
      pass = kFALSE;
  }
  else {
    if (combinedIso>fIsoCuts[3][eb])
      pass = kFALSE;
  }

  return pass;
}


//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassConversionFilter(const Electron *ele, 
                                           const DecayParticleCol *conversions, 
                                           Bool_t WrongHitsRequirement) 
{
  Bool_t isGoodConversion = kFALSE;

  for (UInt_t ifc=0; ifc<conversions->GetEntries(); ifc++) {
    Bool_t ConversionMatchFound = kFALSE;
    for (UInt_t d=0; d<conversions->At(ifc)->NDaughters(); d++) {
      const Track *trk = dynamic_cast<const ChargedParticle*>
        (conversions->At(ifc)->Daughter(d))->Trk();
      if (ele->GsfTrk() == trk) {
        ConversionMatchFound = kTRUE;
        break;
      }
    }

    // if match between the e-track and one of the conversion legs
    if (ConversionMatchFound == kTRUE){
      isGoodConversion =  (conversions->At(ifc)->Prob() > 1e-6) &&
        (conversions->At(ifc)->Lxy() > 0) &&
        (conversions->At(ifc)->Lz() > 0) &&
        (conversions->At(ifc)->Position().Rho() > 2.0);

      if (isGoodConversion == kTRUE) {
        for (UInt_t d=0; d<conversions->At(ifc)->NDaughters(); d++) {
          const Track *trk = dynamic_cast<const ChargedParticle*>
            (conversions->At(ifc)->Daughter(d))->Trk();
          if (trk) {
            // These requirements are not used for the GSF track
            if (!(trk->NHits() >= 3 && trk->Prob() > 1e-6) && trk!=ele->GsfTrk())
              isGoodConversion = kFALSE;
            const StableData *sd = dynamic_cast<const StableData*>
              (conversions->At(ifc)->DaughterDat(d));
            if (WrongHitsRequirement && sd->NWrongHits() != 0)
              isGoodConversion = kFALSE;
          } else {
            isGoodConversion = kFALSE;
          }
        }
      }
    }

    if (isGoodConversion == kTRUE) break;
    
  } // loop over all conversions 

  return !isGoodConversion;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassD0Cut(const Electron *ele, const VertexCol *vertices, Double_t fD0Cut, 
                                Bool_t fReverseD0Cut) 
{
  Bool_t d0cut = kFALSE;
  // d0 cut
  Double_t d0_real = 99999;
  for(UInt_t i0 = 0; i0 < vertices->GetEntries(); i0++) {
    if(vertices->At(i0)->NTracks() > 0){
      Double_t pD0 = ele->GsfTrk()->D0Corrected(*vertices->At(i0));
      d0_real = TMath::Abs(pD0);
      break;
    }
  }
  if(d0_real < fD0Cut) d0cut = kTRUE;
  
  if     (fReverseD0Cut == kTRUE &&
          d0cut == kFALSE && d0_real < 0.05)
    d0cut = kTRUE;
  else if(fReverseD0Cut == kTRUE)
    d0cut = kFALSE;
  
  return d0cut;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassD0Cut(const Electron *ele, const BeamSpotCol *beamspots, Double_t fD0Cut, 
                                Bool_t fReverseD0Cut) 
{
  Bool_t d0cut = kFALSE;
  // d0 cut
  Double_t d0_real = 99999;
  for(UInt_t i0 = 0; i0 < beamspots->GetEntries(); i0++) {
    Double_t pD0 = ele->GsfTrk()->D0Corrected(*beamspots->At(i0));
    if(TMath::Abs(pD0) < TMath::Abs(d0_real)) d0_real = TMath::Abs(pD0);
  }
  if(d0_real < fD0Cut) d0cut = kTRUE;
  
  if     (fReverseD0Cut == kTRUE &&
          d0cut == kFALSE && d0_real < 0.05)
    d0cut = kTRUE;
  else if(fReverseD0Cut == kTRUE)
    d0cut = kFALSE;
  
  return d0cut;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassChargeFilter(const Electron *ele) 
{
  Bool_t passChargeFilter = kTRUE;
  if(ele->TrackerTrk() &&
     ele->TrackerTrk()->Charge() != ele->Charge()) passChargeFilter = kFALSE;

  return passChargeFilter;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassSpikeRemovalFilter(const Electron *ele) 
{
  Bool_t passSpikeRemovalFilter = kTRUE;
  if(ele->SCluster() &&
     ele->SCluster()->Seed()->Energy() > 5.0 && 
     ele->SCluster()->Seed()->EMax() / ele->SCluster()->Seed()->E3x3() > 0.95
    ) {
    passSpikeRemovalFilter = kFALSE;
  }

  // For Now Only use the EMax/E3x3 prescription.
  //   if(ele->SCluster()->Seed()->Energy() > 5.0 && 
  //      (1 - (ele->SCluster()->Seed()->E1x3() + ele->SCluster()->Seed()->E3x1() - 2*ele->SCluster()->Seed()->EMax())) > 0.95
  //     ) {
  //     passSpikeRemovalFilter = kFALSE;
  //   }
    
  return passSpikeRemovalFilter;
}

Bool_t ElectronTools::PassTriggerMatching(const Electron *ele, const TriggerObjectCol *trigobjs)
{
  
  for (UInt_t i=0; i<trigobjs->GetEntries(); ++i) {
    const TriggerObject *trigobj = trigobjs->At(i);
    if (trigobj->TriggerType()==TriggerObject::TriggerCluster || trigobj->TriggerType()==TriggerObject::TriggerElectron || trigobj->TriggerType()==TriggerObject::TriggerPhoton) {
      if (MathUtils::DeltaR(ele->SCluster(),trigobj)<0.3) {
        return kTRUE;
      }
    }
  }
  
  return kFALSE;
  
  
}

//--------------------------------------------------------------------------------------------------
Int_t ElectronTools::Classify(const Electron *ele) {
  
  double eta    = ele->AbsEta();
  double eOverP = ele->ESuperClusterOverP();
  double fBrem  = ele->FBrem();

  int cat = -1;
  if (ele->IsEB()) {
    if ((fBrem >= 0.12) and (eOverP > 0.9) and (eOverP < 1.2))
      cat = 0;
    else if (((eta >  .445   and eta <  .45  ) or
  	      (eta >  .79    and eta <  .81  ) or
  	      (eta > 1.137   and eta < 1.157 ) or
  	      (eta > 1.47285 and eta < 1.4744)))
      cat = 6;
    else if (ele->IsTrackerDriven() and !ele->IsEcalDriven())
      cat = 8;
    else if (fBrem < 0.12)
      cat = 1;
    else
      cat = 2;
  } else {
    if ((fBrem >= 0.2) and (eOverP > 0.82) and (eOverP < 1.22))
      cat = 3;
    else if (eta > 1.5 and eta <  1.58)
      cat = 7;
    else if (ele->IsTrackerDriven() and !ele->IsEcalDriven())
      cat = 8;
    else if (fBrem < 0.2)
      cat = 4;
    else
      cat = 5;
  }

  return cat;
}

//--------------------------------------------------------------------------------------------------
Int_t ElectronTools::PassTightId(const Electron *ele, const VertexCol *vertices, 
                                 const DecayParticleCol *conversions, const Int_t typeCuts){

  Double_t scEt   = ele->SCluster()->Et(); 
  Double_t eOverP = ele->ESuperClusterOverP();

  Double_t fBrem = ele->FBrem();
  Double_t hOverE = ele->HadronicOverEm();
  Double_t sigmaee = ele->CoviEtaiEta();
  Double_t deltaPhiIn = TMath::Abs(ele->DeltaPhiSuperClusterTrackAtVtx());
  Double_t deltaEtaIn = TMath::Abs(ele->DeltaEtaSuperClusterTrackAtVtx());
  Double_t eSeedOverPin = ele->ESeedClusterOverPIn(); 

  Int_t mishits = ele->BestTrk()->NExpectedHitsInner();
  Double_t tkIso   = ele->TrackIsolationDr03();
  Double_t ecalIso = ele->EcalRecHitIsoDr04();
  Double_t ecalIsoPed = (ele->IsEB())?std::max(0.,ecalIso-1.):ecalIso;
  Double_t hcalIso  = ele->HcalTowerSumEtDr04();
  Double_t hcalIso1 = ele->HcalDepth1TowerSumEtDr04();
  Double_t hcalIso2 = ele->HcalDepth2TowerSumEtDr04();
  Double_t e25Max = ele->E25Max();
  Double_t e15 = ele->E15();
  Double_t e55 = ele->E55();
  Double_t e25Maxoe55 = e25Max/e55;
  Double_t e15oe55 = e15/e55;

  int cat = Classify(ele);
  int eb;

  if (ele->IsEB()) 
    eb = 0;
  else 
    eb = 1; 

  Int_t result = 0.;
  
  Int_t bin = 0;
  
  Double_t WantBin = kFALSE;
  if(WantBin) {
    if (scEt < 20.)
      bin = 2;
    else if (scEt > 30.)
      bin = 0;
    else
      bin = 1;
  }

  if (fBrem > 0)
    eSeedOverPin = eSeedOverPin + fBrem;

  if(typeCuts == 3 || typeCuts == 4){
  // Loose robust cuts
    Double_t robustlooseEleIDCuts[52] = {
       0.05,  0.0103,  0.8, 0.00688, -1, -1, 7.33, 4.68, 9999., 9999., 9999.,9999., 9999.,9999., 9999., 9999.,9999., 9999., 0.000, -9999., 9999., 9999., 9999, -1, 0, 0,
       0.0389, 0.0307, 0.7, 0.00944, -1, -1, 7.76, 3.09,  2.23, 9999., 9999.,9999., 9999.,9999., 9999., 9999.,9999., 9999., 0.000, -9999., 9999., 9999., 9999, -1, 0, 0
    };
  // Tight robust cuts
    Double_t robusttightEleIDCuts[52] = {
        0.0201, 0.0102, 0.0211, 0.00606, -1, -1, 2.34, 3.24, 4.51, 9999., 9999.,9999., 9999.,9999., 9999., 9999.,9999., 9999., 0.000, -9999., 9999., 9999., 9999, -1, 0, 0,
        0.00253, 0.0291, 0.022, 0.0032, -1, -1, 0.826, 2.7, 0.255, 9999., 9999.,9999., 9999.,9999., 9999., 9999.,9999., 9999., 0.000, -9999., 9999., 9999., 9999, -1, 0, 0
    };                            
    Double_t cut[52];
    if     (typeCuts == 3) {
      memcpy(cut  ,robustlooseEleIDCuts  ,sizeof(cut));
    }
    else {
      memcpy(cut  ,robusttightEleIDCuts  ,sizeof(cut));
    }

    if ((tkIso > cut[26*eb+6]) || (ecalIso > cut[26*eb+7]) || (hcalIso > cut[26*eb+8]) || (hcalIso1 > cut[26*eb+9]) || (hcalIso2 > cut[26*eb+10]) || 
        (tkIso/ele->Pt() > cut[26*eb+11]) || (ecalIso/ele->Pt() > cut[26*eb+12]) || (hcalIso/ele->Pt() > cut[26*eb+13]) ||
        ((tkIso+ecalIso+hcalIso)>cut[26*eb+14]) || (((tkIso+ecalIso+hcalIso)/ ele->Pt()) > cut[26*eb+15]) || 
        ((tkIso+ecalIsoPed+hcalIso)>cut[26*eb+16]) || (((tkIso+ecalIsoPed+hcalIso)/ ele->Pt()) > cut[26*eb+17])  )
      result = 0.;
    else
      result = 2.;

    if (hOverE > cut[26*eb+0]) 
      return result;    

    if (sigmaee > cut[26*eb+1]) 
      return result;    

    if (fabs(deltaPhiIn) > cut[26*eb+2]) 
      return result;    

    if (fabs(deltaEtaIn) > cut[26*eb+3]) 
      return result;    
    
    if (e25Maxoe55 < cut[26*eb+4] && e15oe55 < cut[26*eb+5])
      return result;
    // some extra electron id cuts
    if (sigmaee < cut[26*eb+18]) // inverted sigmaee cut - spike removal related
      return result;

    if (  eOverP < cut[26*eb+19] ||  eOverP > cut[26*eb+20]) // lower and upper cut in E/P
      return result;
    
    result = result + 1;

    Bool_t passD0cut = PassD0Cut(ele, vertices, cut[26*eb+21], kFALSE);
    if (!passD0cut)
      return result;

    if (mishits > cut[26*eb+22]) // expected missing hits
      return result;

    Bool_t passConvVeto = PassConversionFilter(ele, conversions, kTRUE);
    if (passConvVeto == kFALSE)
      return result;
      
    result += 4;
  }
  else if(typeCuts == 0 || typeCuts == 1 || typeCuts == 2){
    // Loose cuts
  Double_t cutdcotdistLoose[9] = {9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
  };
  Double_t cutdetainLoose[9] = {1.30e-02, 5.95e-03, 3.10e-02, 1.68e-02, 8.44e-03, 1.70e-02, 1.55e-02, 5.13e-02, 1.61e-02
  };
  Double_t cutdphiinLoose[9] = {7.51e-02, 3.30e-01, 4.20e-01, 9.86e-02, 2.84e-01, 3.28e-01, 3.77e-01, 4.32e-01, 3.74e-01
  };
  Double_t cuteseedopcorLoose[9] = {6.31e-01, 3.02e-01, 3.04e-01, 8.10e-01, 2.23e-01, 5.03e-01, 2.78e-01, 3.10e-01, 4.69e-01
  };
  Double_t cutetLoose[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.
  };
  Double_t cutfmishitsLoose[9] = {4.50e+00, 1.50e+00, 1.50e+00, 2.50e+00, 2.50e+00, 1.50e+00, 2.50e+00, 4.50e+00, 5.00e-01
  };
  Double_t cuthoeLoose[9] = {2.47e-01, 7.78e-02, 1.49e-01, 3.82e-01, 4.70e-02, 1.12e-01, 1.16e+00, 5.04e+00, 1.35e+00
  };
  Double_t cutip_gsfLoose[9] = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02
  };
  Double_t cutiso_sumLoose[9] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
  };
  Double_t cutiso_sumoetLoose[9] = {9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
  };
  Double_t cutseeLoose[9] = {1.92e-02, 1.31e-02, 2.53e-02, 5.27e-02, 3.29e-02, 4.19e-02, 2.65e-02, 6.58e-02, 1.38e-01
  };

    // Medium cuts
  Double_t cutdcotdistMedium[9] = {9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
  };
  Double_t cutdetainMedium[9] = {1.19e-02, 4.20e-03, 1.07e-02, 1.49e-02, 6.56e-03, 1.19e-02, 1.16e-02, 5.13e-02, 6.37e-03
  };
  Double_t cutdphiinMedium[9] = {7.51e-02, 2.93e-01, 3.58e-01, 9.53e-02, 1.62e-01, 2.99e-01, 2.76e-01, 4.32e-01, 2.57e-01
  };
  Double_t cuteseedopcorMedium[9] = {6.31e-01, 8.14e-01, 7.60e-01, 8.18e-01, 7.56e-01, 5.35e-01, 6.20e-01, 7.88e-01, 8.85e-01
  };
  Double_t cutetMedium[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.
  };
  Double_t cutfmishitsMedium[9] = {1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 2.50e+00, 1.50e+00, 5.00e-01
  };
  Double_t cuthoeMedium[9] = {2.46e-01, 6.80e-02, 1.35e-01, 3.73e-01, 2.33e-02, 5.58e-02, 8.80e-01, 5.04e+00, 3.78e-02
  };
  Double_t cutip_gsfMedium[9] = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02
  };
  Double_t cutiso_sumMedium[9] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
  };
  Double_t cutiso_sumoetMedium[9] = {9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
  };
  Double_t cutseeMedium[9] = {1.92e-02, 1.13e-02, 1.47e-02, 3.84e-02, 3.05e-02, 3.36e-02, 1.35e-02, 5.05e-02, 2.79e-02
  };

    // Tight cuts
  Double_t cutdcotdistTight[9] = {9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
  };
  Double_t cutdetainTight[9] = {9.28e-03, 3.56e-03, 7.16e-03, 1.31e-02, 5.81e-03, 9.79e-03, 1.15e-02, 1.66e-02, 3.19e-03
  };
  Double_t cutdphiinTight[9] = {4.66e-02, 7.80e-02, 2.64e-01, 4.42e-02, 3.20e-02, 2.37e-01, 8.25e-02, 2.07e-01, 5.39e-02
  };
  Double_t cuteseedopcorTight[9] = {6.48e-01, 8.97e-01, 8.91e-01, 8.39e-01, 8.35e-01, 6.49e-01, 6.76e-01, 8.70e-01, 9.91e-01
  };
  Double_t cutetTight[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.
  };
  Double_t cutfmishitsTight[9] = {1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 5.00e-01, 2.50e+00, 5.00e-01, 5.00e-01
  };
  Double_t cuthoeTight[9] = {9.94e-02, 5.61e-02, 1.05e-01, 9.73e-02, 1.81e-02, 3.06e-02, 5.57e-01, 5.04e+00, 1.06e-03
  };
  Double_t cutip_gsfTight[9] = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02
  };
  Double_t cutiso_sumTight[9] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
  };
  Double_t cutiso_sumoetTight[9] = {9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.
  };
  Double_t cutseeTight[9] = {1.56e-02, 1.07e-02, 1.23e-02, 3.35e-02, 2.98e-02, 3.06e-02, 1.07e-02, 3.79e-02, 1.01e-02
  };

    Double_t cutdcotdist[9];
    Double_t cutdetain[9];
    Double_t cutdphiin[9];
    Double_t cuteseedopcor[9];
    Double_t cutet[9];
    Double_t cutfmishits[9];
    Double_t cuthoe[9];
    Double_t cutip_gsf[9];
    Double_t cutiso_sum[9];
    Double_t cutiso_sumoet[9];
    Double_t cutsee[9];
    if     (typeCuts == 0) {
      memcpy(cutdcotdist  ,cutdcotdistLoose  ,sizeof(cutdcotdistLoose));
      memcpy(cutdetain    ,cutdetainLoose	,sizeof(cutdetainLoose));
      memcpy(cutdphiin    ,cutdphiinLoose	,sizeof(cutdphiinLoose));
      memcpy(cuteseedopcor,cuteseedopcorLoose,sizeof(cuteseedopcorLoose));
      memcpy(cutet        ,cutetLoose	,sizeof(cutetLoose));
      memcpy(cutfmishits  ,cutfmishitsLoose  ,sizeof(cutfmishitsLoose));
      memcpy(cuthoe       ,cuthoeLoose	,sizeof(cuthoeLoose));
      memcpy(cutip_gsf    ,cutip_gsfLoose	,sizeof(cutip_gsfLoose));
      memcpy(cutiso_sum   ,cutiso_sumLoose	,sizeof(cutiso_sumLoose));
      memcpy(cutiso_sumoet,cutiso_sumoetLoose,sizeof(cutiso_sumoetLoose));
      memcpy(cutsee       ,cutseeLoose	,sizeof(cutseeLoose));
    }
    else if(typeCuts == 1) {
      memcpy(cutdcotdist  ,cutdcotdistMedium  ,sizeof(cutdcotdistMedium));
      memcpy(cutdetain    ,cutdetainMedium	,sizeof(cutdetainMedium));
      memcpy(cutdphiin    ,cutdphiinMedium	,sizeof(cutdphiinMedium));
      memcpy(cuteseedopcor,cuteseedopcorMedium,sizeof(cuteseedopcorMedium));
      memcpy(cutet        ,cutetMedium	,sizeof(cutetMedium));
      memcpy(cutfmishits  ,cutfmishitsMedium  ,sizeof(cutfmishitsMedium));
      memcpy(cuthoe       ,cuthoeMedium	,sizeof(cuthoeMedium));
      memcpy(cutip_gsf    ,cutip_gsfMedium	,sizeof(cutip_gsfMedium));
      memcpy(cutiso_sum   ,cutiso_sumMedium	,sizeof(cutiso_sumMedium));
      memcpy(cutiso_sumoet,cutiso_sumoetMedium,sizeof(cutiso_sumoetMedium));
      memcpy(cutsee       ,cutseeMedium	,sizeof(cutseeMedium));
    }
    else {
      memcpy(cutdcotdist  ,cutdcotdistTight  ,sizeof(cutdcotdistTight));
      memcpy(cutdetain    ,cutdetainTight	,sizeof(cutdetainTight));
      memcpy(cutdphiin    ,cutdphiinTight	,sizeof(cutdphiinTight));
      memcpy(cuteseedopcor,cuteseedopcorTight,sizeof(cuteseedopcorTight));
      memcpy(cutet        ,cutetTight	,sizeof(cutetTight));
      memcpy(cutfmishits  ,cutfmishitsTight  ,sizeof(cutfmishitsTight));
      memcpy(cuthoe       ,cuthoeTight	,sizeof(cuthoeTight));
      memcpy(cutip_gsf    ,cutip_gsfTight	,sizeof(cutip_gsfTight));
      memcpy(cutiso_sum   ,cutiso_sumTight	,sizeof(cutiso_sumTight));
      memcpy(cutiso_sumoet,cutiso_sumoetTight,sizeof(cutiso_sumoetTight));
      memcpy(cutsee       ,cutseeTight	,sizeof(cutseeTight));
    }

    // CAREFUL, I HAVE COMMENTED OUT WHAT SANI IS DOING
    //Double_t iso_sum = tkIso + ecalIso + hcalIso;
    Double_t iso_sum = (ele->TrackIsolationDr03() + TMath::Max(ele->EcalRecHitIsoDr03() - 1.0, 0.0) + 
                        ele->HcalTowerSumEtDr03()) / ele->Pt();
    Double_t iso_sum_corrected = iso_sum*pow(40./scEt, 2);
    if ((iso_sum < cutiso_sum[cat+bin*9]) &&
	(iso_sum_corrected < cutiso_sumoet[cat+bin*9]))
      result += 2;

    if (fBrem > -2) {
      if ((hOverE < cuthoe[cat+bin*9]) and
    	  (sigmaee < cutsee[cat+bin*9]) and
    	  (fabs(deltaPhiIn) < cutdphiin[cat+bin*9]) and
    	  (fabs(deltaEtaIn) < cutdetain[cat+bin*9]) and
    	  (eSeedOverPin > cuteseedopcor[cat+bin*9]) and
    	  (scEt > cutet[cat+bin*9]))
	result += 1.;
    }

    Bool_t passD0cut = PassD0Cut(ele, vertices, cutip_gsf[cat+bin*9], kFALSE);
    if (passD0cut)
      result += 4;

    Bool_t passConvVeto = PassConversionFilter(ele, conversions, kTRUE);
    if (mishits < cutfmishits[cat+bin*9] &&
	passConvVeto)
      result += 8;
  } // classbased

  return result;
}
