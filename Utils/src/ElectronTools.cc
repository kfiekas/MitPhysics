// $Id: ElectronTools.cc,v 1.19 2011/02/17 13:44:55 bendavid Exp $

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
                                           const BaseVertex *vtx,
                                           UInt_t nWrongHitsMax,
                                           Double_t probMin,
                                           Double_t lxyMin,
                                           Bool_t matchCkf,
                                           Bool_t requireArbitratedMerged) 
{
  Bool_t isGoodConversion = kFALSE;

  for (UInt_t ifc=0; ifc<conversions->GetEntries(); ifc++) {
    Bool_t ConversionMatchFound = kFALSE;
    for (UInt_t d=0; d<conversions->At(ifc)->NDaughters(); d++) {
      const Track *trk = dynamic_cast<const ChargedParticle*>
        (conversions->At(ifc)->Daughter(d))->Trk();
      if (ele->GsfTrk() == trk || (matchCkf && ele->TrackerTrk()==trk) ) {
        ConversionMatchFound = kTRUE;
        break;
      }
    }

    // if match between the e-track and one of the conversion legs
    if (ConversionMatchFound == kTRUE){
      isGoodConversion =  (conversions->At(ifc)->Prob() > probMin) &&
        (!requireArbitratedMerged || conversions->At(ifc)->Quality().Quality(ConversionQuality::arbitratedMerged)) &&
        (conversions->At(ifc)->LxyCorrected(vtx) > lxyMin);

      if (isGoodConversion == kTRUE) {
        for (UInt_t d=0; d<conversions->At(ifc)->NDaughters(); d++) {
          const Track *trk = dynamic_cast<const ChargedParticle*>
            (conversions->At(ifc)->Daughter(d))->Trk();
          if (trk) {
            const StableData *sd = dynamic_cast<const StableData*>
              (conversions->At(ifc)->DaughterDat(d));
            if (sd->NWrongHits() > nWrongHitsMax)
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
Bool_t ElectronTools::PassD0Cut(const Electron *ele, const VertexCol *vertices, Double_t fD0Cut) 
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
  
  return d0cut;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronTools::PassD0Cut(const Electron *ele, const BeamSpotCol *beamspots, Double_t fD0Cut) 
{
  Bool_t d0cut = kFALSE;
  // d0 cut
  Double_t d0_real = 99999;
  for(UInt_t i0 = 0; i0 < beamspots->GetEntries(); i0++) {
    Double_t pD0 = ele->GsfTrk()->D0Corrected(*beamspots->At(i0));
    if(TMath::Abs(pD0) < TMath::Abs(d0_real)) d0_real = TMath::Abs(pD0);
  }
  if(d0_real < fD0Cut) d0cut = kTRUE;
  
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

// original code on
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/RecoEgamma/ElectronIdentification/src/CutBasedElectronID.cc

  Double_t scEt   = ele->SCluster()->Et();
  Double_t scEta  = ele->SCluster()->Eta();

  Double_t fBrem = ele->FBrem();
  Double_t hOverE = ele->HadronicOverEm();
  Double_t sigmaee = ele->CoviEtaiEta();
  Double_t deltaPhiIn = TMath::Abs(ele->DeltaPhiSuperClusterTrackAtVtx());
  Double_t deltaEtaIn = TMath::Abs(ele->DeltaEtaSuperClusterTrackAtVtx());
  Double_t eSeedOverPin = ele->ESeedClusterOverPIn(); 

  Int_t mishits = ele->BestTrk()->NExpectedHitsInner();
  Double_t tkIso   = ele->TrackIsolationDr03();
  Double_t ecalIso = ele->EcalRecHitIsoDr04();
  Double_t hcalIso  = ele->HcalTowerSumEtDr04();

  int cat = Classify(ele);
  int eb;

  if (ele->IsEB()) 
    eb = 0;
  else 
    eb = 1; 

  // Medium cuts
  Double_t cutdcotdistMedium[9] = {
  3.32e-02, 2.92e-02, 2.49e-02, 3.92e-02, 3.41e-02, 3.96e-02, 2.91e-02, 3.95e-02, 7.71e-03};
  Double_t cutdetainMedium[9] = {
  1.33e-02, 4.48e-03, 9.22e-03, 1.54e-02, 7.26e-03, 1.24e-02, 1.29e-02, 3.84e-02, 1.88e-02};
  Double_t cutdetainlMedium[9] = {
  1.21e-02, 4.22e-03, 9.18e-03, 1.61e-02, 6.45e-03, 1.16e-02, 1.23e-02, 6.20e-02, 2.43e-02};
  Double_t cutdphiinMedium[9] = {
  7.09e-02, 2.43e-01, 2.96e-01, 7.98e-02, 2.35e-01, 2.76e-01, 3.42e-01, 4.04e-01, 2.99e-01};
  Double_t cutdphiinlMedium[9] = {
  7.42e-02, 2.43e-01, 2.97e-01, 9.12e-02, 2.26e-01, 2.76e-01, 3.34e-01, 5.58e-01, 2.91e-01};
  Double_t cuteseedopcorMedium[9] = {
  6.42e-01, 9.44e-01, 4.53e-01, 7.62e-01, 3.67e-01, 5.57e-01, 1.98e-01, 9.15e-01, 6.28e-02};
  Double_t cutfmishitsMedium[9] = {
  4.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01};
  Double_t cuthoeMedium[9] = {
  1.96e-01, 6.30e-02, 1.48e-01, 3.66e-01, 5.66e-02, 1.45e-01, 4.29e-01, 4.28e-01, 3.99e-01};
  Double_t cuthoelMedium[9] = {
  2.19e-01, 6.19e-02, 1.47e-01, 3.58e-01, 4.61e-02, 1.46e-01, 3.26e-01, 3.81e-01, 3.89e-01};
  Double_t cutip_gsfMedium[9] = {
  2.45e-02, 9.74e-02, 1.48e-01, 5.49e-02, 5.65e-01, 3.33e-01, 2.04e-01, 5.41e-01, 1.21e-01};
  Double_t cutip_gsflMedium[9] = {
  1.92e-02, 9.81e-02, 1.33e-01, 4.34e-02, 5.65e-01, 3.24e-01, 2.33e-01, 4.30e-01, 6.44e-02};
  Double_t cutiso_sumMedium[9] = {
  1.44e+01, 1.12e+01, 1.09e+01, 1.08e+01, 6.35e+00, 9.78e+00, 1.30e+01, 1.62e+01, 1.96e+00};
  Double_t cutiso_sumoetMedium[9] = {
  1.01e+01, 6.41e+00, 6.00e+00, 8.14e+00, 3.90e+00, 4.76e+00, 6.86e+00, 6.48e+00, 1.74e+01};
  Double_t cutiso_sumoetlMedium[9] = {
  9.44e+00, 7.67e+00, 7.15e+00, 7.34e+00, 3.35e+00, 4.70e+00, 8.32e+00, 7.55e+00, 6.25e+00};
  Double_t cutseeMedium[9] = {
  1.30e-02, 1.09e-02, 1.18e-02, 3.94e-02, 3.04e-02, 3.28e-02, 1.00e-02, 3.73e-02, 6.69e-02};
  Double_t cutseelMedium[9] = {
  1.42e-02, 1.11e-02, 1.29e-02, 4.32e-02, 2.96e-02, 3.82e-02, 1.01e-02, 4.45e-02, 1.19e-01};

  // Tight cuts
  Double_t cutdcotdistTight[9] = {
  2.68e-02, 2.36e-02, 2.21e-02, 3.72e-02, 3.17e-02, 3.61e-02, 2.55e-02, 3.75e-02, 2.16e-04};
  Double_t cutdetainTight[9] = {
  8.92e-03, 3.96e-03, 8.50e-03, 1.34e-02, 6.27e-03, 1.05e-02, 1.12e-02, 3.09e-02, 1.88e-02};
  Double_t cutdetainlTight[9] = {
  9.23e-03, 3.77e-03, 8.70e-03, 1.39e-02, 5.60e-03, 9.40e-03, 1.07e-02, 6.20e-02, 4.10e-03};
  Double_t cutdphiinTight[9] = {
  6.37e-02, 1.53e-01, 2.90e-01, 7.69e-02, 1.81e-01, 2.34e-01, 3.42e-01, 3.93e-01, 2.84e-01};
  Double_t cutdphiinlTight[9] = {
  6.92e-02, 2.33e-01, 2.96e-01, 8.65e-02, 1.85e-01, 2.76e-01, 3.34e-01, 3.53e-01, 2.90e-01};
  Double_t cuteseedopcorTight[9] = {
  6.52e-01, 9.69e-01, 9.12e-01, 7.79e-01, 3.67e-01, 6.99e-01, 3.28e-01, 9.67e-01, 5.89e-01};
  Double_t cutfmishitsTight[9] = {
  4.50e+00, 1.50e+00, 5.00e-01, 1.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
  Double_t cuthoeTight[9] = {
  1.74e-01, 4.88e-02, 1.46e-01, 3.64e-01, 4.93e-02, 1.45e-01, 4.29e-01, 4.20e-01, 3.99e-01};
  Double_t cuthoelTight[9] = {
  2.19e-01, 5.25e-02, 1.47e-01, 3.57e-01, 4.25e-02, 1.45e-01, 3.26e-01, 3.80e-01, 1.32e-01};
  Double_t cutip_gsfTight[9] = {
  1.58e-02, 8.25e-02, 1.15e-01, 4.05e-02, 5.40e-01, 1.51e-01, 7.74e-02, 4.17e-01, 7.80e-02};
  Double_t cutip_gsflTight[9] = {
  1.27e-02, 6.26e-02, 9.68e-02, 3.02e-02, 5.65e-01, 1.46e-01, 7.90e-02, 4.10e-01, 4.79e-02};
  Double_t cutiso_sumTight[9] = {
  1.23e+01, 9.77e+00, 1.01e+01, 9.77e+00, 6.13e+00, 7.55e+00, 1.30e+01, 1.62e+01, 1.78e+00};
  Double_t cutiso_sumoetTight[9] = {
  7.75e+00, 5.45e+00, 5.67e+00, 5.97e+00, 3.17e+00, 3.86e+00, 6.06e+00, 5.31e+00, 1.05e+01};
  Double_t cutiso_sumoetlTight[9] = {
  7.56e+00, 5.08e+00, 5.77e+00, 5.74e+00, 2.37e+00, 3.32e+00, 4.97e+00, 5.46e+00, 3.82e+00};
  Double_t cutseeTight[9] = {
  1.16e-02, 1.07e-02, 1.08e-02, 3.49e-02, 2.89e-02, 3.08e-02, 9.87e-03, 3.37e-02, 4.40e-02};
  Double_t cutseelTight[9] = {
  1.27e-02, 1.08e-02, 1.13e-02, 4.19e-02, 2.81e-02, 3.02e-02, 9.76e-03, 4.28e-02, 2.98e-02};
  
  Double_t cutdcotdistSuperTight[9] = {
  2.11e-02, 1.86e-02, 1.55e-02, 3.40e-02, 2.85e-02, 3.32e-02, 1.64e-02, 3.75e-02, 1.30e-04};
  Double_t cutdetainSuperTight[9] = {
  7.84e-03, 3.67e-03, 7.00e-03, 1.28e-02, 5.65e-03, 9.53e-03, 1.08e-02, 2.97e-02, 7.24e-03};
  Double_t cutdetainlSuperTight[9] = {
  7.61e-03, 3.28e-03, 6.57e-03, 1.03e-02, 5.05e-03, 8.55e-03, 1.07e-02, 2.94e-02, 4.10e-03};
  Double_t cutdphiinSuperTight[9] = {
  4.83e-02, 7.39e-02, 2.38e-01, 5.74e-02, 1.29e-01, 2.13e-01, 3.31e-01, 3.93e-01, 2.84e-01};
  Double_t cutdphiinlSuperTight[9] = {
  5.79e-02, 7.21e-02, 2.18e-01, 7.70e-02, 1.41e-01, 2.11e-01, 2.43e-01, 3.53e-01, 2.89e-01};
  Double_t cuteseedopcorSuperTight[9] = {
  7.32e-01, 9.77e-01, 9.83e-01, 8.55e-01, 4.31e-01, 7.35e-01, 4.18e-01, 9.99e-01, 5.89e-01};
  Double_t cutfmishitsSuperTight[9] = {
  3.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
  Double_t cuthoeSuperTight[9] = {
  9.19e-02, 4.11e-02, 1.42e-01, 3.35e-01, 3.82e-02, 1.41e-01, 4.29e-01, 4.01e-01, 3.99e-01};
  Double_t cuthoelSuperTight[9] = {
  7.51e-02, 3.81e-02, 1.41e-01, 3.32e-01, 3.10e-02, 1.43e-01, 2.35e-01, 3.80e-01, 1.32e-01};
  Double_t cutip_gsfSuperTight[9] = {
  1.42e-02, 2.66e-02, 1.06e-01, 3.38e-02, 3.23e-01, 1.07e-01, 7.74e-02, 2.32e-01, 7.80e-02};
  Double_t cutip_gsflSuperTight[9] = {
  1.15e-02, 2.72e-02, 8.41e-02, 2.49e-02, 4.17e-01, 1.02e-01, 7.90e-02, 1.69e-01, 4.79e-02};
  Double_t cutiso_sumSuperTight[9] = {
  8.95e+00, 8.18e+00, 8.75e+00, 7.47e+00, 5.43e+00, 5.87e+00, 8.16e+00, 1.02e+01, 1.78e+00};
  Double_t cutiso_sumoetSuperTight[9] = {
  6.45e+00, 5.14e+00, 4.99e+00, 5.21e+00, 2.65e+00, 3.12e+00, 4.52e+00, 4.72e+00, 3.68e+00};
  Double_t cutiso_sumoetlSuperTight[9] = {
  6.02e+00, 3.96e+00, 4.23e+00, 4.73e+00, 1.99e+00, 2.64e+00, 3.72e+00, 3.81e+00, 1.44e+00};
  Double_t cutseeSuperTight[9] = {
  1.09e-02, 1.05e-02, 1.05e-02, 3.24e-02, 2.81e-02, 2.95e-02, 9.77e-03, 2.75e-02, 2.95e-02};
  Double_t cutseelSuperTight[9] = {
  1.12e-02, 1.05e-02, 1.07e-02, 3.51e-02, 2.75e-02, 2.87e-02, 9.59e-03, 2.67e-02, 2.98e-02};

  Double_t cutdcotdist[9];
  Double_t cutdetain[9];
  Double_t cutdetainl[9];
  Double_t cutdphiin[9];
  Double_t cutdphiinl[9];
  Double_t cuteseedopcor[9];
  Double_t cutfmishits[9];
  Double_t cuthoe[9];
  Double_t cuthoel[9];
  Double_t cutip_gsf[9];
  Double_t cutip_gsfl[9];
  Double_t cutiso_sum[9];
  Double_t cutiso_sumoet[9];
  Double_t cutiso_sumoetl[9];
  Double_t cutsee[9];
  Double_t cutseel[9];
  if	 (typeCuts == 0) {
    memcpy(cutdcotdist   ,cutdcotdistMedium   ,sizeof(cutdcotdistMedium));
    memcpy(cutdetain     ,cutdetainMedium     ,sizeof(cutdetainMedium));
    memcpy(cutdetainl    ,cutdetainlMedium    ,sizeof(cutdetainlMedium));
    memcpy(cutdphiin     ,cutdphiinMedium     ,sizeof(cutdphiinMedium));
    memcpy(cutdphiinl    ,cutdphiinlMedium    ,sizeof(cutdphiinlMedium));
    memcpy(cuteseedopcor ,cuteseedopcorMedium ,sizeof(cuteseedopcorMedium));
    memcpy(cutfmishits   ,cutfmishitsMedium   ,sizeof(cutfmishitsMedium));
    memcpy(cuthoe        ,cuthoeMedium	      ,sizeof(cuthoeMedium));
    memcpy(cuthoel       ,cuthoelMedium	      ,sizeof(cuthoelMedium));
    memcpy(cutip_gsf     ,cutip_gsfMedium     ,sizeof(cutip_gsfMedium));
    memcpy(cutip_gsfl    ,cutip_gsflMedium    ,sizeof(cutip_gsflMedium));
    memcpy(cutiso_sum    ,cutiso_sumMedium    ,sizeof(cutiso_sumMedium));
    memcpy(cutiso_sumoet ,cutiso_sumoetMedium ,sizeof(cutiso_sumoetMedium));
    memcpy(cutiso_sumoetl,cutiso_sumoetlMedium,sizeof(cutiso_sumoetlMedium));
    memcpy(cutsee        ,cutseeMedium	      ,sizeof(cutseeMedium));
    memcpy(cutseel       ,cutseelMedium	      ,sizeof(cutseelMedium));
  }
  else if(typeCuts == 1) {
    memcpy(cutdcotdist   ,cutdcotdistTight   ,sizeof(cutdcotdistTight));
    memcpy(cutdetain     ,cutdetainTight     ,sizeof(cutdetainTight));
    memcpy(cutdetainl    ,cutdetainlTight    ,sizeof(cutdetainlTight));
    memcpy(cutdphiin     ,cutdphiinTight     ,sizeof(cutdphiinTight));
    memcpy(cutdphiinl    ,cutdphiinlTight    ,sizeof(cutdphiinlTight));
    memcpy(cuteseedopcor ,cuteseedopcorTight ,sizeof(cuteseedopcorTight));
    memcpy(cutfmishits   ,cutfmishitsTight   ,sizeof(cutfmishitsTight));
    memcpy(cuthoe        ,cuthoeTight	     ,sizeof(cuthoeTight));
    memcpy(cuthoel       ,cuthoelTight	     ,sizeof(cuthoelTight));
    memcpy(cutip_gsf     ,cutip_gsfTight     ,sizeof(cutip_gsfTight));
    memcpy(cutip_gsfl    ,cutip_gsflTight    ,sizeof(cutip_gsflTight));
    memcpy(cutiso_sum    ,cutiso_sumTight    ,sizeof(cutiso_sumTight));
    memcpy(cutiso_sumoet ,cutiso_sumoetTight ,sizeof(cutiso_sumoetTight));
    memcpy(cutiso_sumoetl,cutiso_sumoetlTight,sizeof(cutiso_sumoetlTight));
    memcpy(cutsee        ,cutseeTight	     ,sizeof(cutseeTight));
    memcpy(cutseel       ,cutseelTight	     ,sizeof(cutseelTight));
  }
  else {
    memcpy(cutdcotdist   ,cutdcotdistSuperTight   ,sizeof(cutdcotdistSuperTight));
    memcpy(cutdetain     ,cutdetainSuperTight     ,sizeof(cutdetainSuperTight));
    memcpy(cutdetainl    ,cutdetainlSuperTight    ,sizeof(cutdetainlSuperTight));
    memcpy(cutdphiin     ,cutdphiinSuperTight     ,sizeof(cutdphiinSuperTight));
    memcpy(cutdphiinl    ,cutdphiinlSuperTight    ,sizeof(cutdphiinlSuperTight));
    memcpy(cuteseedopcor ,cuteseedopcorSuperTight ,sizeof(cuteseedopcorSuperTight));
    memcpy(cutfmishits   ,cutfmishitsSuperTight   ,sizeof(cutfmishitsSuperTight));
    memcpy(cuthoe        ,cuthoeSuperTight	  ,sizeof(cuthoeSuperTight));
    memcpy(cuthoel       ,cuthoelSuperTight	  ,sizeof(cuthoelSuperTight));
    memcpy(cutip_gsf     ,cutip_gsfSuperTight     ,sizeof(cutip_gsfSuperTight));
    memcpy(cutip_gsfl    ,cutip_gsflSuperTight    ,sizeof(cutip_gsflSuperTight));
    memcpy(cutiso_sum    ,cutiso_sumSuperTight    ,sizeof(cutiso_sumSuperTight));
    memcpy(cutiso_sumoet ,cutiso_sumoetSuperTight ,sizeof(cutiso_sumoetSuperTight));
    memcpy(cutiso_sumoetl,cutiso_sumoetlSuperTight,sizeof(cutiso_sumoetlSuperTight));
    memcpy(cutsee        ,cutseeSuperTight	  ,sizeof(cutseeSuperTight));
    memcpy(cutseel       ,cutseelSuperTight	  ,sizeof(cutseelSuperTight));
  }

  int result = 0;
  
  const int ncuts = 10;
  std::vector<bool> cut_results(ncuts, false);
  
  float iso_sum = tkIso + ecalIso + hcalIso;
  if(fabs(scEta)>1.5) 
    iso_sum += (fabs(scEta)-1.5)*1.09;
  
  float iso_sumoet = iso_sum*(40./scEt);
  
  float eseedopincor = eSeedOverPin + fBrem;
  if(fBrem < 0)
    eseedopincor = eSeedOverPin;

  float dist = (TMath::Abs(ele->ConvPartnerDist())      == -9999.? 9999:TMath::Abs(ele->ConvPartnerDist()));
  float dcot = (TMath::Abs(ele->ConvPartnerDCotTheta()) == -9999.? 9999:TMath::Abs(ele->ConvPartnerDCotTheta()));

  float dcotdistcomb = ((0.04 - std::max(dist, dcot)) > 0?(0.04 - std::max(dist, dcot)):0);

  Double_t ip = 99999;
  for(UInt_t i0 = 0; i0 < vertices->GetEntries(); i0++) {
    if(vertices->At(i0)->NTracks() > 0){
      Double_t pD0 = ele->GsfTrk()->D0Corrected(*vertices->At(i0));
      ip = TMath::Abs(pD0);
      break;
    }
  }

  for (int cut=0; cut<ncuts; cut++) {
    switch (cut) {
    case 0:
      cut_results[cut] = compute_cut(fabs(deltaEtaIn), scEt, cutdetainl[cat], cutdetain[cat]);
      break;
    case 1:
      cut_results[cut] = compute_cut(fabs(deltaPhiIn), scEt, cutdphiinl[cat], cutdphiin[cat]);
      break;
    case 2:
      cut_results[cut] = (eseedopincor > cuteseedopcor[cat]);
      break;
    case 3:
      cut_results[cut] = compute_cut(hOverE, scEt, cuthoel[cat], cuthoe[cat]);
      break;
    case 4:
      cut_results[cut] = compute_cut(sigmaee, scEt, cutseel[cat], cutsee[cat]);
      break;
    case 5:
      cut_results[cut] = compute_cut(iso_sumoet, scEt, cutiso_sumoetl[cat], cutiso_sumoet[cat]);
      break;
    case 6:
      cut_results[cut] = (iso_sum < cutiso_sum[cat]);
      break;
    case 7:
      cut_results[cut] = compute_cut(fabs(ip), scEt, cutip_gsfl[cat], cutip_gsf[cat]);
      break;
    case 8:
      cut_results[cut] = (mishits < cutfmishits[cat]);
      break;
    case 9:
      cut_results[cut] = (dcotdistcomb < cutdcotdist[cat]);
      break;
    }
  }
  
  // ID part
  if (cut_results[0] & cut_results[1] & cut_results[2] & cut_results[3] & cut_results[4])
    result = result + 1;
  
  // ISO part
  if (cut_results[5] & cut_results[6])
    result = result + 2;
  
  // IP part
  if (cut_results[7])
    result = result + 8;
  
  // Conversion part
  if (cut_results[8] and cut_results[9])
    result = result + 4;
  
  return result;
}

//--------------------------------------------------------------------------------------------------
bool ElectronTools::compute_cut(double x, double et, double cut_min, double cut_max, bool gtn) {

  float et_min = 10;
  float et_max = 40;

  bool accept = false;
  float cut = cut_max; //  the cut at et=40 GeV

  if(et < et_max) {
    cut = cut_min + (1/et_min - 1/et)*(cut_max - cut_min)/(1/et_min - 1/et_max);
  } 
  
  if(et < et_min) {
    cut = cut_min;
  } 

  if(gtn) {   // useful for e/p cut which is gt
    accept = (x >= cut);
  } 
  else {
    accept = (x <= cut);
  }

  return accept;
}
