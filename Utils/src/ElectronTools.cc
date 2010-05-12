// $Id: ElectronTools.cc,v 1.3 2010/05/03 11:35:17 bendavid Exp $

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
    {0.05,  0.05,  0.05,  0.05,  0.04,   0.04,   0.04,   0.04  }, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03  }, //sigmaetaeta
    {0.8,   0.8,   0.8,   0.8,   0.7,    0.7,    0.7,    0.7   }, //deltaphiin
    {0.006, 0.006, 0.006, 0.006, 0.008,  0.008,  0.008,  0.008 }, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0   }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0     }  //extra cuts fbrem and E_Over_P 
  };            


  Double_t VBTFWorkingPoint90[6][8] = {
    {0.05,  0.05,  0.05,  0.05,  0.025,  0.025,  0.025,  0.025 }, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03  }, //sigmaetaeta
    {0.04,  0.04,  0.04,  0.04,  0.025,  0.025,  0.025,  0.025 }, //deltaphiin
    {0.006, 0.006, 0.006, 0.006, 0.008,  0.008,  0.008,  0.008 }, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0   }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0     }  //extra cuts fbrem and E_Over_P 
  };            

  Double_t VBTFWorkingPoint80[6][8] = {
    {0.05,  0.05,  0.05,  0.05,  0.025,  0.025,  0.025,  0.025}, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03 }, //sigmaetaeta
    {0.02,  0.02,  0.02,  0.02,  0.02,   0.02,   0.02,   0.02 }, //deltaphiin
    {0.006, 0.006, 0.006, 0.006, 0.006,  0.006,  0.006,  0.006}, //deltaetain
    {0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,    0.0  }, //eoverp
    {0.0,   0.0,   0,     0,     0,      0,      0,      0    }  //extra cuts fbrem and E_Over_P 
  };            

  Double_t VBTFWorkingPoint70[6][8] = {
    {0.02,  0.02,  0.02,  0.02,  0.025,  0.025,  0.025,  0.025}, //hovere
    {0.01,  0.01,  0.01,  0.01,  0.03,   0.03,   0.03,   0.03 }, //sigmaetaeta
    {0.02,  0.02,  0.02,  0.02,  0.02,   0.02,   0.02,   0.02 }, //deltaphiin
    {0.006, 0.006, 0.006, 0.006, 0.003,  0.003,  0.003,  0.003}, //deltaetain
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

  if ((eOverP < fCuts[5][0]) && (fBrem < fCuts[5][1]))
    return kFALSE;

  if (eOverP < fCuts[5][2]*(1-fBrem))
    return kFALSE;

  Int_t cat = 2;
  if ((ele->IsEB() && fBrem<0.06) || (ele->IsEE() && fBrem<0.1)) 
    cat=1;
  else if (eOverP < 1.2 && eOverP > 0.8) 
    cat=0;

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
Bool_t ElectronTools::PassCustomIso(const Electron *ele, EElIsoType isoType) 
{
  Bool_t pass = kTRUE;
  Double_t fIsoCuts[4][2];          //!custom isolation cuts

  Double_t VBTFWorkingPoint95[4][2] = {
    {7.0 ,  8.0   },   //TrkIso
    {5.0 ,  3.0   },   //ECALIso
    {5.0 ,  2.0   },   //HCALIso
    {9999,  9999  }   //Combined    
  };            

  Double_t VBTFWorkingPoint90[4][2] = {
    {6.0 ,  6.0   },   //TrkIso
    {5.0 ,  2.5   },   //ECALIso
    {5.0 ,  1.5   },   //HCALIso
    {9999,  9999  }   //Combined    
  };            

  Double_t VBTFWorkingPoint80[4][2] = {
    {3.0 ,  1.5   },   //TrkIso
    {4.0 ,  2.5   },   //ECALIso
    {5.0 ,  0.7   },   //HCALIso
    {9999,  9999  }   //Combined    
  };            

  Double_t VBTFWorkingPoint70[4][2] = {
    {2.5 ,  0.8   },   //TrkIso
    {3.0 ,  2.5   },   //ECALIso
    {5.0 ,  0.25  },   //HCALIso
    {9999,  9999  }   //Combined    
  };            

  switch (isoType) {
    case kVBTFWorkingPoint95Iso:
      memcpy(fIsoCuts,VBTFWorkingPoint95,sizeof(fIsoCuts));
      break;
    case kVBTFWorkingPoint90Iso:
      memcpy(fIsoCuts,VBTFWorkingPoint90,sizeof(fIsoCuts));
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

  Double_t trkIso  = ele->TrackIsolationDr03();
  Double_t ecalIso = ele->EcalRecHitIsoDr04();
  Double_t hcalIso = ele->HcalIsolation();
  Double_t combinedIso = ele->TrackIsolationDr03() + ele->EcalRecHitIsoDr04() - 1.5;

  Int_t eb = 1;
  if (ele->IsEB()) 
    eb = 0;
 
  if (trkIso>fIsoCuts[0][eb])
    pass = kFALSE;
  if (ecalIso>fIsoCuts[1][eb])
    pass = kFALSE;
  if (hcalIso>fIsoCuts[2][eb])
    pass = kFALSE;
  if (combinedIso>fIsoCuts[3][eb])
    pass = kFALSE;


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
    Double_t pD0 = ele->GsfTrk()->D0Corrected(*vertices->At(i0));
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
  if(ele->SCluster()->Seed()->Energy() > 5.0 && 
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

//--------------------------------------------------------------------------------------------------
Int_t ElectronTools::Classify(const Electron *ele) {
  
  double eta = ele->AbsEta();
  double eOverP = ele->ESuperClusterOverP();
  double fBrem  = ele->FBrem();

  int cat = -1;
  if (ele->IsEB()) {
    if ((fBrem >= 0.12) && (eOverP > 0.9) && (eOverP < 1.2))
      cat = 0;
    else if ((ele->IsTrackerDriven()) && (!ele->IsEcalDriven()))
      cat = 8;
    else if (fBrem < 0.12)
      cat = 1;
    else if ((eta >  .445   && eta <  .45  ) ||
  	     (eta >  .79    && eta <  .81  ) ||
  	     (eta > 1.137   && eta < 1.157 ) ||
  	     (eta > 1.47285 && eta < 1.4744))
      cat = 6;
    else
      cat = 2;
  } else {
    if ((fBrem >= 0.2) && (eOverP > 0.82) && (eOverP < 1.22))
      cat = 3;
    else if ((ele->IsTrackerDriven()) && (!ele->IsEcalDriven()))
      cat = 8;
    else if (fBrem < 0.2)
      cat = 4;
    else if (eta > 1.5 && eta <  1.58)
      cat = 7;
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
  if (scEt < 20.)
    bin = 2;
  else if (scEt > 30.)
    bin = 0;
  else
    bin = 1;

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

    Double_t d0_real = 99999;
    for(UInt_t i0 = 0; i0 < vertices->GetEntries(); i0++) {
      Double_t pD0 = ele->GsfTrk()->D0Corrected(*vertices->At(i0));
      if(TMath::Abs(pD0) < TMath::Abs(d0_real)) d0_real = TMath::Abs(pD0);
    }
    if (d0_real > cut[26*eb+21])
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
  Double_t cutdcotdistLoose[27] = {
  3.98e-02, 3.97e-02, 4.01e-02, 3.99e-02, 3.97e-02, 3.94e-02, 3.82e-02, 3.93e-02, 3.88e-02,
  3.89e-02, 1.58e-02, 4.18e-03, 3.94e-02, 3.85e-02, 3.86e-02, 2.14e-02, 2.65e-02, 1.33e-02,
  2.86e-02, 3.19e-02, 1.04e-02, 3.86e-02, 3.25e-02, 3.00e-02, 4.77e-04, 3.55e-02, 2.32e-02
  };
  Double_t cutdetainLoose[27] = {
  9.57e-03, 4.87e-03, 1.49e-02, 1.47e-02, 9.13e-03, 1.68e-02, 1.32e-02, 5.18e-02, 3.00e-02,
  9.38e-03, 3.72e-03, 9.06e-03, 1.31e-02, 6.65e-03, 1.23e-02, 1.30e-02, 2.32e-02, 1.18e-02,
  1.07e-02, 3.82e-03, 8.96e-03, 1.40e-02, 6.64e-03, 1.22e-02, 1.20e-02, 2.09e-02, 2.18e-03
  };
  Double_t cutdphiinLoose[27] = {
  4.10e-02, 2.72e-01, 3.69e-01, 4.70e-02, 2.74e-01, 2.91e-01, 3.30e-01, 4.98e-01, 6.38e-01,
  5.95e-02, 9.63e-02, 3.24e-01, 6.87e-02, 6.17e-02, 2.82e-01, 9.06e-02, 2.93e-01, 2.29e-01,
  8.14e-02, 5.19e-02, 2.99e-01, 1.14e-01, 6.05e-02, 3.03e-01, 1.16e-01, 2.56e-01, 4.63e-02
  };
  Double_t cuteseedopcorLoose[27] = {
  7.81e-01, 2.96e-01, 4.81e-01, 9.03e-01, 1.68e-01, 6.44e-01, 1.06e-01, 2.83e-01, 4.65e-01,
  5.99e-01, 3.33e-01, 5.45e-01, 8.17e-01, 7.82e-01, 6.75e-01, 4.17e-01, 8.24e-01, 8.79e-01,
  5.16e-01, 9.38e-01, 8.10e-01, 8.17e-01, 8.50e-01, 5.06e-01, 4.30e-01, 8.33e-01, 7.36e-01
  };
  Double_t cutetLoose[27] = {
  -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
  -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
   1.20e+01,  1.20e+01,  1.20e+01,  1.20e+01,  1.20e+01,  1.20e+01,  1.20e+01,  1.20e+01,  1.42e+01
  };
  Double_t cutfmishitsLoose[27] = {
  4.50e+00, 1.50e+00, 1.50e+00, 2.50e+00, 2.50e+00, 1.50e+00, 2.50e+00, 2.50e+00, 1.50e+00,
  2.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 2.50e+00, 2.50e+00, 5.00e-01,
  2.50e+00, 1.50e+00, 5.00e-01, 1.50e+00, 1.50e+00, 5.00e-01, 2.50e+00, 5.00e-01, 1.50e+00
  };
  Double_t cuthoeLoose[27] = {
  1.65e-01, 7.55e-02, 1.43e-01, 3.70e-01, 5.01e-02, 1.25e-01, 4.06e-01, 2.69e+00, 5.58e-01,
  2.37e-01, 5.43e-02, 1.47e-01, 3.68e-01, 3.12e-02, 1.20e-01, 6.06e-01, 2.81e+00, 1.07e+00,
  1.04e-01, 6.22e-02, 5.73e-02, 2.97e-01, 2.10e-02, 3.08e-02, 5.46e-01, 4.71e+00, 8.00e-01
  };
  Double_t cutip_gsfLoose[27] = {
  4.31e-02, 7.19e-02, 1.40e-01, 8.68e-02, 1.45e-01, 1.57e-01, 9.21e-01, 1.62e-01, 1.86e-01,
  2.09e-02, 4.59e-02, 9.03e-02, 3.65e-02, 1.49e-01, 9.89e-02, 6.24e-02, 1.88e-01, 1.31e-01,
  1.22e-02, 1.23e-02, 7.07e-02, 1.65e-02, 1.23e-01, 6.70e-02, 4.77e-02, 5.92e-02, 1.84e-02
  };
  Double_t cutiso_sumLoose[27] = {
  3.17e+01, 9.85e+00, 8.88e+00, 1.13e+01, 6.05e+00, 7.19e+00, 7.65e+00, 9.84e+00, 3.45e+00,
  1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05,
  1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05
  };
  Double_t cutiso_sumoetLoose[27] = {
  2.58e+01, 1.48e+01, 1.14e+01, 1.93e+01, 7.10e+00, 9.51e+00, 1.09e+01, 9.71e+00, 3.86e+00,
  2.31e+01, 1.58e+01, 1.24e+01, 1.72e+01, 7.37e+00, 8.98e+00, 1.51e+01, 1.31e+01, 6.20e+00,
  2.12e+01, 2.14e+01, 1.75e+01, 1.56e+01, 9.62e+00, 1.06e+01, 2.00e+01, 2.23e+01, 1.37e+01
  };
  Double_t cutseeLoose[27] = {
  1.82e-02, 1.29e-02, 1.68e-02, 3.73e-02, 3.14e-02, 3.30e-02, 1.53e-02, 3.97e-02, 8.41e+00,
  1.62e-02, 1.07e-02, 1.40e-02, 3.62e-02, 3.21e-02, 3.53e-02, 1.19e-02, 3.72e-02, 5.04e+01,
  1.59e-02, 1.12e-02, 1.38e-02, 3.98e-02, 3.27e-02, 4.07e-02, 1.04e-02, 4.42e-02, 1.19e-02
  };

    // Medium cuts
  Double_t cutdcotdistMedium[27] = {
  3.93e-02, 3.01e-02, 3.32e-02, 3.94e-02, 3.91e-02, 3.90e-02, 3.70e-02, 3.88e-02, 3.90e-02,
  3.21e-02, 7.36e-04, 1.95e-04, 3.76e-02, 2.35e-02, 1.13e-02, 1.96e-02, 1.62e-02, 1.04e-03,
  1.34e-03, 2.56e-02, 1.31e-03, 3.26e-02, 8.21e-03, 1.60e-02, 5.17e-03, 2.63e-02, 1.80e-03
  };
  Double_t cutdetainMedium[27] = {
  7.95e-03, 3.91e-03, 9.04e-03, 1.46e-02, 7.02e-03, 1.09e-02, 1.10e-02, 3.87e-02, 1.81e-02,
  9.29e-03, 3.31e-03, 6.94e-03, 1.31e-02, 5.39e-03, 9.18e-03, 1.22e-02, 1.63e-02, 2.82e-03,
  1.11e-02, 4.09e-03, 6.89e-03, 1.08e-02, 5.78e-03, 9.06e-03, 1.19e-02, 1.50e-02, 1.31e-03
  };
  Double_t cutdphiinMedium[27] = {
  4.03e-02, 9.76e-02, 3.24e-01, 4.70e-02, 2.36e-01, 2.73e-01, 3.10e-01, 3.98e-01, 6.38e-01,
  5.53e-02, 5.50e-02, 3.21e-01, 6.57e-02, 2.58e-02, 1.14e-01, 5.92e-02, 2.68e-01, 2.61e-02,
  7.99e-02, 2.47e-02, 1.72e-01, 5.64e-02, 2.45e-02, 8.58e-02, 6.44e-02, 7.18e-02, 3.48e-02
  };
  Double_t cuteseedopcorMedium[27] = {
  7.85e-01, 3.53e-01, 5.77e-01, 9.03e-01, 1.66e-01, 6.43e-01, 1.78e-01, 3.24e-01, 9.02e-01,
  8.15e-01, 9.11e-01, 8.94e-01, 8.48e-01, 8.70e-01, 9.50e-01, 5.34e-01, 9.19e-01, 9.61e-01,
  4.57e-01, 9.51e-01, 9.53e-01, 8.22e-01, 8.74e-01, 9.18e-01, 8.08e-01, 9.19e-01, 1.02e+00
  };
  Double_t cutetMedium[27] = {
  -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
  -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
   1.20e+01,  1.20e+01,  1.20e+01,  1.22e+01,  1.22e+01,  1.22e+01,  1.20e+01,  1.23e+01,  1.84e+01
  };
  Double_t cutfmishitsMedium[27] = {
  2.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 2.50e+00, 2.50e+00, 1.50e+00,
  2.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 2.50e+00, 5.00e-01, 5.00e-01,
  2.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 2.50e+00
  };
  Double_t cuthoeMedium[27] = {
  9.79e-02, 5.42e-02, 1.13e-01, 3.61e-01, 3.73e-02, 7.61e-02, 3.45e-01, 2.69e+00, 4.05e-01,
  1.11e-01, 4.66e-02, 1.39e-01, 3.59e-01, 2.34e-02, 5.75e-02, 6.10e-01, 1.35e+00, 6.35e-01,
  6.07e-02, 5.34e-02, 2.97e-02, 1.16e-01, 1.45e-02, 1.92e-02, 1.30e-01, 5.37e+00, 5.06e-01
  };
  Double_t cutip_gsfMedium[27] = {
  3.95e-02, 5.29e-02, 1.15e-01, 4.56e-02, 1.31e-01, 1.27e-01, 3.13e-01, 9.50e-02, 1.58e-01,
  1.41e-02, 2.00e-02, 5.52e-02, 2.97e-02, 7.65e-02, 4.84e-02, 4.15e-02, 1.73e-01, 1.31e-02,
  1.01e-02, 9.93e-03, 3.60e-02, 1.17e-02, 1.89e-02, 3.43e-02, 1.83e-02, 2.27e-02, 3.30e-03
  };
  Double_t cutiso_sumMedium[27] = {
  1.48e+01, 8.45e+00, 6.96e+00, 8.73e+00, 4.11e+00, 4.97e+00, 5.36e+00, 5.83e+00, 1.89e+00,
  1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05,
  1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05
  };
  Double_t cutiso_sumoetMedium[27] = {
  1.74e+01, 1.13e+01, 8.29e+00, 1.42e+01, 4.72e+00, 6.09e+00, 7.07e+00, 8.32e+00, 2.32e+00,
  1.49e+01, 1.28e+01, 1.01e+01, 1.13e+01, 6.51e+00, 6.93e+00, 1.16e+01, 1.01e+01, 6.16e+00,
  1.40e+01, 1.49e+01, 1.24e+01, 1.02e+01, 6.59e+00, 7.29e+00, 1.51e+01, 1.61e+01, 1.38e+01
  };
  Double_t cutseeMedium[27] = {
  1.53e-02, 1.11e-02, 1.39e-02, 3.31e-02, 3.08e-02, 3.15e-02, 1.26e-02, 3.18e-02, 4.42e-01,
  1.45e-02, 1.04e-02, 1.28e-02, 3.50e-02, 3.02e-02, 3.18e-02, 1.14e-02, 3.36e-02, 3.94e+00,
  1.56e-02, 1.06e-02, 1.28e-02, 3.75e-02, 3.16e-02, 3.65e-02, 9.67e-03, 3.61e-02, 9.56e-03
  };

    // Tight cuts
  Double_t cutdcotdistTight[27] = {
  3.93e-02, 2.56e-02, 8.20e-03, 3.94e-02, 3.85e-02, 3.92e-02, 3.27e-02, 3.71e-02, 3.82e-02,
  2.40e-02, 3.57e-03, 9.45e-03, 3.35e-02, 2.38e-02, 5.29e-04, 1.75e-02, 2.63e-03, 8.06e-05,
  4.87e-04, 2.91e-02, 3.85e-03, 1.94e-02, 6.72e-03, 8.84e-03, 2.41e-04, 2.04e-03, 1.80e-03
  };
  Double_t cutdetainTight[27] = {
  8.57e-03, 3.40e-03, 6.56e-03, 1.05e-02, 6.56e-03, 9.98e-03, 1.04e-02, 1.48e-02, 1.44e-02,
  7.71e-03, 2.59e-03, 5.21e-03, 9.40e-03, 4.25e-03, 8.49e-03, 1.18e-02, 2.29e-02, 1.90e-03,
  8.70e-03, 3.62e-03, 4.92e-03, 9.33e-03, 4.80e-03, 7.72e-03, 1.19e-02, 4.61e-03, 1.34e-03
  };
  Double_t cutdphiinTight[27] = {
  4.05e-02, 5.29e-02, 2.45e-01, 4.05e-02, 4.87e-02, 2.44e-01, 2.85e-01, 2.27e-01, 1.36e-01,
  5.49e-02, 3.41e-02, 1.54e-01, 6.24e-02, 1.87e-02, 3.97e-02, 5.33e-02, 9.76e-02, 7.28e-03,
  4.25e-02, 2.14e-02, 9.01e-02, 5.30e-02, 1.76e-02, 2.48e-02, 6.61e-02, 2.61e-02, 3.48e-02
  };
  Double_t cuteseedopcorTight[27] = {
  7.82e-01, 3.83e-01, 6.00e-01, 9.11e-01, 1.63e-01, 6.39e-01, 5.58e-01, 5.55e-01, 9.62e-01,
  8.35e-01, 9.66e-01, 9.67e-01, 9.23e-01, 8.97e-01, 9.79e-01, 6.25e-01, 9.60e-01, 1.00e+00,
  5.51e-01, 9.67e-01, 9.87e-01, 8.21e-01, 8.56e-01, 1.01e+00, 9.40e-01, 1.01e+00, 1.03e+00
  };
  Double_t cutetTight[27] = {
  -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
  -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05, -1.00e+05,
   1.37e+01,  1.32e+01,  1.36e+01,  1.42e+01,  1.43e+01,  1.41e+01,  1.30e+01,  1.75e+01,  1.83e+01
  };
  Double_t cutfmishitsTight[27] = {
  2.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 1.50e+00, 5.00e-01, 2.50e+00, 5.00e-01, 5.00e-01,
  2.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01,
  2.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 2.50e+00
  };
  Double_t cuthoeTight[27] = {
  7.61e-02, 3.97e-02, 1.07e-01, 1.22e-01, 2.28e-02, 6.10e-02, 1.19e-01, 2.69e+00, 3.14e-01,
  9.01e-02, 4.94e-02, 1.01e-01, 7.59e-02, 1.53e-02, 3.38e-02, 5.41e-01, 1.10e+00, 6.97e-01,
  4.83e-02, 4.88e-02, 2.45e-02, 3.91e-02, 9.84e-03, 1.38e-02, 1.89e-01, 7.10e-01, 4.39e-01
  };
  Double_t cutip_gsfTight[27] = {
  2.45e-02, 4.31e-02, 6.38e-02, 3.33e-02, 7.65e-02, 1.36e-01, 6.79e-01, 9.98e-02, 4.62e-02,
  1.31e-02, 1.45e-02, 4.60e-02, 1.53e-02, 3.43e-02, 2.58e-02, 2.82e-02, 1.29e-01, 3.21e-03,
  8.51e-03, 7.46e-03, 1.77e-02, 1.06e-02, 1.18e-02, 1.13e-02, 1.22e-02, 2.26e-02, 5.29e-03
  };
  Double_t cutiso_sumTight[27] = {
  1.14e+01, 8.08e+00, 6.16e+00, 6.73e+00, 3.18e+00, 4.36e+00, 4.40e+00, 5.30e+00, 1.34e+00,
  1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05,
  1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05, 1.00e+05
  };
  Double_t cutiso_sumoetTight[27] = {
  1.36e+01, 1.02e+01, 7.07e+00, 9.77e+00, 3.57e+00, 4.85e+00, 6.13e+00, 7.89e+00, 2.14e+00,
  1.12e+01, 1.18e+01, 7.73e+00, 8.18e+00, 4.85e+00, 5.05e+00, 1.13e+01, 7.16e+00, 5.76e+00,
  1.04e+01, 1.12e+01, 1.05e+01, 7.53e+00, 5.51e+00, 6.08e+00, 1.21e+01, 1.31e+01, 1.41e+01
  };
  Double_t cutseeTight[27] = {
  1.42e-02, 1.05e-02, 1.25e-02, 3.24e-02, 3.08e-02, 3.01e-02, 1.10e-02, 2.71e-02, 4.69e-02,
  1.33e-02, 1.04e-02, 1.15e-02, 3.31e-02, 2.96e-02, 3.08e-02, 9.90e-03, 2.74e-02, 3.16e-01,
  1.45e-02, 1.05e-02, 1.10e-02, 3.37e-02, 2.97e-02, 3.01e-02, 9.40e-03, 2.77e-02, 9.39e-03
  };

    Double_t cutdcotdist[27];
    Double_t cutdetain[27];
    Double_t cutdphiin[27];
    Double_t cuteseedopcor[27];
    Double_t cutet[27];
    Double_t cutfmishits[27];
    Double_t cuthoe[27];
    Double_t cutip_gsf[27];
    Double_t cutiso_sum[27];
    Double_t cutiso_sumoet[27];
    Double_t cutsee[27];
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

    Double_t iso_sum = tkIso + ecalIso + hcalIso;
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

    Double_t d0_real = 99999;
    for(UInt_t i0 = 0; i0 < vertices->GetEntries(); i0++) {
      Double_t pD0 = ele->GsfTrk()->D0Corrected(*vertices->At(i0));
      if(TMath::Abs(pD0) < TMath::Abs(d0_real)) d0_real = TMath::Abs(pD0);
    }
    if (d0_real < cutip_gsf[cat+bin*9])
      result += 4;

    Bool_t passConvVeto = PassConversionFilter(ele, conversions, kTRUE);
    if (mishits < cutfmishits[cat+bin*9] &&
	passConvVeto)
      result += 8;
  } // classbased

  return result;
}
