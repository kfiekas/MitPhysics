// $Id: ElectronTools.cc,v 1.1 2010/04/10 18:06:51 sixie Exp $

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
