// $Id: ElectronIDMod.cc,v 1.42 2009/10/26 14:33:20 sixie Exp $

#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::ElectronIDMod)

//--------------------------------------------------------------------------------------------------
ElectronIDMod::ElectronIDMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fElectronBranchName(Names::gkElectronBrn),
  fConversionBranchName(Names::gkMvfConversionBrn),
  fGoodElectronsName(ModNames::gkGoodElectronsName),  
  fVertexName(string("PrimaryVertexesBeamSpot").c_str()),
  fElectronIDType("CustomTight"),
  fElectronIsoType("TrackJuraSliding"),
  fElectronPtMin(10),
  fIDLikelihoodCut(0.9),
  fTrackIsolationCut(5.0),
  fCaloIsolationCut(5.0),
  fEcalJuraIsoCut(5.0),
  fHcalIsolationCut(5.0),
  fApplyConvFilter(kTRUE),
  fWrongHitsRequirement(kTRUE),
  fApplyD0Cut(kTRUE),
  fD0Cut(0.025),
  fChargeFilter(kTRUE),
  fReverseIsoCut(kFALSE),
  fReverseD0Cut(kFALSE),
  fElIdType(kIdUndef),
  fElIsoType(kIsoUndef),
  fElectrons(0),
  fConversions(0),
  fVertices(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronIDMod::PassCustomID(const Electron *ele) const
{
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
Bool_t ElectronIDMod::PassIDCut(const Electron *ele, EElIdType idType) const
{

  Bool_t idcut = kFALSE;
  switch (idType) {
    case kTight:
      idcut = ele->PassTightID();
      break;
    case kLoose:
      idcut = ele->PassLooseID();
      break;
    case kLikelihood:
      idcut = (ele->IDLikelihood() > fIDLikelihoodCut);
      break;
    case kNoId:
      idcut = kTRUE;
      break;
    case kCustomIdLoose:
      idcut = ElectronIDMod::PassCustomID(ele);
      break;
    case kCustomIdTight:
      idcut = ElectronIDMod::PassCustomID(ele);
      break;
    case kZeeId:
      if (ele->IsEB()) {
        idcut = (ele->CoviEtaiEta() < 0.01 && ele->DeltaEtaSuperClusterTrackAtVtx() < 0.0071);
      } else {
        idcut = (ele->CoviEtaiEta() < 0.028 && ele->DeltaEtaSuperClusterTrackAtVtx() < 0.0066);
      }
      break;
    default:
      break;
  }
  
  return idcut;
}


//--------------------------------------------------------------------------------------------------
Bool_t ElectronIDMod::PassIsolationCut(const Electron *ele, EElIsoType isoType) const
{

  Bool_t isocut = kFALSE;
  switch (isoType) {
    case kTrackCalo:
      isocut = (ele->TrackIsolationDr03() < fTrackIsolationCut) &&
        (ele->CaloIsolation() < fCaloIsolationCut);
      break;
    case kTrackJura:
      isocut = (ele->TrackIsolationDr03() < fTrackIsolationCut) &&
        (ele->EcalRecHitIsoDr04() < fEcalJuraIsoCut) &&
        (ele->HcalIsolation() < fHcalIsolationCut);
      break;
    case kTrackJuraSliding:
    {
      Double_t totalIso = ele->TrackIsolationDr03() + ele->EcalRecHitIsoDr04() - 1.5;
      if (totalIso < (ele->Pt()-10.0)*4.5/20.0 ||
          totalIso <= 0)
        isocut = kTRUE;
      
      if     (fReverseIsoCut == kTRUE &&
              isocut == kFALSE && totalIso < 10)
        isocut = kTRUE;
      else if(fReverseIsoCut == kTRUE)
        isocut = kFALSE;
    }
    break;
    case kNoIso:
      isocut = kTRUE;
      break;
    case kZeeIso:
      if (ele->IsEB()) {
        isocut = (ele->TrackIsolationDr04() < 7.2 && ele->EcalRecHitIsoDr04() < 5.7 && ele->HcalTowerSumEtDr04() < 8.1);
      } else {
        isocut = (ele->TrackIsolationDr04() < 5.1 && ele->EcalRecHitIsoDr04() < 5.0 && ele->HcalTowerSumEtDr04() < 3.4);
      }      
      break;
    case kCustomIso:
    default:
      break;
  }

  return isocut;
}


//--------------------------------------------------------------------------------------------------
Bool_t ElectronIDMod::PassConversionFilter(const Electron *ele, const DecayParticleCol *conversions) const
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
            if (fWrongHitsRequirement && sd->NWrongHits() != 0)
              isGoodConversion = kFALSE;
            
          } else {
            isGoodConversion = kFALSE;
          }
        }
      }
    }
    
    if (isGoodConversion == kTRUE) break;
    
  } // loop over all conversions 
  
  return isGoodConversion;
}

//--------------------------------------------------------------------------------------------------
Bool_t ElectronIDMod::PassD0Cut(const Electron *ele, const VertexCol *vertices) const
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
Bool_t ElectronIDMod::PassChargeFilter(const Electron *ele) const
{
  Bool_t passChargeFilter = kTRUE;
  if(ele->TrackerTrk() &&
     ele->TrackerTrk()->Charge() != ele->Charge()) passChargeFilter = kFALSE;


  return passChargeFilter;
}


//--------------------------------------------------------------------------------------------------
void ElectronIDMod::Process()
{
  // Process entries of the tree. 

  LoadEventObject(fElectronBranchName, fElectrons);

  ElectronOArr *GoodElectrons = new ElectronOArr;
  GoodElectrons->SetName(fGoodElectronsName);

  for (UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    const Electron *e = fElectrons->At(i);        

    if (e->Pt() <= fElectronPtMin) 
      continue;
    
    //apply id cut
    Bool_t idcut = PassIDCut(e, fElIdType);
    if (!idcut) 
      continue;

    //apply Isolation Cut
    Bool_t isocut = PassIsolationCut(e, fElIsoType);
    if (!isocut)
      continue;

    // apply conversion filter
    Bool_t isGoodConversion = kFALSE;
    if (fApplyConvFilter) {
      LoadEventObject(fConversionBranchName, fConversions);
      isGoodConversion = PassConversionFilter(e, fConversions);      
    }
    if (isGoodConversion) continue;
    
    // apply d0 cut
    if (fApplyD0Cut) {
      LoadEventObject(fVertexName, fVertices);
      Bool_t passD0cut = PassD0Cut(e, fVertices);
      if (!passD0cut)
        continue;
    }

    //apply charge filter
    if(fChargeFilter == kTRUE) {
      Bool_t passChargeFilter = PassChargeFilter(e);
      if (!passChargeFilter) continue;
    } 
    // add good electron
    GoodElectrons->Add(e);
  }

  // sort according to pt
  GoodElectrons->Sort();

  // add to event for other modules to use
  AddObjThisEvt(GoodElectrons);  
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the electron collection branch.

  ReqEventObject(fElectronBranchName, fElectrons, kTRUE);

  if (fApplyConvFilter)
    ReqEventObject(fConversionBranchName, fConversions, kTRUE);

  if (fApplyD0Cut)
    ReqEventObject(fVertexName, fVertices, kTRUE);

  Setup();
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMod::Setup()
{
  // Set all options properly before execution.

  if (fElectronIDType.CompareTo("Tight") == 0) 
    fElIdType = kTight;
  else if (fElectronIDType.CompareTo("Loose") == 0) 
    fElIdType = kLoose;
  else if (fElectronIDType.CompareTo("Likelihood") == 0) 
    fElIdType = kLikelihood;
  else if (fElectronIDType.CompareTo("NoId") == 0) 
    fElIdType = kNoId;
  else if (fElectronIDType.CompareTo("ZeeId") == 0) 
    fElIdType = kZeeId;
  else if (fElectronIDType.CompareTo("CustomLoose") == 0) {
    fElIdType = kCustomIdLoose;
  } else if (fElectronIDType.CompareTo("CustomTight") == 0) {
    fElIdType = kCustomIdTight;
  }
   else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified electron identification %s is not defined.",
              fElectronIDType.Data());
    return;
  }

  SetCustomIDCuts(fElIdType);

  if (fElectronIsoType.CompareTo("TrackCalo") == 0 )
    fElIsoType = kTrackCalo;
  else if (fElectronIsoType.CompareTo("TrackJura") == 0) 
    fElIsoType = kTrackJura;
  else if(fElectronIsoType.CompareTo("TrackJuraSliding") == 0)
    fElIsoType = kTrackJuraSliding;
  else if (fElectronIsoType.CompareTo("NoIso") == 0 )
    fElIsoType = kNoIso;
  else if (fElectronIsoType.CompareTo("ZeeIso") == 0 )
    fElIsoType = kZeeIso;
  else if (fElectronIsoType.CompareTo("Custom") == 0 ) {
    fElIsoType = kCustomIso;
    SendError(kWarning, "SlaveBegin",
              "Custom electron isolation is not yet implemented.");
  } else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified electron isolation %s is not defined.",
              fElectronIsoType.Data());
    return;
  }

}

//--------------------------------------------------------------------------------------------------
void ElectronIDMod::SetCustomIDCuts(EElIdType idt)
{
  // Set cut values based on RecoEgamma/ElectronIdentification/python/electronIdCutBasedExt_cfi.py.
  // The following changes are in sigmaetaeta for endcups and deltaetain.

  Double_t tightcuts[6][8]={
    {0.086, 0.1, 0.052, 0.0, 0.050, 0.059, 0.061, 0.0},   //hovere
    {0.01, 0.01, 0.01, 0.0, 0.033, 0.029, 0.028, 0.0},     //sigmaetaeta
    {0.038, 0.024, 0.038, 0.0, 0.034, 0.011, 0.023, 0.0},   //deltaphiin
    {0.0081, 0.0029, 0.0051, 0.0, 0.0056, 0.0062, 0.0088, 0.0}, //deltaetain
    {0.0, 0.9, 0.0, 0.0, 0.0, 0.78, 0.0, 0.0},             //eoverp
    {0.8,0.2,0.9,0,0,0,0,0}};                              //extra cuts fbrem and E_Over_P 

  Double_t loosecuts[6][8]={
    {0.076, 0.033, 0.07, 0.0, 0.083,0.148, 0.033, 0.0},         //hovere
    {0.0101, 0.0095, 0.0097, 0.0, 0.03, 0.03, 0.03, 0.0},       //sigmaetaeta
    {0.053, 0.0189, 0.059, 0.099, 0.0278,0.0157, 0.042, 0.080}, //deltaphiin
    {0.0078, 0.00259, 0.0062, 0.0, 0.0078,0.0061, 0.0061, 0.0}, //deltaetain
    {0.3, 0.92, 0.211, 0.0, 0.42, 0.88, 0.68, 0.0},             //eoverp
    {0.8,0.2,0,0,0,0,0,0}};                                     //extra cuts fbrem and E_Over_P 


  switch (idt) {
    case kCustomIdTight:    
      memcpy(fCuts,tightcuts,sizeof(fCuts));
      break;
    case kCustomIdLoose:
      memcpy(fCuts,loosecuts,sizeof(fCuts));
      break;
    default:
      memset(fCuts,0,sizeof(fCuts));
      break;
  }
}
