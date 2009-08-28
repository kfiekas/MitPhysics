// $Id: ElectronIDMod.cc,v 1.32 2009/08/12 07:50:31 loizides Exp $

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
  fApplyD0Cut(kTRUE),
  fD0Cut(0.025),
  fReverseIsoCut(kFALSE),
  fReverseD0Cut(kFALSE),
  fElIdType(kIdUndef),
  fElIsoType(kIsoUndef),
  fChargeFilter(kTRUE),
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

  if (eOverP<1.5) {  
    if (deltaPhiIn>fCuts[2][cat+4*eb])
      return kFALSE; 
  } else {
    if(deltaPhiIn>fCuts[2][3+4*eb])
      return kFALSE;
  }

  if(deltaEtaIn>fCuts[3][cat+4*eb])
    return kFALSE; 

  if(eSeedOverPin<fCuts[4][cat+4*eb])
    return kFALSE;

  return kTRUE;
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
    
    Bool_t idcut = kFALSE;
    switch (fElIdType) {
      case kTight:
        idcut = e->PassTightID();
        break;
      case kLoose:
        idcut = e->PassLooseID();
       break;
      case kLikelihood:
        idcut = (e->IDLikelihood() > fIDLikelihoodCut);
        break;
      case kNoId:
        idcut = kTRUE;
        break;
      case kCustomIdLoose:
        idcut = ElectronIDMod::PassCustomID(e);
        break;
      case kCustomIdTight:
        idcut = ElectronIDMod::PassCustomID(e);
        break;
      default:
        break;
    }

    if (!idcut) 
      continue;

    Bool_t isocut = kFALSE;
    switch (fElIsoType) {
      case kTrackCalo:
        isocut = (e->TrackIsolationDr03() < fTrackIsolationCut) &&
                 (e->CaloIsolation() < fCaloIsolationCut);
        break;
      case kTrackJura:
        isocut = (e->TrackIsolationDr03() < fTrackIsolationCut) &&
                 (e->EcalRecHitIsoDr04() < fEcalJuraIsoCut) &&
                 (e->HcalIsolation() < fHcalIsolationCut);
        break;
      case kTrackJuraSliding:
        {
          Double_t totalIso = e->TrackIsolationDr03() + e->EcalRecHitIsoDr04() - 1.5;
          if ((totalIso < (e->Pt()-10.0)*5.0/15.0 && e->Pt() <= 25) ||
              (totalIso < 5.0 && e->Pt() > 25) ||
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
      case kCustomIso:
      default:
        break;
    }

    if (isocut == kFALSE)
      continue;

    // apply conversion filter
    Bool_t isGoodConversion = kFALSE;
    if (fApplyConvFilter) {
      LoadBranch(fConversionBranchName);
      for (UInt_t ifc=0; ifc<fConversions->GetEntries(); ifc++) {
        
	Bool_t ConversionMatchFound = kFALSE;
        for (UInt_t d=0; d<fConversions->At(ifc)->NDaughters(); d++) {
          const Track *trk = dynamic_cast<const ChargedParticle*>
	               (fConversions->At(ifc)->Daughter(d))->Trk();
          if (e->GsfTrk() == trk) {
            ConversionMatchFound = kTRUE;
            break;
          }
        }

        // if match between the e-track and one of the conversion legs
        if (ConversionMatchFound == kTRUE){
          isGoodConversion =  (fConversions->At(ifc)->Prob() > 0.0005) &&
                              (fConversions->At(ifc)->Lxy() > 0) &&
                              (fConversions->At(ifc)->Lz() > 0);

          if (isGoodConversion == kTRUE) {
	    for (UInt_t d=0; d<fConversions->At(ifc)->NDaughters(); d++) {
              const Track *trk = dynamic_cast<const ChargedParticle*>
	                   (fConversions->At(ifc)->Daughter(d))->Trk();
            	      
              if (trk) {
                // These requirements are not used for the GSF track
            	if (!(trk->NHits() > 8 && trk->Prob() > 0.005) && trk!=e->GsfTrk())
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
      
    }
    if (isGoodConversion == kTRUE) continue;

    if (fApplyD0Cut) {
      Bool_t d0cut = kFALSE;
      LoadBranch(fVertexName);
      // d0 cut
      Double_t d0_real = 99999;
      for(UInt_t i0 = 0; i0 < fVertices->GetEntries(); i0++) {
	Double_t pD0 = e->GsfTrk()->D0Corrected(*fVertices->At(i0));
	if(TMath::Abs(pD0) < TMath::Abs(d0_real)) d0_real = TMath::Abs(pD0);
      }
      if(d0_real < fD0Cut) d0cut = kTRUE;

      if     (fReverseD0Cut == kTRUE &&
              d0cut == kFALSE && d0_real < 0.05)
	d0cut = kTRUE;
      else if(fReverseD0Cut == kTRUE)
	d0cut = kFALSE;

      if (d0cut == kFALSE)
        continue;
    }
    if(fChargeFilter == kTRUE &&
       e->TrackerTrk() &&
       e->TrackerTrk()->Charge() != e->Charge()) continue;

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

  if (fElectronIDType.CompareTo("Tight") == 0) 
    fElIdType = kTight;
  else if (fElectronIDType.CompareTo("Loose") == 0) 
    fElIdType = kLoose;
  else if (fElectronIDType.CompareTo("Likelihood") == 0) 
    fElIdType = kLikelihood;
  else if (fElectronIDType.CompareTo("NoId") == 0) 
    fElIdType = kNoId;
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
    {0.05, 0.042, 0.045, 0.0, 0.055, 0.037, 0.05, 0.0},        //hovere
    {0.0125, 0.011, 0.01, 0.0, 0.0295, 0.0292, 0.0283, 0.0},   //sigmaetaeta
    {0.032, 0.016, 0.0525, 0.09, 0.025, 0.035, 0.065, 0.092},  //deltaphiin
    {0.0055, 0.003, 0.0075, 0.0, 0.0070, 0.0055, 0.0085, 0.0}, //deltaetain
    {0.24, 0.94, 0.11, 0.0, 0.32, 0.83, 0.0, 0.0},             //eoverp
    {0.8,0.2,0.9,0,0,0,0,0}};                                  //extra cuts fbrem and E_Over_P 

  Double_t loosecuts[6][8]={
    {0.076, 0.033, 0.07, 0.0, 0.083,0.148, 0.033, 0.0},              //hovere
    {0.0101, 0.0095, 0.0097, 0.0, 0.03, 0.03, 0.03, 0.0},            //sigmaetaeta
    {0.053, 0.0189, 0.059, 0.099, 0.0278,0.0157, 0.042, 0.080},      //deltaphiin
    {0.0078, 0.00259, 0.0062, 0.0, 0.0078,0.0061, 0.0061, 0.0},      //deltaetain
    {0.3, 0.92, 0.211, 0.0, 0.42, 0.88, 0.68, 0.0},                  //eoverp
    {0.8,0.2,0,0,0,0,0,0}};                                          //extra cuts fbrem and E_Over_P 

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

