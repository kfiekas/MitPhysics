// $Id: ElectronIDMod.cc,v 1.28 2009/07/21 16:35:12 bendavid Exp $

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
  fElectronIDType("Tight"),
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
  fElectrons(0),
  fConversions(0),
  fVertices(0)
{
  // Constructor.
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
        idcut = ElectronIDMod::PassCustomID(e,kLooseCustom);
        break;
      case kCustomIdTight:
        idcut = ElectronIDMod::PassCustomID(e,kTightCustom);
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
      double d0_real = 99999;
      for(uint i0 = 0; i0 < fVertices->GetEntries(); i0++) {
	double pD0 = e->GsfTrk()->D0Corrected(*fVertices->At(i0));
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

    // add good electron
    GoodElectrons->Add(e);
  }

  // sort according to pt
  GoodElectrons->Sort();

  // add to event for other modules to use
  AddObjThisEvt(GoodElectrons);  
}

Bool_t ElectronIDMod::PassCustomID(const Electron *ele, EElIdCustomType CustomIDName){

     double tightcuts[6][8]={
         {0.056, 0.0221, 0.037, 0.0, 0.0268, 0.0102, 0.0104, 0.0}, //hovere
	 {0.0095, 0.0094, 0.0094, 0.0, 0.029, 0.029, 0.029, 0.0}, //sigmaetaeta
	 {0.0225, 0.0114, 0.0234, 0.039, 0.0215, 0.0095, 0.0148, 0.0167}, //deltaphiin
	 {0.0043, 0.00282, 0.0036, 0.0, 0.0066, 0.0049, 0.0041, 0.0},  //deltaetain
	 {0.32, 0.94, 0.221, 0.0, 0.74, 0.89, 0.66, 0.0},  //eoverp
         {0.8,0.2,0.9,0,0,0,0,0}}; //extra cuts fbrem and E_Over_P 
      double loosecuts[6][8]={
	 {0.076, 0.033, 0.07, 0.0, 0.083,0.148, 0.033, 0.0}, //hovere
	 {0.0101, 0.0095, 0.0097, 0.0, 0.03, 0.03, 0.03, 0.0}, //sigmaetaeta
	 {0.053, 0.0189, 0.059, 0.099, 0.0278,0.0157, 0.042, 0.080}, //deltaphiin
	 {0.0078, 0.00259, 0.0062, 0.0, 0.0078,0.0061, 0.0061, 0.0}, //deltaetain
	 {0.3, 0.92, 0.211, 0.0, 0.42, 0.88, 0.68, 0.0}, //eoverp
	 {0.8,0.2,0,0,0,0,0,0}}; //extra cuts fbrem and E_Over_P 

  double cuts[6][8];

    switch (CustomIDName) {
      case kTightCustom:    
        for (int i=0;i<6;i++){
         for (int j=0;j<8;j++){
           cuts[i][j]=tightcuts[i][j];
          }
        }
         break;
      case kLooseCustom:
        for (int i=0;i<6;i++){
         for (int j=0;j<8;j++){
           cuts[i][j]=loosecuts[i][j];
          }
        }
       break;
    }

   double eOverP = ele->ESuperClusterOverP();
   double pin  = ele->PIn(); 
   double pout = ele->POut();
   double fBrem = (pin-pout)/pin;
   if ((eOverP < cuts[5][0]) && (fBrem < cuts[5][1])) {
       return false; }
   if (eOverP < cuts[5][2]*(1-fBrem)){
        return false; }
   int cat = ElectronIDMod::Categories(ele);
   double eSeedOverPin = ele->ESeedClusterOverPIn(); 
   double hOverE = ele->HadronicOverEm();
   double sigmaee = ele->CoviEtaiEta();
   double deltaPhiIn = fabs(ele->DeltaPhiSuperClusterTrackAtVtx());
   double deltaEtaIn = fabs(ele->DeltaEtaSuperClusterTrackAtVtx());

  int eb;
  if (ele->IsEB()) 
     eb = 0;
   else {
     eb = 1; 
   }
 
   if(hOverE>cuts[0][cat+4*eb]) {  return false;  }
   if(sigmaee>cuts[1][cat+4*eb]){  return false; }
   if(eOverP<1.5){  
       if(deltaPhiIn>cuts[2][cat+4*eb]){ return false;} }
   else{
       if(deltaPhiIn>cuts[2][3+4*eb]){ return false;}  }
   if(deltaEtaIn>cuts[3][cat+4*eb]){ return false;} 
   if(eSeedOverPin<cuts[4][cat+4*eb]){ return false;} 
   return true;
}


Int_t ElectronIDMod::Categories(const Electron *electron){   
   double eOverP = electron->ESuperClusterOverP();
   double pin  = electron->PIn(); 
   double pout = electron->POut(); 
   double fBrem = (pin-pout)/pin;
   
   int cat;
   
   if( (electron->IsEB() && fBrem<0.06) || (electron->IsEE() && fBrem<0.1)) 
     cat=1;
   else if (eOverP < 1.2 && eOverP > 0.8) 
     cat=0;
   else 
     cat=2;
   
  return cat;
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
  else if (fElectronIDType.CompareTo("CustomLoose") == 0) 
    fElIdType = kCustomIdLoose;
  else if (fElectronIDType.CompareTo("CustomTight") == 0) 
    fElIdType = kCustomIdTight;
   else {
    SendError(kAbortAnalysis, "SlaveBegin",
              "The specified electron identification %s is not defined.",
              fElectronIDType.Data());
    return;
  }

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
