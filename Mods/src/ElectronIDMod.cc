// $Id: ElectronIDMod.cc,v 1.25 2009/06/02 05:30:44 loizides Exp $

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
      case kCustomId:
      default:
        break;
    }

    if (!idcut) 
      continue;

    Bool_t isocut = kFALSE;
    switch (fElIsoType) {
      case kTrackCalo:
        isocut = (e->TrackIsolation() < fTrackIsolationCut) &&
                 (e->CaloIsolation() < fCaloIsolationCut);
        break;
      case kTrackJura:
        isocut = (e->TrackIsolation() < fTrackIsolationCut) &&
                 (e->EcalJurassicIsolation() < fEcalJuraIsoCut) &&
                 (e->HcalIsolation() < fHcalIsolationCut);
        break;
      case kTrackJuraSliding:
        {
          Double_t totalIso = e->TrackIsolation() + e->EcalJurassicIsolation() - 1.5;
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
  else if (fElectronIDType.CompareTo("Custom") == 0) {
    fElIdType = kCustomId;
    SendError(kWarning, "SlaveBegin",
              "Custom electron identification is not yet implemented.");
  } else {
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
