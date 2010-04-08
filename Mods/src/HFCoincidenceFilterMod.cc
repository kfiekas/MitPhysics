// $Id: HFCoincidenceFilterMod.cc,v 1.2 2010/01/18 17:24:53 bendavid Exp $

#include "MitPhysics/Mods/interface/HFCoincidenceFilterMod.h"
#include <TFile.h>
#include <TTree.h>
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/Vertex.h"
#include "MitAna/DataTree/interface/CaloTower.h"

using namespace mithep;

ClassImp(mithep::HFCoincidenceFilterMod)

//--------------------------------------------------------------------------------------------------
HFCoincidenceFilterMod::HFCoincidenceFilterMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fAbort(kTRUE),
  fMinNCoincidentCaloTowers(1),
  fCaloTowerEnergyThreshold(4.0),
  fCaloTowersName(Names::gkCaloTowerBrn),
  fNEvents(0),
  fNAcceped(0),
  fNFailed(0),
  fCaloTowers(0)
{
  // Constructor. 
}

//--------------------------------------------------------------------------------------------------
HFCoincidenceFilterMod::~HFCoincidenceFilterMod() 
{
  // Destructor.
}


//--------------------------------------------------------------------------------------------------
void HFCoincidenceFilterMod::BeginRun()
{
  
}

//--------------------------------------------------------------------------------------------------
void HFCoincidenceFilterMod::Process()
{
  
  LoadBranch(fCaloTowersName);
  
  // Increment counters and stop further processing of an event if current run is excluded

  ++fNEvents; 
  Bool_t PassHFCoincidence = kFALSE;
  Int_t NPositiveEtaHFCaloTowers = 0;
  Int_t NNegativeEtaHFCaloTowers = 0;

  for (UInt_t i=0; i<fCaloTowers->GetEntries(); ++i) {
    
    const CaloTower *tower = fCaloTowers->At(i);
    
    if (tower->Eta() > 3.0 && tower->Eta() < 5.0 && tower->E() > fCaloTowerEnergyThreshold) {
      NPositiveEtaHFCaloTowers++;
    }

    if (tower->Eta() < -3.0 && tower->Eta() > -5.0 && tower->E() > fCaloTowerEnergyThreshold) {
      NNegativeEtaHFCaloTowers++;
    }

    if (NPositiveEtaHFCaloTowers >= fMinNCoincidentCaloTowers && 
        NNegativeEtaHFCaloTowers >= fMinNCoincidentCaloTowers) {
      PassHFCoincidence = kTRUE;
      break;
    }    
  }
  
  // take action if failed
  if (!PassHFCoincidence) {
    ++fNFailed;
    OnFailed();
    if (fAbort) {
      SkipEvent(); // abort processing of this event by sub-modules
    }
    return;
  } 

  // take action if accepted
  ++fNAcceped;
  IncNEventsProcessed();
  OnAccepted();
}

//--------------------------------------------------------------------------------------------------
void HFCoincidenceFilterMod::SlaveBegin()
{

  ReqBranch(fCaloTowersName, fCaloTowers);
    
}

//--------------------------------------------------------------------------------------------------
void HFCoincidenceFilterMod::SlaveTerminate()
{
  // Save number of accepted events.

  SaveNEventsProcessed();
}
