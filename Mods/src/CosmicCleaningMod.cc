// $Id: CosmicCleaningMod.cc,v 1.0

#include "MitPhysics/Mods/interface/CosmicCleaningMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/Track.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::CosmicCleaningMod)

//--------------------------------------------------------------------------------------------------
CosmicCleaningMod::CosmicCleaningMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCosmicsName("random"),        
  fCleanMuonsName(ModNames::gkCleanMuonsName),        
  fCleanCosmicsName(ModNames::gkCleanCosmicsName),
  fDeltaR(0.1),
  fCosmics(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void CosmicCleaningMod::Process()
{
  // Process entries of the tree. 

  // get input collection
  const MuonCol     *CleanMuons    = GetObjThisEvt<MuonCol>(fCleanMuonsName);
  LoadEventObject(fCosmicsName, fCosmics);

  // Go through all cosmics and remove cosmic overlaps with collision muons and duplicates.
  std::vector<const Muon*> CleanCosTemp;

  if (fCosmics) {
    for (UInt_t i=0; i<fCosmics->GetEntries(); ++i) {    
      const Muon *u = fCosmics->At(i);   

      FourVectorM mom(u->Mom());

      // Check whether it overlaps with a good muon: If the muon and cosmic both have
      // tracker tracks then compare the tracks, otherwise
      Bool_t isMuonOverlap = kFALSE;
      if (CleanMuons) {
        UInt_t n = CleanMuons->GetEntries();
        for (UInt_t j=0; j<n; ++j) {
          Double_t deltaR = MathUtils::DeltaR(CleanMuons->At(j)->Mom(), mom);     
          if (deltaR < fDeltaR) {
            isMuonOverlap = kTRUE;
            break;         
          }
        }
      } else {
        cout << "Warning: GoodMuons collection " << fCleanMuonsName << " was not found.\n";
      }
        
      if (isMuonOverlap)
        continue;

      // if no overlaps then add to clean cosmics
      CleanCosTemp.push_back(fCosmics->At(i));   
    } 
  } else {
    cout << "Warning: fCosmics collection " << fCosmicsName << " was not found.\n";
  }

  // Fill the muon array with the contents of the vector:
  MuonOArr *CleanCosmics = new MuonOArr;
  CleanCosmics->SetName(fCleanCosmicsName);

  for (UInt_t j=0; j<CleanCosTemp.size(); ++j) CleanCosmics->Add(CleanCosTemp[j]);
  CleanCosmics->Sort();
       
  // add to event for other modules to use
  AddObjThisEvt(CleanCosmics);
}

//--------------------------------------------------------------------------------------------------
void CosmicCleaningMod::SlaveBegin()
{
  ReqEventObject(fCosmicsName, fCosmics, kTRUE);

}
