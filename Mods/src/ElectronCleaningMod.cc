// $Id: ElectronCleaningMod.cc,v 1.4 2008/11/28 09:13:50 loizides Exp $

#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::ElectronCleaningMod)

//--------------------------------------------------------------------------------------------------
ElectronCleaningMod::ElectronCleaningMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fGoodElectronsName(ModNames::gkGoodElectronsName),        
  fCleanMuonsName(ModNames::gkCleanMuonsName),        
  fCleanElectronsName(ModNames::gkCleanElectronsName)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void ElectronCleaningMod::Process()
{
  // Process entries of the tree. 

  // get input collection
  MuonOArr *CleanMuons = GetObjThisEvt<MuonOArr>(fCleanMuonsName);
  ElectronOArr *GoodElectrons =  GetObjThisEvt<ElectronOArr>(fGoodElectronsName);

  // Go through all electrons and remove electron overlaps with muons and duplicates.
  std::vector<const Electron*> CleanElTemp;
  for (UInt_t i=0; i<GoodElectrons->GetEntries(); ++i) {    
    const Electron *e = GoodElectrons->At(i);   

    FourVector mom(e->Mom());
    const Track *trtrack = e->TrackerTrk();

    // Check whether it overlaps with a good muon: If the muon and electron both have
    // tracker tracks then compare the tracks, otherwise
    Bool_t isMuonOverlap = kFALSE;
    for (UInt_t j=0; j<CleanMuons->GetEntries(); ++j) {
      
      if (trtrack && CleanMuons->At(j)->TrackerTrk()) {
        if (CleanMuons->At(j)->TrackerTrk() == trtrack) {
          isMuonOverlap = kTRUE;
          break;	 
        }
      } else {      
        Double_t deltaR = MathUtils::DeltaR(CleanMuons->At(j)->Mom(), mom);     
        if (deltaR < 0.1) {
          isMuonOverlap = kTRUE;
          break;	 
        }
      }
    }
    
    if (isMuonOverlap)
      continue;

    // Check whether it overlaps with another electron candidate:
    // Here I check whether we have two electron candidates with the same super cluster
    // or two electron candidates with the same track. At the end we also check deltaR
    // to be sure. If there is a duplicate we swap the old one with the new one if the new
    // one has E/P closer to 1.0.
    bool isElectronOverlap = kFALSE;
    for (UInt_t j=0; j<CleanElTemp.size(); ++j) {

      if (e->SCluster() == CleanElTemp[j]->SCluster() ||
          e->Trk() == CleanElTemp[j]->Trk())
        isElectronOverlap = kTRUE;

      if (!isElectronOverlap) {
        Double_t deltaR = MathUtils::DeltaR(CleanElTemp[j]->Mom(), mom);
        if (deltaR < 0.1)
          isElectronOverlap = kTRUE;        
      }

      if (isElectronOverlap) {
        if (TMath::Abs(CleanElTemp[j]->ESuperClusterOverP() - 1) > 
            TMath::Abs(e->ESuperClusterOverP() - 1)) {	   
          CleanElTemp[j] = e;
        }        
        break;
      }
    }

    if (isElectronOverlap)
      continue;

    // if no overlaps then add to clean electrons
    CleanElTemp.push_back(GoodElectrons->At(i));   
  } 

  // Fill the electron array with the contents of the vector:
  ElectronOArr *CleanElectrons = new ElectronOArr;
  CleanElectrons->SetName(fCleanElectronsName);

  for (UInt_t j=0; j<CleanElTemp.size(); ++j) 
    CleanElectrons->Add(CleanElTemp[j]);
       
  // add to event for other modules to use
  AddObjThisEvt(CleanElectrons);  
}
