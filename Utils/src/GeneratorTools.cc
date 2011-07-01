// $Id  $

#include "MitPhysics/Utils/interface/GeneratorTools.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataUtil/interface/Debug.h"

using namespace mithep;
using namespace std;

//--------------------------------------------------------------------------------------------------
void GeneratorTools::PrintHepMCTable(const mithep::Collection<mithep::MCParticle> *particles, 
                                     Bool_t showDaughters = true, int suppressEntriesAfterThisIndex = -1) {
   
  //Print entire HepMC Decay Table

  for (UInt_t d=0; d<particles->GetEntries(); d++) {
    if (suppressEntriesAfterThisIndex == -1 || int(d) <= suppressEntriesAfterThisIndex) {
      cout << "Particle " << d << " : " << particles->At(d)->PdgId() << "    " 
           << particles->At(d)->Status() << " " 
           << particles->At(d)->Pt() << " " 
           << particles->At(d)->Eta() << " " << particles->At(d)->Phi() 
           << " IsParton:" << particles->At(d)->IsParton() 
           << " Decays at " << particles->At(d)->DecayVertex().Rho() << " ";
      if (particles->At(d)->HasMother()) {
        cout << " Mother: " << particles->At(d)->Mother()->PdgId() << " "
             << particles->At(d)->Mother()->Status() << " " 
             << particles->At(d)->Mother()->Pt() << " "
             << particles->At(d)->Mother()->Eta() << " "
             << particles->At(d)->Mother()->Phi() << " ";
      }
      cout << " " << particles->At(d)->NDaughters() << " Daughters : ";
      if (showDaughters) {
        for (UInt_t e=0;e<particles->At(d)->NDaughters();e++) {
          if (e <= 4 || particles->At(d)->Daughter(e)->Pt() > 5.0 ) {
            cout << " | " 
                 << particles->At(d)->Daughter(e)->PdgId() << " "
                 << particles->At(d)->Daughter(e)->Status() << " "
                 << particles->At(d)->Daughter(e)->Pt() << " "
                 << particles->At(d)->Daughter(e)->Eta() << " "
                 << particles->At(d)->Daughter(e)->Phi() << " "
                 << " ";
          } else {
            cout << " ... ";
          }
        }
      }
      cout << endl;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------------------
void GeneratorTools::PrintNearbyParticles(const mithep::Collection<mithep::MCParticle> *particles, 
                                     Double_t eta, Double_t phi, Double_t deltaR = 0.3) {
   
  //Print Particles in the HepMC table which are within the given DeltaR distance

  mithep::ObjArray<mithep::MCParticle> *descendants = new mithep::ObjArray<mithep::MCParticle>;

  //Print all particles and their descendants that are within the dR Cone.
  for (UInt_t d=0; d<particles->GetEntries(); d++) {
    Double_t dR = MathUtils::DeltaR(phi,eta,particles->At(d)->Phi(),particles->At(d)->Eta());

    //print out the particle if it is inside the dR cone, or if it is a descendant
    if (dR < deltaR || descendants->HasObject(particles->At(d))) {
      cout << "Particle " << d << " : " << particles->At(d)->PdgId() << "    " 
           << particles->At(d)->Status() << " " 
           << particles->At(d)->Pt() << " " 
           << particles->At(d)->Eta() << " " << particles->At(d)->Phi() 
           << " IsParton:" << particles->At(d)->IsParton() 
           << " Decays at " << particles->At(d)->DecayVertex().Rho() 
           << " " << particles->At(d)->NDaughters() << " Daughters : ";

      for (UInt_t e=0;e<particles->At(d)->NDaughters();e++) {
        descendants->Add(particles->At(d)->Daughter(e));
        if (e <= 4 || particles->At(d)->Daughter(e)->Pt() > 5.0) {
          cout << " | " 
               << particles->At(d)->Daughter(e)->PdgId() << " "
               << particles->At(d)->Daughter(e)->Status() << " "
               << particles->At(d)->Daughter(e)->Pt() << " "
               << particles->At(d)->Daughter(e)->Eta() << " "
               << particles->At(d)->Daughter(e)->Phi() << " "
               << " ";
        }
      }
      cout << " ... ";
      cout << endl;
    }
  }
  return;
}

//***********************************************************************************************
//Match Electron to a SimTrack
//***********************************************************************************************  
const MCParticle* GeneratorTools::MatchElectronToSimParticle(const mithep::Collection<mithep::MCParticle> *particles,
                               const mithep::Track *eleTrack, Bool_t isGsfTrack, Int_t printDebugLevel, 
                                                             Bool_t printHepMCTable) {
  
  //Match Clean Electron to a SimTrack and then trace back up the mother chain.
  const MCParticle *match = NULL;

  Double_t matchDeltaPtOverPt = 0.0;
  Double_t secondaryMatchDeltaPtOverPt = 0.0;
  Double_t minProductionVertexRho = 10000.0;
  Double_t secondaryMatchMinDR = 5000.0;
  Double_t minDR = 5000.0;
  const MCParticle *matchedSimChargedParticle = NULL;
  const MCParticle *secondaryMatchedSimChargedParticle = NULL;
  
  if (printDebugLevel >= 5) cout << "\nSimulated particles near Electron\n";
  for(UInt_t j=0;j<particles->GetEntries();j++) {
    Bool_t isStable = false;
    Bool_t isTrackable = false;
    if(particles->At(j)->NDaughters() == 0 || particles->At(j)->DecayVertex().Rho() > 10.0)
      isStable = true;
    
    if (particles->At(j)->Mother() && particles->At(j)->Mother()->DecayVertex().Rho() <= 40)
      isTrackable = true;
    
    if (particles->At(j)->IsSimulated()) {
      Double_t DR = MathUtils::DeltaR(eleTrack->Phi(), eleTrack->Eta(), particles->At(j)->Phi(),
                                      particles->At(j)->Eta());
      
      Double_t dPtOverPt = (eleTrack->Pt() - particles->At(j)->Pt())/eleTrack->Pt();
      
      //remove Neutrals
      if (particles->At(j)->Charge() != 0) { 
        if( DR < 0.3) {           
          if (printDebugLevel >= 5) {
            cout << "ChargedSimulatedParticle: " << j << " : " << fabs(dPtOverPt) << " : " << DR 
                 << "   " << particles->At(j)->PdgId()  
                 << " " << particles->At(j)->Status() << " " << particles->At(j)->Pt()
                 << " " << particles->At(j)->Eta() << " " << particles->At(j)->Phi()
                 << " produced at ";
          
            if (particles->At(j)->Mother())
              cout << particles->At(j)->Mother()->DecayVertex().Rho();
            else 
              cout << "0.0";
            cout << " decays at ";
            if (particles->At(j)->NDaughters() > 0)
              cout << particles->At(j)->DecayVertex().Rho();
            else
              cout << " Stable ";
            cout << endl;
          }
        }
        
        //regular match
        if (DR < minDR && fabs(dPtOverPt) < 0.5 && DR < 0.3 && isStable && isTrackable ) {
          minDR = DR;
          matchedSimChargedParticle = particles->At(j); 
        }   
        
        Double_t productionVertexRho = 20000.0;
        if (particles->At(j)->Mother())
          productionVertexRho = particles->At(j)->Mother()->DecayVertex().Rho();
        
        //perform secondary match for GSF track only electrons
        if ( DR < 0.3 && particles->At(j)->AbsPdgId() == 11 && particles->At(j)->Pt() > 3.0
             && productionVertexRho < minProductionVertexRho ) {
          minProductionVertexRho = productionVertexRho;
          secondaryMatchedSimChargedParticle = particles->At(j);
        }                        
      } //is charged particle
    } //is simulated
  }// for all particles

  if (printDebugLevel >= 5) {
    if (matchedSimChargedParticle && minDR < 0.3) {      
      cout << "Match : " << matchDeltaPtOverPt << " : " << minDR 
           << "   " << matchedSimChargedParticle->PdgId()  
           << " " << matchedSimChargedParticle->Status() << " " << matchedSimChargedParticle->Pt()
           << " " << matchedSimChargedParticle->Eta() << " " << matchedSimChargedParticle->Phi()
           << " produced at ";
      if (matchedSimChargedParticle->Mother())
        cout << matchedSimChargedParticle->Mother()->DecayVertex().Rho();
      else 
        cout << "0.0";
      cout << " decays at ";
      if (matchedSimChargedParticle->NDaughters() > 0)
        cout << matchedSimChargedParticle->DecayVertex().Rho();
        else
          cout << " Stable ";
      cout << endl;
    }
  }

  if (!(matchedSimChargedParticle && minDR < 0.3)) {
    if (secondaryMatchedSimChargedParticle && isGsfTrack) {
      matchedSimChargedParticle = secondaryMatchedSimChargedParticle;
      matchDeltaPtOverPt = secondaryMatchDeltaPtOverPt;
      minDR = secondaryMatchMinDR;
    }

    if (printDebugLevel >= 5) {
      if (matchedSimChargedParticle && minDR < 0.3) {     
        cout << "secondaryMatch : " << matchDeltaPtOverPt << " : " << minDR 
             << "   " << matchedSimChargedParticle->PdgId()  
             << " " << matchedSimChargedParticle->Status() << " " << matchedSimChargedParticle->Pt()
             << " " << matchedSimChargedParticle->Eta() << " " << matchedSimChargedParticle->Phi()
             << " produced at ";
        if (matchedSimChargedParticle->Mother())
          cout << matchedSimChargedParticle->Mother()->DecayVertex().Rho();
        else 
          cout << "0.0";
        cout << " decays at ";
        if (matchedSimChargedParticle->NDaughters() > 0)
          cout << matchedSimChargedParticle->DecayVertex().Rho();
        else
          cout << " Stable ";
        cout << endl;
      }
    }    
  }
  
  if (matchedSimChargedParticle && minDR < 0.3)
    match = matchedSimChargedParticle;
  
  //MDB(kAnalysis, 6) PrintNearbyParticles(particles, ele->Eta(), ele->Phi(), 0.3);
  MDB(kAnalysis, 7) if (printHepMCTable) PrintHepMCTable(particles);

  return match;
}

//***********************************************************************************************
//Return the meson ancestor given the input matched simulated particle. Intended to be the
//particle corresponding to the large x fragmentation resulting in an electron fake.
//***********************************************************************************************  
const MCParticle* GeneratorTools::FindElectronFakeAncestor(const mithep::MCParticle *matchedSimChargedParticle) {

  if (!matchedSimChargedParticle)
    return 0;

  //If it's an electron, try 
  //1) to match to a conversion, then find ancestor pi0
  //2) find B / D meson if it's semi-leptonic decay
  //3) if it's a prompt photon, return null.
  if (matchedSimChargedParticle->AbsPdgId() == 11) {

    //if no mother then return null
    if (!matchedSimChargedParticle->Mother()) {
      MDB(kAnalysis, 5) cout << "no mother\n";
      return 0;
    }


    //while look to look for photon mother
    Bool_t FoundPhotonMother = false;
    const MCParticle *tmpParticle = matchedSimChargedParticle;
    while (tmpParticle->Mother() && tmpParticle->Mother()->PdgId() != 22 ) {     
      tmpParticle = tmpParticle->Mother();
    }
    if (tmpParticle->Mother() && tmpParticle->Mother()->PdgId() == 22) {
      FoundPhotonMother = true;
    }

    //conversions
    if (FoundPhotonMother) {
      //find the meson ancestor.
      const MCParticle *tmpParticle = matchedSimChargedParticle;
      while (tmpParticle->Mother() 
             && 
             tmpParticle->Mother()->PdgId() != 111 
             && 
             tmpParticle->Mother()->PdgId() != 221 
             && 
             tmpParticle->Mother()->PdgId() != 223 
             && 
             !(tmpParticle->Mother()->PdgId() == 22 && tmpParticle->Mother()->Mother() && 
               tmpParticle->Mother()->Mother()->IsParton())             
             &&
             tmpParticle->Mother()->AbsPdgId() != 24
             && 
             tmpParticle->Mother()->AbsPdgId() != 23
             && 
             tmpParticle->Mother()->AbsPdgId() != 92
             && 
             !tmpParticle->Mother()->IsParton()
        ) {           
        MDB(kAnalysis, 5) {
          cout << "find : " << tmpParticle->Mother()->PdgId() << " " << tmpParticle->Mother()->Status() 
               << " " << tmpParticle->Mother()->Pt() << " " << tmpParticle->Mother()->Eta() << " " 
               << tmpParticle->Mother()->Phi() << endl;  
        }
        tmpParticle = tmpParticle->Mother();
      }
      if (tmpParticle->Mother()->PdgId() == 111 
          ||tmpParticle->Mother()->PdgId() == 221 
          ||tmpParticle->Mother()->PdgId() == 223
        ) {
        return tmpParticle->Mother();
      } else if (tmpParticle->Mother()->PdgId() == 22 
                 || tmpParticle->Mother()->PdgId() == 23
                 || tmpParticle->Mother()->PdgId() == 24
        ) {
        //prompt photon or FSR photon or Brem photon
        return 0;
      } else {
        return 0;
      }
      //end if sim electron mother was a photon
    } else {
      //now look for semi-leptonic decays 
      const MCParticle *tmpParticle = matchedSimChargedParticle;
      while (tmpParticle->Mother() 
             && 
             !((tmpParticle->Mother()->AbsPdgId() >= 400 && 
                tmpParticle->Mother()->AbsPdgId() < 600) ||
               (tmpParticle->Mother()->AbsPdgId() >= 4000 && 
                tmpParticle->Mother()->AbsPdgId() < 6000)
               )
        ) {
        MDB(kAnalysis, 5) {
          cout << "find : " << tmpParticle->Mother()->PdgId() << " " << tmpParticle->Mother()->Status() 
               << " " << tmpParticle->Mother()->Pt() << " " << tmpParticle->Mother()->Eta() << " " 
               << tmpParticle->Mother()->Phi() << endl;  
        }
        tmpParticle = tmpParticle->Mother();
      }
      if (tmpParticle->Mother() && tmpParticle->Mother()->HasDaughter(11)) {
        return tmpParticle->Mother();
      }   
    } 
    //end if sim particle was an electron
  } else if (matchedSimChargedParticle->AbsPdgId() == 211 ) {
    return matchedSimChargedParticle;
  } else if (matchedSimChargedParticle->AbsPdgId() == 321) {
    return matchedSimChargedParticle;
  } else {            
    return 0;
  }
  return 0;
}


//***********************************************************************************************
//Categorizes Electron Fakes based on the input matched simulated particle
//***********************************************************************************************  
Int_t GeneratorTools::CategorizeFakeElectron(const mithep::MCParticle *matchedSimChargedParticle) {

  //Categorizes Electron Fakes based on the input matched simulated particle
  Int_t FakeCategory = 0;
  
  if (!matchedSimChargedParticle)
    return FakeCategory;

  //Try to match to conversion
  if (matchedSimChargedParticle->AbsPdgId() == 11) {
    const MCParticle *mother = matchedSimChargedParticle->Mother();
    
    if (!mother) {
      MDB(kAnalysis, 5) cout << "no mother\n";
      FakeCategory = 7;
    }  else {
      //try to find an ancestor that was a pi0 or a heavy meson
      MDB(kAnalysis, 5) cout << "find pi0 or heavy meson ancestor\n";
      const MCParticle *tmpParticle = matchedSimChargedParticle;
      while (tmpParticle->Mother() 
             && 
             tmpParticle->Mother()->PdgId() != 111 
             && 
             tmpParticle->Mother()->PdgId() != 221 
             && 
             tmpParticle->Mother()->PdgId() != 223 
             && 
             !(tmpParticle->Mother()->PdgId() == 22 && tmpParticle->Mother()->Mother() && 
               tmpParticle->Mother()->Mother()->IsParton())
             &&
             !((tmpParticle->Mother()->AbsPdgId() >= 400 && 
                tmpParticle->Mother()->AbsPdgId() < 600) ||
                       (tmpParticle->Mother()->AbsPdgId() >= 4000 && 
                        tmpParticle->Mother()->AbsPdgId() < 6000)
               ) 
             && 
             tmpParticle->Mother()->AbsPdgId() != 24
             && 
             tmpParticle->Mother()->AbsPdgId() != 23
             && 
             tmpParticle->Mother()->AbsPdgId() != 92
             && 
             !tmpParticle->Mother()->IsParton()
        ) {           
        MDB(kAnalysis, 5) {
          cout << "find : " << tmpParticle->Mother()->PdgId() << " " << tmpParticle->Mother()->Status() 
               << " " << tmpParticle->Mother()->Pt() << " " << tmpParticle->Mother()->Eta() << " " 
               << tmpParticle->Mother()->Phi() << endl;  
        }
        tmpParticle = tmpParticle->Mother();
      }
      const MCParticle *Ancestor = tmpParticle->Mother();
      if (Ancestor) {

        MDB(kAnalysis, 5) {
          cout << "Found ancestor: " << Ancestor->PdgId() << " " 
               << Ancestor->Status() << " " << Ancestor->Pt()
               << " " << Ancestor->Eta() << " " 
               << Ancestor->Phi() << "  HasElectronDaughter: " << Ancestor->HasDaughter(11)
               << endl;
        }
        
        if (Ancestor->PdgId() == 111 || Ancestor->PdgId() == 221 || Ancestor->PdgId() == 223) {
          FakeCategory = 1; //pi0->photon->conversion
        } else if (Ancestor->PdgId() == 22) {
          FakeCategory = 6; //prompt photon -> conversion
        } else if ((Ancestor->AbsPdgId() >= 400 && 
                    Ancestor->AbsPdgId() < 600) ||
                   ((Ancestor->AbsPdgId() >= 4000 && 
                    Ancestor->AbsPdgId() < 6000) &&
                   Ancestor->HasDaughter(11))) {
          FakeCategory = 4; //heavy flavor 
        } else if (Ancestor->AbsPdgId() == 24 || Ancestor->AbsPdgId() == 23) {
          FakeCategory = 5;
        } else {
          FakeCategory = 7;
        } 
      } else {
        FakeCategory = 7;
      }
    }
  } else if (matchedSimChargedParticle->AbsPdgId() == 211 ) {
    FakeCategory = 2;
  } else if (matchedSimChargedParticle->AbsPdgId() == 321) {
    FakeCategory = 3;
  } else {            
    FakeCategory = 7;
  }
  return FakeCategory;
}


//***********************************************************************************************
//Match Muon to a SimTrack 
//***********************************************************************************************  
const MCParticle* GeneratorTools::MatchMuonToSimParticle(
  const mithep::Collection<mithep::MCParticle> *particles,const mithep::Track *muonTrack, 
  Bool_t isTrackerTrack, Int_t printDebugLevel, Bool_t printHepMCTable) {
  
  //Match Muon to a SimTrack
  const MCParticle *match = NULL;

  Double_t minDR = 5000.0;
  Double_t matchDeltaPtOverPt = 0.0;
  Double_t secondaryMatchMinDR = 5000.0;
  Double_t secondaryMatchDeltaPtOverPt = 0.0;

  const MCParticle *matchedSimChargedParticle = NULL;
  const MCParticle *secondaryMatchedSimChargedParticle = NULL;
  
  if (printDebugLevel >= 5) cout << "\nSimulated particles near Electron\n";
  for(UInt_t j=0;j<particles->GetEntries();j++) {
    Bool_t isStable = false;
    Bool_t isTrackable = false;

    if (particles->At(j)->IsSimulated()) {
      Double_t DR = MathUtils::DeltaR(muonTrack->Phi(), muonTrack->Eta(), particles->At(j)->Phi(),
                                      particles->At(j)->Eta());
      
      Double_t dPtOverPt = (muonTrack->Pt() - particles->At(j)->Pt())/muonTrack->Pt();
      
      //remove Neutrals
      if (particles->At(j)->Charge() != 0) { 
        if( DR < 0.3) {           
          if (printDebugLevel >= 5) {
            cout << "ChargedSimulatedParticle: " << j << " : " << fabs(dPtOverPt) << " : " << DR 
                 << "   " << particles->At(j)->PdgId()  
                 << " " << particles->At(j)->Status() << " " << particles->At(j)->Pt()
                 << " " << particles->At(j)->Eta() << " " << particles->At(j)->Phi()
                 << " produced at ";
          
            if (particles->At(j)->Mother())
              cout << particles->At(j)->Mother()->DecayVertex().Rho();
            else 
              cout << "0.0";
            cout << " decays at ";
            if (particles->At(j)->NDaughters() > 0)
              cout << particles->At(j)->DecayVertex().Rho();
            else
              cout << " Stable ";
            cout << endl;
          }
        }
        
        //regular match. match to a muon
        if (DR < minDR && DR < 0.3 && particles->At(j)->AbsPdgId() == 13 
            && particles->At(j)->Pt() > 5.0
            && (particles->At(j)->Mother() && 
                particles->At(j)->Mother()->DecayVertex().Rho() <= 150
              )
            ) {
          minDR = DR;
          matchedSimChargedParticle = particles->At(j); 
          matchDeltaPtOverPt = dPtOverPt;
        }
                
        //perform secondary match for non muons
        if (DR < secondaryMatchMinDR && DR < 0.3 && fabs(dPtOverPt) < 0.5 ) {
          secondaryMatchMinDR = DR;
          secondaryMatchedSimChargedParticle = particles->At(j);
          secondaryMatchDeltaPtOverPt = dPtOverPt;
        }   
        
      } //is charged particle
    } //is simulated
  }// for all particles

  if (printDebugLevel >= 5) {
    if (matchedSimChargedParticle && minDR < 0.3) {      
      cout << "Match : " << matchDeltaPtOverPt << " : " << minDR 
           << "   " << matchedSimChargedParticle->PdgId()  
           << " " << matchedSimChargedParticle->Status() << " " << matchedSimChargedParticle->Pt()
           << " " << matchedSimChargedParticle->Eta() << " " << matchedSimChargedParticle->Phi()
           << " produced at ";
      if (matchedSimChargedParticle->Mother())
        cout << matchedSimChargedParticle->Mother()->DecayVertex().Rho();
      else 
        cout << "0.0";
      cout << " decays at ";
      if (matchedSimChargedParticle->NDaughters() > 0)
        cout << matchedSimChargedParticle->DecayVertex().Rho();
        else
          cout << " Stable ";
      cout << endl;
    }
  }

  if (!matchedSimChargedParticle) {
    if (secondaryMatchedSimChargedParticle) {
      matchedSimChargedParticle = secondaryMatchedSimChargedParticle;
      minDR = secondaryMatchMinDR;
      matchDeltaPtOverPt = secondaryMatchDeltaPtOverPt;
    }

    if (printDebugLevel >= 5) {
      if (matchedSimChargedParticle && minDR < 0.3) {     
        cout << "secondaryMatch : " << matchDeltaPtOverPt << " : " << minDR 
             << "   " << matchedSimChargedParticle->PdgId()  
             << " " << matchedSimChargedParticle->Status() << " " << matchedSimChargedParticle->Pt()
             << " " << matchedSimChargedParticle->Eta() << " " << matchedSimChargedParticle->Phi()
             << " produced at ";
        if (matchedSimChargedParticle->Mother())
          cout << matchedSimChargedParticle->Mother()->DecayVertex().Rho();
        else 
          cout << "0.0";
        cout << " decays at ";
        if (matchedSimChargedParticle->NDaughters() > 0)
          cout << matchedSimChargedParticle->DecayVertex().Rho();
        else
          cout << " Stable ";
        cout << endl;
      }
    }    
  }
  
  if (matchedSimChargedParticle && minDR < 0.3)
    match = matchedSimChargedParticle;
  
  return match;
}


//***********************************************************************************************
//Categorizes Electron Fakes based on the input matched simulated particle
//***********************************************************************************************  
Int_t GeneratorTools::CategorizeFakeMuon(const mithep::MCParticle *matchedSimChargedParticle) {

  //Categorizes Electron Fakes based on the input matched simulated particle
  Int_t FakeCategory = 0;
  
  if (!matchedSimChargedParticle)
    return FakeCategory;

  //look for a muon first
  if (matchedSimChargedParticle->AbsPdgId() == 13) {
    const MCParticle *mother = matchedSimChargedParticle->Mother();
    
    if (!mother) {
      MDB(kAnalysis, 5) cout << "no mother\n";
      FakeCategory = 7;
    }  else {
      //try to find an ancestor that was a pi+,K+ or a heavy meson
      const MCParticle *tmpParticle = matchedSimChargedParticle;
      while (tmpParticle->Mother() 
             &&
             tmpParticle->Mother()->AbsPdgId() != 24 
             &&
             tmpParticle->Mother()->AbsPdgId() != 23 
             &&
             tmpParticle->Mother()->AbsPdgId() != 211 
             &&
             tmpParticle->Mother()->AbsPdgId() != 321 
             &&
             tmpParticle->Mother()->AbsPdgId() != 113 
             &&
             !((tmpParticle->Mother()->AbsPdgId() >= 400 && 
                tmpParticle->Mother()->AbsPdgId() < 600) ||
                       (tmpParticle->Mother()->AbsPdgId() >= 4000 && 
                        tmpParticle->Mother()->AbsPdgId() < 6000)
               )
        ) {           
        MDB(kAnalysis, 5) {
          cout << "find : " << tmpParticle->Mother()->PdgId() << " " << tmpParticle->Mother()->Status() 
               << " " << tmpParticle->Mother()->Pt() << " " << tmpParticle->Mother()->Eta() << " " 
               << tmpParticle->Mother()->Phi() << endl;  
        }
        tmpParticle = tmpParticle->Mother();
      }
      const MCParticle *Ancestor = tmpParticle->Mother();

      if (Ancestor) {
        MDB(kAnalysis, 5) {
          cout << "Found ancestor: " << Ancestor->PdgId() << " " 
               << Ancestor->Status() << " " << Ancestor->Pt()
               << " " << Ancestor->Eta() << " " 
               << Ancestor->Phi() << "  HasElectronDaughter: " << Ancestor->HasDaughter(13)
               << endl;        
        }
        
        if (Ancestor->AbsPdgId() == 211) {
          FakeCategory = 1; //pi+ decay in flight
        } else if (Ancestor->AbsPdgId() == 321) {
          FakeCategory = 2; //K+ decay in flight
        } else if ((Ancestor->AbsPdgId() >= 400 && 
                    Ancestor->AbsPdgId() < 600) ||
                  ((Ancestor->AbsPdgId() >= 4000 && 
                    Ancestor->AbsPdgId() < 6000) &&
                    Ancestor->HasDaughter(13))) {
          FakeCategory = 4; //heavy flavor semileptonic decay
        } else if ( Ancestor->AbsPdgId() == 24 || Ancestor->AbsPdgId() == 23 ) {
          FakeCategory = 7; //comes from real muon. bad global fit direction
        } else {
          FakeCategory = 3; //muon from something else
        }
      } else {
        FakeCategory = -1; //no ancestor found
      }
    }
  } else if (matchedSimChargedParticle->AbsPdgId() == 211 
             || matchedSimChargedParticle->AbsPdgId() == 321
             || matchedSimChargedParticle->AbsPdgId() == 2212 ) {
    FakeCategory = 5; //pion punchthrough
  } else {            
    FakeCategory = 6; //other punchthrough?
  }
  return FakeCategory;
}



//--------------------------------------------------------------------------------------------------
TString GeneratorTools::ConvertPdgIdToName(Int_t pdgId) {

  //Convert the pdgID code into a particle name

  TString answer = "";

  switch ( pdgId )
  {
    case 11:
      answer = "e-";
      break;
    case -11:
      answer = "e+"; 
      break;
    case 12:
      answer = "nu_e";
      break;
    case -12:
      answer = "nu_e bar"; 
      break;
    case 13:
      answer = "mu-";
      break;
    case -13:
      answer = "mu+"; 
      break;
    case 15:
      answer = "tau-";
      break;
    case -15:
      answer = "tau+"; 
      break;
    case 21:
      answer = "gluon"; 
      break;
    case 22:
      answer = "photon";
      break;
    case 23:
      answer = "Z"; 
      break;
    case 24:
      answer = "W+";
      break;
    case -24:
      answer = "W-";
      break;
    case 25:
      answer = "Higgs";
      break;
    case 111:
      answer = "Pi0";
      break;
    case 113:
      answer = "Rho0(770)";
      break;
    case 115:
      answer = "a_2^0(1320)";
      break;
    case 117:
      answer = "rho_3^0(1690)";
      break;
    case 119:
      answer = "a_4^0(2040)";
      break;
    case 130:
      answer = "K_L^0";
      break;



    case 211:
      answer = "Pi+";
      break;
    case -211:
      answer = "Pi-";
      break;
    case 213:
      answer = "Rho+";
      break;
    case -213:
      answer = "Rho-";
      break;
    case 221:
      answer = "Eta";
      break;
    case 223:
      answer = "omega(782)";
      break;
    case 310:
      answer = "K_S^0";
      break;
    case 311:
      answer = "K0";
      break;
    case 321:
      answer = "K+";
      break;
    case -321:
      answer = "K-";
      break;
    case 313:
      answer = "K0*(892)";
      break;
    case 323:
      answer = "K+*(892)";
      break;
    case -323:
      answer = "K-*(892)";
      break;
    case 315:
      answer = "K_2^0*(1430)";
      break;
    case 325:
      answer = "K_2^+*(1430)";
      break;
    case -325:
      answer = "K_2^-*(1430)";
      break;
    case 317:
      answer = "K_3^0*(1780)";
      break;
    case 327:
      answer = "K_3^+*(1780)";
      break;
    case -327:
      answer = "K_3^-*(1780)";
      break;
    case 319:
      answer = "K_4^0*(2045)";
      break;
    case 329:
      answer = "K_4^+*(2045)";
      break;
    case -329:
      answer = "K_4^-*(2045)";
      break;
    case 331:
      answer = "Eta'(958)";
      break;
    case 333:
      answer = "phi'(1020)";
      break;


    case 411:
      answer = "D+";
      break;
    case -411:
      answer = "D-";
      break;
    case 421:
      answer = "D0";
      break;
    case -421:
      answer = "D0 bar";
      break;
    case 413:
      answer = "D^+*(2010)";
      break;
    case -413:
      answer = "D^-*(2010)";
      break;
    case 423:
      answer = "D^0*(2010)";
      break;
    case -423:
      answer = "D^0*(2010) bar";
      break;
    case 415:
      answer = "D_2^+*(2460)";
      break;
    case -415:
      answer = "D_2^-*(2460)";
      break;
    case 425:
      answer = "D_2^0*(2460)";
      break;
    case -425:
      answer = "D_2^0*(2460) bar";
      break;
    case 431:
      answer = "D_s^+";
      break;
    case -431:
      answer = "D_s^-";
      break;
    case 433:
      answer = "D_s^+*";
      break;
    case -433:
      answer = "D_s^-*";
      break;
    case 435:
      answer = "D_s2^+*";
      break;
    case -435:
      answer = "D_s2^-*";
      break;

    case 511:
      answer = "B0";
      break;
    case -511:
      answer = "B0 bar";
      break;
    case -521:
      answer = "B-";
      break;
    case 521:
      answer = "B+";
      break;
    case 513:
      answer = "B0*";
      break;
    case -513:
      answer = "B0* bar";
      break;
    case -523:
      answer = "B-*";
      break;
    case 523:
      answer = "B+*";
      break;
    case 515:
      answer = "B_2^0*";
      break;
    case -515:
      answer = "B_2^0* bar";
      break;
    case -525:
      answer = "B_2^-*";
      break;
    case 525:
      answer = "B_2^+*";
      break;
    case 531:
      answer = "B_s^0";
      break;
    case -531:
      answer = "B_s^0 bar";
      break;
    case 533:
      answer = "B_s^0*";
      break;
    case -533:
      answer = "B_s^0* bar";
      break;
    case 535:
      answer = "B_s2^0*";
      break;
    case -535:
      answer = "B_s2^0* bar";
      break;
    case 541:
      answer = "B_c^+";
      break;
    case -541:
      answer = "B_c^-";
      break;
    case 543:
      answer = "B_c^+*";
      break;
    case -543:
      answer = "B_c^-*";
      break;
     case 545:
      answer = "B_c2^+*";
      break;
    case -545:
      answer = "B_c2^-*";
      break;

     case 441:
      answer = "eta_c(1S)";
      break;
     case 443:
      answer = "JPsi(1S)";
      break;
     case 445:
      answer = "chi_c2(1P)";
      break;
     case 551:
      answer = "eta_b(1S)";
      break;
     case 553:
      answer = "Upsilon(1S)";
      break;
     case 555:
      answer = "chi_b2(1P)";
      break;
     case 557:
      answer = "Upsilon(1D)";
      break;

     case 2212:
      answer = "proton";
      break;
     case 2112:
      answer = "neutron";
      break;
     case 2224:
      answer = "Delta++";
      break;
     case 2214:
      answer = "Delta+";
      break;
     case 2114:
      answer = "Delta0";
      break;
     case 1114:
      answer = "Delta-";
      break;
     case 3122:
      answer = "Lambda";
      break;
     case 3222:
      answer = "Sigma+";
      break;
     case 3212:
      answer = "Sigma0";
      break;
     case 3112:
      answer = "Sigma-";
      break;
     case 3224:
      answer = "Sigma+*";
      break;
     case 3214:
      answer = "Sigma0*";
      break;
     case 3114:
      answer = "Sigma-*";
      break;
     case 3322:
      answer = "Xi0";
      break;
     case 3312:
      answer = "Xi-";
      break;
     case 3324:
      answer = "Xi0*";
      break;
     case 3314:
      answer = "Xi-*";
      break;
     case 3334:
      answer = "Omega-";
      break;

    default:
      answer = "";
  }  
  return answer;
}
