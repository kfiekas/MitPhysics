// $Id: GenFakesMod.cc,v 1.1 2009/06/30 10:47:17 loizides Exp $

#include "MitPhysics/FakeMods/interface/GenFakesMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/FakeMods/interface/FakeEventHeader.h"

using namespace mithep;

ClassImp(mithep::GenFakesMod)

//--------------------------------------------------------------------------------------------------
GenFakesMod::GenFakesMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fElectronFRFilename("InputRequired"),
  fMuonFRFilename("InputRequired"),
  fUse2DFakeRate(false),
  fUseFitFunction(false),
  fElectronFRFunctionName("InputRequired"),
  fMuonFRFunctionName("InputRequired"),
  fElectronFRHistName("InputRequired"),
  fMuonFRHistName("InputRequired"),
  fCleanElectronsName(ModNames::gkCleanElectronsName),        
  fCleanMuonsName(ModNames::gkCleanMuonsName),        
  fCleanPhotonsName(ModNames::gkCleanPhotonsName),        
  fCleanJetsName(ModNames::gkCleanJetsName),
  fElFakeableObjsName(ModNames::gkElFakeableObjsName),
  fMuFakeableObjsName(ModNames::gkMuFakeableObjsName),
  fMCLeptonsName(ModNames::gkMCLeptonsName),
  fMCTausName(ModNames::gkMCTausName),
  fFakeEventHeadersName(ModNames::gkFakeEventHeadersName)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void GenFakesMod::LoadFakeRate()
{ 
  //Load FakeRate Probabilities.
  fFakeRate = new FakeRate(fElectronFRFilename,fMuonFRFilename,fElectronFRFunctionName,
                           fMuonFRFunctionName,fElectronFRHistName,
                           fMuonFRHistName,fUse2DFakeRate, fUseFitFunction );

}


//--------------------------------------------------------------------------------------------------
void GenFakesMod::Process()
{
  // Process entries of the tree.

  // get input Fakeable object collections
   const ElectronCol *ElFakeableObjs = 0;
   if (!fElFakeableObjsName.IsNull())
     ElFakeableObjs = GetObjThisEvt<ElectronCol>(fElFakeableObjsName);
   const MuonCol *MuFakeableObjs = 0;
   if (!fMuFakeableObjsName.IsNull())
     MuFakeableObjs = GetObjThisEvt<MuonCol>(fMuFakeableObjsName);

  // get input clean object collections
  const ElectronCol *CleanElectrons = 0;
  if (!fCleanElectronsName.IsNull())
    CleanElectrons = GetObjThisEvt<ElectronCol>(fCleanElectronsName);
  const MuonCol *CleanMuons = 0;
  if (!fCleanMuonsName.IsNull())
    CleanMuons = GetObjThisEvt<MuonCol>(fCleanMuonsName);
  const PhotonCol   *CleanPhotons   = 0;
  if (!fCleanPhotonsName.IsNull())
  CleanPhotons    = GetObjThisEvt<PhotonCol>(fCleanPhotonsName);
  const JetCol      *CleanJets       = 0;
  if (!fCleanJetsName.IsNull())
    CleanJets = GetObjThisEvt<JetCol>(fCleanJetsName);
  mithep::ParticleOArr *CleanLeptons = dynamic_cast<mithep::ParticleOArr*>
    (FindObjThisEvt(ModNames::gkMergedLeptonsName));

  //get monte carlo collections
  const MCParticleCol *GenLeptons = 0;
  if (!fMCLeptonsName.IsNull())
    GenLeptons = GetObjThisEvt<MCParticleCol>(fMCLeptonsName);
  const MCParticleCol *GenTaus = 0;
  if (!fMCTausName.IsNull())
    GenTaus = GetObjThisEvt<MCParticleCol>(fMCTausName);
  ObjArray<MCParticle> *GenLeptonsAndTaus = new ObjArray<MCParticle>;
  if (GenLeptons) {
    for (UInt_t i=0; i<GenLeptons->GetEntries(); i++)
      GenLeptonsAndTaus->Add(GenLeptons->At(i));
  }
  if (GenTaus) {
    for (UInt_t i=0; i<GenTaus->GetEntries(); i++)
      GenLeptonsAndTaus->Add(GenTaus->At(i));
  }

  // create final output collection
  ObjArray <FakeEventHeader> *FakeEventHeaders = new  ObjArray <FakeEventHeader> ;
  FakeEventHeaders->SetOwner(kTRUE);

  //initialize with one fake event containing no fake objects and all jets.
  FakeEventHeader *initialFakeEvent = new FakeEventHeader();
  for (UInt_t j=0;j<CleanJets->GetEntries();j++)
    initialFakeEvent->AddJets(CleanJets->At(j));

  FakeEventHeaders->AddOwned(initialFakeEvent);
  
  // *****************************************************************************************
  // Fake into Muons
  // Loop through all Muon Fakeable objects and consider the fake possibility.
  // *****************************************************************************************
  for (UInt_t n = 0; n < MuFakeableObjs->GetEntries();n++) {

    //make temporary fake event headers array
    ObjArray <FakeEventHeader> *tempFakeEventHeaders = new ObjArray <FakeEventHeader> ;
    tempFakeEventHeaders->SetOwner(kTRUE);

    //loop over all fake events generated so far - and perform an additional fake if necessary
    for (UInt_t i=0; i<FakeEventHeaders->GetEntries();i++) {           

      // *****************************************************************************************
      // Determine if the fakeable object was a clean lepton
      // *****************************************************************************************
      Bool_t isCleanLepton = false;
      for (UInt_t j = 0; j < CleanLeptons->GetEntries() ; j++) {
        Double_t deltaR = MathUtils::DeltaR(MuFakeableObjs->At(n)->Phi(),
                                            MuFakeableObjs->At(n)->Eta(),
                                            CleanLeptons->At(j)->Phi(), CleanLeptons->At(j)->Eta());

        if (deltaR < 0.3) {
          isCleanLepton = true;
          break;
        }
      }

      // *****************************************************************************************
      // Determine if the fakeable object was a real muon from Monte Carlo        
      // *****************************************************************************************
      Bool_t isGenMuon = false;
      for (UInt_t l=0; l<GenLeptonsAndTaus->GetEntries(); l++) {
        if (MathUtils::DeltaR(MuFakeableObjs->At(n)->Phi(), MuFakeableObjs->At(n)->Eta(),
                              GenLeptonsAndTaus->At(l)->Phi(), 
                              GenLeptonsAndTaus->At(l)->Eta()) < 0.3) {
          isGenMuon = true;
          break;
        }
      }

      //this is used to determine the weight of the unfaked event.
      double totalCumulativeFakeProbability = 0.0;

      // *****************************************************************************************
      // Perform Muon Fake
      // *****************************************************************************************
      
      //match fake to one of the jets
      int fakeToJetMatch = -1; //index of the jet that matches to the fake
      double minDR = 5000;
      for (UInt_t jj=0;jj<FakeEventHeaders->At(i)->NJets();jj++) {
        Double_t deltaR = MathUtils::DeltaR(FakeEventHeaders->At(i)->UnfakedJet(jj)->Mom(),
                                            MuFakeableObjs->At(n)->Mom());
        if (deltaR < minDR) {
          minDR = deltaR;
          fakeToJetMatch = jj;
        }
      }
      if (!(minDR < 0.5)) {
        fakeToJetMatch = -1;
      }
      
      //Obtain the muon FakeRate 
      Double_t muonFakeProb = 0.0;
      Double_t muonFakeProbLowError = 0.0;
      Double_t muonFakeProbHighError = 0.0;
      if(fFakeRate) {
        muonFakeProb = fFakeRate->MuonFakeRate(MuFakeableObjs->At(n)->Et(),
                                               MuFakeableObjs->At(n)->Eta(),
                                               MuFakeableObjs->At(n)->Phi());
        muonFakeProbLowError = fFakeRate->MuonFakeRateError(MuFakeableObjs->At(n)->Et(),
                                                            MuFakeableObjs->At(n)->Eta(),
                                                            MuFakeableObjs->At(n)->Phi());
        muonFakeProbHighError = fFakeRate->MuonFakeRateError(MuFakeableObjs->At(n)->Et(),
                                                             MuFakeableObjs->At(n)->Eta(),
                                                             MuFakeableObjs->At(n)->Phi());                
      } else {
        cerr << "Error: fFakeRate is a NULL pointer.\n";
        assert(false);
      }

      //only fake into a muon if the fakeable object did not match to a clean lepton
      if (!isCleanLepton) {

        //create new fake event header
        FakeEventHeader *fakeMuonEvent = new FakeEventHeader();
        fakeMuonEvent->SetWeight(FakeEventHeaders->At(i)->Weight() * muonFakeProb);
        fakeMuonEvent->SetWeightLowError(FakeEventHeaders->At(i)->Weight()*muonFakeProb*
                                         TMath::Sqrt((FakeEventHeaders->At(i)->WeightLowError()/
                                                      FakeEventHeaders->At(i)->Weight())*
                                                     (FakeEventHeaders->At(i)->WeightLowError()/
                                                      FakeEventHeaders->At(i)->Weight()) +
                                                     (muonFakeProbLowError/muonFakeProb)*
                                                     (muonFakeProbLowError/muonFakeProb)                                                      
                                           ));
        fakeMuonEvent->SetWeightHighError(FakeEventHeaders->At(i)->Weight()*muonFakeProb*
                                         TMath::Sqrt((FakeEventHeaders->At(i)->WeightHighError()/
                                                      FakeEventHeaders->At(i)->Weight())*
                                                     (FakeEventHeaders->At(i)->WeightHighError()/
                                                      FakeEventHeaders->At(i)->Weight()) +
                                                     (muonFakeProbHighError/muonFakeProb)*
                                                     (muonFakeProbHighError/muonFakeProb)                                                      
                                           ));

        //add all previous fakes
        for (UInt_t f=0;f<FakeEventHeaders->At(i)->FakeObjsSize();f++) {
          fakeMuonEvent->AddFakeObject(FakeEventHeaders->At(i)->FakeObj(f));        
        }
        //add new fake
        fakeMuonEvent->AddFakeObject(MuFakeableObjs->At(n),kMuon,isCleanLepton,isGenMuon);

        //add all previous jets except the one matching to the new fake
        for (UInt_t jj=0;jj<FakeEventHeaders->At(i)->NJets();jj++) {
          if (jj != fakeToJetMatch)
            fakeMuonEvent->AddJets(FakeEventHeaders->At(i)->UnfakedJet(jj));
        }
        
        //add fake event to the temporary fake event header array
        tempFakeEventHeaders->AddOwned(fakeMuonEvent);
        
        //increase cumulative fake probability
        totalCumulativeFakeProbability += muonFakeProb;
      }

      // *****************************************************************************************
      // Do not fake into Muon
      // *****************************************************************************************
      //create new fake event header
      FakeEventHeader *notFakeEvent = new FakeEventHeader();
      notFakeEvent->SetWeight(FakeEventHeaders->At(i)->Weight() * 
                              (1-totalCumulativeFakeProbability));
      //add previous fakes
      for (UInt_t f=0;f<FakeEventHeaders->At(i)->FakeObjsSize();f++) {
        notFakeEvent->AddFakeObject(FakeEventHeaders->At(i)->FakeObj(f));        
      }
      //add previous jets
      for (UInt_t jj=0;jj<FakeEventHeaders->At(i)->NJets();jj++) {
          notFakeEvent->AddJets(FakeEventHeaders->At(i)->UnfakedJet(jj));
      }
      tempFakeEventHeaders->AddOwned(notFakeEvent);   

    } //loop over all current fake event headers

    //replace current fake event headers with the new temporary ones.
    delete FakeEventHeaders;
    FakeEventHeaders = tempFakeEventHeaders;
  } //loop over all muon fakeable objects


  // *****************************************************************************************
  // Fake into Electrons
  // Loop through all Electron Fakeable objects and consider the fake possibility.
  // *****************************************************************************************
  for (UInt_t n = 0; n < ElFakeableObjs->GetEntries();n++) {

    //make temporary fake event headers array
    ObjArray <FakeEventHeader> *tempFakeEventHeaders = new ObjArray <FakeEventHeader> ;
    tempFakeEventHeaders->SetOwner(kTRUE);

    //loop over all fake events generated so far - and perform an additional fake if necessary
    for (UInt_t i=0; i<FakeEventHeaders->GetEntries();i++) {      
      
      // *****************************************************************************************
      // Determine if the fakeable object was a clean lepton
      // *****************************************************************************************
      Bool_t isCleanLepton = false;
      for (UInt_t j = 0; j < CleanLeptons->GetEntries() ; j++) {
        Double_t deltaR = MathUtils::DeltaR(ElFakeableObjs->At(n)->Phi(),
                                            ElFakeableObjs->At(n)->Eta(),
                                            CleanLeptons->At(j)->Phi(), CleanLeptons->At(j)->Eta());

        if (deltaR < 0.3) {
          isCleanLepton = true;
          break;
        }
      }

      // *****************************************************************************************
      // Determine if the fakeable object was a real electron from Monte Carlo        
      // *****************************************************************************************
      Bool_t isGenElectron = false;
      for (UInt_t l=0; l<GenLeptonsAndTaus->GetEntries(); l++) {
        if (MathUtils::DeltaR(ElFakeableObjs->At(n)->Phi(), 
                              ElFakeableObjs->At(n)->Eta(),
                              GenLeptonsAndTaus->At(l)->Phi(), 
                              GenLeptonsAndTaus->At(l)->Eta()) < 0.3) {
          isGenElectron = true;
        }
      }


      //this is used to determine the weight of the unfaked event.
      double totalCumulativeFakeProbability = 0.0;
      
      // *****************************************************************************************
      // Determine if the current electron fakeable object already corresponds to one of the
      // fake muons already in the FakeEventHeader, determined based on deltaR proximity.
      // If the current electron fakeable object corresponds to one of the fake muon, then
      // we do not allow it to fake an electron, since one denominator cannot fake two leptons.
      // *****************************************************************************************
      Bool_t alreadyFaked = false;
      for (UInt_t f = 0; f < FakeEventHeaders->At(i)->FakeObjsSize() ; f++) {
        double deltaR = MathUtils::DeltaR(FakeEventHeaders->At(i)->FakeObj(f)->Mom(),
                                          ElFakeableObjs->At(n)->Mom());
        if (deltaR < 0.3)
          alreadyFaked = true;
      }
      if (!alreadyFaked) {

        // *****************************************************************************************
        // Perform Electron Fake
        // *****************************************************************************************
  
        //match fake to one of the jets
        int fakeToJetMatch = -1; //index of the jet that matches to the fake
        double minDR = 5000;
        for (UInt_t jj=0;jj<FakeEventHeaders->At(i)->NJets();jj++) {          
          Double_t deltaR = MathUtils::DeltaR(FakeEventHeaders->At(i)->UnfakedJet(jj)->Mom(),
                                              ElFakeableObjs->At(n)->Mom());
          if (deltaR < minDR) {
            minDR = deltaR;
            fakeToJetMatch = jj;
          }
        }
        if (!(minDR < 0.5)) { 
          fakeToJetMatch = -1;
        }

        //Obtain the electron FakeRate 
        Double_t electronFakeProb = 0.0;
        Double_t electronFakeProbLowError = 0.0;
        Double_t electronFakeProbHighError = 0.0;
        if(fFakeRate) {
          electronFakeProb = fFakeRate->ElectronFakeRate(ElFakeableObjs->At(n)->Et(),
                                                         ElFakeableObjs->At(n)->Eta(),
                                                         ElFakeableObjs->At(n)->Phi());
          electronFakeProbLowError = fFakeRate->ElectronFakeRateError(
            ElFakeableObjs->At(n)->Et(),
            ElFakeableObjs->At(n)->Eta(),
            ElFakeableObjs->At(n)->Phi());
          electronFakeProbHighError = fFakeRate->ElectronFakeRateError(
            ElFakeableObjs->At(n)->Et(),
            ElFakeableObjs->At(n)->Eta(),
            ElFakeableObjs->At(n)->Phi());                          
        } else {
          cerr << "Error: fFakeRate is a NULL pointer.\n";
          assert(false);
        }

        //only fake into a muon if the fakeable object did not match to a clean lepton
        if (!isCleanLepton) {
          //create new fake event header
          FakeEventHeader *fakeElectronEvent = new FakeEventHeader();
          fakeElectronEvent->SetWeight(FakeEventHeaders->At(i)->Weight() * electronFakeProb);
          fakeElectronEvent->SetWeightLowError(FakeEventHeaders->At(i)->Weight()*electronFakeProb*
                                         TMath::Sqrt((FakeEventHeaders->At(i)->WeightLowError()/
                                                      FakeEventHeaders->At(i)->Weight())*
                                                     (FakeEventHeaders->At(i)->WeightLowError()/
                                                      FakeEventHeaders->At(i)->Weight()) +
                                                     (electronFakeProbLowError/electronFakeProb)*
                                                     (electronFakeProbLowError/electronFakeProb)                                                      
                                           ));
          fakeElectronEvent->SetWeightHighError(FakeEventHeaders->At(i)->Weight()*electronFakeProb*
                                         TMath::Sqrt((FakeEventHeaders->At(i)->WeightHighError()/
                                                      FakeEventHeaders->At(i)->Weight())*
                                                     (FakeEventHeaders->At(i)->WeightHighError()/
                                                      FakeEventHeaders->At(i)->Weight()) +
                                                     (electronFakeProbHighError/electronFakeProb)*
                                                     (electronFakeProbHighError/electronFakeProb)                                                      
                                           ));


          //add previous fakes
          for (UInt_t f=0;f<FakeEventHeaders->At(i)->FakeObjsSize();f++) {
            fakeElectronEvent->AddFakeObject(FakeEventHeaders->At(i)->FakeObj(f));
          }
          //add the new fake
          fakeElectronEvent->AddFakeObject(ElFakeableObjs->At(n),
                                           kElectron,isCleanLepton,isGenElectron);
          //add previous jets that do not match to the new fake
          for (UInt_t jj=0;jj<FakeEventHeaders->At(i)->NJets();jj++) {
            if (jj != fakeToJetMatch)
              fakeElectronEvent->AddJets(FakeEventHeaders->At(i)->UnfakedJet(jj));
          }
          
          //add fake event to the temporary fake event header array
          tempFakeEventHeaders->AddOwned(fakeElectronEvent);
          //increase cumulative fake probability
          totalCumulativeFakeProbability += electronFakeProb;
        }
      }
      
      // *****************************************************************************************
      // Do not fake into anything
      // *****************************************************************************************
      //create new fake event header
      FakeEventHeader *notFakeEvent = new FakeEventHeader();
      notFakeEvent->SetWeight(FakeEventHeaders->At(i)->Weight() * 
                              (1-totalCumulativeFakeProbability));
      //add previous fakes
      for (UInt_t f=0;f<FakeEventHeaders->At(i)->FakeObjsSize();f++) {
        notFakeEvent->AddFakeObject(FakeEventHeaders->At(i)->FakeObj(f));        
      }
      //add previous jets
      for (UInt_t jj=0;jj<FakeEventHeaders->At(i)->NJets();jj++) {
          notFakeEvent->AddJets(FakeEventHeaders->At(i)->UnfakedJet(jj));
      }
      tempFakeEventHeaders->AddOwned(notFakeEvent);      

    } //for all current fake event headers

    //replace current fake event headers with the new temporary ones.
    delete FakeEventHeaders;
    FakeEventHeaders = tempFakeEventHeaders;
  } //loop over all fakeable objects

  // export FakeEventHeaders for other modules to use
  FakeEventHeaders->SetName(fFakeEventHeadersName);
  AddObjThisEvt(FakeEventHeaders);

  //delete temporary collections
  delete GenLeptonsAndTaus;
}
