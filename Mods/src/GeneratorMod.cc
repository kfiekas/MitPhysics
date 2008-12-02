// $Id: GeneratorMod.cc,v 1.12 2008/12/01 11:44:34 ceballos Exp $

#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include <TH1D.h>
#include <TH2D.h>

using namespace mithep;

ClassImp(mithep::GeneratorMod)

//--------------------------------------------------------------------------------------------------
GeneratorMod::GeneratorMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fFillHist(kFALSE),
  fMCPartName(Names::gkMCPartBrn),
  fMCLeptonsName(ModNames::gkMCLeptonsName),
  fMCAllLeptonsName(ModNames::gkMCAllLeptonsName),
  fMCTausName(ModNames::gkMCTausName),
  fMCNeutrinosName(ModNames::gkMCNeutrinosName),
  fMCQuarksName(ModNames::gkMCQuarksName),
  fMCqqHsName(ModNames::gkMCqqHsName),
  fMCBosonsName(ModNames::gkMCBosonsName),
  fMCPhotonsName(ModNames::gkMCPhotonsName),
  fParticles(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void GeneratorMod::Process()
{
  // Process entries of the tree.

  // these arrays will be filled in the loop of particles
  MCParticleOArr *GenLeptons    = new MCParticleOArr;
  GenLeptons->SetName(fMCLeptonsName);
  MCParticleOArr *GenAllLeptons = new MCParticleOArr;
  GenAllLeptons->SetName(fMCAllLeptonsName);
  MCParticleOArr *GenTaus       = new MCParticleOArr; 
  GenTaus->SetName(fMCTausName);
  GenTaus->SetOwner(kTRUE);
  MCParticleOArr *GenNeutrinos  = new MCParticleOArr;
  GenNeutrinos->SetName(fMCNeutrinosName);
  MCParticleOArr *GenQuarks     = new MCParticleOArr;
  GenQuarks->SetName(fMCQuarksName);
  MCParticleOArr *GenqqHs       = new MCParticleOArr;
  GenqqHs->SetName(fMCqqHsName);
  MCParticleOArr *GenBosons     = new MCParticleOArr;
  GenBosons->SetName(fMCBosonsName);
  MCParticleOArr *GenPhotons    = new MCParticleOArr;
  GenPhotons->SetName(fMCPhotonsName);

  // load MCParticle branch
  LoadBranch(fMCPartName);

  Bool_t isqqH = kFALSE;
  for (UInt_t i=0; i<fParticles->GetEntries(); ++i) {
    const MCParticle *p = fParticles->At(i);

    if (!p->IsGenerated()) continue;

    // muons/electrons from W/Z decays
    if ((p->Is(MCParticle::kEl) || p->Is(MCParticle::kMu)) && p->Status() == 1) {
      if (p->Pt() > 3.0 && TMath::Abs(p->Eta()) < 3.0) {
        GenAllLeptons->Add(p);
      }
      Bool_t isGoodLepton = kFALSE;
      const MCParticle *pm = p;
      while (pm->HasMother() && isGoodLepton == kFALSE) {
        if (pm->Mother()->Is(MCParticle::kZ) || pm->Mother()->Is(MCParticle::kW)) {
          GenLeptons->Add(p);
          isGoodLepton = kTRUE;
          break;
        } else if (pm->Mother()->Is(MCParticle::kPi0) || pm->Mother()->Is(MCParticle::kEta)) {
          // this is fake, but it is a trick to get rid of these cases and abort the loop
          isGoodLepton = kTRUE;
          break;
        } 
        pm = pm->Mother();
      }
    }

    // hadronic taus
    else if (p->Is(MCParticle::kTau) && p->Status() == 2) {
      if (!p->HasDaughter(MCParticle::kEl) && !p->HasDaughter(MCParticle::kMu)) {
        const MCParticle *tv = p->FindDaughter(MCParticle::kTauNu);
        if (tv) {
          MCParticle *pm_f = new MCParticle(*p);
          pm_f->SetMom(p->Px()-tv->Px(), p->Py()-tv->Py(),
                       p->Pz()-tv->Pz(), p->E()-tv->E());
          GenTaus->AddOwned(pm_f);
        } else {
          SendError(kWarning, "Process", "Could not find tau neutrino!");
        }
      }
    }

    // neutrinos
    else if (p->Status() == 1 && p->IsNeutrino()) {
      GenNeutrinos->Add(p);
    }

    // quarks from W/Z decays or top particles
    else if (p->IsQuark() && p->HasMother()) {
      if (p->Mother()->Is(MCParticle::kZ) || p->Mother()->Is(MCParticle::kW) ||
          p->Is(MCParticle::kTop)         || p->Mother()->Is(MCParticle::kTop)) {
        GenQuarks->Add(p);
      }
    }

    // qqH, information about the forward jets
    else if (isqqH == kFALSE && p->Is(MCParticle::kH)) {
      isqqH = kTRUE;
      MCParticle *pq1 = fParticles->At(i-1);
      MCParticle *pq2 = fParticles->At(i-2);
      if (!pq1 || !pq2) {
          SendError(kWarning, "Process", "Could not find quark pair!");
      } else if (pq1->IsQuark()   && pq2->IsQuark()   && 
                 pq1->HasMother() && pq2->HasMother() &&
                 pq1->Mother()->PdgId() == p->Mother()->PdgId() &&
                 pq2->Mother()->PdgId() == p->Mother()->PdgId()) {
        GenqqHs->Add(pq1);
        GenqqHs->Add(pq2);
      }
      if (p->Status() == 3)  
        GenBosons->Add(p); // take higgs boson in  account here rather in next else if 
    }

    // information about bosons: W, Z, h, Z', W', H0, A0, H+
    else if (p->Status() == 3 &&
             (p->Is(MCParticle::kZ)  || p->Is(MCParticle::kW)   || p->Is(MCParticle::kH) ||
              p->Is(MCParticle::kZp) || p->Is(MCParticle::kZpp) ||
              p->Is(MCParticle::kH0) || p->Is(MCParticle::kA0)  || p->Is(MCParticle::kHp))) {
      GenBosons->Add(p);
    }

    // photons, only with Pt > 15
    else if (p->Status() == 1 && p->Is(MCParticle::kGamma) &&
             p->Pt() > 15) {
      GenPhotons->Add(p);
    }

  }

  // add objects to this event for other modules to use
  AddObjThisEvt(GenLeptons);  
  AddObjThisEvt(GenAllLeptons);
  AddObjThisEvt(GenTaus);
  AddObjThisEvt(GenNeutrinos);
  AddObjThisEvt(GenQuarks);
  AddObjThisEvt(GenqqHs);
  AddObjThisEvt(GenBosons);
  AddObjThisEvt(GenPhotons);

  // fill histograms if requested
  if (fFillHist) {

    // leptons
    hDGenLeptons[0]->Fill(GenLeptons->GetEntries());
    Int_t   idxMaxLep[2] = {-1, -1};
    Double_t ptMaxLep[2] = {-1, -1};
    for(UInt_t i=0; i<GenLeptons->GetEntries(); i++) {
      hDGenLeptons[1]->Fill(GenLeptons->At(i)->Pt());
      hDGenLeptons[2]->Fill(TMath::Abs(GenLeptons->At(i)->Eta()));
      hDGenLeptons[3]->Fill(GenLeptons->At(i)->PhiDeg());
      for(UInt_t j=i+1; j<GenLeptons->GetEntries(); j++) {
        CompositeParticle *dilepton = new CompositeParticle();
        dilepton->AddDaughter(GenLeptons->At(i));
        dilepton->AddDaughter(GenLeptons->At(j));
        hDGenLeptons[4]->Fill(dilepton->Mass());
	delete dilepton;
      }
      // selecting the two highest pt leptons
      if     (GenLeptons->At(i)->Pt() > ptMaxLep[0]) {
        ptMaxLep[1]  = ptMaxLep[0];
	idxMaxLep[1] = idxMaxLep[0];
        ptMaxLep[0]  = GenLeptons->At(i)->Pt();
	idxMaxLep[0] = i;
      }
      else if (GenLeptons->At(i)->Pt() > ptMaxLep[1]) {
        ptMaxLep[1]  = GenLeptons->At(i)->Pt();
	idxMaxLep[1] = i;
      }
    }
    // looking at events with at least two leptons
    if (ptMaxLep[0] > 0 && ptMaxLep[1] > 0) {
      hDGenLeptons[5]->Fill(TMath::Min(TMath::Max(TMath::Abs(GenLeptons->At(idxMaxLep[0])->Eta()),
                                                  TMath::Abs(GenLeptons->At(idxMaxLep[1])->Eta())),
                                       4.999));
      hDGenLeptons[6]->Fill(TMath::Min(TMath::Min(TMath::Abs(GenLeptons->At(idxMaxLep[0])->Eta()),
                                                  TMath::Abs(GenLeptons->At(idxMaxLep[1])->Eta())),
                                       4.999));
      if (TMath::Abs(GenLeptons->At(idxMaxLep[0])->Eta()) < 2.5 &&
          TMath::Abs(GenLeptons->At(idxMaxLep[1])->Eta()) < 2.5) {
        hDGenLeptons[7]->Fill(TMath::Min(GenLeptons->At(idxMaxLep[0])->Pt(),199.999));
        if (GenLeptons->At(idxMaxLep[0])->Pt() > 20.0) {
          hDGenLeptons[8]->Fill(TMath::Min(GenLeptons->At(idxMaxLep[1])->Pt(),199.999));
          if (GenLeptons->At(idxMaxLep[1])->Pt() > 10.0) {
            CompositeParticle *dilepton = new CompositeParticle();
            dilepton->AddDaughter(GenLeptons->At(idxMaxLep[0]));
            dilepton->AddDaughter(GenLeptons->At(idxMaxLep[1]));
            hDGenLeptons[9]->Fill(TMath::Min(dilepton->Mass(),999.999));
	    hDGenLeptons[10]->Fill(MathUtils::DeltaPhi(GenLeptons->At(idxMaxLep[0])->Phi(),
	                                               GenLeptons->At(idxMaxLep[1])->Phi())
						       * 180./ TMath::Pi());
	    delete dilepton;
	  }
	}
      }
    }

    // all leptons
    hDGenAllLeptons[0]->Fill(GenAllLeptons->GetEntries());
    for(UInt_t i=0; i<GenAllLeptons->GetEntries(); i++) {
      hDGenAllLeptons[1]->Fill(GenAllLeptons->At(i)->Pt());
      hDGenAllLeptons[2]->Fill(GenAllLeptons->At(i)->Eta());
      hDGenAllLeptons[3]->Fill(GenAllLeptons->At(i)->PhiDeg());
    }

    // taus
    hDGenTaus[0]->Fill(GenTaus->GetEntries());
    for(UInt_t i=0; i<GenTaus->GetEntries(); i++) {
      hDGenTaus[1]->Fill(GenTaus->At(i)->Pt());
      hDGenTaus[2]->Fill(GenTaus->At(i)->Eta());
      hDGenTaus[3]->Fill(GenTaus->At(i)->PhiDeg());
    }

    // neutrinos
    hDGenNeutrinos[0]->Fill(GenNeutrinos->GetEntries());
    CompositeParticle *neutrinoTotal = new CompositeParticle();
    for(UInt_t i=0; i<GenNeutrinos->GetEntries(); i++) {
      if (GenNeutrinos->At(i)->HasMother())
        neutrinoTotal->AddDaughter(GenNeutrinos->At(i));
    }
    if (GenNeutrinos->GetEntries() > 0) {
      hDGenNeutrinos[1]->Fill(neutrinoTotal->Pt());
      hDGenNeutrinos[2]->Fill(neutrinoTotal->Eta());
      hDGenNeutrinos[3]->Fill(neutrinoTotal->PhiDeg());    
    }
    delete neutrinoTotal;
    
    // quarks
    hDGenQuarks[0]->Fill(GenQuarks->GetEntries());
    for(UInt_t i=0; i<GenQuarks->GetEntries(); i++) {
      for(UInt_t j=i+1; j<GenQuarks->GetEntries(); j++) {
        CompositeParticle *dijet = new CompositeParticle();
        dijet->AddDaughter(GenQuarks->At(i));
        dijet->AddDaughter(GenQuarks->At(j));
        hDGenQuarks[1]->Fill(dijet->Pt());
        hDGenQuarks[2]->Fill(dijet->Mass());
	if (TMath::Abs(GenQuarks->At(i)->Eta()) < 2.5 && 
            TMath::Abs(GenQuarks->At(j)->Eta()) < 2.5) {
          hDGenQuarks[3]->Fill(dijet->Pt());
          hDGenQuarks[4]->Fill(dijet->Mass());
	}
	delete dijet;
      }
      // b quark info
      if    (GenQuarks->At(i)->AbsPdgId() == 5) {
        hDGenQuarks[5]->Fill(GenQuarks->At(i)->Pt());
        hDGenQuarks[6]->Fill(GenQuarks->At(i)->Eta());      
        hDGenQuarks[7]->Fill(GenQuarks->At(i)->Phi());      
        if (ptMaxLep[0] > 20 && ptMaxLep[1] > 15) {
          if (TMath::Abs(GenLeptons->At(idxMaxLep[0])->Eta()) < 2.5 &&
              TMath::Abs(GenLeptons->At(idxMaxLep[1])->Eta()) < 2.5) {
            hDGenQuarks[8]->Fill(GenQuarks->At(i)->Pt());	
            hDGenQuarks[9]->Fill(GenQuarks->At(i)->Eta());      
            hDGenQuarks[10]->Fill(GenQuarks->At(i)->Phi());	
	  }
	}
      }
      // t quark info
      else if (GenQuarks->At(i)->AbsPdgId() == 6) {
        hDGenQuarks[11]->Fill(GenQuarks->At(i)->Pt());
        hDGenQuarks[12]->Fill(GenQuarks->At(i)->Eta());      
        hDGenQuarks[13]->Fill(GenQuarks->At(i)->Phi());      
      }
      // light quark info
      else {
        hDGenQuarks[14]->Fill(GenQuarks->At(i)->Pt());
        hDGenQuarks[15]->Fill(GenQuarks->At(i)->Eta());      
        hDGenQuarks[16]->Fill(GenQuarks->At(i)->Phi());      
      }
    }

    // wbf
    if (GenqqHs->GetEntries() == 2) {
      hDGenWBF[0]->Fill(MathUtils::DeltaPhi(GenqqHs->At(0)->Phi(),
    					    GenqqHs->At(1)->Phi()) * 180./ TMath::Pi());
      hDGenWBF[1]->Fill(TMath::Abs(GenqqHs->At(0)->Eta()-GenqqHs->At(1)->Eta()));
      hDGenWBF[2]->Fill(TMath::Max(GenqqHs->At(0)->Pt(),GenqqHs->At(1)->Pt()));
      hDGenWBF[3]->Fill(TMath::Min(GenqqHs->At(0)->Pt(),GenqqHs->At(1)->Pt()));
      CompositeParticle *diqq = new CompositeParticle();
      diqq->AddDaughter(GenqqHs->At(0));
      diqq->AddDaughter(GenqqHs->At(1));
      hDGenWBF[4]->Fill(diqq->Mass());
      delete diqq;
    }

    // bosons
    hDGenBosons[0]->Fill(GenBosons->GetEntries());
    for(UInt_t i=0; i<GenBosons->GetEntries(); i++) {
      hDGenBosons[1]->Fill(GenBosons->At(i)->Pt());
      hDGenBosons[2]->Fill(GenBosons->At(i)->Eta());
      hDGenBosons[3]->Fill(GenBosons->At(i)->Mass());
    } 

    // photons
    hDGenPhotons[0]->Fill(GenPhotons->GetEntries());
    for(UInt_t i=0; i<GenPhotons->GetEntries(); i++) {
      hDGenPhotons[1]->Fill(GenPhotons->At(i)->Pt());
      hDGenPhotons[2]->Fill(GenPhotons->At(i)->Eta());
    } 
  }
}

//--------------------------------------------------------------------------------------------------
void GeneratorMod::SlaveBegin()
{
  // Book branch and histograms if wanted.

  ReqBranch(fMCPartName, fParticles);

  // fill histograms
  if (fFillHist == kTRUE) {
    char sb[1024];

    // leptons from W
    sprintf(sb,"hDGenLeptons_%d", 0);  hDGenLeptons[0]  = new TH1D(sb,sb,10,-0.5,9.5); 
    sprintf(sb,"hDGenLeptons_%d", 1);  hDGenLeptons[1]  = new TH1D(sb,sb,100,0.0,200.0); 
    sprintf(sb,"hDGenLeptons_%d", 2);  hDGenLeptons[2]  = new TH1D(sb,sb,50,0.0,5.0); 
    sprintf(sb,"hDGenLeptons_%d", 3);  hDGenLeptons[3]  = new TH1D(sb,sb,90,0.0,180.0); 
    sprintf(sb,"hDGenLeptons_%d", 4);  hDGenLeptons[4]  = new TH1D(sb,sb,1000,0.0,1000.0); 
    sprintf(sb,"hDGenLeptons_%d", 5);  hDGenLeptons[5]  = new TH1D(sb,sb,50,0.0,5.0); 
    sprintf(sb,"hDGenLeptons_%d", 6);  hDGenLeptons[6]  = new TH1D(sb,sb,50,0.0,5.0); 
    sprintf(sb,"hDGenLeptons_%d", 7);  hDGenLeptons[7]  = new TH1D(sb,sb,100,0.0,200.0); 
    sprintf(sb,"hDGenLeptons_%d", 8);  hDGenLeptons[8]  = new TH1D(sb,sb,100,0.0,200.0); 
    sprintf(sb,"hDGenLeptons_%d", 9);  hDGenLeptons[9]  = new TH1D(sb,sb,1000,0.0,1000.0); 
    sprintf(sb,"hDGenLeptons_%d",10);  hDGenLeptons[10] = new TH1D(sb,sb,90,0.0,180.0); 
    for(Int_t i=0; i<11; i++) AddOutput(hDGenLeptons[i]);

    // all leptons
    sprintf(sb,"hDGenAllLeptons_%d", 0);  hDGenAllLeptons[0]  = new TH1D(sb,sb,10,-0.5,9.5); 
    sprintf(sb,"hDGenAllLeptons_%d", 1);  hDGenAllLeptons[1]  = new TH1D(sb,sb,100,0.0,200.0); 
    sprintf(sb,"hDGenAllLeptons_%d", 2);  hDGenAllLeptons[2]  = new TH1D(sb,sb,60,-3.0,3.0); 
    sprintf(sb,"hDGenAllLeptons_%d", 3);  hDGenAllLeptons[3]  = new TH1D(sb,sb,90,0.0,180.0); 
    for(Int_t i=0; i<4; i++) AddOutput(hDGenAllLeptons[i]);

    // taus
    sprintf(sb,"hDGenTaus_%d", 0);  hDGenTaus[0]  = new TH1D(sb,sb,10,-0.5,9.5); 
    sprintf(sb,"hDGenTaus_%d", 1);  hDGenTaus[1]  = new TH1D(sb,sb,100,0.0,200.0); 
    sprintf(sb,"hDGenTaus_%d", 2);  hDGenTaus[2]  = new TH1D(sb,sb,100,-5.0,5.0); 
    sprintf(sb,"hDGenTaus_%d", 3);  hDGenTaus[3]  = new TH1D(sb,sb,90,0.0,180.0); 
    for(Int_t i=0; i<4; i++) AddOutput(hDGenTaus[i]);

    // neutrinos
    sprintf(sb,"hDGenNeutrinos_%d", 0);  hDGenNeutrinos[0]  = new TH1D(sb,sb,10,-0.5,9.5); 
    sprintf(sb,"hDGenNeutrinos_%d", 1);  hDGenNeutrinos[1]  = new TH1D(sb,sb,100,0.0,200.0); 
    sprintf(sb,"hDGenNeutrinos_%d", 2);  hDGenNeutrinos[2]  = new TH1D(sb,sb,100,-5.0,5.0); 
    sprintf(sb,"hDGenNeutrinos_%d", 3);  hDGenNeutrinos[3]  = new TH1D(sb,sb,90,0.0,180.0); 
    for(Int_t i=0; i<4; i++) AddOutput(hDGenNeutrinos[i]);

    // quarks
    sprintf(sb,"hDGenQuarks_%d", 0);  hDGenQuarks[0]  = new TH1D(sb,sb,10,-0.5,9.5); 
    sprintf(sb,"hDGenQuarks_%d", 1);  hDGenQuarks[1]  = new TH1D(sb,sb,200,0.0,400.); 
    sprintf(sb,"hDGenQuarks_%d", 2);  hDGenQuarks[2]  = new TH1D(sb,sb,2000,0.0,2000.);
    sprintf(sb,"hDGenQuarks_%d", 3);  hDGenQuarks[3]  = new TH1D(sb,sb,200,0.0,400.); 
    sprintf(sb,"hDGenQuarks_%d", 4);  hDGenQuarks[4]  = new TH1D(sb,sb,2000,0.0,2000.);
    sprintf(sb,"hDGenQuarks_%d", 5);  hDGenQuarks[5]  = new TH1D(sb,sb,200,0.0,400.); 
    sprintf(sb,"hDGenQuarks_%d", 6);  hDGenQuarks[6]  = new TH1D(sb,sb,200,-10.0,10.); 
    sprintf(sb,"hDGenQuarks_%d", 7);  hDGenQuarks[7]  = new TH1D(sb,sb,200,-TMath::Pi(),
                                                                 TMath::Pi()); 
    sprintf(sb,"hDGenQuarks_%d", 8);  hDGenQuarks[8]  = new TH1D(sb,sb,200,0.0,400.); 
    sprintf(sb,"hDGenQuarks_%d", 9);  hDGenQuarks[9]  = new TH1D(sb,sb,200,-10.0,10.); 
    sprintf(sb,"hDGenQuarks_%d",10);  hDGenQuarks[10] = new TH1D(sb,sb,200,-TMath::Pi(),
                                                                 TMath::Pi()); 
    sprintf(sb,"hDGenQuarks_%d",11);  hDGenQuarks[11] = new TH1D(sb,sb,200,0.0,400.); 
    sprintf(sb,"hDGenQuarks_%d",12);  hDGenQuarks[12] = new TH1D(sb,sb,200,-10.0,10.); 
    sprintf(sb,"hDGenQuarks_%d",13);  hDGenQuarks[13] = new TH1D(sb,sb,200,-TMath::Pi(),
                                                                 TMath::Pi()); 
    sprintf(sb,"hDGenQuarks_%d",14);  hDGenQuarks[14] = new TH1D(sb,sb,200,0.0,400.); 
    sprintf(sb,"hDGenQuarks_%d",15);  hDGenQuarks[15] = new TH1D(sb,sb,200,-10.0,10.); 
    sprintf(sb,"hDGenQuarks_%d",16);  hDGenQuarks[16] = new TH1D(sb,sb,200,-TMath::Pi(),
                                                                 TMath::Pi()); 
    for(Int_t i=0; i<17; i++) AddOutput(hDGenQuarks[i]);

    // qqH
    sprintf(sb,"hDGenWBF_%d", 0);  hDGenWBF[0]  = new TH1D(sb,sb,90,0.0,180.);
    sprintf(sb,"hDGenWBF_%d", 1);  hDGenWBF[1]  = new TH1D(sb,sb,100,0.0,10.);   
    sprintf(sb,"hDGenWBF_%d", 2);  hDGenWBF[2]  = new TH1D(sb,sb,200,0.0,400.); 
    sprintf(sb,"hDGenWBF_%d", 3);  hDGenWBF[3]  = new TH1D(sb,sb,200,0.0,400.);
    sprintf(sb,"hDGenWBF_%d", 4);  hDGenWBF[4]  = new TH1D(sb,sb,200,0.0,4000.);
    for(Int_t i=0; i<5; i++) AddOutput(hDGenWBF[i]);

    // bosons
    sprintf(sb,"hDGenBosons_%d", 0);  hDGenBosons[0]  = new TH1D(sb,sb,10,-0.5,9.5); 
    sprintf(sb,"hDGenBosons_%d", 1);  hDGenBosons[1]  = new TH1D(sb,sb,200,0.0,400.0); 
    sprintf(sb,"hDGenBosons_%d", 2);  hDGenBosons[2]  = new TH1D(sb,sb,100,-5.0,5.0); 
    sprintf(sb,"hDGenBosons_%d", 3);  hDGenBosons[3]  = new TH1D(sb,sb,2000,0.0,2000.0); 
    for(Int_t i=0; i<4; i++) AddOutput(hDGenBosons[i]);

    // photons
    sprintf(sb,"hDGenPhotons_%d", 0);  hDGenPhotons[0]  = new TH1D(sb,sb,10,-0.5,9.5); 
    sprintf(sb,"hDGenPhotons_%d", 1);  hDGenPhotons[1]  = new TH1D(sb,sb,200,0.0,400.0); 
    sprintf(sb,"hDGenPhotons_%d", 2);  hDGenPhotons[2]  = new TH1D(sb,sb,50,-2.5,2.5); 
    for(Int_t i=0; i<3; i++) AddOutput(hDGenPhotons[i]);
  }
}
