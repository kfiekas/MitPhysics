// $Id: GeneratorMod.cc,v 1.38 2009/04/22 11:02:55 ceballos Exp $

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
  fPrintDebug(kFALSE),
  fMCPartName(Names::gkMCPartBrn),
  fMCMETName(ModNames::gkMCMETName),
  fMCLeptonsName(ModNames::gkMCLeptonsName),
  fMCAllLeptonsName(ModNames::gkMCAllLeptonsName),
  fMCTausName(ModNames::gkMCTausName),
  fMCNeutrinosName(ModNames::gkMCNeutrinosName),
  fMCQuarksName(ModNames::gkMCQuarksName),
  fMCqqHsName(ModNames::gkMCqqHsName),
  fMCBosonsName(ModNames::gkMCBosonsName),
  fMCPhotonsName(ModNames::gkMCPhotonsName),
  fMCRadPhotonsName(ModNames::gkMCRadPhotonsName),
  fMCISRPhotonsName(ModNames::gkMCISRPhotonsName),
  fPtLeptonMin(0.0),
  fEtaLeptonMax(5.0),
  fPtPhotonMin(0.0),
  fEtaPhotonMax(5.0),
  fPtRadPhotonMin(0.0),
  fEtaRadPhotonMax(5.0),
  fPdgIdCut(0),
  fMassMinCut(-FLT_MAX),
  fMassMaxCut(FLT_MAX),
  fApplyISRFilter(kFALSE),
  fParticles(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void GeneratorMod::Process()
{
  // Process entries of the tree.

  // these arrays will be filled in the loop of particles
  MetOArr *GenMet               = new MetOArr;
  GenMet->SetName(fMCMETName);
  GenMet->SetOwner(kTRUE);
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
  MCParticleOArr *GenRadPhotons = new MCParticleOArr;
  GenRadPhotons->SetName(fMCRadPhotonsName);
  MCParticleOArr *GenISRPhotons = new MCParticleOArr;
  GenISRPhotons->SetName(fMCISRPhotonsName);

  if(fPrintDebug) printf("\n************ Next Event ************\n\n");

  // load MCParticle branch
  LoadBranch(fMCPartName);

  Double_t totalMET[3] = {0.0, 0.0, 0.0};
  Bool_t isqqH = kFALSE;
  for (UInt_t i=0; i<fParticles->GetEntries(); ++i) {
    const MCParticle *p = fParticles->At(i);

    if(fPrintDebug) p->Print("l");

    // Rad photons
    if(p->Is(MCParticle::kGamma) && p->HasMother() && p->Mother()->Status() == 3 &&
       (p->Mother()->Is(MCParticle::kEl) || p->Mother()->Is(MCParticle::kMu) ||
        p->Mother()->Is(MCParticle::kTau)) &&
       p->Pt() > fPtRadPhotonMin && p->AbsEta() < fEtaRadPhotonMax) {
      GenRadPhotons->Add(p);
    }

    // ISR photons
    if(p->Is(MCParticle::kGamma) && p->HasMother() && p->Mother()->IsQuark()) {
      GenISRPhotons->Add(p);
    }

    // MET computation at generation level
    if (p->Status() == 1 && !p->IsNeutrino()) {
      totalMET[0] = totalMET[0] + p->Px();
      totalMET[1] = totalMET[1] + p->Py();
      totalMET[2] = totalMET[2] + p->Pz();
    }

    if (!p->IsGenerated()) continue;

    // muons/electrons from W/Z decays
    if ((p->Is(MCParticle::kEl) || p->Is(MCParticle::kMu)) && p->Status() == 1) {
      if (p->Pt() > fPtLeptonMin && p->AbsEta() < fEtaLeptonMax) {
        GenAllLeptons->Add(p);
      }
      Bool_t isGoodLepton = kFALSE;
      const MCParticle *pm = p;
      while (pm->HasMother() && isGoodLepton == kFALSE) {
        if (pm->PdgId() == 92) // string reached, terminate loop
          break;
        if (pm->Mother()->Is(MCParticle::kZ)  || pm->Mother()->Is(MCParticle::kW)  ||
            pm->Mother()->Is(MCParticle::kZp) || pm->Mother()->Is(MCParticle::kWp) ||
            pm->Mother()->Is(MCParticle::kH)) {
          GenLeptons->Add(p);
          isGoodLepton = kTRUE;
          break;
        } else if (pm->Mother()->Is(MCParticle::kPi0) || pm->Mother()->Is(MCParticle::kEta)) {
          // this is fake, but it is a trick to get rid of these cases and abort the loop
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
          SendError(kWarning, "Process", "Could not find a tau neutrino!");
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
      const MCParticle *pq1 = fParticles->At(i-1);
      const MCParticle *pq2 = fParticles->At(i-2);
      if (!pq1 || !pq2) {
          SendError(kWarning, "Process", "Could not find quark pair!");
      } else if (pq1->IsQuark()   && pq2->IsQuark()   && 
                 pq1->HasMother() && pq2->HasMother() &&
                 pq1->Mother() == pq2->Mother()) {
        GenqqHs->Add(pq1);
        GenqqHs->Add(pq2);
      }

      if (p->Status() == 3)  
        GenBosons->Add(p); // take higgs boson in account here rather in next else if 
    }

    // information about bosons: W, Z, h, Z', W', H0, A0, H+
    else if (p->Status() == 3 &&
             (p->Is(MCParticle::kZ)  || p->Is(MCParticle::kW)   || p->Is(MCParticle::kH) ||
              p->Is(MCParticle::kZp) || p->Is(MCParticle::kZpp) ||
              p->Is(MCParticle::kH0) || p->Is(MCParticle::kA0)  || p->Is(MCParticle::kHp))) {
      GenBosons->Add(p);
    }

    // photons
    else if (p->Status() == 1 && p->Is(MCParticle::kGamma) &&
             p->Pt() > fPtPhotonMin && p->AbsEta() < fEtaPhotonMax) {
      GenPhotons->Add(p);
    }

    // W/Z -> lnu for Madgraph
    if (p->IsParton() && p->NDaughters() >= 2) {
      CompositeParticle *diBoson = new CompositeParticle();
      if (p->HasDaughter(MCParticle::kMu) && p->HasDaughter(MCParticle::kMuNu)) {
        diBoson->AddDaughter(p->FindDaughter(MCParticle::kMu));
        diBoson->AddDaughter(p->FindDaughter(MCParticle::kMuNu));
        if (GetFillHist()) 
          hDVMass[0]->Fill(TMath::Min(diBoson->Mass(),199.999));
        const MCParticle *tmp_mu = p->FindDaughter(MCParticle::kMu);
        while (tmp_mu->HasDaughter(MCParticle::kMu) && 
               tmp_mu->FindDaughter(MCParticle::kMu)->IsGenerated())
          tmp_mu = tmp_mu->FindDaughter(MCParticle::kMu);       
        GenLeptons->Add(tmp_mu);
      }
      else if (p->HasDaughter(MCParticle::kEl) && p->HasDaughter(MCParticle::kElNu)) {
        diBoson->AddDaughter(p->FindDaughter(MCParticle::kEl));
        diBoson->AddDaughter(p->FindDaughter(MCParticle::kElNu));
        if (GetFillHist()) 
          hDVMass[1]->Fill(TMath::Min(diBoson->Mass(),199.999));
        const MCParticle *tmp_e = p->FindDaughter(MCParticle::kEl);
        while (tmp_e->HasDaughter(MCParticle::kEl) && 
               tmp_e->FindDaughter(MCParticle::kEl)->IsGenerated())
          tmp_e = tmp_e->FindDaughter(MCParticle::kEl);       
        GenLeptons->Add(tmp_e);
      }
      else if (p->HasDaughter(MCParticle::kTau) && p->HasDaughter(MCParticle::kTauNu)) {
        diBoson->AddDaughter(p->FindDaughter(MCParticle::kTau));
        diBoson->AddDaughter(p->FindDaughter(MCParticle::kTauNu));
        if (GetFillHist()) 
          hDVMass[2]->Fill(TMath::Min(diBoson->Mass(),199.999));
        const MCParticle *tau = p->FindDaughter(MCParticle::kTau);
        if (tau->HasDaughter(MCParticle::kMu)) 
          GenLeptons->Add(tau->FindDaughter(MCParticle::kMu));
        if (tau->HasDaughter(MCParticle::kEl)) 
          GenLeptons->Add(tau->FindDaughter(MCParticle::kEl));
        if (tau->HasDaughter(MCParticle::kTau)) {
          const MCParticle *tau_second = tau->FindDaughter(MCParticle::kTau);
          if (tau_second->HasDaughter(MCParticle::kMu)) 
            GenLeptons->Add(tau_second->FindDaughter(MCParticle::kMu));
          if (tau_second->HasDaughter(MCParticle::kEl)) 
            GenLeptons->Add(tau_second->FindDaughter(MCParticle::kEl));
        }
      }
      else if (p->HasDaughter(MCParticle::kMu,kTRUE) && p->HasDaughter(-1*MCParticle::kMu,kTRUE)) {
        diBoson->AddDaughter(p->FindDaughter(MCParticle::kMu,kTRUE));
        diBoson->AddDaughter(p->FindDaughter(-1*MCParticle::kMu,kTRUE));
        if (GetFillHist()) 
          hDVMass[3]->Fill(TMath::Min(diBoson->Mass(),199.999));
        const MCParticle *tmp_mu0 = p->FindDaughter(MCParticle::kMu,kTRUE);
        while (tmp_mu0->HasDaughter(MCParticle::kMu) && 
               tmp_mu0->FindDaughter(MCParticle::kMu)->IsGenerated())
          tmp_mu0 = tmp_mu0->FindDaughter(MCParticle::kMu);       
        const MCParticle *tmp_mu1 = p->FindDaughter(-1*MCParticle::kMu,kTRUE);
        while (tmp_mu1->HasDaughter(MCParticle::kMu) && 
               tmp_mu1->FindDaughter(MCParticle::kMu)->IsGenerated())
          tmp_mu1 = tmp_mu1->FindDaughter(MCParticle::kMu);       
        GenLeptons->Add(tmp_mu0);
        GenLeptons->Add(tmp_mu1);
      }
      else if (p->HasDaughter(MCParticle::kEl,kTRUE) && p->HasDaughter(-1*MCParticle::kEl,kTRUE)) {
        diBoson->AddDaughter(p->FindDaughter(MCParticle::kEl,kTRUE));
        diBoson->AddDaughter(p->FindDaughter(-1*MCParticle::kEl,kTRUE));
        if (GetFillHist()) 
          hDVMass[4]->Fill(TMath::Min(diBoson->Mass(),199.999));
        const MCParticle *tmp_e0 = p->Daughter(0);
        while (tmp_e0->HasDaughter(MCParticle::kEl) && 
               tmp_e0->FindDaughter(MCParticle::kEl)->IsGenerated())
          tmp_e0 = tmp_e0->FindDaughter(MCParticle::kEl);       
        const MCParticle *tmp_e1 = p->Daughter(1);
        while (tmp_e1->HasDaughter(MCParticle::kEl) && 
               tmp_e1->FindDaughter(MCParticle::kEl)->IsGenerated())
          tmp_e1 = tmp_e1->FindDaughter(MCParticle::kEl);       
       GenLeptons->Add(tmp_e0);
       GenLeptons->Add(tmp_e1);
      }
      else if (p->HasDaughter(MCParticle::kTau,kTRUE) && p->HasDaughter(-1*MCParticle::kTau,kTRUE)) {
        diBoson->AddDaughter(p->FindDaughter(MCParticle::kTau,kTRUE));
        diBoson->AddDaughter(p->FindDaughter(-1*MCParticle::kTau,kTRUE));
        if (GetFillHist()) 
          hDVMass[5]->Fill(TMath::Min(diBoson->Mass(),199.999));
        const MCParticle *tau0 = p->Daughter(0);
        if (tau0->HasDaughter(MCParticle::kMu)) 
          GenLeptons->Add(tau0->FindDaughter(MCParticle::kMu));
        if (tau0->HasDaughter(MCParticle::kEl)) 
          GenLeptons->Add(tau0->FindDaughter(MCParticle::kEl));
        const MCParticle *tau1 = p->Daughter(1);
        if (tau1->HasDaughter(MCParticle::kMu)) 
          GenLeptons->Add(tau1->FindDaughter(MCParticle::kMu));
        if (tau1->HasDaughter(MCParticle::kEl)) 
          GenLeptons->Add(tau1->FindDaughter(MCParticle::kEl));
        if (tau0->HasDaughter(MCParticle::kTau)) {
          const MCParticle *tau0_second = tau0->FindDaughter(MCParticle::kTau);
          if (tau0_second->HasDaughter(MCParticle::kMu)) 
            GenLeptons->Add(tau0_second->FindDaughter(MCParticle::kMu));
          if (tau0_second->HasDaughter(MCParticle::kEl)) 
            GenLeptons->Add(tau0_second->FindDaughter(MCParticle::kEl));
        }
        if (tau1->HasDaughter(MCParticle::kTau)) {
          const MCParticle *tau1_second = tau1->FindDaughter(MCParticle::kTau);
          if (tau1_second->HasDaughter(MCParticle::kMu)) 
            GenLeptons->Add(tau1_second->FindDaughter(MCParticle::kMu));
          if (tau1_second->HasDaughter(MCParticle::kEl)) 
            GenLeptons->Add(tau1_second->FindDaughter(MCParticle::kEl));
        }        
      }
      delete diBoson;
    }

    // t -> lnu for Madgraph
    if (p->Is(MCParticle::kTop)) {
      CompositeParticle *diBoson = new CompositeParticle();
      if (p->HasDaughter(MCParticle::kMu) && p->HasDaughter(MCParticle::kMuNu)) {
        diBoson->AddDaughter(p->FindDaughter(MCParticle::kMu));
        diBoson->AddDaughter(p->FindDaughter(MCParticle::kMuNu));
        if (GetFillHist()) 
          hDVMass[6]->Fill(TMath::Min(diBoson->Mass(),199.999));
        GenLeptons->Add(p->FindDaughter(MCParticle::kMu));
      }    
      else if (p->HasDaughter(MCParticle::kEl) && p->HasDaughter(MCParticle::kElNu)) {
        diBoson->AddDaughter(p->FindDaughter(MCParticle::kEl));
        diBoson->AddDaughter(p->FindDaughter(MCParticle::kElNu));
        if (GetFillHist()) 
          hDVMass[7]->Fill(TMath::Min(diBoson->Mass(),199.999));
        GenLeptons->Add(p->FindDaughter(MCParticle::kEl));
      }    
      else if (p->HasDaughter(MCParticle::kTau) && p->HasDaughter(MCParticle::kTauNu)) {
        diBoson->AddDaughter(p->FindDaughter(MCParticle::kTau));
        diBoson->AddDaughter(p->FindDaughter(MCParticle::kTauNu));
        if (GetFillHist()) 
          hDVMass[8]->Fill(TMath::Min(diBoson->Mass(),199.999));
        const MCParticle *tau = p->FindDaughter(MCParticle::kTau);
        if (tau->HasDaughter(MCParticle::kMu)) 
          GenLeptons->Add(tau->FindDaughter(MCParticle::kMu));
        if (tau->HasDaughter(MCParticle::kEl)) 
          GenLeptons->Add(tau->FindDaughter(MCParticle::kEl));
      }
      else if (!p->HasDaughter(MCParticle::kW)) {
        for(UInt_t nd=0; nd<p->NDaughters(); ++nd)
          if (p->Daughter(nd)->IsNot(MCParticle::kBottom) &&
              p->Daughter(nd)->IsNot(MCParticle::kGamma)) 
            diBoson->AddDaughter(p->Daughter(nd));
        if (GetFillHist()) 
          hDVMass[9]->Fill(TMath::Min(diBoson->Mass(),199.999));
      }
      delete diBoson;
    }
    
    // mass cut for given pid
    if(fPdgIdCut && p->Is(fPdgIdCut) && 
       (p->Mass() < fMassMinCut || p->Mass() > fMassMaxCut)) {
      SkipEvent();
      return;
    }   
  } // end loop of particles

  Met *theMET = new Met(totalMET[0], totalMET[1]);
  theMET->SetElongitudinal(totalMET[2]);
  GenMet->AddOwned(theMET);

  // sort according to pt
  GenLeptons->Sort();
  GenAllLeptons->Sort();
  GenTaus->Sort();
  GenNeutrinos->Sort();
  GenQuarks->Sort();
  GenqqHs->Sort();
  GenBosons->Sort();
  GenPhotons->Sort();
  GenRadPhotons->Sort();
  GenISRPhotons->Sort();

  // add objects to this event for other modules to use
  AddObjThisEvt(GenMet);  
  AddObjThisEvt(GenLeptons);  
  AddObjThisEvt(GenAllLeptons);
  AddObjThisEvt(GenTaus);
  AddObjThisEvt(GenNeutrinos);
  AddObjThisEvt(GenQuarks);
  AddObjThisEvt(GenqqHs);
  AddObjThisEvt(GenBosons);
  AddObjThisEvt(GenPhotons);
  AddObjThisEvt(GenRadPhotons);
  AddObjThisEvt(GenISRPhotons);

  // fill histograms if requested
  if (GetFillHist()) {

    // MET
    hDGenMet[0]->Fill(GenMet->At(0)->Pt());
    hDGenMet[1]->Fill(GenMet->At(0)->Px());
    hDGenMet[2]->Fill(GenMet->At(0)->Py());
    hDGenMet[3]->Fill(GenMet->At(0)->Elongitudinal());

    // leptons
    hDGenLeptons[0]->Fill(GenLeptons->GetEntries());
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
    }
    // looking at events with two leptons
    if (GenLeptons->GetEntries() == 2) {
      hDGenLeptons[5]->Fill(TMath::Min(TMath::Max(TMath::Abs(GenLeptons->At(0)->Eta()),
                                                  TMath::Abs(GenLeptons->At(1)->Eta())),
                                       4.999));
      hDGenLeptons[6]->Fill(TMath::Min(TMath::Min(TMath::Abs(GenLeptons->At(0)->Eta()),
                                                  TMath::Abs(GenLeptons->At(1)->Eta())),
                                       4.999));
      if (TMath::Abs(GenLeptons->At(0)->Eta()) < 2.5 &&
          TMath::Abs(GenLeptons->At(1)->Eta()) < 2.5) {
        hDGenLeptons[7]->Fill(TMath::Min(GenLeptons->At(0)->Pt(),199.999));
        if (GenLeptons->At(0)->Pt() > 20.0) {
          hDGenLeptons[8]->Fill(TMath::Min(GenLeptons->At(1)->Pt(),199.999));
          if (GenLeptons->At(1)->Pt() > 10.0) {
            CompositeParticle *dilepton = new CompositeParticle();
            dilepton->AddDaughter(GenLeptons->At(0));
            dilepton->AddDaughter(GenLeptons->At(1));
            hDGenLeptons[9]->Fill(TMath::Min(dilepton->Mass(),999.999));
	    if(dilepton->Mass() > 12.0){
	      hDGenLeptons[10]->Fill(MathUtils::DeltaPhi(GenLeptons->At(0)->Phi(),
	                                                 GenLeptons->At(1)->Phi())
						         * 180./ TMath::Pi());
	      hDGenLeptons[11]->Fill(MathUtils::DeltaR(GenLeptons->At(0)->Eta(), GenLeptons->At(0)->Phi(),
	                                               GenLeptons->At(1)->Eta(), GenLeptons->At(1)->Phi()));
	    }
	    delete dilepton;
	  }
	}
      }
    }
    // looking at events with three leptons
    if (GenLeptons->GetEntries() == 3) {
      if (TMath::Abs(GenLeptons->At(0)->Eta()) < 2.5 &&
          TMath::Abs(GenLeptons->At(1)->Eta()) < 2.5 &&
          TMath::Abs(GenLeptons->At(2)->Eta()) < 2.5) {
        hDGenLeptons[12]->Fill(TMath::Min(GenLeptons->At(0)->Pt(),199.999));
        if (GenLeptons->At(0)->Pt() > 20.0) {
          hDGenLeptons[13]->Fill(TMath::Min(GenLeptons->At(1)->Pt(),199.999));
          hDGenLeptons[14]->Fill(TMath::Min(GenLeptons->At(2)->Pt(),199.999));
          if (GenLeptons->At(1)->Pt() > 10.0 && GenLeptons->At(2)->Pt() > 10.0) {
            CompositeParticle *dilepton01 = new CompositeParticle();
            dilepton01->AddDaughter(GenLeptons->At(0));
            dilepton01->AddDaughter(GenLeptons->At(1));
            CompositeParticle *dilepton02 = new CompositeParticle();
            dilepton02->AddDaughter(GenLeptons->At(0));
            dilepton02->AddDaughter(GenLeptons->At(2));
            CompositeParticle *dilepton12 = new CompositeParticle();
            dilepton12->AddDaughter(GenLeptons->At(1));
            dilepton12->AddDaughter(GenLeptons->At(2));
            hDGenLeptons[15]->Fill(TMath::Min(dilepton01->Mass(),999.999));
            hDGenLeptons[15]->Fill(TMath::Min(dilepton02->Mass(),999.999));
            hDGenLeptons[15]->Fill(TMath::Min(dilepton12->Mass(),999.999));
            CompositeParticle *trilepton = new CompositeParticle();
            trilepton->AddDaughter(GenLeptons->At(0));
            trilepton->AddDaughter(GenLeptons->At(1));
            trilepton->AddDaughter(GenLeptons->At(2));
            hDGenLeptons[16]->Fill(TMath::Min(trilepton->Mass(),999.999));
	    Double_t deltaR[3] = {MathUtils::DeltaR(GenLeptons->At(0)->Eta(), 
                                                    GenLeptons->At(0)->Phi(),
                                                    GenLeptons->At(1)->Eta(), 
                                                    GenLeptons->At(1)->Phi()),
                                  MathUtils::DeltaR(GenLeptons->At(0)->Eta(), 
                                                    GenLeptons->At(0)->Phi(),
                                                    GenLeptons->At(2)->Eta(), 
                                                    GenLeptons->At(2)->Phi()),
                                  MathUtils::DeltaR(GenLeptons->At(1)->Eta(), 
                                                    GenLeptons->At(1)->Phi(),
                                                    GenLeptons->At(2)->Eta(), 
                                                    GenLeptons->At(2)->Phi())};
	    Double_t deltaRMin = deltaR[0];
            for(Int_t i=1; i<3; i++) if(deltaRMin > deltaR[i]) deltaRMin = deltaR[i];
            hDGenLeptons[17]->Fill(deltaRMin);

	    delete dilepton01;
	    delete dilepton02;
	    delete dilepton12;
	    delete trilepton;
	  }
	}
      }
    }
    // looking at events with four leptons
    if (GenLeptons->GetEntries() == 4) {
      if (TMath::Abs(GenLeptons->At(0)->Eta()) < 2.5 &&
          TMath::Abs(GenLeptons->At(1)->Eta()) < 2.5 &&
          TMath::Abs(GenLeptons->At(2)->Eta()) < 2.5 &&
          TMath::Abs(GenLeptons->At(3)->Eta()) < 2.5) {
        hDGenLeptons[18]->Fill(TMath::Min(GenLeptons->At(0)->Pt(),199.999));
        if (GenLeptons->At(0)->Pt() > 20.0) {
          hDGenLeptons[19]->Fill(TMath::Min(GenLeptons->At(1)->Pt(),199.999));
          hDGenLeptons[20]->Fill(TMath::Min(GenLeptons->At(2)->Pt(),199.999));
          hDGenLeptons[21]->Fill(TMath::Min(GenLeptons->At(3)->Pt(),199.999));
          if (GenLeptons->At(1)->Pt() > 10.0 && GenLeptons->At(2)->Pt() > 10.0 &&
	      GenLeptons->At(3)->Pt() > 10.0) {
            CompositeParticle *dilepton01 = new CompositeParticle();
            dilepton01->AddDaughter(GenLeptons->At(0));
            dilepton01->AddDaughter(GenLeptons->At(1));
            CompositeParticle *dilepton02 = new CompositeParticle();
            dilepton02->AddDaughter(GenLeptons->At(0));
            dilepton02->AddDaughter(GenLeptons->At(2));
            CompositeParticle *dilepton03 = new CompositeParticle();
            dilepton03->AddDaughter(GenLeptons->At(0));
            dilepton03->AddDaughter(GenLeptons->At(3));
            CompositeParticle *dilepton12 = new CompositeParticle();
            dilepton12->AddDaughter(GenLeptons->At(1));
            dilepton12->AddDaughter(GenLeptons->At(2));
            CompositeParticle *dilepton13 = new CompositeParticle();
            dilepton13->AddDaughter(GenLeptons->At(1));
            dilepton13->AddDaughter(GenLeptons->At(3));
            CompositeParticle *dilepton23 = new CompositeParticle();
            dilepton23->AddDaughter(GenLeptons->At(2));
            dilepton23->AddDaughter(GenLeptons->At(3));
            hDGenLeptons[22]->Fill(TMath::Min(dilepton01->Mass(),999.999));
            hDGenLeptons[22]->Fill(TMath::Min(dilepton02->Mass(),999.999));
            hDGenLeptons[22]->Fill(TMath::Min(dilepton03->Mass(),999.999));
            hDGenLeptons[22]->Fill(TMath::Min(dilepton12->Mass(),999.999));
            hDGenLeptons[22]->Fill(TMath::Min(dilepton13->Mass(),999.999));
            hDGenLeptons[22]->Fill(TMath::Min(dilepton23->Mass(),999.999));
            CompositeParticle *fourlepton = new CompositeParticle();
            fourlepton->AddDaughter(GenLeptons->At(0));
            fourlepton->AddDaughter(GenLeptons->At(1));
            fourlepton->AddDaughter(GenLeptons->At(2));
            fourlepton->AddDaughter(GenLeptons->At(3));
            hDGenLeptons[23]->Fill(TMath::Min(fourlepton->Mass(),999.999));
	    Double_t deltaR[6] = {MathUtils::DeltaR(GenLeptons->At(0)->Eta(),
                                                    GenLeptons->At(0)->Phi(),
                                                    GenLeptons->At(1)->Eta(), 
                                                    GenLeptons->At(1)->Phi()),
                                  MathUtils::DeltaR(GenLeptons->At(0)->Eta(), 
                                                    GenLeptons->At(0)->Phi(),
                                                    GenLeptons->At(2)->Eta(), 
                                                    GenLeptons->At(2)->Phi()),
                                  MathUtils::DeltaR(GenLeptons->At(0)->Eta(), 
                                                    GenLeptons->At(0)->Phi(),
                                                    GenLeptons->At(3)->Eta(), 
                                                    GenLeptons->At(3)->Phi()),
                                  MathUtils::DeltaR(GenLeptons->At(1)->Eta(), 
                                                    GenLeptons->At(1)->Phi(),
                                                    GenLeptons->At(2)->Eta(), 
                                                    GenLeptons->At(2)->Phi()),
                                  MathUtils::DeltaR(GenLeptons->At(1)->Eta(), 
                                                    GenLeptons->At(1)->Phi(),
                                                    GenLeptons->At(3)->Eta(), 
                                                    GenLeptons->At(3)->Phi()),
                                  MathUtils::DeltaR(GenLeptons->At(2)->Eta(), 
                                                    GenLeptons->At(2)->Phi(),
                                                    GenLeptons->At(3)->Eta(), 
                                                    GenLeptons->At(3)->Phi())};
	    Double_t deltaRMin = deltaR[0];
            for(Int_t i=1; i<6; i++) if(deltaRMin > deltaR[i]) deltaRMin = deltaR[i];
            hDGenLeptons[24]->Fill(deltaRMin);

	    delete dilepton01;
	    delete dilepton02;
	    delete dilepton03;
	    delete dilepton12;
	    delete dilepton13;
	    delete dilepton23;
	    delete fourlepton;
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
        if (GenLeptons->GetEntries() >= 2 && 
	    GenLeptons->At(0)->Pt() > 20  &&
	    GenLeptons->At(1)->Pt() > 15) {
          if (TMath::Abs(GenLeptons->At(0)->Eta()) < 2.5 &&
              TMath::Abs(GenLeptons->At(1)->Eta()) < 2.5) {
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

    // Rad photons
    hDGenRadPhotons[0]->Fill(GenRadPhotons->GetEntries());
    for(UInt_t i=0; i<GenRadPhotons->GetEntries(); i++) {
      hDGenRadPhotons[1]->Fill(TMath::Min(GenRadPhotons->At(i)->Pt(),199.999));
      hDGenRadPhotons[2]->Fill(TMath::Min(GenRadPhotons->At(i)->AbsEta(),4.999));
      hDGenRadPhotons[3]->Fill(TMath::Min((double)GenRadPhotons->At(i)->Mother()->Status(),19.499));
      hDGenRadPhotons[4]->Fill(GenRadPhotons->At(i)->IsGenerated()+2*GenRadPhotons->At(i)->IsSimulated());
      hDGenRadPhotons[5]->Fill(TMath::Min(
                               MathUtils::DeltaR(GenRadPhotons->At(i)->Eta(), GenRadPhotons->At(i)->Phi(),
	                                         GenRadPhotons->At(i)->Mother()->Eta(), GenRadPhotons->At(i)->Mother()->Phi()),4.999));
      Int_t Mother = 0;
      if(GenRadPhotons->At(i)->Mother()->Is(MCParticle::kMu)) Mother = 1;
      hDGenRadPhotons[6]->Fill(Mother);
    }

    // ISR photons
    hDGenISRPhotons[0]->Fill(GenISRPhotons->GetEntries());
    for(UInt_t i=0; i<GenISRPhotons->GetEntries(); i++) {
      hDGenISRPhotons[1]->Fill(TMath::Min(GenISRPhotons->At(i)->Pt(),199.999));
      hDGenISRPhotons[2]->Fill(TMath::Min(GenISRPhotons->At(i)->AbsEta(),4.999));
      hDGenISRPhotons[3]->Fill(TMath::Min((double)GenISRPhotons->At(i)->Mother()->Status(),19.499));
      hDGenISRPhotons[4]->Fill(GenISRPhotons->At(i)->IsGenerated()+2*GenISRPhotons->At(i)->IsSimulated());
      hDGenISRPhotons[5]->Fill(TMath::Min(
                               MathUtils::DeltaR(GenISRPhotons->At(i)->Eta(), GenISRPhotons->At(i)->Phi(),
	                                         GenISRPhotons->At(i)->Mother()->Eta(), GenISRPhotons->At(i)->Mother()->Phi()),4.999));
    }
  }

  // Apply ISR filter (but filling all histograms)
  if(fApplyISRFilter == kTRUE && GenISRPhotons->GetEntries() > 0 &&
     GenISRPhotons->At(0)->Pt() > 15.0){
    SkipEvent();
  }
}

//--------------------------------------------------------------------------------------------------
void GeneratorMod::SlaveBegin()
{
  // Book branch and histograms if wanted.

  ReqBranch(fMCPartName, fParticles);

  // fill histograms
  if (GetFillHist()) {
    char sb[1024];
    // MET
    sprintf(sb,"hDGenMet_%d", 0);  hDGenMet[0]  = new TH1D(sb,sb,200,0,200); 
    sprintf(sb,"hDGenMet_%d", 1);  hDGenMet[1]  = new TH1D(sb,sb,400,-200,200); 
    sprintf(sb,"hDGenMet_%d", 2);  hDGenMet[2]  = new TH1D(sb,sb,400,-200,200); 
    sprintf(sb,"hDGenMet_%d", 3);  hDGenMet[3]  = new TH1D(sb,sb,400,-1000,1000); 
    for(Int_t i=0; i<4; i++) AddOutput(hDGenMet[i]);

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
    sprintf(sb,"hDGenLeptons_%d",11);  hDGenLeptons[11] = new TH1D(sb,sb,100,0.0,5.0); 
    sprintf(sb,"hDGenLeptons_%d",12);  hDGenLeptons[12] = new TH1D(sb,sb,100,0.0,200.0); 
    sprintf(sb,"hDGenLeptons_%d",13);  hDGenLeptons[13] = new TH1D(sb,sb,100,0.0,200.0); 
    sprintf(sb,"hDGenLeptons_%d",14);  hDGenLeptons[14] = new TH1D(sb,sb,100,0.0,200.0); 
    sprintf(sb,"hDGenLeptons_%d",15);  hDGenLeptons[15] = new TH1D(sb,sb,1000,0.0,1000.0); 
    sprintf(sb,"hDGenLeptons_%d",16);  hDGenLeptons[16] = new TH1D(sb,sb,1000,0.0,1000.0); 
    sprintf(sb,"hDGenLeptons_%d",17);  hDGenLeptons[17] = new TH1D(sb,sb,100,0.0,5.0); 
    sprintf(sb,"hDGenLeptons_%d",18);  hDGenLeptons[18] = new TH1D(sb,sb,100,0.0,200.0); 
    sprintf(sb,"hDGenLeptons_%d",19);  hDGenLeptons[19] = new TH1D(sb,sb,100,0.0,200.0); 
    sprintf(sb,"hDGenLeptons_%d",20);  hDGenLeptons[20] = new TH1D(sb,sb,100,0.0,200.0); 
    sprintf(sb,"hDGenLeptons_%d",21);  hDGenLeptons[21] = new TH1D(sb,sb,100,0.0,200.0); 
    sprintf(sb,"hDGenLeptons_%d",22);  hDGenLeptons[22] = new TH1D(sb,sb,1000,0.0,1000.0); 
    sprintf(sb,"hDGenLeptons_%d",23);  hDGenLeptons[23] = new TH1D(sb,sb,1000,0.0,1000.0); 
    sprintf(sb,"hDGenLeptons_%d",24);  hDGenLeptons[24] = new TH1D(sb,sb,100,0.0,5.0); 
    for(Int_t i=0; i<25; i++) AddOutput(hDGenLeptons[i]);

    // all leptons
    sprintf(sb,"hDGenAllLeptons_%d", 0);  hDGenAllLeptons[0]  = new TH1D(sb,sb,10,-0.5,9.5); 
    sprintf(sb,"hDGenAllLeptons_%d", 1);  hDGenAllLeptons[1]  = new TH1D(sb,sb,100,0.0,200.0); 
    sprintf(sb,"hDGenAllLeptons_%d", 2);  hDGenAllLeptons[2]  = new TH1D(sb,sb,100,-5.0,5.0); 
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
    sprintf(sb,"hDGenPhotons_%d", 2);  hDGenPhotons[2]  = new TH1D(sb,sb,100,-5.0,5.0); 
    for(Int_t i=0; i<3; i++) AddOutput(hDGenPhotons[i]);

    //  Rad photons
    sprintf(sb,"hDGenRadPhotons_%d", 0);  hDGenRadPhotons[0]  = new TH1D(sb,sb,10,-0.5,9.5); 
    sprintf(sb,"hDGenRadPhotons_%d", 1);  hDGenRadPhotons[1]  = new TH1D(sb,sb,400,0.0,200.0); 
    sprintf(sb,"hDGenRadPhotons_%d", 2);  hDGenRadPhotons[2]  = new TH1D(sb,sb,100,0.0,5.0); 
    sprintf(sb,"hDGenRadPhotons_%d", 3);  hDGenRadPhotons[3]  = new TH1D(sb,sb,20,-0.5,19.5); 
    sprintf(sb,"hDGenRadPhotons_%d", 4);  hDGenRadPhotons[4]  = new TH1D(sb,sb,4,-0.5,3.5); 
    sprintf(sb,"hDGenRadPhotons_%d", 5);  hDGenRadPhotons[5]  = new TH1D(sb,sb,500,0.0,5.0); 
    sprintf(sb,"hDGenRadPhotons_%d", 6);  hDGenRadPhotons[6]  = new TH1D(sb,sb,2,-0.5,1.5); 
    for(Int_t i=0; i<7; i++) AddOutput(hDGenRadPhotons[i]);

    //  ISR photons
    sprintf(sb,"hDGenISRPhotons_%d", 0);  hDGenISRPhotons[0]  = new TH1D(sb,sb,10,-0.5,9.5); 
    sprintf(sb,"hDGenISRPhotons_%d", 1);  hDGenISRPhotons[1]  = new TH1D(sb,sb,400,0.0,200.0); 
    sprintf(sb,"hDGenISRPhotons_%d", 2);  hDGenISRPhotons[2]  = new TH1D(sb,sb,100,0.0,5.0); 
    sprintf(sb,"hDGenISRPhotons_%d", 3);  hDGenISRPhotons[3]  = new TH1D(sb,sb,20,-0.5,19.5); 
    sprintf(sb,"hDGenISRPhotons_%d", 4);  hDGenISRPhotons[4]  = new TH1D(sb,sb,4,-0.5,3.5); 
    sprintf(sb,"hDGenISRPhotons_%d", 5);  hDGenISRPhotons[5]  = new TH1D(sb,sb,500,0.0,5.0); 
    for(Int_t i=0; i<6; i++) AddOutput(hDGenISRPhotons[i]);

    // auxiliar
    sprintf(sb,"hDVMass_%d", 0);  hDVMass[0]  = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDVMass_%d", 1);  hDVMass[1]  = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDVMass_%d", 2);  hDVMass[2]  = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDVMass_%d", 3);  hDVMass[3]  = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDVMass_%d", 4);  hDVMass[4]  = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDVMass_%d", 5);  hDVMass[5]  = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDVMass_%d", 6);  hDVMass[6]  = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDVMass_%d", 7);  hDVMass[7]  = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDVMass_%d", 8);  hDVMass[8]  = new TH1D(sb,sb,200,0.,200.); 
    sprintf(sb,"hDVMass_%d", 9);  hDVMass[9]  = new TH1D(sb,sb,200,0.,200.);
    for(Int_t i=0; i<10; i++) AddOutput(hDVMass[i]);
  }
}
