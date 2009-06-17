// $Id: DilepSelMod.cc,v 1.2 2009/06/15 15:00:22 loizides Exp $

#include "MitPhysics/SelMods/interface/DilepSelMod.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/CompositeParticle.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include <TH1D.h>

using namespace mithep;

ClassImp(mithep::DilepSelMod)

//--------------------------------------------------------------------------------------------------
DilepSelMod::DilepSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCleanLeptonsName(ModNames::gkMergedLeptonsName),
  fMinPt(10),
  fDilMinMass(12),
  fMinZMass(70),
  fMaxZMass(110),
  fIgnoreElCharge(kTRUE),
  fNAccCounters(0),
  fAllDiLepMass(0),
  fDiElMass(0),
  fDiMuMass(0),
  fElMuMass(0),
  fAllDiLepMassAcc(0),
  fDiElMassAcc(0),
  fDiMuMassAcc(0),
  fElMuMassAcc(0),
  fNLeptons(0),
  fNGPairs(0),
  fNZPairs(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void DilepSelMod::Process()
{
  // Process entries of the tree.

  fNAccCounters->Fill(0);

  const ParticleCol *leptons = GetObjThisEvt<ParticleCol>(fCleanLeptonsName);
  if (!leptons) {
    SkipEvent();
    return;
  }

  fNAccCounters->Fill(1);

  // make sure have found at least 2 leptons
  if (leptons->GetEntries()<2) {
    SkipEvent();
    return;
  }

  UInt_t nLeps = leptons->GetEntries();
  UInt_t nZPairs    = 0;
  UInt_t nGoodPairs = 0;

  for (UInt_t i=0; i<nLeps; ++i) {
    const Particle *li = leptons->At(i);

    if (li->Pt()<fMinPt)
      continue;

    for (UInt_t j=0; j<i; ++j) {
      const Particle *lj = leptons->At(j);

      if (lj->Pt()<fMinPt)
        continue;

      CompositeParticle dil;
      dil.AddDaughter(li);
      dil.AddDaughter(lj);
      Double_t mass = dil.Mass();
      if (mass<fDilMinMass)
        continue;

      fAllDiLepMass->Fill(mass);

      if (li->ObjType()!=lj->ObjType()) {
        fElMuMass->Fill(mass);
        ++nGoodPairs;
        continue;
      }

      if (li->Is(kMuon)) {
        if (li->Charge()!=lj->Charge()) {
          fDiMuMass->Fill(mass);
          if ((mass>fMinZMass) && (mass<fMaxZMass)) {
            ++nZPairs;
            continue;
          }
        }
        ++nGoodPairs;
        continue;
      }

      if (li->Is(kElectron)) {
        if (fIgnoreElCharge || (li->Charge()!=lj->Charge())) {
          fDiElMass->Fill(mass);
          if ((mass>fMinZMass) && (mass<fMaxZMass)) {
            ++nZPairs;
            continue;
          }
        }
        ++nGoodPairs;
        continue;
      }
    }
  }

  fNLeptons->Fill(nLeps);
  fNGPairs->Fill(nGoodPairs);
  fNZPairs->Fill(nZPairs);
  fNAccCounters->Fill(2);

  // cut on number of Z pairs
  if (nZPairs>=1) {
    SkipEvent();
    return;
  }

  fNAccCounters->Fill(3);

  // cut on number of good pairs
  if (nGoodPairs<1) {
    SkipEvent();
    return;
  }

  fNAccCounters->Fill(4);
  for (UInt_t i=0; i<nLeps; ++i) {
    const Particle *li = leptons->At(i);

    if (li->Pt()<fMinPt)
      continue;

    for (UInt_t j=0; j<i; ++j) {
      const Particle *lj = leptons->At(j);

      if (lj->Pt()<fMinPt)
        continue;

      CompositeParticle dil;
      dil.AddDaughter(li);
      dil.AddDaughter(lj);
      Double_t mass = dil.Mass();
      if (mass<fDilMinMass)
        continue;

      fAllDiLepMassAcc->Fill(mass);

      if (li->ObjType()!=lj->ObjType()) {
        fElMuMassAcc->Fill(mass);
        continue;
      }

      if (li->Is(kMuon)) {
        if (li->Charge()!=lj->Charge()) {
          fDiMuMassAcc->Fill(mass);
        continue;
        }
      }

      if (li->Is(kElectron)) {
        if (fIgnoreElCharge || (li->Charge()!=lj->Charge())) {
          fDiElMassAcc->Fill(mass);
        }
        continue;
      }
    }
  }
}

//--------------------------------------------------------------------------------------------------
void DilepSelMod::SlaveBegin()
{
  // Create and add histograms to the output list.

  AddTH1(fNAccCounters,"hNAccCounters",";cut;#",6,-0.5,5.5);
  if (1) {
    TAxis *xa = fNAccCounters->GetXaxis();
    for(Int_t i=1;i<=fNAccCounters->GetNbinsX();++i)
      xa->SetBinLabel(i,"unused");
    xa->SetBinLabel(1,"Enter");
    xa->SetBinLabel(2,"Objs");
    xa->SetBinLabel(3,"2Lep");
    xa->SetBinLabel(4,"ZPair");
    xa->SetBinLabel(5,"GPair");
    xa->SetRangeUser(0,4);
  }
  AddTH1(fAllDiLepMass,"hAllDiLepMass",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fDiElMass,"hDiElMass",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fDiMuMass,"hDiMuMass",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fElMuMass,"hElMuMass",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fAllDiLepMassAcc,"hAllDiLepMassAcc",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fDiElMassAcc,"hDiElMassAcc",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fDiMuMassAcc,"hDiMuMassAcc",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fElMuMassAcc,"hElMuMassAcc",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fNLeptons,"hNLeptons",";leptons;#",10,-0.5,9.5);
  AddTH1(fNGPairs,"hNGoodPairs",";leptons;#",10,-0.5,9.5);
  AddTH1(fNZPairs,"hNZPairs",";leptons;#",10,-0.5,9.5);
}
