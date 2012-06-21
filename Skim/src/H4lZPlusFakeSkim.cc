#include "MitPhysics/Skim/interface/H4lZPlusFakeSkim.h"

using namespace mithep;

ClassImp(mithep::H4lZPlusFakeSkim)

H4lZPlusFakeSkim::H4lZPlusFakeSkim(const char *name, const char *title):
   BaseH4lSkim(name,title)
{
  // Constructor
}

//--------------------------------------------------------------------------------------------------
H4lZPlusFakeSkim::~H4lZPlusFakeSkim()
{
  // Destructor
}	

//--------------------------------------------------------------------------------------------------      
void H4lZPlusFakeSkim::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void H4lZPlusFakeSkim::SlaveBegin()
{

  ReqBranch(fMuonName,            fMuons);
  ReqBranch(fElectronName,        fElectrons);
  ReqBranch(fPrimVtxName,         fPrimVerts);
  ReqBranch(fPfCandidateName,     fPfCandidates);  
  
}

//--------------------------------------------------------------------------------------------------
void H4lZPlusFakeSkim::SlaveTerminate()
{
  BaseH4lSkim::SlaveTerminate();
}

//--------------------------------------------------------------------------------------------------
void H4lZPlusFakeSkim::Terminate()
{
}

//--------------------------------------------------------------------------------------------------
void H4lZPlusFakeSkim::BeginRun()
{
}

//--------------------------------------------------------------------------------------------------
void H4lZPlusFakeSkim::EndRun()
{
}

//--------------------------------------------------------------------------------------------------
void H4lZPlusFakeSkim::Process()
{
  fTotal++;

  LoadBranch(fMuonName);
  LoadBranch(fElectronName);
  LoadBranch(fPrimVtxName);
  LoadBranch(fPfCandidateName);

  SetBestPv();

  // select good leptons
  vector<const mithep::Muon*> goodMuons;
  vector<const mithep::Electron*> goodElectrons;
  for(UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i);
    // doing okay with no selection here at the moment
    goodMuons.push_back(mu);
  }
  for(UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    const Electron *ele = fElectrons->At(i);
    // doing okay with no selection here at the moment
    goodElectrons.push_back(ele);
  }

  // Look for z pairs
  int NZCandidates=0;
  const mithep::Muon *zmu1=0,*zmu2=0;
  const mithep::Electron *zEle1=0,*zEle2=0;
  for(unsigned imu=0; imu<goodMuons.size(); imu++) {
    TLorentzVector mu1;
    mu1.SetPtEtaPhiM(goodMuons[imu]->Pt(),goodMuons[imu]->Eta(),goodMuons[imu]->Phi(),105.658369e-3);
    for(unsigned jmu=imu+1; jmu<goodMuons.size(); jmu++) {
      TLorentzVector mu2;
      mu2.SetPtEtaPhiM(goodMuons[jmu]->Pt(),goodMuons[jmu]->Eta(),goodMuons[jmu]->Phi(),105.658369e-3);
      TLorentzVector dimu(mu1+mu2);
      if(dimu.M() < 60 || dimu.M() > 120)
	continue;
    if(goodMuons[imu]->Charge() == goodMuons[jmu]->Charge())
      continue;

    zmu1 = goodMuons[imu];
    zmu2 = goodMuons[jmu];
    NZCandidates++;
    }
  }
  for(unsigned iele=0; iele<goodElectrons.size(); iele++) {
    TLorentzVector ele1;
    ele1.SetPtEtaPhiM(goodElectrons[iele]->Pt(),goodElectrons[iele]->Eta(),goodElectrons[iele]->Phi(),105.658369e-3);
    for(unsigned jele=iele+1; jele<goodElectrons.size(); jele++) {
      TLorentzVector ele2;
      ele2.SetPtEtaPhiM(goodElectrons[jele]->Pt(),goodElectrons[jele]->Eta(),goodElectrons[jele]->Phi(),105.658369e-3);
      TLorentzVector diele(ele1+ele2);
      if(diele.M() < 60 || diele.M() > 120)
	continue;
      if(goodElectrons[iele]->Charge() == goodElectrons[jele]->Charge())
	continue;

      zEle1 = goodElectrons[iele];
      zEle2 = goodElectrons[jele];
      NZCandidates++;
    }
  }

  // Look for additional fakes
  int nFakes=0;
  for(UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i);
    if(mu==zmu1 || mu==zmu2)
      continue;
    if(mu->Pt() < 5)
      continue;

    nFakes++;
  }
  for(UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    const Electron *ele = fElectrons->At(i);
    if(ele==zEle1 || ele==zEle2)
      continue;
    if(ele->SCluster()->Et() < 5)
      continue;

    nFakes++;
  }

  if(NZCandidates==1 && nFakes>0)
    fSelected++;
  else
    SkipEvent();

}
